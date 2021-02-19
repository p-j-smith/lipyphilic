# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# lipyphilic --- lipyphilic.readthedocs.io
#
# Released under the GNU Public Licence, v2 or any higher version
#

"""Neighbours --- :mod:`lipyphilic.lib.neighbours`
================================================

:Author: Paul Smith
:Year: 2021
:Copyright: GNU Public License v2

This module provides methods for finding neighbouring lipids in a bilayer.

Two lipids are considered neighbours if they have any atoms within a given
cutoff distance of one another.

Input
------

Required:
  - *universe* : an MDAnalysis Universe object.
  - *lipid_sel* : atom selection for lipids in the bilayer

Optional:
  - *cutoff* : lipids are considered to be neighbouring if they have at least one pair of atoms less than this distance apart (in Å)
  

Output
------

  - *neighbouring* : a binary variable, equal to 1 if two lipids are in contact, and 0 otherwise

For efficient use of memory, an adjacency matrix of neighbouring lipids is stored
in a :class:`scipy.sparse.csc_matrix` sparse matrix for each frame of the analysis. The data
are stored in the :attr:`neighbours.neighbours` attribute. The matrix has shape
(n_residues, n_residues * n_frames), and the data for frame *n* is accessed via::

  neighbours.neighbours[:, (n * n_residues):((n + 1) * n_residues)]

A matrix element *[i, j]* is equal to 1 if lipid *i* neighbours lipid *j* and equal to 0 otherwise.
The matrix is symmetric: if lipid *i* neighbours lipid *j* then *j* must neighbour *i*.


Example usage of :class:`Neighbours`
--------------------------------------

An MDAnalysis Universe must first be created before using :class:`Neighbours`::

  import MDAnalysis as mda
  from lipyphilic.lib.neighbours import Neighbours

  u = mda.Universe(tpr, trajectory)

We can now create our :class:`Neighbours` object::

  neighbours = Neighbours(
      universe=u,
      lipid_sel="name GL1 GL2 ROH",  # assuming we're using the MARTINI forcefield
      cutoff=12.0
  )
  
A lipid will be considered to be neighbouring a cholesterol molecule if either its *GL1* or *GL2* bead
is within *12* Å of the ROH bead of the cholesterol. For neighbouring lipids, the distances
between there respective *GL1* and "GL2* beads will be considered.
  
We then select which frames of the trajectory to analyse (`None` will use every
frame) and select to display a progress bar (`verbose=True`)::
  
  neighbours.run(
    start=None,
    stop=None,
    step=None,
    verbose=True
  )
  
The results are then available in the :attr:`neighbours.Neighbours` attribute as a
:class:`scipy.sparse.csc_matrix`.

The class and its methods
-------------------------

.. autoclass:: Neighbours
    :members:

"""
from tqdm.auto import tqdm
import numpy as np
import scipy.sparse
import pandas as pd
from MDAnalysis.lib.distances import capped_distance

from lipyphilic.lib import base


class Neighbours(base.AnalysisBase):
    """Find neighbouring lipids in a bilayer.
    """

    def __init__(self, universe,
                 lipid_sel,
                 cutoff=10.0
                 ):
        """Set up parameters for finding neighbouring lipids.
        
        Parameters
        ----------
        universe : Universe
            MDAnalysis Universe object
        lipid_sel : str
            Selection string for lipids in the bilayer.
        cutoff : float
            To be considered neighbours, two lipids must have at least one pair of atoms within
            this cutoff distance (in Å)
            
        Tip
        ---
        
        The resultant sparse matrix can be used to calculate the number of each lipid species
        neighbouring each individual lipid at each frame using :func:`lipyphilic.lib.neighbours.count_neighbours`
        
        """
        self.u = universe
        self._trajectory = self.u.trajectory
        self.membrane = self.u.select_atoms(lipid_sel, updating=False)
        
        if cutoff <= 0:
            raise ValueError("'cutoff' must be greater than 0")
        
        self.cutoff = cutoff
           
        self.neighbours = None
        
    def _prepare(self):
        
        # Output array
        self.neighbours = scipy.sparse.lil_matrix(
            (self.membrane.n_residues, self.membrane.n_residues * self.n_frames),
            dtype=np.int8
        )
        
    def _single_frame(self):
        
        pairs = capped_distance(
            self.membrane.positions,
            self.membrane.positions,
            max_cutoff=self.cutoff,
            box=self._ts.dimensions,
            return_distances=False
        )
        self.pairs = pairs
        
        # Find unique pairs of residues interacting
        # Currently we have pairs of atoms
        ref, neigh = np.unique(self.membrane.resindices[pairs], axis=0).T
        
        # Dont keep self-interactions between lipids
        different = ref != neigh
        ref = ref[different]
        neigh = neigh[different]
        
        # store neighbours for this frame
        frame_start = self._frame_index * self.membrane.n_residues
        self.neighbours[ref, neigh + frame_start] = 1
        self.neighbours[neigh, ref + frame_start] = 1

    def _conclude(self):
        
        # csc sparse matrices are faster for slicing columns
        # we'll want this for selecting neighbours at given frames
        self.neighbours = self.neighbours.tocsc()

    def count_neighbours(self, fractional_enrichment=False, count_by=None):
        """Count the number of each neighbour type at each frame.

        Parameters
        ----------
        fractional_enrichment : bool
            If true, will also calculate the fractional enrichment of each lipid at each frame,
            defined as [A_local] / [A_bulk], where [A_local] is the number of lipids of a given
            species neighbouring a reference lipid and [A_bulk] is the molar concentration
            of this lipid species in the bilayer. Defaults to None, in which case the fractional
            enrichment is not calculated.
        count_by : numpy.ndarray
            An array containing ordinal data describing each lipid at each frame. For example,
            it may be an array containing information on the ordered state or each lipid.
            Defaults to None, in which case the lipid species
            (resnames) are used for counting neighbours.
        count_by_labels : list
            A dictionary of labels describing what each unique value in `count_by` refers to, e.g
            if `count_by` contains information on the ordered state of each lipid at each frame, whereby
            0 corresponds to disordered and 1 corresponds to disordered, then
            count_by_labels = {0: 'Ld', 1: 'Lo'}. There must be precisely one label for each unique
            value in 'count_by'.
        
        Returns
        -------
        
        counts : pandas.DataFrame
            A DataFrame containing the following data for each lipid at each frame: lipid identifier,
            frame number, number of neighbour of each species (or of each type in :attr:'count_by' if)
            this is provided, as well as the total number of neighbours. If :attr:'fractional_enrichment'
            is true
            
        Note
        ----
        
        Neighbours must be found by using `Neighbours.run()` before calling this method.
        """
        
        # create output array
        all_resnames = self.membrane.resnames
        unique_resnames = np.unique(all_resnames)
        all_counts = np.full(
            (self.membrane.n_residues, self.n_frames, unique_resnames.size),
            fill_value=0,
            dtype=np.uint8  # count can't be negative, and no lipid will have more than 255 neighbours
        )
        
        if fractional_enrichment:
            all_fe = np.full(
                (self.membrane.n_residues, self.n_frames, unique_resnames.size),
                fill_value=0,
                dtype=float
            )
        
        # For counts we need to know which column of the output array to add counts to for each lipid type
        type_index = {resname: i for i, resname in enumerate(unique_resnames)}
        
        # Get counts at each frame
        n_residues = self.membrane.n_residues
        for frame_index in tqdm(np.arange(self.n_frames)):
        
            ref, neigh = self.neighbours[:, frame_index * n_residues:(frame_index + 1) * n_residues].nonzero()
            
            unique, counts = np.unique([ref, [type_index[t] for t in all_resnames[neigh]]], axis=1, return_counts=True)
            r, t = unique  # reference index (r) and type index (t)
            all_counts[r, frame_index, t] = counts

        # Assemble data for the DataFrame
        all_counts = all_counts.reshape(n_residues * self.n_frames, unique_resnames.size)
        total_counts = np.sum(all_counts, axis=1)

        frames = np.full((n_residues, self.n_frames), fill_value=self.frames)
        frames = frames.reshape(n_residues * self.n_frames)

        resnames = np.full((n_residues, self.n_frames), fill_value=self.membrane.residues.resnames[:, np.newaxis])
        resnames = resnames.reshape(n_residues * self.n_frames)
        
        data = np.concatenate((resnames[:, np.newaxis], frames[:, np.newaxis], all_counts, total_counts[:, np.newaxis]), axis=1)
        
        # Create DataFrame
        columns = ["Resname", "Frame"] + ["n" + lipid for lipid in unique_resnames] + ["Total"]
        df = pd.DataFrame(
            data=data,
            columns=columns
        )
        
        return df
