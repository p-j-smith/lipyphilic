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
cutoff of one another.

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

An MDAnalysis Universe must first be created before using :class:`FlipFlop`::

  import MDAnalysis as mda
  from lipyphilic.lib.neighbours import Neighbours

  u = mda.Universe(tpr, trajectory)

We can now create our :class:`Neighbours` object::

  neighbours = Neighbours(
      universe=u,
      lipid_sel="name GL1 GL2 ROH",  # assuming we're using the MARTINI forcefield
      cutoff=12.0
  )
  
A lipid will be considered a cholesterol molecule if either its *GL1* or *GL2* bead
is within *12* Å of the ROH bead of the cholesterol. For neighbouring lipids, the distances
between there respective *GL1* and "GL2* beads will be considered.
  
We then select which frames of the trajectory to analyse (`None` will use every
frame)::
  
  neighbours.run(
    start=None,
    stop=None,
    step=None
  )
  
The results are then available in the :attr:`neighbours.Neighbours` attribute as a
:class:`scipy.sparse.csc_matrix`.

The class and its methods
-------------------------

.. autoclass:: Neighbours
    :members:

"""
import numpy as np
import scipy.sparse
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
        neighbouring each individual lipid at each frame using :class:`lipyphilic.lib.neighbours.NumNeighbours`
        
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
        
        # Dont keep self-interactions between lipid
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
