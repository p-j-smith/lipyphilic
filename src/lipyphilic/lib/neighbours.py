# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# lipyphilic --- lipyphilic.readthedocs.io
#
# Released under the GNU Public Licence, v2 or any higher version
#

"""Neighbours --- :mod:`lipyphilic.lib.neighbours`
==================================================

:Author: Paul Smith
:Year: 2021
:Copyright: GNU Public License v2

This module provides methods for finding neighbouring lipids in a bilayer,
calculating local lipid compositions and finding the largest cluster of
lipids over time.

Two lipids are considered neighbours if they have any atoms within a given
cutoff distance of one another.

Input
-----

Required:
  - *universe* : an MDAnalysis Universe object.
  - *lipid_sel* : atom selection for lipids in the bilayer

Optional:
  - *cutoff* : lipids are considered to be neighbouring if they have at least one pair of atoms less than this distance apart (in Å)
  

Output
------

  - *neighbours* : a sparse matrix of binary variables, equal to 1 if two lipids are in contact, and 0 otherwise

For efficient use of memory, an adjacency matrix of neighbouring lipids is stored
in a :class:`scipy.sparse.csc_matrix` sparse matrix for each frame of the analysis. The data
are stored in the :attr:`neighbours.neighbours` attribute. The matrix has shape
(n_residues, n_residues * n_frames), and the data for frame *n* is accessed via::

  neighbours.neighbours[:, (n * n_residues):((n + 1) * n_residues)]

A matrix element *[i, j]* is equal to 1 if lipid *i* neighbours lipid *j* and equal to 0 otherwise.
The matrix is symmetric: if lipid *i* neighbours lipid *j* then *j* must neighbour *i*.


Tip
---

The resultant sparse matrix can be used to calculate the loca lipid composition of each individual lipid
at each frame using :func:`lipyphilic.lib.neighbours.count_neighbours`, or to find the largest cluster of
lipids at each frame using :func:`lipyphilic.lib.neighbours.largest_cluster`.


Example usage of :class:`Neighbours`
------------------------------------

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

Counting the number of neighbours: by lipid species
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to compute the number of each lipid species around each lipid at each frame,
after generating the neighbour matrix we can use the :func:`count_neighbours`
method::

  counts = neighbours.count_neighbours()

*Counts* is a :class:`pandas.DataFrame` in which each row contains the following
information (if there are N distinct species in the membrane)::

    [
        <lipid identifier>,  # by default, the lipid resname
        <lipid resindex>,
        <frame>,
        <num species_1 neighbours>,
        ...
        <num species_N neighbours>,
        <total num neighbours>
    ]

Counting the number of neighbours: by user-defined labels
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Instead of using the lipid resname to identify neighbouring lipids, any ordinal data may
be used for counting lipid neighbours. This is done through use of the :attr:`count_by` and
:attr:`count_by_labels` parameters::

  counts = neighbours.count_neighbours(
    count_by=lipid_order_data,
    count_by_labels={'Ld': 0, 'Lo': 1}
  )

Here we assume that 'lipid_order_data' contains information on whether each lipid is in
the liquid-disordered phase or the liquid-ordered phase at each frame. It must take
the shape '(n_residues, n_frames)', and in this example 'lipid_order_data[i, j]' would
be equal to '0' if lipid 'i' is liquid-disordered at frame 'j' and equal to '1' if it is
liquid-ordered. 'count_by_labels' is used to signify that the value '0' corresponds to
the liquid-disordered (Ld) phase and the value '1' to the liquid-ordered  (Lo) phase. In
this example, the returned :class:`pandas.DataFrame` would contain the following information
in each row::
  
    [
        <Ld or Lo>,
        <lipid resindex>,
        <frame>,
        <num Ld neighbours>,
        <num Lo neighbours>,
        <total num neighbours>
    ]


Find the largest cluster
^^^^^^^^^^^^^^^^^^^^^^^^

To find the largest cluster of a set of lipid species we can use the :func:`largest_cluster`
method::

  largest_cluster = neighbours.largest_cluster(
    cluster_sel="resname CHOL DPPC"
  )
  
The results are returned in a :class:`numpy.ndarray` and contain the number of lipids in the largest
cluster at each frame.

  
Find the largest cluster in a given leaflet
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The previous example will compute the largest cluster formed by cholesterol and DPPC molecules at each
frame. In large coarse-grained systems where there is substantial flip-flop of sterols, this cluster may
span both leaflets. In order to find the largest cluster at each frame within a given leaflet, we can
tell :func:`largest_cluster` to consider only lipids in the upper leaflet by using the
`filter_by` parameter.

First, though, we need to know which leaflet each lipid is in at each frame. This may be done using
:class:`lipyphilic.lib.assign_leaflets.AssignLeaflets`::

  leaflets = AssignLeaflets(
    universe=u,
    lipid_sel="name GL1 GL2 ROH"  # pass the same selection that was passed to Neighbours
  )
  leaflets.run()  # run the analysis on the same frames as Neighbours.run()
  
The leaflets data are stored in the :attr:`leaflets.leaflets` attribute, will be equal to '1' if the
lipid is in the upper leaflet at a given frame and equal to '-1' if it is in the lower leaflet. See
:class:`lipyphilic.lib.assign_leaflets.AssignLeaflets` for more information. We can now find the
largest cluster over time in the upper (1) leaflet.

The :attr:`filter_by` parameter takes as input a 2D :class:`numpy.ndarray` of shape
(n_residues, n_frames). The array should be a `boolean mask
<https://docs.scipy.org/doc/numpy-1.15.0/user/basics.indexing.html#boolean-or-mask-index-arrays>`__,
where `True` indicates that we should include this lipid in the neighbour calculation::

  upper_leaflet_mask = leaflet.leaflets == 1

  largest_cluster_upper_leaflet = neighbours.largest_cluster(
    cluster_sel="resname CHOL DPPC",
    leaflets=upper_leaflet_mask
  )

Now, lipids either in the lower leaflet (-1) or the midplane (0) will not be included when determining
the largest cluster.

Get residue indices of lipids in the largest cluster
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If we want to know not just the cluster size but also which lipids are in the largest cluster at each
frame, we can set the :attr:`return_indices` parameter to `True`::

  largest_cluster, largest_cluster_indices = neighbours.largest_cluster(
    cluster_sel="resname CHOL DPPC",
    return_indices=True
  )

The residue indices will be returned as list of `numpy.ndarray` arrays - one per frame of the analysis. Each
array contains the residue indices of the lipids in the largest cluster at that frame


The class and its methods
-------------------------

.. autoclass:: Neighbours
    :members:

"""
from tqdm.auto import tqdm
import numpy as np
import scipy.stats
import scipy.sparse
import pandas as pd
from MDAnalysis.lib.distances import capped_distance
from MDAnalysis.exceptions import NoDataError

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
        cutoff : float, optional
            To be considered neighbours, two lipids must have at least one pair of atoms within
            this cutoff distance (in Å). The default is `10.0`.
        
        """
        super(Neighbours, self).__init__(universe.trajectory)
        
        self.u = universe
        self.membrane = self.u.select_atoms(lipid_sel, updating=False)
        
        # to allow for non-sequential resindices
        self._sorted_membrane_resindices = scipy.stats.rankdata(
            self.membrane.resindices,
            method="dense"
        ) - 1
        
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
        ref, neigh = np.unique(self._sorted_membrane_resindices[pairs], axis=0).T
        
        # Dont keep self-interactions between lipids
        different = ref != neigh
        ref = ref[different]
        neigh = neigh[different]
        
        # store neighbours for this frame
        frame_start = self._frame_index * self.membrane.n_residues
        self.neighbours[ref, neigh + frame_start] = 1
        self.neighbours[neigh, ref + frame_start] = 1  # the neighbour matrix must be symmetric

    def _conclude(self):
        
        # csc sparse matrices are faster for slicing columns
        # we'll want this for selecting neighbours at given frames
        self.neighbours = self.neighbours.tocsc()

    def count_neighbours(self, count_by=None, count_by_labels=None):
        """Count the number of each neighbour type at each frame.

        Parameters
        ----------
        count_by : numpy.ndarray, optional
            An array containing ordinal data describing each lipid at each frame. For example,
            it may be an array containing information on the ordered state or each lipid.
            Defaults to None, in which case the lipid species (resnames) are used for counting neighbours.
        count_by_labels : dict, optional
            A dictionary of labels describing what each unique value in `count_by` refers to, e.g
            if `count_by` contains information on the ordered state of each lipid at each frame, whereby
            0 corresponds to disordered and 1 corresponds to ordered, then
            `count_by_labels = {'Ld': 0, 'Lo': 1}`. There **must** be precisely one label for each unique
            value in 'count_by'. If `count_by` is given but `count_by_labels` is left as `None`, the values
            in `count_by` will be used as the labels.
        
        Returns
        -------
        
        counts : pandas.DataFrame
            A DataFrame containing the following data for each lipid at each frame: lipid identifier
            (default is resname), lipid residue index, frame number, number of neighbours of each species
            (or of each type in 'count_by' if this is provided), as well as the total number of neighbours.
        
        """
        
        if self.neighbours is None:
            raise NoDataError(".neighbours attribute is None: use .run() before calling .count_neighbours()")
        
        # create output array
        if count_by is None:
            
            # Use lipid resnames to distinguish lipids
            count_by = np.full(
                (self.membrane.n_residues, self.n_frames),
                fill_value=self.membrane.residues.resnames[:, np.newaxis],
            )
            count_by_labels = {label: index for index, label in enumerate(np.unique(self.membrane.resnames))}
        
        elif count_by_labels is None:
            
            # Use values in 'count_by' as the labels
            count_by_labels = {label: index for index, label in enumerate(np.unique(count_by))}
            
        else:
            
            # the ordinal values in 'count_by' now take on the string labels supplied
            max_label_size = max([len(label) for label in count_by_labels])
            new_count_by = np.full_like(count_by, dtype=f'<U{max_label_size}', fill_value="")
            for label in count_by_labels:
                new_count_by[count_by == count_by_labels[label]] = label
            count_by = new_count_by
            del new_count_by
        
        # create output array
        all_counts = np.full(
            (self.membrane.n_residues, self.n_frames, len(count_by_labels)),
            fill_value=0,
            dtype=np.uint8  # count can't be negative, and no lipid will have more than 255 neighbours
        )
        
        # For counts we need to know which column of the output array to add counts to for each lipid type
        type_index = {value: index for index, value in enumerate(count_by_labels)}
        
        # Get counts at each frame
        n_residues = self.membrane.n_residues
        for frame_index in tqdm(np.arange(self.n_frames)):
        
            ref, neigh = self.neighbours[:, frame_index * n_residues:(frame_index + 1) * n_residues].nonzero()
            unique, counts = np.unique([ref, [type_index[t] for t in count_by[neigh, frame_index]]], axis=1, return_counts=True)
            
            r, t = unique  # reference index (r) and type index (t)
            all_counts[r, frame_index, t] = counts

        # Assemble data for the DataFrame
        labels = np.array([list(count_by_labels)[type_index[frame_index]] for lipid in count_by for frame_index in lipid])
        
        resindices = np.full((n_residues, self.n_frames), fill_value=self.membrane.residues.resindices[:, np.newaxis])
        resindices = resindices.reshape(n_residues * self.n_frames)
        
        frames = np.full((n_residues, self.n_frames), fill_value=self.frames)
        frames = frames.reshape(n_residues * self.n_frames)

        all_counts = all_counts.reshape(n_residues * self.n_frames, len(count_by_labels))
        total_counts = np.sum(all_counts, axis=1)
        
        data = np.concatenate(
            (labels[:, np.newaxis], resindices[:, np.newaxis], frames[:, np.newaxis], all_counts, total_counts[:, np.newaxis]),
            axis=1
            )
        
        # Create DataFrame
        columns = ["Label", "Resindex", "Frame"] + [f"n{label}" for label in count_by_labels] + ["Total"]
        df = pd.DataFrame(
            data=data,
            columns=columns
        )
        
        # make every column except the label take on integer values
        for column in df.columns[1:]:
            df[column] = pd.to_numeric(df[column])
        
        return df
    
    def largest_cluster(self, cluster_sel=None, filter_by=None, return_indices=False):
        """Find the largest cluster of lipids at each frame.
        
        Parameters
        ----------
        cluster_sel : str, optional
            Selection string for lipids to include in the cluster analysis. The default is `None`, in
            which case all lipid used in identiying neighbouring lipids will be used for finding
            the largest cluster.
        filter_by : numpy.ndarray, optional
            A boolean array indicating whether or not to include each lipid in the cluster analysis. If
            the array is 1D and of shape (n_lipids), the same lipids will be used in the cluster
            analysis at every frame. If the array is 2D and of shape (n_lipids, n_frames), the boolean
            value of each lipid at each frame will be taken into account. The default is `None`, in which
            case all lipids used in identiying neighbours will be used for finding
            the largest cluster.
        return_indices : bool, optional
            If `True`, a list of NumPy arrays will also be returned, on for each frame. Each NumPy array
            will contain the residue indices of the lipids in the largest cluster at that frame. Note, if
            there are two largest clusters of equal size, only the residue indices of lipids in one
            cluster will be returned (the cluster that has the lipid with the smallest residue index). The
            default is `False`, in which case no reidue indices are returned.
            
        
        Returns
        -------
        
        largest_cluster : numpy.ndarray
            An array containing the number of lipids in the largest cluster at each frame.
            
        Note
        ----
        
        Neighbours must be found by using `Neighbours.run()` before calling either
        `Neighbours.count_neighbours()` or `Neighbours.largest_cluster()`.
        
        """

        if self.neighbours is None:
            raise NoDataError(".neighbours attribute is None: use .run() before calling .largest_cluster()")
        
        if filter_by is not None and np.array(filter_by).ndim not in [1, 2]:
            raise ValueError("'filter_by' must either be a 1D array containing non-changing boolean"
                             "values for each lipid, or a 2D array of shape (n_residues, n_frames)"
                             " containing a boolean value for each lipid at each frame."
                             )

        elif filter_by is not None and len(filter_by) != self.membrane.n_residues:
            raise ValueError("The shape of 'filter_by' must be (n_residues,)")
        
        # determine which lipids to use in the analysis at each frame
        if filter_by is None:
            
            filter_by = np.full(
                (self.membrane.n_residues, self.n_frames),
                fill_value=True,
                dtype=bool
            )
        elif filter_by.ndim == 1:
            
            filter_by = np.full(
                (self.membrane.n_residues, self.n_frames),
                fill_value=filter_by[:, np.newaxis],
                dtype=bool
            )
            
        # also create mask based on `cluster_sel`
        if cluster_sel is None:
            
            filter_lipids = np.full(
                self.membrane.n_residues,
                fill_value=True,
                dtype=bool
            )
        else:
            
            lipids = self.u.select_atoms(cluster_sel).residues
            
            if lipids.n_residues == 0:
                raise ValueError(
                    "'cluster_sel' produces atom empty AtomGroup. Please check the selection string."
                )
            
            filter_lipids = np.in1d(
                self.membrane.residues.resindices,
                lipids.resindices
            )
            
        # combine the masks
        filter_by[filter_lipids == False] = False  # noqa: E712
                
        # output arrays
        largest_cluster = np.zeros(self.n_frames, dtype=int)
        largest_cluster_resindices = np.full(self.n_frames, fill_value=0, dtype=object)

        n_residues = self.membrane.n_residues
        for frame_index, frame in enumerate(self.frames):
            
            frame_filter = filter_by[:, frame_index]
            
            frame_neighbours = self.neighbours[:, frame_index * n_residues:(frame_index + 1) * n_residues]
            frame_neighbours = frame_neighbours[frame_filter][:, frame_filter].toarray()
            
            # find all connected components
            _, com_labels = scipy.sparse.csgraph.connected_components(frame_neighbours)
            
            unique_com_labels, counts = np.unique(com_labels, return_counts=True)
            largest_label = unique_com_labels[np.argmax(counts)]
            
            # largest cluster and resindices of lipids in the cluster
            largest_cluster[frame_index] = max(counts)
            
            frame_resindices = self.membrane.residues.resindices[frame_filter]
            largest_cluster_resindices[frame_index] = frame_resindices[com_labels == largest_label]
            
        if return_indices is True:
            return largest_cluster, largest_cluster_resindices
        else:
            return largest_cluster
