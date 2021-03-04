# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# lipyphilic --- lipyphilic.readthedocs.io
#
# Released under the GNU Public Licence, v2 or any higher version
#

r"""Lipid order parameter --- :mod:`lipyphilic.lib.order_parameter`
===================================================================

:Author: Paul Smith
:Year: 2021
:Copyright: GNU Public License v2

This module provides methods for calculating the orientational order parameter
of lipid tails in a bilayer.

Coarse-grained order parameter
------------------------------

The class :class:`liyphilic.lib.order_parameter.SCC` calculates the coarse=grained
order parameter, as defined in
`Seo et al. (2020) <https://pubs.acs.org/doi/full/10.1021/acs.jpclett.0c01317>`__.
The coarse-grained order parameter, :math:`S_{CC}`, is defined as:

.. math::

  S_{CC} = \displaystyle \frac{\big \langle 3 \cos^2 \theta - 1 \big \rangle}{2}

where :math:`\theta` is the angle between the membrane normal and the vector connecting
two consecutive tail beads. Angular brackets denote averages over all beads
in an acyl tail.

See `Piggot et al. (2017) <https://pubs.acs.org/doi/full/10.1021/acs.jctc.7b00643>`__
for an excellent discussion on calculating acyl tali order parameters in molecular
dynamics simulations.

Input
-----

Required:
  - *universe* : an MDAnalysis Universe object
  - *tail_sel* : atom selection for beads in the acyl tail

Options:
  - *normals* : local membrane normals for each tail at each frame

Output
------

  - *SCC* : order parameter of each tail at each frame

  
The :math:`S_{CC}` order parameter data are returned in a :class:`numpy.ndarray`, where each
row corresponds to an individual lipid and each column corresponds to an individual frame.

Warning
-------

`tail_sel` should select beads in **either** the *sn1* **or** *sn2* tails, not both tails.


Example usage of :class:`Scc`
-----------------------------

An MDAnalysis Universe must first be created before using `SCC`::

  import MDAnalysis as mda
  from lipyphilic.lib.order_parameter import SCC

  u = mda.Universe(tpr, trajectory)

If we have used the MARTINI forcefield to study phospholipid/cholesterol mixture,
we can calculate the order parameter of the *sn1* of tails of DPPC and DOPC as follows::

  scc_sn1 = SCC(
    universe=u,
    tail_sel="name ??A"  # selects C1A, C2A, D2A, C3A, and C4A
  )
  
This will calculate :math:`S_{CC}` of each DOPC and DPPC sn1 tail.

We then select which frames of the trajectory to analyse (`None` will use every
frame) and choose to display a progress bar (`verbose=True`)::
  
  scc_sn1.run(
    start=None,
    stop=None,
    step=None,
    verbose=True
  )
  
The results are then available in the :attr:`scc_sn1.SCC` attribute as a
:class:`numpy.ndarray`. The array has the shape (n_residues, n_frames). Each row
corresponds to an individual lipid and each column to an individual frame.

Likewise, to calculate the :math:`S_{CC}` of the *sn2* tails, we can do::

  scc_sn2 = SCC(
    universe=u,
    tail_sel="name ??B"  # selects C1B, C2B, D2B, C3B, and C4B
  )
  scc_sn2.run(verbose=True)

And the get a weighted-average :math:`S_{CC}` we can do::

  SCC.weighted_average(scc_sn1, scc_sn2)
  
which will take into account the number of beads in each tail and return
an weighted-average :math:`S_{CC}` for each lipid at each frame.

Local membrane normals
----------------------

By default, the :math:`S_{CC}` is calculated as the angle between the positive
:math:`z` axis and the vector between two consecutive beads in an acyl tail.
However, it is also possible to pass to :class:`SCC` local membrane normals
to use instead of the positive :math:`z` axis.

A simple example would involve first determinig which leaflet each lipid is in
--- *-1* to indicate the lower leaflet and *1* to indicate the upper leaflet ---
then passing this information to :class:`SCC`::

  import MDAndalysis as mda
  from lipyphilic.lib.assign_leaflets import AssignLeaflets
  from lipyphilic.lib.order_parameter import SCC
  
  u = mda.Universe(tpr, trajectory)
  
  leaflets = AssignLeaflets(
    universe=u,
    lipid_sel="name GL1 GL2 ROH"  # assign all lipids to leaflets, incluing cholesterol (ROH)
  )
  leaflets.run(verbose=True)
  
  scc_sn1 = SCC(
    universe=u,
    tail_sel="name ??A",
    normals=leaflets.filter_leaflets(lipid_sel="name ??A")  # pass only DPPC/DOPC leaflet info
  )
  scc_sn1.run(verbose=True)
  
In this case, `leaflets.filter_leaflets("name ??A")` is a :class:`numpy.ndarray`
that contains the leaflet information for each DPPC and DOPC molecule at each frame.
It has a shape of (n_residues, n_frames). :class:`SCC` assumes that the data therefore
correspond to the sign of the local :math:`z`-axis --- *-1* for lipids in the lower
leaflet and *1* for lipids in the upper leaflet.

You can also calculate local membrane normals using, for example, `MemSurfer
<https://pubs.acs.org/doi/abs/10.1021/acs.jctc.9b00453>`__. If you store the local
membrane normals in a :class:`numpy.ndarray` called *normals*, with shape
(n_residues, n_frames, 3), then you can simply pass these normals to `SCC`::

  scc_sn1 = SCC(
    universe=u,
    tail_sel="name ??A",
    normals=normals
  )
  scc_sn1.run(verbose=True)

The class and its methods
-------------------------

.. autoclass:: SCC
    :members:

"""

import numpy as np

from lipyphilic.lib import base


class SCC(base.AnalysisBase):
    """Calculate coarse-grained acyl tail ordere parameter.
    """

    def __init__(self, universe,
                 tail_sel,
                 normals=None):
        """Set up parameters for calculating the SCC.

        Parameters
        ----------
        universe : Universe
            MDAnalysis Universe object
        tail_sel : str
            Selection string for atoms in either the sn1 **or** sn2 tail of lipids in the
            membrane
        normals : numpy.ndarray, optional
            Local membrane normals. If the array is 2D and of shape (n_residues, n_frames),
            the values in the array are taken to correspond to the sign of the 'z'-axis:
            '-1' for the lower leaflet and '1' for the upper leaflet. If the array is 3D
            and of shape (n_residues, n_frames, 3) then the values are taken to
            correspond to the vector components of local membrane normals.

        Tip
        ----

        Data for 'normals' can be generated using `lipyphilic.lib.assign_leaflets.AssignLeaflets`.
        """
        super(SCC, self).__init__(universe.trajectory)

        self.u = universe
        self.tails = self.u.select_atoms(tail_sel, updating=False)
        
        # For fancy slicing of atoms for each species
        self.tail_atom_mask = {species: self.tails.resnames == species for species in np.unique(self.tails.resnames)}
        
        # For assigning Scc to correct lipids
        self.tail_residue_mask = {species: self.tails.residues.resnames == species for species in np.unique(self.tails.resnames)}
        
        normals = np.array(normals)
        if normals.ndim not in [0, 2, 3]:
            raise ValueError("'normals' must either be a 2D array containing leaflet ids "
                             "of each lipid, or a 3D array containing local membrane normals."
                             )

        if normals.ndim > 0 and len(normals) != self.tails.n_residues:
            raise ValueError("The shape of 'normals' must be (n_residues,)")
        
        if normals.ndim == 2:
            
            normals = np.concatenate(
                (np.zeros_like(normals)[:, :, np.newaxis], np.zeros_like(normals)[:, :, np.newaxis], np.sign(normals)[:, :, np.newaxis]),
                axis=2
            )
            
        self.normals = normals
        self.SCC = None
        
    def _prepare(self):
        
        if self.normals.ndim == 0:
            
            self.normals = np.zeros((self.tails.n_residues, self.n_frames, 3))
            self.normals[:, :, 2] = 1
        
        if self.normals.shape[1] != self.n_frames:
            raise ValueError("The frames to analyse must be identical to those used "
                             "for calculating local membrane normals."
                             )
        
        # Output array
        self.SCC = np.full(
            (self.tails.n_residues, self.n_frames),
            fill_value=np.NaN
        )

    def _single_frame(self):
        
        for species in np.unique(self.tails.resnames):
                
            # Reshape positions so the first axis is per residue
            species_atoms = self.tails[self.tail_atom_mask[species]]
            tail_pos = species_atoms.positions.reshape(species_atoms.n_residues, -1, 3)
            
            cc_vectors = np.diff(tail_pos, axis=1)
            
            # Account for PBC
            for dim in range(3):

                cc_vectors[:, :, dim][cc_vectors[:, :, dim] > self._ts.dimensions[dim] / 2.0] -= self._ts.dimensions[dim]
                cc_vectors[:, :, dim][cc_vectors[:, :, dim] < -self._ts.dimensions[dim] / 2.0] += self._ts.dimensions[dim]
            
            # Calclate SCC using membrane normals
            tails_normals = self.normals[self.tail_residue_mask[species], self._frame_index][:, np.newaxis, :]
            cos_theta = np.sum(cc_vectors * tails_normals, axis=2) / (np.linalg.norm(cc_vectors, axis=2) * np.linalg.norm(tails_normals, axis=2))
            s = (3 * cos_theta**2 - 1) * 0.5
            
            self.SCC[self.tail_residue_mask[species], self._frame_index] = np.mean(s, axis=1)

    @staticmethod
    def weighted_average(sn1_scc, sn2_scc, return_indices=False):
        """Calculate the weighted average Scc of two tails.

        Given two SCC objects, a weighted average of the Scc of each lipid is calculated.

        Parameters
        ----------
        sn1_scc : SCC
            An SCC object for which the order parameters have been calculated.
        sn2_scc : SCC
            An SCC object for which the order parameters have been calculated.
        return_indices : bool, optional
            Whether to return the residue indices of each lipid. This is useful if there are
            some residues present in 'sn1_scc' that are not in 'sn2_scc', or vice-versa.
            
        Returns
        -------
        scc : numpy.ndarray
            An array containing the weighted average Scc of each lipid at each frame.
            
        Warning
        -------
        The frames used in analysing 'sn1_scc' and 'sn2_scc' must be the same - i.e. the 'start',
        'stop', and 'step' parameters passed to the '.run()' methods must be identical.
        
        """
        
        if not ((sn1_scc.n_frames == sn2_scc.n_frames) and (sn1_scc.frames == sn2_scc.frames).all()):
            raise ValueError("sn1_scc and sn2_scc must have been run with the same frames")
        
        sn1_resindices = sn1_scc.tails.residues.resindices
        sn2_resindices = sn2_scc.tails.residues.resindices
        combined_resindices = np.unique([sn1_resindices, sn2_resindices])
        n_residues = combined_resindices.size
        
        scc = np.zeros((n_residues, sn1_scc.n_frames))
        
        for species in np.unique([sn1_scc.tails.resnames, sn2_scc.tails.resnames]):
            
            if species not in sn1_scc.tails.resnames:
                
                # Use sn2 tail only
                species_scc = sn2_scc.SCC[sn2_scc.tail_residue_mask[species]]
                species_resindices = np.in1d(combined_resindices, sn2_resindices[sn2_scc.tail_residue_mask[species]])
                scc[species_resindices] = species_scc
            
            elif species not in sn2_scc.tails.resnames:
                
                # Use sn1 tail only
                species_scc = sn1_scc.SCC[sn1_scc.tail_residue_mask[species]]
                species_resindices = np.in1d(combined_resindices, sn1_resindices[sn1_scc.tail_residue_mask[species]])
                scc[species_resindices] = species_scc
                
            else:
                
                # Calculate mean SCC for the lipid based on the number of beads in each tail
                sn1_species_scc = sn1_scc.SCC[sn1_scc.tail_residue_mask[species]]
                sn1_n_atoms_per_lipid = sn1_scc.tails[sn1_scc.tail_atom_mask[species]].n_atoms / len(sn1_species_scc)
                
                sn2_species_scc = sn2_scc.SCC[sn2_scc.tail_residue_mask[species]]
                sn2_n_atoms_per_lipid = sn2_scc.tails[sn2_scc.tail_atom_mask[species]].n_atoms / len(sn2_species_scc)
                
                species_scc = np.average(
                    np.array([sn1_species_scc, sn2_species_scc]),
                    axis=0,
                    weights=[sn1_n_atoms_per_lipid - 1, sn2_n_atoms_per_lipid - 1]  # - 1 to obain the number of C-C bonds
                )
                
                species_resindices = np.in1d(combined_resindices, sn1_resindices[sn1_scc.tail_residue_mask[species]])
                scc[species_resindices] = species_scc
                
        if return_indices:
            return scc, combined_resindices
        
        else:
            return scc
