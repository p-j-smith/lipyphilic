# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# lipyphilic --- lipyphilic.readthedocs.io
#
# Released under the GNU Public Licence, v2 or any higher version
#

"""Lipid `z` thickness --- :mod:`lipyphilic.lib.z_thickness`
============================================================

:Author: Paul Smith
:Year: 2021
:Copyright: GNU Public License v2

This module provides methods for calculating the thickness in :math:`z` of lipids
or lipid tails.

The thickness of lipid tails is a useful input feature for creating Hidden Markov
Models (HMM) to detect phase separation in lipid bilayers. See `Park and Im (2019)
<https://pubs.acs.org/doi/abs/10.1021/acs.jctc.8b00828>`__ for a description of
using HMMs in lipid membrane analysis.

Input
------

Required:
  - *universe* : an MDAnalysis Universe object
  - *lipid_sel* : atom selection for the atoms to be used in calculating the thickness of a lipid

Output
------

  - *z_thickness* : thickness in :math:`z` of each lipid in the bilayer
  
The :math:`z` thickness data are returned in a :class:`numpy.ndarray`, where each row corresponds
to an individual lipid and each column corresponds to an individual frame.


Example usage of :class:`ZThickness`
------------------------------------

An MDAnalysis Universe must first be created before using ZThickness::

  import MDAnalysis as mda
  from lipyphilic.lib.z_thickness import ZThickness

  u = mda.Universe(tpr, trajectory)


If we have used the MARTINI forcefield to study a phospholipid/cholesterol mixture,
we can calculate the thickness of cholesterol and *sn1* tails in the bilayer as follows::

  z_thickness_sn1 = ZThickness(
    universe=u,
    lipid_sel="(name ??1 ??A) or (resname CHOL and not name ROH)"
  )


Above, our *lipid_sel* selection will select sn1 beads and cholesterol beads in the MARTINI forcefield,
making use of the powerful MDAnalysis atom selection language.

We then select which frames of the trajectory to analyse (`None` will use every
frame) and choose to display a progress bar (`verbose=True`)::
  
  z_thickness_sn1.run(
    start=None,
    stop=None,
    step=None,
    verbose=True
  )


The results are then available in the :attr:`z_thickness_sn1.z_thickness` attribute as a
:class:`numpy.ndarray`. The array has the shape (n_residues, n_frames). Each row
corresponds to an individual lipid and each column to an individual frame.

Averaging the thickness of two tails
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Above we saw how to calculate the thickness of the *sn1* tail of lipids along with cholesterol.
Similarly, we can calculate the thickness of the *sn2* tails::

  z_thickness_sn2 = ZThickness(
    universe=u,
    lipid_sel="(name ??1 ??A) or (resname CHOL and not name ROH)"
  )
  z_thickness_sn2.run(verbose=True)

Now, if we would like to know the mean thickness of acyl tails across both *sn1* and *sn2* tails,
we can use the :func:`average` method of :class:`ZThickness`::

  z_thickness = ZThickness.average(
      z_thickness_sn1,
      z_thickness_sn2
  )

This will average the thickness of the two tails, leaving the cholesterol thicknesses (from
*z_thickness_sn1*) unchanged, and return a new :class:`ZThickness` object containing the averaged data
in its :attr:`z_thickness` attribute.

The class and its methods
-------------------------

.. autoclass:: ZThickness
    :members:

"""

import numpy as np

from lipyphilic.lib import base


class ZThickness(base.AnalysisBase):
    """Calculate the thickness in z of lipids in a bilayer.
    """

    def __init__(self, universe,
                 lipid_sel):
        """Set up parameters for calculating lipid thicknesses.

        Parameters
        ----------
        universe : Universe
            MDAnalysis Universe object
        lipid_sel : str
            Selection string for atoms to use in calculating lipid thickneses
        
        """
        super(ZThickness, self).__init__(universe.trajectory)

        self.u = universe
        self.lipids = self.u.select_atoms(lipid_sel, updating=False)
        
        # For fancy slicing of atoms for each species
        self.lipid_atom_mask = {species: self.lipids.resnames == species for species in np.unique(self.lipids.resnames)}
        
        # For assigning thickness to correct lipids
        self.lipid_residue_mask = {species: self.lipids.residues.resnames == species for species in np.unique(self.lipids.resnames)}
        
        self.z_thickness = None
        
    def _prepare(self):
        
        # Output array
        self.z_thickness = np.full(
            (self.lipids.n_residues, self.n_frames),
            fill_value=np.NaN
        )

    def _single_frame(self):

        for species in np.unique(self.lipids.resnames):
            
            # Reshape positions so the first axis is per residue
            species_atoms = self.lipids[self.lipid_atom_mask[species]]
            lipid_zpos = species_atoms.positions.reshape(species_atoms.n_residues, -1, 3)[:, :, 2]
            thicknesses = np.max(lipid_zpos, axis=1) - np.min(lipid_zpos, axis=1)
            
            # Check no lipids are split across PBC
            z_dim = self._ts.dimensions[2]
            thicknesses[thicknesses > z_dim / 2] = z_dim - thicknesses[thicknesses > z_dim / 2]
            
            self.z_thickness[self.lipid_residue_mask[species], self._frame_index] = thicknesses

    @staticmethod
    def average(sn1_thickness, sn2_thickness):
        """Calculate the average thickness of two tails.

        Given two ZThickness objects, typically each representing either the sn1 or sn2 tails of the lipids,
        an averagte thickness of each lipid is calculated.

        Parameters
        ----------
        sn1_thickness : ZThickness
            A ZThickness object for which the thicknesses have been calculated.
        sn2_thickness : ZThickness
            A ZThickness object for which the thicknesses have been calculated.
            
        Returns
        -------
        z_thickness : ZThickness
            A new `ZThickness` object containing the averaged data in its `z_thickness` attribute.
        
        Warning
        -------
        The frames used in analysing 'sn1_thickness' and 'sn2_thickness' must be the same - i.e. the 'start',
        'stop', and 'step' parameters passed to the '.run()' methods must be identical.
        
        """
        
        if not ((sn1_thickness.n_frames == sn2_thickness.n_frames) and (sn1_thickness.frames == sn2_thickness.frames).all()):
            raise ValueError("sn1_thickness and sn2_thickness must have been run with the same frames")
        
        sn1_resindices = sn1_thickness.lipids.residues.resindices
        sn2_resindices = sn2_thickness.lipids.residues.resindices
        combined_resindices = np.unique(np.hstack([sn1_resindices, sn2_resindices]))
        n_residues = combined_resindices.size
        
        z_thickness = np.zeros((n_residues, sn1_thickness.n_frames))
        
        resnames = np.unique(np.hstack([sn1_thickness.lipids.residues.resnames, sn2_thickness.lipids.residues.resnames]))
        for species in resnames:
            
            if species not in sn1_thickness.lipids.resnames:
                
                # Use sn2 tail only
                species_thickness = sn2_thickness.z_thickness[sn2_thickness.lipid_residue_mask[species]]
                species_resindices = np.in1d(combined_resindices, sn2_resindices[sn2_thickness.lipid_residue_mask[species]])
                z_thickness[species_resindices] = species_thickness
            
            elif species not in sn2_thickness.lipids.resnames:
                
                # Use sn1 tail only
                species_thickness = sn1_thickness.z_thickness[sn1_thickness.lipid_residue_mask[species]]
                species_resindices = np.in1d(combined_resindices, sn1_resindices[sn1_thickness.lipid_residue_mask[species]])
                z_thickness[species_resindices] = species_thickness
                
            else:
                
                # Calculate mean thickness for the lipid based on the number of atoms in both tails
                sn1_species_thickness = sn1_thickness.z_thickness[sn1_thickness.lipid_residue_mask[species]]
                sn2_species_thickness = sn2_thickness.z_thickness[sn2_thickness.lipid_residue_mask[species]]
                species_thickness = (sn1_species_thickness + sn2_species_thickness) / 2
                                
                species_resindices = np.in1d(combined_resindices, sn1_resindices[sn1_thickness.lipid_residue_mask[species]])
                z_thickness[species_resindices] = species_thickness
        
        # Create a new ZThickness object
        sn1_atom_indices = sn1_thickness.lipids.indices
        sn2_atom_indices = sn2_thickness.lipids.indices
        combined_atom_indices = np.unique(np.hstack([sn1_atom_indices, sn2_atom_indices]))
        
        new_thickness = ZThickness(
          universe=sn1_thickness.u,
          lipid_sel=f"index {' '.join(combined_atom_indices.astype(str))}",
          
        )
        
        new_thickness.start, new_thickness.stop, new_thickness.step = sn1_thickness.start, sn1_thickness.stop, sn1_thickness.step
        new_thickness.frames = np.arange(new_thickness.start, new_thickness.stop, new_thickness.step)
        new_thickness.n_frames = new_thickness.frames.size
        new_thickness.times = sn1_thickness.times
        new_thickness._trajectory = sn1_thickness._trajectory
        new_thickness.z_thickness = z_thickness
        
        return new_thickness
