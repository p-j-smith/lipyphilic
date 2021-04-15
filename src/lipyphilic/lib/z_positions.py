# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# lipyphilic --- lipyphilic.readthedocs.io
#
# Released under the GNU Public Licence, v2 or any higher version
#

"""Lipid `z` positions --- :mod:`lipyphilic.lib.z_positions`
============================================================

:Author: Paul Smith
:Year: 2021
:Copyright: GNU Public License v2

This module provides methods for calculating the distance in :math:`z` of lipids
to the bilayer center.

The class :class:`lipyphilic.lib.z_position.ZPositions` assigns the membrane
midpoint to be at :math:`z = 0` Lipids in the upper leaflet will have positive
:math:`z` values and those in the lower leaflet will have negative :math:`z` values.

Input
------

Required:
  - *universe* : an MDAnalysis Universe object
  - *lipid_sel* : atom selection for all lipids in the bilayer
  - *height_sel* : atom selection for the molecules for which the :math:`z` position will be calculated

Options:
  - *n_bins* : split the membrane into *n_bins \\* n_bins* patches, and calculate local membrane midpoints for each patch

Output
------

  - *z_position* : height in :math:`z` of each selected molecule in the bilayer
  
The :math:`z` positions data are returned in a :class:`numpy.ndarray`, where each row corresponds
to an individual molecule and each column corresponds to an individual frame.


Example usage of :class:`ZPositions`
------------------------------------

An MDAnalysis Universe must first be created before using ZPositions::

  import MDAnalysis as mda
  from lipyphilic.lib.z_positions import ZPositions

  u = mda.Universe(tpr, trajectory)

If we have used the MARTINI forcefield to study a phospholipid/cholesterol mixture,
we can calculate the height of cholesterol in the bilayer as follows::

  z_positions = ZPositions(
    universe=u,
    lipid_sel="name GL1 GL2 ROH",
    height_sel="name ROH"
  )
  
:attr:`lipid_sel` is an atom selection that covers all lipids in the bilayer. This
is used for calculating the membrane midpoint. :attr:`height_sel` selects which
atoms to use for caclulating the height of each each molecule.

Note
----

In the above example we are calculating the height of cholesterol in the bilayer, although
the height of any molecule - even those not in the bilayer, such as peptides - can be
calculated instead.


We then select which frames of the trajectory to analyse (`None` will use every
frame) and choose to display a progress bar (`verbose=True`)::
  
  z_positions.run(
    start=None,
    stop=None,
    step=None,
    verbose=True
  )
  
The results are then available in the :attr:`z_positions.z_positions` attribute as a
:class:`numpy.ndarray`. The array has the shape (n_residues, n_frames). Each row
corresponds to an individual molecule and each column to an individual frame.
The height is signed (not absolute) --- positive and negative values correspond to
the molecule being in the upper of lower leaflet respecitvely.

:math:`z` positions based on local membrane midpoints
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first example computes a global membrane midpoint based on all the atoms
of the lipids in the membrane. :math:`z` positions are then calculated as the distance
to this midpoint. This is okay for planar bilayers, but can lead to inaccurate
results in membranes with undulations. If your bilayer has
undulations, `ZPositions` can account for this by creating a grid in :math:`xy`
of your membrane, calculating the local membrane midpoint in each patch,
then find the distance of each molecule to its local midpoint. This is done through
use of `n_bins`::

  z_positions = ZPositions(
    universe=u,
    lipid_sel="name GL1 GL2 ROH",
    height_sel="name ROH"
    n_bins=10
  )
  
In this example, the membrane will be split into a *10 x 10* grid and a lipid
:math:`z` positions calculated based on the distance to the midpoint of the patch the
molecule is in.

Warning
-------

Using `n_bins` can account for small undulations. However, if you have large unulations in
your bilayer the calculated height will be inaccurate.


The class and its methods
-------------------------

.. autoclass:: ZPositions
    :members:

"""

import numpy as np
import scipy.stats

from lipyphilic.lib import base


class ZPositions(base.AnalysisBase):
    """Calculate the :math:`z` position of molecules in a bilayer.
    """

    def __init__(self, universe,
                 lipid_sel,
                 height_sel,
                 n_bins=1):
        """Set up parameters for calculating :math:`z` positions.

        Parameters
        ----------
        universe : Universe
            MDAnalysis Universe object
        lipid_sel : str
            Selection string for the lipids in a membrane. Atoms in this selection are used
            for calculating membrane midpoints.
        height_sel :  str
            Selection string for molecules for which the height in :math:`z` will be calculated.
            Any residues not in this selection will not have their :math:`z` positions calculated.
        n_bins : int, optional
            Number of bins in *x* and *y* to use to create a grid of membrane patches.
            Local membrane midpoints are computed for each patch, and lipid :math:`z`
            positions calculated based on the distance to their local membrane midpoint. The
            default is `1`, which is equivalent to computing a single global
            midpoint.
            
        Note
        ----

        :attr:`height_sel` must be a subset of :attr:`lipid_sel`
        """
        super(ZPositions, self).__init__(universe.trajectory)

        self.u = universe
        self.membrane = self.u.select_atoms(lipid_sel, updating=False)
        self._height_atoms = self.u.select_atoms(height_sel, updating=False)
            
        # lipid species for which the height in z will be calculated
        self._height_species = np.unique(self._height_atoms.resnames)
        # number of each lipid species
        num_lipids = {lipid: sum(self._height_atoms.residues.resnames == lipid) for lipid in self._height_species}
        # number of atoms (seeds) used in the Voronoi tessellation per molecule for each species
        self._n_atoms_per_lipid = {
            lipid: sum(self._height_atoms.resnames == lipid) // num_lipids[lipid] for lipid in self._height_species
        }
        
        self.n_bins = n_bins
        self.z_positions = None
        
    def _prepare(self):
        
        # Output array
        self.z_positions = np.full(
            (self._height_atoms.n_residues, self.n_frames),
            fill_value=np.NaN
        )

    def _single_frame(self):
        
        # Atoms must be wrapped before creating a lateral grid of the membrane
        self.membrane.wrap(inplace=True)
        self._height_atoms.wrap(inplace=True)

        # Find the midpoint of the bilayer as a function of (x,y), using
        # `n_bins` grid points in each dimensions
        # Use all atoms in the membrane to get better statistics
        if self.n_bins > 1:
            bins = np.linspace(0.0, self._ts.dimensions[0], self.n_bins + 1)
        else:
            # scipy.stats.binned_statistics raises Value error if there is only one bin
            bins = [0.0, self._ts.dimensions[0] + 1, self._ts.dimensions[0] + 2]
        
        memb_midpoint_xy = scipy.stats.binned_statistic_2d(
            x=self.membrane.positions[:, 0],
            y=self.membrane.positions[:, 1],
            values=self.membrane.positions[:, 2],
            statistic="mean",
            bins=bins,
            expand_binnumbers=True
        )
        
        # The height in z of each lipid is calculated as the mean heigh
        # of its selected atoms
        for species in self._height_species:
                
            species_indices = self._height_atoms.resnames == species
            species_atoms = self._height_atoms[species_indices]
            
            # get the binnumbers for each lipid
            species_x_bins, species_y_bins = scipy.stats.binned_statistic_2d(
                x=species_atoms.positions[:, 0],
                y=species_atoms.positions[:, 1],
                values=species_atoms.positions[:, 2],
                statistic="mean",
                bins=bins,
                expand_binnumbers=True
            ).binnumber -1  # These were bin numbers, now bin indices  # noqa: E225
            
            # find the mean height in z of the atoms for each individual lipid
            species_zpos = species_atoms.positions[:, 2] - memb_midpoint_xy.statistic[species_x_bins, species_y_bins]
            
            if self._n_atoms_per_lipid[species] > 1:
                species_zpos = species_zpos.reshape((species_atoms.n_residues, self._n_atoms_per_lipid[species]))
                species_zpos = np.mean(species_zpos, axis=1)
            
            # store z position for current lipid species
            species_resindices = np.in1d(
                self._height_atoms.residues.resindices,
                species_atoms.residues.resindices,
                assume_unique=True
            )
            
            self.z_positions[species_resindices, self._frame_index] = species_zpos
