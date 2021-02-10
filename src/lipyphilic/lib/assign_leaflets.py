# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# lipyphilic --- lipyphilic.readthedocs.io
#
# Released under the GNU Public Licence, v2 or any higher version
#

"""Assign leaflets --- :mod:`lipyphilic.lib.assign_leaflets`
============================================================

:Author: Paul Smith
:Year: 2021
:Copyright: GNU Public License v2

This module provides methods for assigning lipids to leaflets in a bilayer.

The :class:`AssignLeaflets` calculates local membrane midpoints then assigns
each lipid to a leaflet based on the distance in *z* to its local midpoint.

Input
------

Required:
  - *universe* : an MDAnalysis Universe object
  - *lipid_sel*: an MDAnalysis AtomGroup selection.

Options:
  - *midplane_sel* :
  - *midplane_cutoff* :
  - *n_bins* : 


Output
------

  - *leaflet* : leaflet to which each lipid is assigned at each frame
  
Leaflet data are returned in a :class:`numpy.ndarray`, where each column corresponds
to an individual lipid and each columns corresponds to an individual frame.

Example use of :class:`AssignLeaflets`
--------------------------------------


The class and its methods
-------------------------

.. autoclass:: AssignLeaflets
   :members:
"""
import logging

import numpy as np
import scipy.stats

from MDAnalysis.analysis.base import AnalysisBase


class AssignLeaflets(AnalysisBase):
    """
    Assign lipids to leaflets in an MDAnalysis trajectory.
    """

    def __init__(self, universe,
                 lipid_sel=None,
                 midplane_sel=None, midplane_cutoff=0.0,
                 n_bins=1):
        """Set up parameters for finding lipids in a leaflet.

        Parameters
        ----------
        universe : Universe
            MDAnalysis Universe object
        lipid_sel : str
            Selection string for the lipids in a membrane. The selection
            should cover **all** residues in the membrane, including cholesterol
        midplane_sel :  str
            Selection string for residues that may be midplane. Any residues not
            in this selection will be assigned to a leaflet regardless of its
            proximity to the midplane.
            The default is `None`, in which case all lipids will be assigned to
            either the upper or lower leaflet.
        midplane_cutoff : float
            Minimum distance in *z* an atom must be from the midplane to be assigned
            to a leaflet rather than the midplane. The default is `0`, in which case
            all lipids will be assigned to either the upper or lower leaflet.
        n_bins : int
            Number

        Note
        ----

        Typically, `midplane_sel` should select only sterols as other lipids have
        flip-flop rates that are currently unaccessible with MD simulations.

        """
        
        self.u = universe
        self._trajectory = self.u.trajectory
        self.membrane = self.u.select_atoms(lipid_sel, updating=False)
        self.potential_midplane = self.u.select_atoms(midplane_sel, updating=False) if midplane_sel else None
        self.midplane_cutoff = midplane_cutoff
        self.n_bins = n_bins
        self.leaflets = None

    def _prepare(self):
        
        # Output array
        self.leaflets = np.full(
            (self.membrane.n_residues, self.n_frames),
            fill_value=0,
            dtype=np.int8  # smallest sized `np.int` is 1 byte, still 8 times smaller than using `int`
        )

    def _single_frame(self):
        
        # Atoms must be wrapped before creating a lateral grid of the membrane
        self.membrane.wrap(inplace=True)

        # Find the midpoint of the bilayer as a function of (x,y), using
        # `n_bins` grid points in each dimensions
        # Use all atoms in the membrane to get better statistics      
        bins = np.linspace(0.0, self._ts.dimensions[0], self.n_bins+1)
        
        memb_midpoint_xy = scipy.stats.binned_statistic_2d(
            x=self.membrane.residues.atoms.positions[:, 0],
            y=self.membrane.residues.atoms.positions[:, 1],
            values=self.membrane.residues.atoms.positions[:, 2],
            statistic="mean",
            bins=bins,
            expand_binnumbers=True
        )
        
        # get the binnumbers for each lipid
        lipid_x_bins, lipid_y_bins = scipy.stats.binned_statistic_2d(
            x=self.membrane.positions[:, 0],
            y=self.membrane.positions[:, 1],
            values=self.membrane.positions[:, 2],
            statistic="mean",
            bins=bins,
            expand_binnumbers=True
        ).binnumber -1  # These were bin numbers, now bin indices
        
        # Assign leaflets
        upper_leaflet = self.membrane[
            self.membrane.positions[:, 2] > (memb_midpoint_xy.statistic[lipid_x_bins, lipid_y_bins] +
                                             self.midplane_cutoff)
        ]
        self.leaflets[
            np.in1d(self.membrane.residues.resindices, upper_leaflet.residues.resindices),
            self._frame_index
        ] = 1
        
        lower_leaflet = self.membrane[
            self.membrane.positions[:, 2] < (memb_midpoint_xy.statistic[lipid_x_bins, lipid_y_bins] -
                                             self.midplane_cutoff)
        ]
        self.leaflets[
            np.in1d(self.membrane.residues.resindices, lower_leaflet.residues.resindices),
            self._frame_index
        ] = -1
        
        # if necessary, check for midplane residues
        if (self.potential_midplane is not None) and self.midplane_cutoff > 0.0:
            
            # TODO: Move this code block into a method
            self.potential_midplane.wrap(inplace=True)

            midplane_x_bins, midplane_y_bins = scipy.stats.binned_statistic_2d(
                x=self.potential_midplane.positions[:, 0],
                y=self.potential_midplane.positions[:, 1],
                values=self.potential_midplane.positions[:, 2],
                statistic="mean",
                bins=bins,
                expand_binnumbers=True
            ).binnumber -1  # These were bin numbers, now bin indices
            
            # First assume they're all midplane
            # Then find residues that have at least one atom further than
            # `midplane_cutoff` from the local midplane
            midplane_mask = np.full(self.potential_midplane.n_residues, fill_value=True, dtype=bool)

            not_midplane = np.abs(
                self.potential_midplane.positions[:, 2] - memb_midpoint_xy.statistic[midplane_x_bins, midplane_y_bins]
            ) > self.midplane_cutoff

            # These residues have at least one atom in `potential_midplane`
            # that is more the `midplane_cutoff` from the local midplane
            midplane_mask[
                np.in1d(self.potential_midplane.residues.resindices, self.potential_midplane[not_midplane].resindices),
            ] = False

            midplane_residues = self.potential_midplane.residues[midplane_mask]

            # Assign midplane
            self.leaflets[
                np.in1d(self.membrane.residues.resindices, midplane_residues.resindices),
                self._frame_index
            ] = 0
