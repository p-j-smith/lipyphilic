# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# lipyphilic --- lipyphilic.readthedocs.io
#
# Released under the GNU Public Licence, v2 or any higher version
#

"""Membrane thickness --- :mod:`lipyphilic.lib.memb_thickness`
==============================================================

:Author: Paul Smith
:Year: 2021
:Copyright: GNU Public License v2

This module provides methods for calculating the membrane thickness over time.

The membrane thickness is defined as the mean distance in :math:`z` between lipids
in the upper and lower leaflets. A discrete intrinsic surface is constructed for each
leaflet based on user-defined lipid headgroup atoms, and the mean distance in :math:`z`
between the two surfaces defines the membrane thickness.

Input
------

Required:
  - *universe* : an MDAnalysis Universe object
  - *leaflets* : a NumPy array containing the leaflet membership of each lipid at each frame
  - *lipid_sel* : atom selection for lipids in the upper leaflet to used in the thickness calculation

Options:
- *n_bins* : a discrete intrinsic surface of each leaflet is created with *n_bins \\* n_bins* patches

Output
------

  - *thickness* : the mean membrane thickness at each frame
  
Thickness data are returned in a :class:`numpy.ndarray`, containing the mean membrane thickness at each
frame.


Example usage of :class:`MembThickness`
---------------------------------------

An MDAnalysis Universe must first be created before using :class:`MembThickness`::

  import MDAnalysis as mda
  from lipyphilic.lib.memb_thickness import MembThickness

  u = mda.Universe(tpr, trajectory)


Then we need to know which leaflet each lipid is in at each frame. This may be done using
:class:`lipyphilic.lib.assign_leaflets.AssignLeaflets`::

  leaflets = AssignLeaflets(
    universe=u,
    lipid_sel="name GL1 GL2 ROH" # assuming we are using the MARTINI forcefield
  )
  leaflets.run()


The leaflets data are stored in the :attr:`leaflets.leaflets` attribute. We can now create our
MembThickness object by passing the results of :class:`lipyphilic.lib.assign_leaflets.AssignLeaflets`
:class:`MembThickness` along with an atom selection for the lipids::

  memb_thickness = MembThickness(
    universe=u,
    leaflets=leaflets.filter_leaflets("resname DOPC and DPPC"),  # exclude cholesterol from thickness calculation
    lipid_sel="resname DPPC DOPC and name PO4"
  )
  
To calculate the membrane thickness, based on interleaflet PO4 to PO4 distances, we need to use the
:func:`.run()` method. We select which frames of the trajectory to analyse (`None` will use every
frame) and choose to display a progress bar (`verbose=True`)::

  memb_thickness.run(
    start=None,
    stop=None,
    step=None,
    verbose=True
  )

The results are then available in the :attr:`memb_thickness.memb_thickness` attribute as a
:class:`numpy.ndarray`.


Changing the resolution of the 2D grid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, the lipid positions of each leaflet are binned into a two-dimensional
histogram using 1 bins in each dimension. This is equivalent to calculating the mean
height of all headgroup atoms in the bilayer, without discretising the surface.

It is also possible to specify the number of bins to use for creating the 2D
histogram of the surface::

  memb_thickness = MembThickness(
    universe=u,
    leaflets=leaflets.filter_leaflets("resname DOPC and DPPC"),  # exclude cholesterol from thickness calculation
    lipid_sel="resname DPPC DOPC and name PO4",
    n_bins=10
  )

This will use *10* bins in each dimension for creating the two-dimensional histogram.

The class and its methods
-------------------------

.. autoclass:: MembThickness
    :members:

"""


import numpy as np
import scipy.stats

from lipyphilic.lib import base


class MembThickness(base.AnalysisBase):
    """Calculate the bilayer thickness.
    """

    def __init__(self, universe,
                 lipid_sel,
                 leaflets,
                 n_bins=1,):
        """Set up parameters for calculating membrane thickness.

        Parameters
        ----------
        universe : Universe
            MDAnalysis Universe object
        leaflets : numpy.ndarray (n_lipids,)
            An array of leaflet membership in which: -1 corresponds to the lower leaflet;
            1 corresponds to the upper leaflet; and 0 corresponds to the midplane.
            If the array is 1D and of shape (n_lipids), each lipid is taken to
            remain in the same leaflet over the trajectory. If the array is 2D and
            of shape (n_lipids, n_frames), the leaflet to which each lipid is
            assisgned at each frame will be taken into account when calculating
            the area per lipid.
        lipid_sel : str, optional
            Selection string for lipid atoms to be used in the thickness calculation.
            The default is `None`, in which case all atoms of the lipids will be used.
        n_bins : int, optional
            Number of bins in *x* and *y* to use to create a grid of membrane patches.
            The intrinsic surface of a leaflet is constructed via the height in `z`
            of each patch. The default is `1`, which is equivalent to computing a
            single global leaflet height.
        
        Tip
        ---
        
        Leaflet membership can be determined using :class:`lipyphilic.lib.assign_leaflets.AssignLeaflets`.
        
        """
        super(MembThickness, self).__init__(universe.trajectory)

        self.u = universe
        self.lipid_sel = lipid_sel if lipid_sel is not None else "all"
        self.membrane = self.u.select_atoms(self.lipid_sel, updating=False)
        
        if np.array(leaflets).ndim not in [1, 2]:
            raise ValueError("'leaflets' must either be a 1D array containing non-changing "
                             "leaflet ids of each lipid, or a 2D array of shape (n_residues, n_frames)"
                             " containing the leaflet id of each lipid at each frame."
                             )

        if len(leaflets) != self.membrane.n_residues:
            raise ValueError("The shape of 'leaflets' must be (n_residues,), but 'lipid_sel' "
                             f"generates an AtomGroup with {self.membrane.n_residues} residues"
                             f" and 'leaflets' has shape {leaflets.shape}."
                             )
        
        self.leaflets = np.array(leaflets)
        self.n_bins = n_bins
        
        self.memb_thickness = None
        
    def _prepare(self):
        
        # Output array
        self.memb_thickness = np.full(self.n_frames, fill_value=np.NaN)

    def _single_frame(self):
        
        # Atoms must be wrapped before creating a lateral grid of the membrane
        self.membrane.wrap(inplace=True)

        # Find the height of each leaflet as a function of (x,y), using
        # `n_bins` grid points in each dimensions
        # Use all atoms in the membrane to get better statistics
        # get the bins for the 2d histograms
        x_length = int(np.ceil(self._ts.dimensions)[0])
        if self.n_bins > 1:
            bins = np.linspace(0.0, x_length, self.n_bins + 1)
        else:
            # scipy.stats.binned_statistics raises Value error if there is only one bin
            bins = [0.0, x_length + 1, x_length + 2]
        
        # Upper leaflet 2d histogram
        upper_res = self.membrane.residues[self.leaflets[:, self._frame_index] == 1]
        upper_atoms = upper_res.atoms.select_atoms(self.lipid_sel)
        
        upper_surface = scipy.stats.binned_statistic_2d(
            x=upper_atoms.positions[:, 0],
            y=upper_atoms.positions[:, 1],
            values=upper_atoms.positions[:, 2],
            statistic="mean",
            bins=bins
        ).statistic
        
        # Lower leaflet 2d histogram
        lower_res = self.membrane.residues[self.leaflets[:, self._frame_index] == -1]
        lower_atoms = lower_res.atoms.select_atoms(self.lipid_sel)
        
        lower_surface = scipy.stats.binned_statistic_2d(
            x=lower_atoms.positions[:, 0],
            y=lower_atoms.positions[:, 1],
            values=lower_atoms.positions[:, 2],
            statistic="mean",
            bins=bins
        ).statistic
        
        thickness = upper_surface - lower_surface
        self.memb_thickness[self._frame_index] = np.mean(thickness)
