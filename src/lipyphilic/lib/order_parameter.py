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

The class :class:`liyphilic.lib.order_parameter.SCC` calculates the coarse-grained
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

  
The order parameter data are returned in a :class:`numpy.ndarray`, where each
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

If we have used the MARTINI forcefield to study a DPPC/DOPC/cholesterol mixture,
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

And then get a weighted-average :math:`S_{CC}` we can do::

  SCC.weighted_average(scc_sn1, scc_sn2)
  
which will take into account the number of beads in each tail and return a new :class:`SCC`
object whose :attr:`.SCC` attribute contains the weighted-average :math:`S_{CC}` for
each lipid at each frame.

Local membrane normals
----------------------

By default, the :math:`S_{CC}` is calculated as the angle between the positive
:math:`z` axis and the vector between two consecutive beads in an acyl tail.
However, it is also possible to pass to :class:`SCC` local membrane normals
to use instead of the positive :math:`z` axis.

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

  
:math:`S_{CC}` projected onto the membrane plane
------------------------------------------------

Once the :math:`S_{CC}` has been calculated, it is possible to create a 2D plot of time-averaged
:math:`S_{CC}` values projected onto the membrane plane. This can be done using the
:func:`liypphilic.lib.SCC.project_SCC` method, which is a wrapper around the more general
:class:`liypphilic.lib.plotting.ProjectionPlot` class.

If the lipids have been assigned to leaflets, and the weighted average of the sn1 and sn2 tails
stored in an :class:`SCC` object named `scc`, we can plot the projection of the coarse-grained
order parameter onto the membrane plane as follows::
  
  scc_projection = scc.project_SCC(
    lipid_sel="name ??A ??B",
    start=-100,
    stop=None,
    step=None,
    filter_by=leaflets.filter_by("name ??A ??B") == -1
  )
  
The order parameter of each lipid parameter will be averaged over the final 100 frames, as
specified by the :attr:`start` argument. The frame in the middle of the selected frames will be used
for determining lipid positions. In the above case, the lipid positions at frame :math:`-50` will be used.
The :attr:`lipid_sel` specifies that the center of mass of the sn1 ("??A") and sn2 ("??B") atoms will
be used for projecting lipid positions onto the membrane plane. And the `filter_by` argument is used here
to specificy that only lipids in the lower (-1) leaflet should be used for plotting the projected
:math:`S_{CC}` values.
  

The class and its methods
-------------------------

.. autoclass:: SCC
    :members:

"""

import numpy as np

from lipyphilic.lib import base
from lipyphilic.lib.plotting import ProjectionPlot


class SCC(base.AnalysisBase):
    """Calculate coarse-grained acyl tail order parameter.
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
            Local membrane normals, a 3D array of shape (n_residues, n_frames, 3), containing
            `x`, `y` and `z` vector components of the local membrane normals.
        
        """
        
        super(SCC, self).__init__(universe.trajectory)

        self.u = universe
        self.tails = self.u.select_atoms(tail_sel, updating=False)
        
        # For fancy slicing of atoms for each species
        self.tail_atom_mask = {species: self.tails.resnames == species for species in np.unique(self.tails.resnames)}
        
        # For assigning Scc to correct lipids
        self.tail_residue_mask = {species: self.tails.residues.resnames == species for species in np.unique(self.tails.resnames)}
        
        normals = np.array(normals)
        if normals.ndim not in [0, 3]:
            raise ValueError("'normals' must be a 3D array containing local membrane normals of each lipi at each frame.")

        if normals.ndim > 0 and len(normals) != self.tails.n_residues:
            raise ValueError("The shape of 'normals' must be (n_residues, n_frames, 3)")
            
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
    def weighted_average(sn1_scc, sn2_scc):
        """Calculate the weighted average Scc of two tails.

        Given two SCC objects, a weighted average of the Scc of each lipid is calculated.

        Parameters
        ----------
        sn1_scc : SCC
            An SCC object for which the order parameters have been calculated.
        sn2_scc : SCC
            An SCC object for which the order parameters have been calculated.
                
        Returns
        -------
        scc : SCC
            An SCC object with the weighted average Scc of each lipid at each frame stored in the
            scc.SCC attirbute
            
        Warning
        -------
        The frames used in analysing 'sn1_scc' and 'sn2_scc' must be the same - i.e. the 'start',
        'stop', and 'step' parameters passed to the '.run()' methods must be identical.
        
        """
        
        if not ((sn1_scc.n_frames == sn2_scc.n_frames) and (sn1_scc.frames == sn2_scc.frames).all()):
            raise ValueError("sn1_scc and sn2_scc must have been run with the same frames")
        
        sn1_resindices = sn1_scc.tails.residues.resindices
        sn2_resindices = sn2_scc.tails.residues.resindices
        combined_resindices = np.unique(np.hstack([sn1_resindices, sn2_resindices]))
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
                
        # Create a new SCC object
        sn1_atom_indices = sn1_scc.tails.indices
        sn2_atom_indices = sn2_scc.tails.indices
        combined_atom_indices = np.unique(np.hstack([sn1_atom_indices, sn2_atom_indices]))
        
        new_scc = SCC(
          universe=sn1_scc.u,
          tail_sel=f"index {' '.join(combined_atom_indices.astype(str))}",
          
        )
        
        new_scc.start, new_scc.stop, new_scc.step = sn1_scc.start, sn1_scc.stop, sn1_scc.step
        new_scc.frames = np.arange(new_scc.start, new_scc.stop, new_scc.step)
        new_scc.n_frames = new_scc.frames.size
        new_scc.times = sn1_scc.times
        new_scc._trajectory = sn1_scc._trajectory
        new_scc.SCC = scc
                
        return new_scc

    def project_SCC(
        self,
        lipid_sel=None,
        start=None, stop=None, step=None,
        filter_by=None,
        bins=None,
        ax=None,
        cmap=None,
        vmin=None, vmax=None,
        cbar=True,
        cbar_kws={},
        imshow_kws={}
    ):
        """Project the SCC values onto the xy plane of the membrane.
        
        The SCC values, averaged over a selected range of frames, are projected onto the xy
        plane based on the center of mass of each lipid. The atoms to be used in calculating
        the center of mass of the lipids can be specified using the `lipid_sel` arugment.
        
        This method creates an instance of `lipyphilic.lib.plotting.ProjectionPlot` with
        the projected :math:`S_{CC}` interpolated across periodic boundaries.
        The plot is returned so further modification can be performed if needed.
        
        Note
        ----
        The lipid positions are taken from the middle frame of the selected range.
        
        Parameters
        ----------
        lipid_sel: MDAnalysis atom selection, optional
            The center of mass of each lipid will be determined via this selection.
            The default is `None`, in which case every atom of a lipid is used to
            determine its center of mass.
        start: int, optional
            Start frame for averaging the SCC results.
        stop: int, optional
            Final frame for averaging the SCC results.
        step: int, optional
            Number of frames to skip
        filter_by: array-like, optional
            A Boolean mask for selecting a subset of lipids.
            It may take the following shapes:
            
            ``(n_lipids)``
            The mask is used to select a subset of lipids for projecting the SCC
            onto the membrane plane.
            
            ``(n_lipids, n_frames)``
            This is the same shape as the NumPy array created by the
            `lipyphilic.lib.SCC.run()` method. Boolean values are used only from the column
            corresponding to the middle frame of the range selected by `start`, `stop`, and
            `step`.
            
            The default is `None`, in which case no filtering is applied.
        
        bins: int or array_like or [int, int] or [array, array]
            The bin specification:
            
            ``int``
              If int, the number of bins for the two dimensions (nx=ny=bins).
              
            ``array-like``
              If array_like, the bin edges for the two dimensions (x_edges=y_edges=bins).
              
            ``[int, int]``
              If [int, int], the number of bins in each dimension (nx, ny = bins).
              
            ``[array, array]``
              If [array, array], the bin edges in each dimension (x_edges, y_edges = bins).
              
            ``combination``
              A combination [int, array] or [array, int], where int is the number of bins and array is the bin edges.
              
              The default is `None`, in which case a grid with 1 x 1 Angstrom resolution is created.
        
        ax: Axes, optional
            Matplotlib Axes on which to plot the projection. The default is `None`,
            in which case a new figure and axes will be created.
        cmap : str or `~matplotlib.colors.Colormap`, optional
            The Colormap instance or registered colormap name used to map
            scalar data to colors.
        vmin, vmax : float, optional
            Define the data range that the colormap covers. By default,
            the colormap covers the complete value range of the supplied
            data.
        cbar : bool, optional
            Whether or not to add a colorbar to the plot.
        cbar_kws : dict, optional
            A dictionary of keyword options to pass to matplotlib.pyplot.colorbar.
        imshow_kws : dict, optional
            A dictionary of keyword options to pass to matplotlib.pyplot.imshow, which
            is used to plot the 2D density map.
        
        
        Returns
        -------
        scc_projection: ProjectionPlot
            The ProjectionPlot object containing the SCC data and the matplotlob.pyplot.imshow
            plot of the projection.
        
        """
        
        if filter_by is not None:
            filter_by = np.array(filter_by)
            
            if not ((self.SCC.shape == filter_by.shape) or (self.SCC.shape[:1] == filter_by.shape)):
                raise ValueError("The shape of `filter_by` must either be (n_lipids, n_frames) or (n_lipids)")
        
        # Check which lipids to use
        lipid_sel = "all" if lipid_sel is None else lipid_sel
        lipids = self.tails.residues.atoms.select_atoms(lipid_sel)
        keep_lipids = np.in1d(self.tails.residues.resindices, lipids.residues.resindices)

        # Check which frames to use
        start, stop, step = self.u.trajectory.check_slice_indices(start, stop, step)
        frames = np.arange(start, stop, step)
        keep_frames = np.in1d(self.frames, frames)
        frames = self.frames[keep_frames]
        
        # Data fro projecting and frame from which to extract lipid positions
        scc = self.SCC[keep_lipids][:, keep_frames]
        mid_frame = frames[frames.size // 2]
        
        # Check whether we need to filter the lipids
        if filter_by is None:
            filter_by = np.full(scc.shape[0], fill_value=True)
            
        elif filter_by.shape == self.SCC.shape[:1]:
            filter_by = filter_by[keep_lipids]
            
        else:
            mid_frame_index = np.min(np.where(self.frames == mid_frame))
            filter_by = filter_by[keep_lipids][:, mid_frame_index]
        
        # get x and y positions, and make sure the COM is in the unit cell, otherwise is will not be included in the plot
        self.u.trajectory[mid_frame]
        residues = lipids.groupby("resindices")
        lipid_com = np.array([residues[res].center_of_mass(unwrap=True) for res in residues])
        for dim in range(3):
            lipid_com[:, dim][lipid_com[:, dim] > self.u.dimensions[dim]] -= self.u.dimensions[dim]
            lipid_com[:, dim][lipid_com[:, dim] < 0.0] += self.u.dimensions[dim]
        lipids_xpos, lipids_ypos, _ = lipid_com.T
        
        # now we can filter lipids and their values if necessary
        lipids_xpos = lipids_xpos[filter_by]
        lipids_ypos = lipids_ypos[filter_by]
        values = np.mean(scc, axis=1)[filter_by]
        
        # And finally we can create our ProjectionPlot
        scc_projection = ProjectionPlot(lipids_xpos, lipids_ypos, values)
        
        # create grid of values
        if bins is None:
              
            x_dim = self.u.dimensions[0]
            x_bins = np.linspace(0.0, np.ceil(x_dim), int(np.ceil(x_dim)) + 1)
            
            y_dim = self.u.dimensions[1]
            y_bins = np.linspace(0.0, np.ceil(y_dim), int(np.ceil(y_dim)) + 1)
              
            bins = (x_bins, y_bins)
        
        scc_projection.project_values(bins=bins)
        scc_projection.interpolate()
        scc_projection.plot_projection(ax=ax, cmap=cmap, vmin=vmin, vmax=vmax, cbar=cbar, cbar_kws=cbar_kws, imshow_kws=imshow_kws)
        
        return scc_projection
