# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# lipyphilic --- lipyphilic.readthedocs.io
#
# Released under the GNU Public Licence, v2 or any higher version
#

"""Area per lipid --- :mod:`lipyphilic.analysis.area_per_lipid`
===============================================================

:Author: Paul Smith
:Year: 2021
:Copyright: GNU Public License v2

This module provides methods for calculating the area per lipid in a bilayer.

The class :class:`lipyphilic.analysis.area_per_lipid.AreaPerLipid` calculates the
area of each lipid via a 2D Voronoi tessellation of atomic positions. See
`Lukat et al. (2013) <https://pubs.acs.org/doi/full/10.1021/ci400172g>`_ for
a description of calculating the area per lipid via Voronoi tessellations.

This class uses `Freud <https://freud.readthedocs.io/en/stable/index.html#>`_
for performing the Voronoi tessellations from which the area per lipid is
calculated.

Input
------

Required:
  - *universe* : an MDAnalysis Universe object.
  - *lipid_sel* : atom selection for lipids in the bilayer. These atoms will be used to perform the Voronoi tessellation.
  - *leaflets* : leaflet membership (-1: lower leaflet, 0: midplane, 1: upper leaflet) of each lipid in the membrane.


Output
------

  - *area* : area per lipid of each lipid as each frame


Area data are returned in a :class:`numpy.ndarray`, where each row corresponds
to an individual lipid and each column corresponds to an individual frame, i.e.
*areas[i, j]* refers to the area of lipid *i* at frame *j*. The results are
accessible via the :attr:`AreaPerLipid.areas` attribute.

Note
----

No area can be calculated for molecules that are in the midplane,
i.e. those for which `leaflets==0`. These molecules will have `NaN` values
in the results array for the frames at which they are in the midplane.


Example usage of :class:`AreaPerLipid`
--------------------------------------

An MDAnalysis Universe must first be created before using AreaPerLipid::

  import MDAnalysis as mda
  from lipyphilic.leaflets.assign_leaflets import AssignLeaflets

  u = mda.Universe(tpr, trajectory)

Then we need to know which leaflet each lipid is in at each frame. This may be done using
:class:`lipyphilic.leaflets.assign_leaflets.AssignLeaflets`::

  leaflets = AssignLeaflets(
    universe=u,
    lipid_sel="name GL1 GL2 ROH" # assuming we are using the MARTINI forcefield
  )
  leaflets.run()

The leaflet data are stored in the :attr:`leaflets.leaflets` attribute. We can now create our
AreaPerLipid object::

  areas = AreaPerLipid(
    universe=u,
    lipid_sel="name GL1 GL2 ROH",
    leaflets=leaflets.leaflets
  )

The above will use GL1 and GL2 beads to calculate the area of each phospholipid, and the
ROH bead to calculate the area of each sterol. Two Voronoi tessellations will be performed at each
frame --- one for the upper leaflet and one for the lower leaflet.

We then select which frames of the trajectory to analyse (`None` will use every
frame) and choose to display a progress bar (`verbose=True`)::

  areas.run(
    start=None,
    stop=None,
    step=None,
    verbose=True
  )

The results are then available in the :attr:`areas.areas` attribute as a
:class:`numpy.ndarray`. Each row corresponds to an individual lipid and each column
to an individual frame, i.e `areas.areas[i, j]` contains the area of lipid *i*
at frame *j*.

Excluding protein and solvent
------------------------------

By default, the Voronoi tessellation only considers lipid atoms when calculating area per lipid.
However, in systems containing proteins or solvent molecules, these can be included in the
tessellation to account for the space they occupy, effectively reducing the area of the lipids.
Use the `exclude_sel` and `exclude_cutoff` parameters to include and control the influence
of non-lipid atoms::

  # Include protein atoms within 10 Angstroms of lipids in the tessellation
  areas = AreaPerLipid(
    universe=u,
    lipid_sel="name GL1 GL2 ROH",
    leaflets=leaflets.leaflets,
    exclude_sel="protein",
    exclude_cutoff=10.0
  )
  areas.run()

Setting `exclude_cutoff` ensures that only excluded atoms in close proximity to the leaflet
are included, which prevents hovering atoms of the protein far from the membrane to
artificially reduce the calculated area per lipid when they are projected onto the leaflet.

Warning
-------
If your membrane is highly curved the calculated area per lipid will be inaccurate.
In this case we recommend you use either `FATSlim <https://pythonhosted.org/fatslim/>`__,
`MemSurfer <https://github.com/LLNL/MemSurfer>`__ or
`ML-LPA <https://vivien-walter.github.io/mllpa/>`__.


The class and its methods
-------------------------

.. autoclass:: AreaPerLipid
    :members:

"""
import freud.locality
import MDAnalysis
from MDAnalysis.analysis.base import AnalysisBase
import numpy as np

from lipyphilic.plotting import ProjectionPlot

__all__ = [
    "AreaPerLipid",
]


class AreaPerLipid(AnalysisBase):
    """Calculate the area of lipids in each leaflet of a bilayer."""

    def __init__(
        self,
        universe: "MDAnalysis.Universe",
        lipid_sel: str,
        leaflets: np.ndarray,
        exclude_sel: str  = "",
        exclude_cutoff: float = 0,
        exclude_dim: int = 3,
    ):
        """Set up parameters for calculating areas.

        Parameters
        ----------
        universe : Universe
            MDAnalysis Universe object
        lipid_sel : str
            Selection string for lipids in the bilayer. Typically, in all-atom
            simulations, one atom per sterol and three atoms per non-sterol lipid
            would be used. In coarse-grained simulations, one bead per sterol and
            two beads per non-sterol lipid would typically be used.
        leaflets : numpy.ndarray (n_lipids,)
            An array of leaflet membership in which: -1 corresponds to the lower leaflet;
            1 corresponds to the upper leaflet; and 0 corresponds to the midplane.
            If the array is 1D and of shape (n_lipids), each lipid is taken to
            remain in the same leaflet over the trajectory. If the array is 2D and
            of shape (n_lipids, n_frames), the leaflet to which each lipid is
            assigned at each frame will be taken into account when calculating
            the area per lipid.
        exclude_sel : str, optional
            Atom selection to be excluded from contributing to lipid area in the Voronoi tessellation.
            This may be used to account for protein or solvent atoms that occupy space in the membrane.
        exclude_cutoff : float, optional
            Cutoff distance from the leaflet in Angstroms for the excluded atoms.
            This is necessary to ensure that excluded atoms are in close proximity to the leaflet.
            By default, all excluded atoms are used for both leaflets (no cutoff).
        exclude_dim: int, optional
            Dimensions along which to apply the `exclude_cutoff`. Must be either 1 or 3, corresponding
            to z or 3D dimensions, respectively. By default, the cutoff is applied in all three dimensions.
            For one dimension, the distance is between the excluded atom and themean z position of the leaflet.
            For three dimensions, the distance is between each atom of the leaflet and the excluded atoms.
        ---

        Leaflet membership can be determined using :class:`lipyphilic.leaflets.assign_leaflets.AssignLeaflets`.

        """
        super().__init__(universe.trajectory)

        self.u = universe
        self.membrane = self.u.select_atoms(lipid_sel, updating=False)
        self.excluded = self.u.select_atoms(exclude_sel)
        self.exclude_cutoff = exclude_cutoff
        self.exclude_dim = exclude_dim

        if not np.allclose(self.u.dimensions[3:], 90.0):
            _msg = "AreaPerLipid requires an orthorhombic box - triclinic systems are not supported."
            raise ValueError(_msg)

        if np.array(leaflets).ndim not in [1, 2]:
            _msg = (
                "'leaflets' must either be a 1D array containing non-changing "
                "leaflet ids of each lipid, or a 2D array of shape (n_residues, n_frames)"
                " containing the leaflet id of each lipid at each frame.",
            )
            raise ValueError(_msg)

        if len(leaflets) != self.membrane.n_residues:
            _msg = (
                "The shape of 'leaflets' must be (n_residues,), but 'lipid_sel' "
                f"generates an AtomGroup with {self.membrane.n_residues} residues"
                f" and 'leaflets' has shape {leaflets.shape}.",
            )
            raise ValueError(_msg)

        if self.exclude_dim not in [1, 3]:
            _msg = "'exclude_dim' must be either 1 or 3."
            raise ValueError(_msg)

        self.leaflets = np.array(leaflets)

        # lipid species in the membrane
        self._lipid_species = np.unique(self.membrane.resnames)
        # number of each lipid species in the membrane
        num_lipids = {lipid: sum(self.membrane.residues.resnames == lipid) for lipid in self._lipid_species}
        # number of atoms (seeds) used in the Voronoi tessellation per molecule for each species
        self._num_seeds = {
            lipid: sum(self.membrane.resnames == lipid) // num_lipids[lipid] for lipid in self._lipid_species
        }

        self.results.areas = None

    @property
    def areas(self):
        return self.results.areas

    def _prepare(self):
        if (self.leaflets.ndim == 2) and (self.leaflets.shape[1] != self.n_frames):
            _msg = ("The frames to analyse must be identical to those used in assigning lipids to leaflets.",)
            raise ValueError(_msg)

        # Output array
        self.results.areas = np.full(
            (self.membrane.n_residues, self.n_frames),
            fill_value=np.nan,
            dtype=float,
        )

    def _single_frame(self):
        # Atoms must be wrapped before creating a lateral grid of the membrane
        self.membrane.wrap(inplace=True)
        if self.excluded:
            self.excluded.wrap(inplace=True)
        frame_leaflets = self.leaflets[:, self._frame_index] if self.leaflets.ndim == 2 else self.leaflets

        # Calculate area per lipid for the lower (-1) and upper (1) leaflets
        # Areas cannot be calculated for midplane (0) molecules.
        for leaflet_sign in [-1, 1]:
            # freud.order.Voronoi requires z positions set to 0
            leaflet = self.membrane.residues[frame_leaflets == leaflet_sign].atoms
            atoms = leaflet.atoms.intersection(self.membrane)
            if self.excluded:
                filter_excluded = self._filter_exclude(atoms)
                atoms = atoms.union(filter_excluded)
            pos = atoms.positions
            pos[:, 2] = 0

            # Check whether any atoms are overlapping in the xy-plane
            self._remove_overlapping(positions=pos)

            # Voronoi tessellation to get area per atom
            areas = self._get_atom_areas(positions=pos)

            # Calculae area per lipid in the current leaflet
            # by considering the contribution of each
            # atom of a given lipid
            self._get_area_per_lipid(
                atoms=atoms,
                atom_areas=areas,
            )

    def _remove_overlapping(self, positions: np.ndarray) -> None:
        """Ensure no two atoms are overlapping in the xy plane.

        Given an Nx3 array of atomic positions, make minor adjustments to xy positions
        if any pair of xy coordinates are identical.

        If atoms are overlapping in xy, Freud will complain when attempting to perform the
        Voronoi tessellation.

        Parameters
        ----------
        positions : numpy ndarray
            Array of shape (n_atoms, 3) containing atomic coordinates.

        Returns
        -------
        None
            The positions are modified in place.

        """

        # Check whether any atoms are overlapping in the xy-plane
        # This may be an issue in CG sims with cholesteorl flip-flop
        # but is unlikely to be so in all-atom sims
        _, indices, counts = np.unique(
            positions,
            return_index=True,
            return_counts=True,
            axis=0,
        )

        # If so, add a small distance between the two atoms (1e-3 A)
        # in the x dimension
        if max(counts > 1):
            for duplicate_index in indices[counts > 1]:
                positions[duplicate_index, 0] += 0.001


    def _get_atom_areas(self, positions: np.ndarray) -> np.ndarray:
        """Calculate area per atom.

        Given xy coordinates of atomic positions, perform a Voronoi
        tessellation and return the area in xy occupied by each Voronoi cell.

        Parameters
        ----------
        positions : numpy ndarray
            Array of shape (n_atoms, 3) containing atomic coordinates.

        Returns
        -------
        areas : numpy ndarray
            Array of shape (n_atoms) containing the lateral area per atom.

        """

        voro = freud.locality.Voronoi()
        areas = voro.compute(
            system=(
                {
                    "Lx": self._ts.dimensions[0],
                    "Ly": self._ts.dimensions[1],
                    "dimensions": 2,
                },
                positions,
            ),
        ).volumes

        return areas

    def _get_area_per_lipid(self, atoms: "MDAnalysis.AtomGroup", atom_areas: np.ndarray) -> None:
        """Calculate the area per lipid given the areas of every Voronoi cell in a tessellation.

        This involves summing contributions from each atom of a given lipid.

        Parameters
        ----------
        atoms : MDAnalysis AtomGroup
            AtomGroup for which a 2D Voronoi tessellation was performed
        atom_areas : numpy ndarray
            Array of areas of each atom in the 2D Voronoi tessellation

        Returns
        -------
        None
            The lipid areas are modified in place.
        """

        for species in self._lipid_species:
            species_indices = atoms.resnames == species

            # We need to sum the area contribution of each cell for a given lipid
            species_apl = atom_areas[species_indices]
            species_atoms = atoms[species_indices]
            species_apl = np.sum(
                species_apl.reshape(species_atoms.n_residues, self._num_seeds[species]),
                axis=1,
            )

            # store apl for current lipid species
            species_resindices = np.isin(
                self.membrane.residues.resindices,
                species_atoms.residues.resindices,
                assume_unique=True,
            )

            self.results.areas[species_resindices, self._frame_index] = species_apl

    def _filter_exclude(self, lipids: "MDAnalysis.AtomGroup") -> "MDAnalysis.AtomGroup":
        """Filter excluded atoms by distance from the lipids."""
        if self.exclude_cutoff > 0:
            if self.exclude_dim == 3:
                filtered_exclude = self.excluded.select_atoms(
                    f"around {self.exclude_cutoff} global group leaflet",
                    leaflet=lipids,
                )
            else:  # exclude_dim == 1
                leaflet_mean_z = lipids.positions[:, 2].mean()
                excluded_distances = np.abs(self.excluded.positions[:, 2] - leaflet_mean_z)
                filtered_exclude = self.excluded[excluded_distances <= self.exclude_cutoff]
            return filtered_exclude
        return self.excluded


    def project_area(
        self,
        lipid_sel=None,
        start=None,
        stop=None,
        step=None,
        filter_by=None,
        bins=None,
        ax=None,
        cmap=None,
        vmin=None,
        vmax=None,
        cbar=True,
        interpolate_kws=None,
        cbar_kws=None,
        imshow_kws=None,
    ):
        """Project the area per lipid onto the xy plane of the membrane.

        The areas per lipid, averaged over a selected range of frames, are projected onto the xy
        plane based on the center of mass of each lipid. The atoms to be used in calculating
        the center of mass of the lipids can be specified using the `lipid_sel` arugment.

        This method creates an instance of `lipyphilic.plotting.ProjectionPlot` with
        the projected areas interpolated across periodic boundaries.
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
            Start frame for averaging the area per lipid results.
        stop: int, optional
            Final frame for averaging the area per lipid results.
        step: int, optional
            Number of frames to skip
        filter_by: array-like, optional
            A Boolean mask for selecting a subset of lipids.
            It may take the following shapes:

            ``(n_lipids)``
            The mask is used to select a subset of lipids for projecting the areas
            onto the membrane plane.

            ``(n_lipids, n_frames)``
            This is the same shape as the NumPy array created by the
            `lipyphilic.analysis.AreaPerLipid.run()` method. Boolean values are used only from the column
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
        interpolate_kws : dict, optional
            A dictionary of keyword options to pass to `scipy.interpolate.griddata`, which is used
            to interpolate the projected areas across periodic boundaries.
        cbar_kws : dict, optional
            A dictionary of keyword options to pass to matplotlib.pyplot.colorbar.
        imshow_kws : dict, optional
            A dictionary of keyword options to pass to matplotlib.pyplot.imshow, which
            is used to plot the 2D density map.


        Returns
        -------
        area_projection: ProjectionPlot
            The ProjectionPlot object containing the area per lipid data and the matplotlob.pyplot.imshow
            plot of the projection.

        """

        if filter_by is not None:
            filter_by = np.array(filter_by)

            if not (
                (self.results.areas.shape == filter_by.shape)
                or (self.results.areas.shape[:1] == filter_by.shape)
            ):
                _msg = "The shape of `filter_by` must either be (n_lipids, n_frames) or (n_lipids)"
                raise ValueError(_msg)

        # Check which lipids to use
        lipid_sel = "all" if lipid_sel is None else lipid_sel
        lipids = self.membrane.residues.atoms.select_atoms(lipid_sel)
        keep_lipids = np.isin(self.membrane.residues.resindices, lipids.residues.resindices)

        # Check which frames to use
        start, stop, step = self.u.trajectory.check_slice_indices(start, stop, step)
        frames = np.arange(start, stop, step)
        keep_frames = np.isin(self.frames, frames)
        frames = self.frames[keep_frames]

        # Data for projecting and frame from which to extract lipid positions
        areas = self.results.areas[keep_lipids][:, keep_frames]
        mid_frame = frames[frames.size // 2]

        # Check whether we need to filter the lipids
        if filter_by is None:
            filter_by = np.full(areas.shape[0], fill_value=True)

        elif filter_by.shape == self.results.areas.shape[:1]:
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
        values = np.nanmean(areas, axis=1)[
            filter_by
        ]  # some molecules may be midplane during the period considered

        # If excluded atoms were provided, include their positions with zero value
        # so that they appear as zero-area seeds in the projection.
        if self.excluded:
            filter_excluded = self._filter_exclude(lipids)
            n_x, n_y, _ = filter_excluded.positions.T
            # Append excluded positions with zero values
            lipids_xpos = np.concatenate((lipids_xpos, n_x))
            lipids_ypos = np.concatenate((lipids_ypos, n_y))
            values = np.concatenate((values, np.zeros(n_x.shape[0], dtype=float)))
        # And finally we can create our ProjectionPlot
        area_projection = ProjectionPlot(lipids_xpos, lipids_ypos, values)

        # create grid of values
        if bins is None:
            x_dim = self.u.dimensions[0]
            x_bins = np.linspace(0.0, np.ceil(x_dim), int(np.ceil(x_dim)) + 1)

            y_dim = self.u.dimensions[1]
            y_bins = np.linspace(0.0, np.ceil(y_dim), int(np.ceil(y_dim)) + 1)

            bins = (x_bins, y_bins)

        area_projection.project_values(bins=bins)

        interpolate_kws = interpolate_kws or {}
        area_projection.interpolate(**interpolate_kws)
        area_projection.plot_projection(
            ax=ax,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            cbar=cbar,
            cbar_kws=cbar_kws,
            imshow_kws=imshow_kws,
        )

        return area_projection
