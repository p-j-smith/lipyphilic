# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# lipyphilic --- lipyphilic.readthedocs.io
#
# Released under the GNU Public Licence, v2 or any higher version
#

"""Area per lipid --- :mod:`lipyphilic.lib.area_per_lipid`
==========================================================

:Author: Paul Smith
:Year: 2021
:Copyright: GNU Public License v2

This module provides methods for calculating the area per lipid in a bilayer.

The class :class:`lipyphilic.lib.area_per_lipid.AreaPerLipid` calculates the
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
to an individual lipid and each columns corresponds to an individual frame. The
results are accessible via the :attr:`AssignLeaflets.leaflets` attribute.

Area data are returned in a :class:`numpy.ndarray`, where each row corresponds
to an individual lipid and each column corresponds to an individual frame, i.e.
areas[i, j] refers to the area of lipid *i* at frame *j*. The results are
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
  from lipyphilic.lib.assign_leaflets import AssignLeaflets

  u = mda.Universe(tpr, trajectory)

Then we need to know which leaflet each lipid is in at each frame. This may be done using
:class:`lipyphilic.lib.assign_leaflets.AssignLeaflets`::

  leaflets = AssignLeaflets(
    universe=u,
    lipid_sel="name GL1 GL2 ROH" # assuming we are using the MARTINI forcefield
  )
  leaflets.run()

The leaflets data are stored in the :attr:`leaflets.leaflets` attribute. We can now create our
AreaPerLipid object::

  areas = AreaPerLipid(
      universe=u,
      lipid_sel="name GL1 GL2 ROH",
      leaflets=leaflets.leaflets
  )
  
The above will use GL1 and GL2 beads to calculate the area of each phospholipid, and the
ROH bead to calculate the area of each sterol.

We then select which frames of the trajectory to analyse (`None` will use every
frame) and choose to display a progress bar (`verbose=True`)::
  
  areas.run(
    start=None,
    stop=None,
    step=None,
    verbose=True
  )
  
Warning
-------
    
The frames used in calculating the area per lipid must be the same as those used for
assigning lipids to leaflets, i.e. the :attr:`start`, :attr:`stop` and :attr:`step` parameters must
be identical for :func:`AssignLeaflets.run()` and :func:`AreaPerLipid.run()`.


The results are then available in the :attr:`areas.areas` attribute as a
:class:`numpy.ndarray`. Each row corresponds to an individual lipid and each column
to an individual frame, i.e `areas.areas[i, j]` contains the area of lipid *i*
at frame *j*.


The class and its methods
-------------------------

.. autoclass:: AreaPerLipid
    :members:

"""
import numpy as np
import freud.locality

from lipyphilic.lib import base


class AreaPerLipid(base.AnalysisBase):
    """Calculate the area of lipids in each leaflet of a bilayer.
    """

    def __init__(self, universe,
                 lipid_sel,
                 leaflets
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
            assisgned at each frame will be taken into account when calculating
            the area per lipid.
            
        Tip
        ---
        
        Leaflet membership can be determined using :class:`lipyphilic.lib.assign_leaflets.AssignLeaflets`.
        
        """
        super(AreaPerLipid, self).__init__(universe.trajectory)

        self.u = universe
        self.membrane = self.u.select_atoms(lipid_sel, updating=False)
        
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
        
        # lipid species in the membrane
        self._lipid_species = np.unique(self.membrane.resnames)
        # number of each lipid species in the membrane
        num_lipids = {lipid: sum(self.membrane.residues.resnames == lipid) for lipid in self._lipid_species}
        # number of atoms (seeds) used in the Voronoi tessellation per molecule for each species
        self._num_seeds = {
            lipid: sum(self.membrane.resnames == lipid) // num_lipids[lipid] for lipid in self._lipid_species
        }
        
        self.areas = None
          
    def _prepare(self):
        
        if (self.leaflets.ndim == 2) and (self.leaflets.shape[1] != self.n_frames):
            raise ValueError("The frames to analyse must be identical to those used "
                             "in assigning lipids to leaflets."
                             )
        
        # Output array
        self.areas = np.full(
            (self.membrane.n_residues, self.n_frames),
            fill_value=np.NaN,
            dtype=float
        )
        
    def _single_frame(self):
        
        # Atoms must be wrapped before creating a lateral grid of the membrane
        self.membrane.wrap(inplace=True)
        frame_leaflets = self.leaflets[:, self._frame_index] if self.leaflets.ndim == 2 else self.leaflets
        
        # Calculate area per lipid for the lower (-1) and upper (1) leaflets
        # Areas cannot be calculated for midplane (0) molecules.
        for leaflet_sign in [-1, 1]:
            
            # freud.order.Voronoi requires z positions set to 0
            leaflet = self.membrane.residues[frame_leaflets == leaflet_sign].atoms
            atoms = leaflet.atoms.intersection(self.membrane)
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
                atom_areas=areas
            )
            
    def _remove_overlapping(self, positions):
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
            positions, return_index=True, return_counts=True, axis=0
        )
        
        # If so, add a small distance between the two atoms (1e-3 A)
        # in the x dimension
        if max(counts > 1):
            for duplicate_index in indices[counts > 1]:
                positions[duplicate_index, 0] += 0.001
                
        return None

    def _get_atom_areas(self, positions):
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
                    "dimensions": 2
                },
                positions
            )
        ).volumes
        
        return areas
    
    def _get_area_per_lipid(self, atoms, atom_areas):
        """Calclate the area per lipid given the areas of every Voronoi cell in a tessellation.
        
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
                axis=1
            )

            # store apl for current lipid species
            species_resindices = np.in1d(
                self.membrane.residues.resindices,
                species_atoms.residues.resindices,
                assume_unique=True
            )
            
            self.areas[species_resindices, self._frame_index] = species_apl
            
        return None
