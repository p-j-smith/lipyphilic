# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# lipyphilic --- lipyphilic.readthedocs.io
#
# Released under the GNU Public Licence, v2 or any higher version
#

"""
Trajectory transformations --- :mod:`lipyphilic.transformations.transformations`
================================================================================

This module contains methods for applying on-the-fly trajectory transformations
with MDAnalysis.


Fix membranes broken across periodic boundaries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The callable class :class:`center_membrane` can be used to fix a membrane
split across periodic boundaries and then center it in the unit cell. The membrane is iteratively
shifted along a dimension until it is no longer split across periodic boundaries. It is then
moved it to the center of the box in this dimension.

The on-the-fly transformation can be added to your trajectory after loading it with
MDAnalysis:

.. code:: python

  import MDAnalysis as mda
  import lipyphilic as lpp

  u = mda.Universe("production.tpr", "production.xtc")

  ag = u.select_atoms("resname DPPC DOPC CHOL")

  u.trajectory.add_transformations(lpp.transformations.center_membrane(ag))

This will center a DPPC/DOPC/cholesterol membrane in :math:`z` every time a new frame is loaded
into memory by MDAnalysis, such as when you iterate over the trajectory:

.. code:: python

  for ts in u.trajectory:

      # do some nice analysis with your centered membrane

Note
----

`ag` should be an AtomGroup that contains *all* atoms in the membrane.


.. autoclass:: nojump

"""

import pathlib
import warnings

import MDAnalysis as mda  # noqa: N813
import numpy as np
from tqdm.auto import tqdm

__all__ = [
    "center_membrane",
]


class triclinic_to_orthorhombic:  # noqa: N801
    """The class is deprecated and has been removed as its usage leads to incorrect results."""

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        _msg = (
            "`lipyphilic.transformations.triclinic_to_orthorhombic` is deprecated and has been removed "
            "as its usage leads to incorrect results."
        )
        raise NotImplementedError(_msg)


class nojump:  # noqa: N801
    """Prevent atoms jumping across periodic boundaries.

    This is useful if you would like to calculate the diffusion coefficient
    of lipids in your membrane.

    This transformation does an initial pass over the trajectory to determine at which frames
    each atom crosses a boundary, keeping a record of the net movement across each boundary.
    Then, as a frame is loaded into memory, atom positions are translated according to their
    total displacement, taking into account crossing of boundaries as well box fluctuations
    in the box volume.

    By default, atoms are only unwrapped in :math:`xy`, as it is assumed the membrane
    is a bilayer. To unwrap in all dimensions, :attr:`center_z` must also be set to `True`.

     .. deprecated:: 0.11.0
        Please use `MDAnalysis.transformations.NoJump` instead.

    """

    def __init__(self, ag, nojump_x=True, nojump_y=True, nojump_z=False, filename=None):
        """

        Parameters
        ----------
        ag : AtomGroup
            MDAnalysis AtomGroup to which to apply the transformation
        nojump_x : bool, optional
            If true, atoms will be prevented from jumping across periodic boundaries
            in the x dimension.
        nojump_y : bool, optional
            If true, atoms will be prevented from jumping across periodic boundaries
            in the y dimension.
        nojump_z : bool, optional
            If true, atoms will be prevented from jumping across periodic boundaries
            in the z dimension.
        filename : str, optional
            File in which to write the unwrapped, nojump trajectory. The default is `None`,
            in which case the transformation will be applied on-the-fly.py

        Returns
        -------
        :class:`MDAnalysis.coordinates.base.Timestep` object, or `None` if a filename is provided.

        Note
        ----
        The `nojump` transformation is memory intensive to perform on-the-fly. If you have a long
        trajectory or a large number of atoms to be unwrapped, you can write the unwrapped coordinates
        to a new file by providing a :attr:`filename` to :class:`nojump`.

        Warning
        -------
        The current implementation of `nojump` can only unwrap coordinates in orthorhombic systems.

        """
        _msg = (
            "`lipyphilic.transformations.nojump` is deprecated and will be removed "
            "in a later version. Please instead use `MDAnalysis.transformations.NoJump`."
        )
        warnings.warn(
            _msg,
            DeprecationWarning,
            stacklevel=2,
        )

        self.ag = ag
        self.nojump_xyz = np.array([nojump_x, nojump_y, nojump_z], dtype=bool)
        self._nojump_indices = self.nojump_xyz.nonzero()[0]

        if not np.allclose(self.ag.universe.dimensions[3:], 90.0):
            _msg = "nojump requires an orthorhombic box - triclinic systems are not supported."
            raise ValueError(_msg)

        self.ref_pos = ag.positions

        self.filename = pathlib.Path(filename) if filename is not None else filename

        if filename is None:
            self.translate = np.zeros(
                (self.ag.n_atoms, self.ag.universe.trajectory.n_frames, self.nojump_xyz.sum()),
                dtype=np.float64,
            )

            self._on_the_fly()

        else:
            # make the output directory if required
            self.filename.parent.resolve().mkdir(exist_ok=True, parents=True)

            # And we only need to know by the translation vectors for a given frame
            self.translate = np.zeros(
                (self.ag.n_atoms, 3),
                dtype=np.float64,
            )

            self._static_transformation()

    def _on_the_fly(self):
        """Apply the nojump transformation on-the-fly.

        This requires that the translation vector for each atom at each frame can be stored
        in memory, i.e n_atoms_to_unwrap * n_nojump_dimesions * n_frames * 64 bits.

        """

        # First, wrap all atoms into the unit cell and check if that causes them to cross periodic boundaries.
        # Some atoms may be outside of the unit cell because they were made whole.
        # To make performing nojump easier, we will wrap every atom at every frame, so we need to check
        # that wrapping the atom doesn't move it across boundaries at this first step.
        self.ag.universe.trajectory[0]
        self.ref_pos = self.ag.positions  # previous frame minus current frame
        self.ag.wrap(inplace=True)
        diff = self.ref_pos - self.ag.positions

        for index, dim in enumerate(self._nojump_indices):
            # Atoms that moved across the negative direction will have a large positive diff
            crossed_pbc = np.nonzero(diff[:, dim] > self.ag.universe.dimensions[dim] / 2)[0]
            self.translate[crossed_pbc, :, index] += self.ag.universe.dimensions[dim]

            # Atoms that moved across the positive direction will have a large negative diff
            crossed_pbc = np.nonzero(diff[:, dim] < -self.ag.universe.dimensions[dim] / 2)[0]
            self.translate[crossed_pbc, :, index] -= self.ag.universe.dimensions[dim]

        self.ref_pos = self.ag.positions

        for ts in tqdm(self.ag.universe.trajectory[1:], desc="Calculating nojump translation vectors"):
            self.ag.wrap(inplace=True)
            diff = self.ref_pos - self.ag.positions
            for index, dim in enumerate(self._nojump_indices):
                # Atoms that moved across the negative direction will have a large positive diff
                crossed_pbc = np.nonzero(diff[:, dim] > ts.dimensions[dim] / 2)[0]
                self.translate[crossed_pbc, ts.frame :, index] += ts.dimensions[dim]

                # Atoms that moved across the positive direction will have a large negative diff
                crossed_pbc = np.nonzero(diff[:, dim] < -ts.dimensions[dim] / 2)[0]
                self.translate[crossed_pbc, ts.frame :, index] -= ts.dimensions[dim]

            self.ref_pos = self.ag.positions

    def _static_transformation(self):
        """Apply the transformation to one frame at a time, and write the unwrapped coordinates to a file."""

        with mda.Writer(self.filename.as_posix(), "w") as W:
            self.ag.universe.trajectory[0]
            self.ref_pos = self.ag.positions  # previous frame minus current frame
            self.ag.wrap(inplace=True)
            diff = self.ref_pos - self.ag.positions

            for index, dim in enumerate(self._nojump_indices):
                # Atoms that moved across the negative direction will have a large positive diff
                crossed_pbc = np.nonzero(diff[:, dim] > self.ag.universe.dimensions[dim] / 2)[0]
                self.translate[crossed_pbc, index] += self.ag.universe.dimensions[dim]

                # Atoms that moved across the positive direction will have a large negative diff
                crossed_pbc = np.nonzero(diff[:, dim] < -self.ag.universe.dimensions[dim] / 2)[0]
                self.translate[crossed_pbc, index] -= self.ag.universe.dimensions[dim]

            self.ref_pos = self.ag.positions
            self.ag.translate(self.translate)
            W.write(self.ag)

            for ts in tqdm(self.ag.universe.trajectory[1:], desc="Writing NoJump trajectory"):
                self.ag.wrap(inplace=True)
                diff = self.ref_pos - self.ag.positions
                for index, dim in enumerate(self._nojump_indices):
                    # Atoms that moved across the negative direction will have a large positive diff
                    crossed_pbc = np.nonzero(diff[:, dim] > ts.dimensions[dim] / 2)[0]
                    self.translate[crossed_pbc, index] += ts.dimensions[dim]

                    # Atoms that moved across the positive direction will have a large negative diff
                    crossed_pbc = np.nonzero(diff[:, dim] < -ts.dimensions[dim] / 2)[0]
                    self.translate[crossed_pbc, index] -= ts.dimensions[dim]

                self.ref_pos = self.ag.positions
                self.ag.translate(self.translate)
                W.write(self.ag)

    def __call__(self, ts):
        """Unwrap atom coordinates."""

        # Do nothing if it was a static transformation
        if self.filename is not None:
            return ts

        self.ag.wrap(inplace=True)
        translate = np.zeros((self.ag.n_atoms, 3), dtype=np.float64)
        for index, dim in enumerate(self._nojump_indices):
            translate[:, dim] = self.translate[:, ts.frame, index]

        self.ag.translate(translate)

        return ts


class center_membrane:  # noqa: N801
    """Fix a membrane split across periodic boundaries and center it in the primary unit cell.

    If, for example, the bilayer is split across :math:`z`, it will be iteratively
    translated in :math:`z` until it is no longer broken. Then it will
    be moved to the center of the box.

    A membrane with a maximum extent almost the same size as the box length in a given dimension
    will be considered to be split across that dimension.

    By default, the membrane is only centered in :math:`z`, as it is assumed the membrane
    is a bilayer. To center a micelle, :attr:`center_x` and :attr:`center_y` must also be set to `True`.

    Note
    ____

    This transformation only works with orthorhombic boxes.

    """

    def __init__(self, ag, shift=20, center_x=False, center_y=False, center_z=True, min_diff=10):
        """

        Parameters
        ----------
        ag : AtomGroup
            MDAnalysis AtomGroup containing *all* atoms in the membrane.
        shift : float, optional
            The distance by which a bilayer will be iteratively translated. This
            *must* be smaller than the thickness of your bilayer or the diameter
            of your micelle.
        min_diff : float, optional
            Minimum difference between the box size and the maximum extent of the
            membrane in order for the membrane to be considered unwrapped.
        center_x : bool, optional
            If true, the membrane will be iteratively shifted in x until it is
            not longer split across periodic boundaries.
        center_y : bool, optional
            If true, the membrane will be iteratively shifted in y until it is
            not longer split across periodic boundaries.
        center_z : bool, optional
            If true, the membrane will be iteratively shifted in z until it is
            not longer split across periodic boundaries.

        Returns
        -------
        :class:`MDAnalysis.coordinates.base.Timestep` object

        """

        self.membrane = ag
        self.shift = shift
        self.center_xyz = np.array([center_x, center_y, center_z], dtype=bool)
        self.min_diff = min_diff

        if not np.allclose(self.membrane.universe.dimensions[3:], 90.0):
            _msg = "center_membrane requires an orthorhombic box - triclinic systems are not supported."
            raise ValueError(_msg)

    def __call__(self, ts):
        """Fix a membrane split across periodic boundaries."""

        self.membrane.universe.atoms.wrap(inplace=True)

        for dim in range(3):
            if self.center_xyz[dim] == False:  # noqa: E712
                continue

            # Get the maximum membrane thickness in this dimension at the current frame
            # If it's almost the same size as the the box size, then the membrane is split across
            # a periodic boundary
            dim_pos = self.membrane.positions[:, dim]
            max_thickness = np.ptp(dim_pos)

            # If the mmebrane is split, iteratively shift the membrane until it is whole
            # Note: the below conditional will cause an infinite loop if the water/vacuum
            # is less than min_diff Angstrom thick in this dimension
            while max_thickness > (ts.dimensions[dim] - self.min_diff):
                # shift the membrane
                translate_atoms = np.array([0, 0, 0])
                translate_atoms[dim] = self.shift
                self.membrane.universe.atoms.translate(translate_atoms)
                self.membrane.universe.atoms.wrap()

                # check if it is still broken
                dim_pos = self.membrane.positions[:, dim]
                max_thickness = np.ptp(dim_pos)

            # now shift the bilayer to the centre of the box in this dimension
            midpoint = np.mean(self.membrane.positions[:, dim])
            move_to_center = np.array([0, 0, 0])
            move_to_center[dim] = -midpoint + self.membrane.universe.dimensions[dim] / 2
            self.membrane.universe.atoms.translate(move_to_center)

        return ts
