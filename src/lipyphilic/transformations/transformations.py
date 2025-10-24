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


import numpy as np

__all__ = [
    "center_membrane",
]


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
