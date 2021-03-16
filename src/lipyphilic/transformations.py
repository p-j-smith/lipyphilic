# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# lipyphilic --- lipyphilic.readthedocs.io
#
# Released under the GNU Public Licence, v2 or any higher version
#

"""
Trajectory transformations --- :mod:`lipyphilic.transformations`
================================================================

This module contains methods for applying on-the-fly trajectory transformations
with MDAnalysis.

Fix membranes broken across periodic boundaries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The callable class :class:`lipyphilic.transformations.center_membrane` can be used to center
a bilayer or a micelle in a box. The membrane is iteratively shifted along a dimension until
it is no longer split across periodic boundaries. It is then moved it to the center of the
box in this dimension.

The on-the-fly transformation can be added to your trajectory after loading it with
MDAnalysis:

.. code:: python

  import MDAnalysis as mds
  from lipyphilic.transformations import center_membrane

  u = mda.Universe("production.tpr", "production.xtc")
  
  ag = u.select_atoms("resname DPPC DOPC CHOL")
  
  u.trajectory.add_transformations(center_membrane(ag))

This will center a DPPC/DOPC/cholesterol membrane in :math:`z` every time a new frame is loaded
into memory by MDAnalysis, such as when you iterate over the trajectory:

.. code:: python

  for ts in u.trajectory:
      
      # do some nice analysis with your centered membrane

Note
----

`ag` should be an AtomGroup that contains *all* atoms in the membrane.


.. autoclass:: center_membrane

"""

import numpy as np


class center_membrane:
    """Center a membrane in the primary unit cell.
    
    This is useful if your membrane is split across periodi
    boundaries and you would like to make it whole.
    
    If the bilayer is split across z, it will be iteratively
    translated in z until it is no longer broken. Then it will
    be moved to the center of the box.
    
    By default, the membrane is only centered in z, as it is assumed the membrane
    is a bilayer. To center a micelle, `center_x` and `center_y` must also be set to `True`.
        
    """
    
    def __init__(self, ag, shift=20, center_x=False, center_y=False, center_z=True):
        """
        
        Parameters
        ----------
        ag : AtomGroup
            MDAnalysis AtomGroup containing *all* atoms in the membrane.
        shift : float, optional
            The distance by which a bilayer will be iteratively translated. This
            *must* be smaller than the thickness of your bilayer or the diameter
            of your micelle.
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
        
    def __call__(self, ts):
        """Center the membrane in the box at a single timestep.
        """
        
        self.membrane.universe.atoms.wrap(inplace=True)

        for dim in range(3):
            
            if self.center_xyz[dim] == False:  # noqa: E712
                continue
        
            # get the maximum membrane thickness in this dim at this frame
            # If it's huge, then the membrane is split across PBC
            dim_pos = self.membrane.positions[:, dim]
            max_thickness = np.ptp(dim_pos)

            # If necessary, shift the membrane
            # Note: the below conditional will cause an infinite loop if the water/vacuum
            # is less than 10 Angstrom thick in this dim
            while max_thickness > (ts.dimensions[dim] - 10):
                
                # shift the membrane
                translate_atoms = np.array([0, 0, 0])
                translate_atoms[dim] = self.shift
                self.membrane.universe.atoms.translate(translate_atoms)
                self.membrane.universe.atoms.wrap()

                # check if it is still broekn
                dim_pos = self.membrane.positions[:, dim]
                max_thickness = np.ptp(dim_pos)

            # now shift the bilayer to the centre of the box in z
            midpoint = np.mean(self.membrane.positions[:, dim])
            move_to_center = np.array([0, 0, 0])
            move_to_center[dim] = -midpoint + self.membrane.universe.dimensions[dim] / 2
            self.membrane.universe.atoms.translate(move_to_center)
        
        return ts
