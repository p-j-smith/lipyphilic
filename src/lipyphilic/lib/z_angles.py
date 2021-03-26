# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# lipyphilic --- lipyphilic.readthedocs.io
#
# Released under the GNU Public Licence, v2 or any higher version
#

r"""Lipid `z` angles --- :mod:`lipyphilic.lib.z_angles`
=======================================================

:Author: Paul Smith
:Year: 2021
:Copyright: GNU Public License v2

This module provides methods for calculating the angle lipids make with the
positive :math:`z` axis.

Two atoms must be selected per lipid, and the angle between the :math:`z` axis
and the vector joining the two atoms will be calculated for each lipid. The
vector will always point from atom B to atom A, even for lipids in the lower leaflet.
This means the angle :math:`\theta_{ABz}` will be in the range
:math:`-180° < \theta < 180°`.

Input
-----

Required:
  - *universe* : an MDAnalysis Universe object
  - *atom_A_sel* : atom selection for atom A in each lipid
  - *atom_B_sel* : atom selection for atom B in each lipid

Options:
  - rad : boolean variable specifying whether to return the angle in radians

Output
------

  - *z_angles* : angle made between the :math:`z`-axis and the vector from :math:`B` to :math:`A`

  
The :math:`z` angles data are returned in a :class:`numpy.ndarray`, where each row corresponds
to an individual lipid and each column corresponds to an individual frame.


Example usage of :class:`ZAngles`
---------------------------------

An MDAnalysis Universe must first be created before using ZAngles::

  import MDAnalysis as mda
  from lipyphilic.lib.z_angles import ZAngles

  u = mda.Universe(tpr, trajectory)

If we have used the MARTINI forcefield to study a phospholipid/cholesterol mixture,
we can calculate the orientation of cholesterol in the bilayer as follows::

  z_angles = ZAngles(
    universe=u,
    atom_A_sel="name ROH",
    atom_B_sel="name R5"
  )
  
This will calculate the angle between the :math:`z`-axis and the vector from the
`R5` bead to the `ROH` bead of cholesterol.

We then select which frames of the trajectory to analyse (`None` will use every
frame) and choose to display a progress bar (`verbose=True`)::
  
  z_angles.run(
    start=None,
    stop=None,
    step=None,
    verbose=True
  )
  
The results are then available in the :attr:`z_angles.z_angles` attribute as a
:class:`numpy.ndarray`. The array has the shape (n_residues, n_frames). Each row
corresponds to an individual lipid and each column to an individual frame.

Calculate the angle in radians
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, the results are returned in degrees. We can also specify that the
results should be returned in radians::

  z_angles = ZAngles(
    universe=u,
    atom_A_sel="name ROH",
    atom_B_sel="name R5",
    rad=True
  )


The class and its methods
-------------------------

.. autoclass:: ZAngles
    :members:

"""

import numpy as np

from lipyphilic.lib import base


class ZAngles(base.AnalysisBase):
    """Calculate the orientation of lipids in a bilayer.
    """

    def __init__(self, universe,
                 atom_A_sel,
                 atom_B_sel,
                 rad=False):
        """Set up parameters for calculating the orientations.

        Parameters
        ----------
        universe : Universe
            MDAnalysis Universe object
        atom_A_sel : str
            Selection string for atom `A` of lipids in the membrane.
        atom_B_sel : str
            Selection string for atom `B` of lipids in the membrane.
        rad : bool, optional
            Whether to return the angles in radians. The default is `False`, in which
            case the results are returned in degrees.
            
        Note
        ----

        The orientation is defined as the angle between :math:`z` and the vector from
        atom `B` to atom `A`.
        """
        super(ZAngles, self).__init__(universe.trajectory)

        self.u = universe
        self.atom_A = self.u.select_atoms(atom_A_sel, updating=False)
        self.atom_B = self.u.select_atoms(atom_B_sel, updating=False)
        
        if self.atom_A.n_atoms != self.atom_B.n_atoms:
            raise ValueError("atom_A_sel and atom_B_sel must select the same number of atoms")
        
        self.rad = rad
        self.z_angles = None
        
    def _prepare(self):
        
        # Output array
        self.z_angles = np.full(
            (self.atom_A.n_residues, self.n_frames),
            fill_value=np.NaN
        )

    def _single_frame(self):
            
        v = self.atom_A.positions - self.atom_B.positions

        # Fix PBC
        for i, dim in enumerate(self._ts.dimensions[:2]):
            v[:, i][v[:, i] > (dim / 2.0)] -= dim
            v[:, i][v[:, i] < (-dim / 2.0)] += dim
        
        angles = np.arccos(np.dot(v, [0, 0, 1]) /
                           (np.linalg.norm(v, axis=1) * np.linalg.norm([0, 0, 1]))
                           )
            
        if self.rad is False:
            angles = np.rad2deg(angles)
        
        self.z_angles[:, self._frame_index] = angles
