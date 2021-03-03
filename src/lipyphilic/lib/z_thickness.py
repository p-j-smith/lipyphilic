# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# lipyphilic --- lipyphilic.readthedocs.io
#
# Released under the GNU Public Licence, v2 or any higher version
#

"""Lipid :math:`z` thickness --- :mod:`lipyphilic.lib.z_thickness`
==================================================================

:Author: Paul Smith
:Year: 2021
:Copyright: GNU Public License v2

This module provides methods for calculating the thickness in :math:`z` of lipids
or lipid tails.

The thickness of lipid tails is a useful input feature for creating Hidden Markov
Models (HMM) to detect phase separation in lipid bilayers. See `Park and Im (2019)
<https://pubs.acs.org/doi/abs/10.1021/acs.jctc.8b00828>`__ for a description of
using HMMs in lipid membrane analysis.

Input
------

Required:
  - *universe* : an MDAnalysis Universe object
  - *lipid_sel* : atom selection for the atoms to be used in calculating the thickness of a lipid

Output
------

  - *z_thickness* : thickness in :math:`z` of each lipid in the bilayer
  
The :math:`z` thickness data are returned in a :class:`numpy.ndarray`, where each row corresponds
to an individual lipid and each column corresponds to an individual frame.


Example usage of :class:`ZThickness`
------------------------------------

An MDAnalysis Universe must first be created before using ZThickness::

  import MDAnalysis as mda
  from lipyphilic.lib.z_thickness import ZThickness

  u = mda.Universe(tpr, trajectory)

If we have used the MARTINI forcefield to study phospholipid/cholesterol mixture,
we can calculate the thickness of cholesterol and *sn1* tails in the bilayer as follows::

  z_thickness_sn1 = ZThickness(
    universe=u,
    lipid_sel="(name ??1 ??A) or (resname CHOL and not name ROH)"
    height_sel="name ROH"
  )
  
Above, our *lipid_sel* selection will select sn1 beads and cholesterol beads in the MARTINI forcefield,
making use of the powerful MDAnalysis atom selection language.

We then select which frames of the trajectory to analyse (`None` will use every
frame) and choose to display a progress bar (`verbose=True`)::
  
  z_thickness_sn1.run(
    start=None,
    stop=None,
    step=None,
    verbose=True
  )
  
The results are then available in the :attr:`z_thickness_sn1.z_thickness` attribute as a
:class:`numpy.ndarray`. The array has the shape (n_residues, n_frames). Each row
corresponds to an individual lipid and each column to an individual frame.


Averaging the thickness of two tails
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Above we saw how to calculate the thickness of the *sn1* tail of lipids along with cholesterol.
Similarly, we can calculate the thickness of the *sn2* tails:

  z_thickness_sn2 = ZThickness(
    universe=u,
    lipid_sel="(name ??1 ??A) or (resname CHOL and not name ROH)"
    height_sel="name ROH"
  )
  z_thickness_sn2.run(verbose=True)

Now, if we would like to know the mean thickness of acyl tails across both *sn1* and *sn2* tails,
we can use the :func:`average` method of :class:`ZThickness`::

  z_thickness = ZThickness.average(
      z_thickness_sn1,
      z_thickness_sn2
  )

This will average the thickness of the two tails, and leave the cholesterol thicknesses (from
*z_thickness_sn1*) unchanges. If *z_thickness_sn1* and *z_thickness)_sn2* contain different lipids,
as they do in our case with cholesterol, it may also be useul to know which row in the returned
array corresponds to which lipid. We can get this by setting :attr:`return_indices` to true::

 z_thickness, indices = ZThickness.average(
      z_thickness_sn1,
      z_thickness_sn2,
      return_indices=True
  )

which will return mean :math:`z` thicknesses of the lipids along with the residue indices of
the lipid in each row.

The class and its methods
-------------------------

.. autoclass:: ZThickness
    :members:

"""

import numpy as np

from lipyphilic.lib import base

