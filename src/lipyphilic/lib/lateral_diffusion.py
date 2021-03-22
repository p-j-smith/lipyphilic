# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# lipyphilic --- lipyphilic.readthedocs.io
#
# Released under the GNU Public Licence, v2 or any higher version
#

r"""
Lateral diffusion --- :mod:`lipyphilic.lib.lateral_diffusion`
=============================================================

This module contains methods for calculating the lateral diffusion coefficient
of lipids in a bilayer.

The class :class:`lipyphilic.lib.lateral_diffusion.MSD` calculates the two-dimensional
mean squared displacent (MSD) of lipids in a bilayer. The `Fast Correlation Algorithm
<https://www.sciencedirect.com/science/article/pii/001046559500048K>`__, implemented
in `tidynamics <http://lab.pdebuyl.be/tidynamics/>`__ is used to calculate the MSD of
each lipid, with optional removal of the center of mass motion of the bilayer.

:class:`lipyphilic.lib.lateral_diffusion.MSD` also contains a method for calculating the
lateral diffusion coefficient, :math:`D_{xy}`, via the Einstein relation:

.. math::

    D_{xy} = \frac{1}{4} \lim_{t\to\infty} \frac{d}{dt} \displaystyle \Bigg\langle \frac{1}{N} \sum_{i=1}^{N} \left  | r_i(t_0 + \Delta t) - r_i(t_0) \right |^2 \displaystyle \Bigg\rangle_{t_0}

where :math:`N` is the number of lipids, :math:`r_i(t0)` is the center of mass in :math:`xy`
of lipids :math:`i` at a time origin `t_0`, :math:`r_i(t0 + \Delta t)` is the same lipid's
center of mass at a lagtime :math:`\Delta t`, and the angular brackets denote an average
over all time origins, :math:`t_0`.

Typically, the MSD is averaged over all molecules. However, :class:`lipyphilic.lib.lateral_diffusion.MSD`
will return the MSD for each individual lipid. This makes it simple to later calculate the diffusion
coefficient using a subset of the lipids, such as a specific lipid species or lipids near a protein.

Input
-----

Required:
  - *universe* : an MDAnalysis Universe object
  - *lipid_sel* : atom selection for calculating the MSD

Optional:
  - *com_removal_sel* : atom selection for center of mass removal from the MSD
  - *dt* : time period betwen consecutive frames in the MSD analysis
  
Output
------

  - *msd* : the mean squared displacement of each lipid at each lagtime, :math:`\Delta t`, in :math:`nm^2`:
  - *lagtimes* : a NumPy array of lagtimes (in :math:`ns`)

The data are stored in the :attr:`MSD.msd` and :attr:`MSD.lagtimes` attribute.


Warning
-------

Before using `lipyphilic.lib.lateral_diffusion.MSD` you *must* ensure that the coordinates have
been unwrapped using, for example, :class:`lipyphilic.transformations.no_jump`.


Example usage of :class:`MSD`
-----------------------------

To calculate the MSD of each lipid in a bilayer we must first load a trajectory using
MDAnalysis::

  import MDAnalysis as mda
  from lipyphilic.lib.lateral_diffusion import MSD

  u = mda.Universe(tpr, trajectory)

If we have used the MARTINI forcefield to study a phospholipid/cholesterol mixture,
we can calculate the MSD of each lipid as follows::

  msd = MSD(
    universe=u,
    lipid_sel="resname DPPC DOPC CHOL"
  )
  
We then select which frames of the trajectory to analyse (`None` will use every
frame) and choose to display a progress bar (`verbose=True`)::
  
  leaflets.run(
    start=None,
    stop=None,
    step=None,
    verbose=True
  )
    
The results are then available in the :attr:`msd.MSD` attribute as a
:class:`numpy.ndarray`. Each row corresponds to an individual lipid and each column
to a different lagtime`.

Center of mass removal
----------------------

During your simulation, it is likely that you removed the center of mass motion of your
bilayer in the :math:`z` direction. However, it is not possible to remove the :math:`x`
and :math:`y` center of mass motions until you have unwrapped your lipid positions.

You may select which lipids to use for the center of mass motion removal using the
:attr:`com_removal_sel` keyword::

  msd = MSD(
    universe=u,
    lipid_sel="resname DPPC",
    com_removla_sel="resname DPPC DOPC CHOL"
  )

In this case, the MSD of DPPC will be calculated with the center of mass motion of the
bilayer will be subtracted from it.

Plotting the MSD of each species
--------------------------------

If you have calculated the MSD of DPPC, DOPC and cholesterol as in the first example, you
can plot the MSD of each species as follows::

  for species in ["DPPC", "DOPC", "CHOL"]:
    
    plt.loglog(
      msd.lagtimes,
      np.mean(msd.msd[msd.membrane.residues.resnames == species], axis=0),
      label=species
    )
  
  plt.legend()

The linear part of the log-log plot can be used for fitting a line and calculating the
diffusion coefficient.

Calculating the lateral diffusion coefficient
---------------------------------------------

After calculating the MSD and identifying the linear portion of the plot, the
:func: `diffusion_coefficient` method of :class:`lipyphilic.lib.lateral_diffusion.MSD`
can be used to calculate :math:`D_{xy}`. We need to pass the time at which to start and
stop the linear fit::

  d, sem = msd.diffusion_coefficient(
    start_fit=400,
    end_fit=600
  )

This will calculate a diffusion coefficient for each individual lipid and return the mean
and standard error of the distribution of coefficients.

To calculate the diffusion coefficient of a subset of lipids we can use the :attr:lipid_sel``
keyword::

  d, sem = msd.diffusion_coefficient(
    start_fit=400,
    end_fit=600,
    lipid_sel="resname CHOL"
  )

which will calculate the lateral diffusion coefficient for cholesterol, using a fit to the MSD
curve from lagtime :math:`\Delta t = 400` to lagtime :math:`\Delta t = 600`.

.. autoclass:: MSD
    :members:

"""

import numpy as np
import scipy.stats
import tidynamics

from lipyphilic.lib import base


class MSD(base.AnalysisBase):
    """Calculate the mean-squared lateral displacement of lipids in a bilayer.
    
    The MSD is returned in units of :math:`nm^2/ns`.
    
    """

    def __init__(self, universe,
                 lipid_sel,
                 com_removal_sel=None,
                 dt=None):
        """Set up parameters for assigning lipids to a leaflet.

        Parameters
        ----------
        universe : Universe
            MDAnalysis Universe object
        lipid_sel : str
            Selection string for calculating the mean-squared displacemnt. IF multiple atoms
            per lipid are selected, the center-of-mass of these atoms will be used for
            calculating the MSD.
        com_removal_sel : str, optional
            The MSD of the center of mass of atoms in this selection will be subtracted from
            all individual lipid MSDs. The default is `None`, in which case no center of mass
            motion removal is performed.
        dt : float, optional
            The time, in nanoseconds, between consecutive frames in `universe.trajectory`.
            The defualt is `None`, in which case `dt` is taken to be `universe.trajectory.dt`
            divided by 1000.
        
        """
        super(MSD, self).__init__(universe.trajectory)

        self.u = universe
        self.membrane = self.u.select_atoms(lipid_sel, updating=False)
        self.com_removal = self.u.select_atoms(com_removal_sel)
        self.dt = dt if dt is not None else self.u.trajectory.dt

        self.msd = None
        self.lagtimes = None
        
        return None
        
    def _prepare(self):
        
        # Output arrays
        self.lipid_com_pos = np.full(
            (self.membrane.n_residues, self.n_frames, 2),
            fill_value=np.NaN,
            dtype=np.float64
        )
        
        self.msd = np.full(
            (self.membrane.n_residues, self.n_frames),
            fill_value=np.NaN
        )
        
        return None

    def _single_frame(self):
        
        self.lipid_com_pos[:, self._frame_index] = self.membrane.center_of_mass(compound="residues")[:, :2]
    
        # Remove COM motion if necessary
        if self.com_removal.n_atoms > 0:
            self.lipid_com_pos[:, self._frame_index] -= self.com_removal.center_of_mass()[:2]
        
        return None
    
    def _conclude(self):
        
        self.lagtimes = self.frames * self.dt
        
        for lipid_index in range(self.membrane.n_residues):
            
            self.msd[lipid_index] = tidynamics.msd(self.lipid_com_pos[lipid_index])

        # MSD must start at 0
        self.msd[:, 0] = 0.0

        # Convert A^2 to nm^2
        self.msd /= 100
        
        return None
    
    def diffusion_coefficient(self, start_fit=None, stop_fit=None, lipid_sel=None):
        """Calculate the lateral diffusion coefficient via the Einstein relation.
        
        A diffusion is calculated for each lipid through a linear fit to its MSD curve.
        The mean and standard error of the diffusion coefficient is returned.

        Parameters
        ----------
        
        start_fit : float, optional
            The time at which to start the linear fit to the MSD curve. The default is
            `None`, in which case the fit will exclude the first 20% of the MSD data.
        stop_fit : float, optional
            The time at which to stop the linear fit to the MSD curve. The default is
            `None`, in which case the fit will exclude the final 20% of the MSD data.
        lipid_sel : str, optional
            Selection string for lipids to include in calculating the diffusion
            coefficient.
            
        Returns
        -------
        
        d : float
            The mean lateral diffusion coefficient, in
            :math:`cm^2/s`., averaged over all lipids in `lipid_sel`.
        sem : float
            The standard error of the diffusion coefficients.
        """
        
        if start_fit is None:
            start_fit_index = self.lagtimes.size * 20 // 100
        else:
            start_fit_index = np.searchsorted(self.lagtimes, start_fit)
            
        if stop_fit is None:
            stop_fit_index = self.lagtimes.size * 80 // 100
        else:
            stop_fit_index = np.searchsorted(self.lagtimes, stop_fit)
            
        if lipid_sel is None:
            mask = np.full(self.membrane.n_residues, fill_value=True, dtype=bool)
        else:
            keep_lipids = self.u.select_atoms(lipid_sel)
            mask = np.in1d(self.membrane.residues.resindices, keep_lipids.residues.resindices)
        
        all_coeffs = np.full(sum(mask), fill_value=np.NaN)
        for index, msd in enumerate(self.msd[mask, start_fit_index:stop_fit_index]):
            
            linear_fit = scipy.stats.linregress(self.lagtimes[start_fit_index:stop_fit_index], msd)
            slope = linear_fit.slope
            d = slope * 1 / 4 * 1e-5  # 1e-5 converts nm^2/ns to cm^2/s
            all_coeffs[index] = d
        
        return np.mean(all_coeffs), scipy.stats.sem(all_coeffs)
