# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# lipyphilic --- lipyphilic.readthedocs.io
#
# Released under the GNU Public Licence, v2 or any higher version
#

r"""Registration --- :mod:`lipyphilic.lib.registration`
======================================================

:Author: Paul Smith
:Year: 2021
:Copyright: GNU Public License v2

This module provides methods for determining registration of leaflets in a bilayer.

The degree of registration is calculated as the pearson correlation coefficient of
densities in the upper and lower leaflets. First, the 2D density of each leaflet,
:math:`L`, is calculated:

.. math::

  \rho(x, y)_{L} = \displaystyle \int\limits_{-\infty}^{\infty} \frac{1}{\sqrt{2\pi}\sigma} \exp \Bigg({-}\frac{1}{2} \bigg(\frac{x' - x}{\sigma} \bigg)^2 \Bigg) \,dx' dy'


where the :math:`(x, y)` positions of lipid atoms in leaflet :math:`L` are binned
into two-dimensional histograms, then convolved with a circular Gaussian density
of standard deviation :math:`\sigma`. :math:`L` is either the upper (:math:`u`) or
lower (:math:`l`) leaflet.

The correlation between the two leaflets, :math:`r_{u/l}`, is then calculated as
the pearson correlation coefficient between :math:`\rho(x, y)_{u}` and
:math:`\rho(x, y)_{l}`, where values of:

* :math:`1` correspond to perfectly registered
* :math:`-1` correspond to perfectly anti-registered

For more information on interleaflet registration in bilayers
see `Thallmair et al. (2018) <https://pubs.acs.org/doi/abs/10.1021/acs.jpclett.8b01877>`__.

Input
-----

Required:
  - *universe* : an MDAnalysis Universe object.
  - *upper_sel* : atom selection for lipids in the upper leaflet to use in the registration calculation
  - *lower_sel* : atom selection for lipids in the lower leaflet to use in the registration calculation
  - *leaflets* : leaflet membership (-1: lower leaflet, 0: midplane, 1: upper leaflet) of each lipid in the membrane.

Optional:
  - *filter_by* : boolean mask for determining which lipids to include in the registration calculation
  - *n_bins* : the number of bins to use in *x* and *y* for the 2D histogram
  - *gaussian_sd* : the standard deviation of the circular Gaussian to convole with the grid densities
  
Output
------

  - *registration* : the degree of interleaglet registration at each frame

The data are stored in the :attr:`registration.registration` attribute, containing the pearson
correlation coefficient of the two-dimensional leaflet densities at each frame.

Example usage of :class:`Registration`
--------------------------------------

An MDAnalysis Universe must first be created before using :class:`Registration`::

  import MDAnalysis as mda
  from lipyphilic.lib.registration import Registration

  u = mda.Universe(tpr, trajectory)


Then we need to know which leaflet each lipid is in at each frame. This may be done using
:class:`lipyphilic.lib.assign_leaflets.AssignLeaflets`::

  leaflets = AssignLeaflets(
    universe=u,
    lipid_sel="name GL1 GL2 ROH" # assuming we are using the MARTINI forcefield
  )
  leaflets.run()


The leaflets data are stored in the :attr:`leaflets.leaflets` attribute. We can now create our
Registration object by passing our :class:`lipyphilic.lib.assign_leaflets.AssignLeaflets`
object to :class:`Registration` along with atom selections for the lipids::

  registration = Registration(
    leaflets=leaflets.filter_leaflets("resname CHOL and name ROH"),
    upper_sel="resname CHOL and name ROH",
    lower_sel="resname CHOL and name ROH",
  )
  
To calculate the interleaflet correlation of cholesterol molecules using their ROH
beads we then need to use the :func:`run()` method. We select which frames of the trajectory
to analyse (`None` will use every frame) and choose to display a progress bar (`verbose=True`)::
  
  registration.run(
    start=None,
    stop=None,
    step=None,
    verbose=True
  )

The results are then available in the :attr:`registration.registration` attribute as a
:class:`numpy.ndarray`. Again, :math:`1` corresponds to the leaflets being perfectly in register
and :math:`-1` corresponds to the leaflets being perfectly anti-registered.

Selecting a subset of lipids for the registration analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The previous example will compute the registration of cholesterol across the upper
and lower leaflets. In, for example, simulations of phase-separation domains, it is useful
to know the registration of liquid-ordered domains (regardless of the species in the domain)
rather than the registrtion of specific lipid species.

If we have a 2D array, 'lipid_order_data', that contains information on whether each lipid is in
the liquid-disordered phase or the liquid-ordered phase at each frame, we can used this to
calculate the registration of ordered domains. The array must take the shape
'(n_residues, n_frames)', and in the below example 'lipid_order_data[i, j]' would be equal to `-1`
if lipid 'i' is liquid-disordered at frame 'j' and equal to `1` if it is liquid-ordered::

  registration = Registration(
    leaflets=leaflets,
    upper_sel="name PO4 ROH",
    lower_sel="name PO4 ROH",
    filter_by=lipid_order_data == 1
  )


If we have a ternary mixture of DPPC/DOPC/Cholesterol, we can also specifcy that we wish to
consider only DPPC and cholesterol in the liquid-ordered phase::

  registration = Registration(
    leaflets=leaflets,
    upper_sel="(resname CHOL and name ROH) or (resname DPPC and name PO4)",
    lower_sel="(resname CHOL and name ROH) or (resname DPPC and name PO4)",
    filter_by=lipid_order_data == 1
  )
  

Changing the resolution of the 2D grid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, the lipid positions of each leaflet are binned into a two-dimensional
histogram using :math:`n\_bins_x = \lceil{x}\rceil`, where :math:`n\_bins_x` is the
numer of bins in :math:`x` and :math:`\lceil{x}\rceil` is the size of system in :math:`x`
rounded up to the nearest integer. This gives a grid resolution of  *1* Å.

It is also possible to specify the number of bins to use for binning the data::

  registration = Registration(
    leaflets=leaflets,
    upper_sel="resname CHOL and name ROH",
    lower_sel="resname CHOL and name ROH",
    n_bins=100
  )


This will use *100* bins for creating the two-dimensional histogram. Fewer bins
will result in a performance increase but at the cost of spatial resolution. For
all but the largest systems, the default of *1* Å is appropriate. If your system
is larger than a few hundred nm in one dimension, you will likely want to set
:attr:`n_bins` to 2000 or less.

Changing the standard deviation of the circular Gaussian density
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The defualt value of :math:`\sigma` is *15*, which is the value used by
`Thallmair et al. (2018) <https://pubs.acs.org/doi/abs/10.1021/acs.jpclett.8b01877>`__
for determining interleaflet cholesterol correlations. This deault value can be
changed using the :attr:`gaussian_sd` parameter::

  registration = Registration(
    leaflets=leaflets,
    upper_sel="resname CHOL and name ROH",
    lower_sel="resname CHOL and name ROH",
    gaussian_sd=12
  )


Figure 2d of `Thallmair et al. (2018)
<https://pubs.acs.org/doi/abs/10.1021/acs.jpclett.8b01877>`__ shows how correlation
tends to increase with increasing :attr:`gaussian_sd`. This is because the density of
atomic positions is more diffuse and thus more likely to overlap between the two
leaflets. Increasing :attr:`gaussian_sd` also incurs a performance cost.

The class and its methods
-------------------------

.. autoclass:: Registration
    :members:

"""

import numpy as np
import scipy.stats
import scipy.ndimage

from lipyphilic.lib import base


class Registration(base.AnalysisBase):
    """Calculate interleaflet registration in a bilayer.
    """

    def __init__(self,
                 universe,
                 upper_sel,
                 lower_sel,
                 leaflets,
                 filter_by=None,
                 n_bins=None,
                 gaussian_sd=15
                 ):
        """Set up parameters for the registration calculation.
        
        Parameters
        ----------
        universe : Universe
            MDAnalysis Universe object
        upper_sel : str
            Selection string for lipids in the upper leaflet of the bilayer to be used
            for determining registration.
        lower_sel : str
            Selection string for lipids in the lower leaflet of the bilayer to be used
            for determining registration.
        leaflets : numpy.ndarray
            An array of leaflet membership in which: -1 corresponds to the lower leaflet;
            1 corresponds to the upper leaflet; and 0 corresponds to the midplane.
            If the array is 1D and of shape (n_lipids), each lipid is taken to
            remain in the same leaflet over the trajectory. If the array is 2D and
            of shape (n_lipids, n_frames), the leaflet to which each lipid is
            assisgned at each frame will be taken into account when calculating
            the area per lipid.
        filter_by : numpy.ndarray, optional
            A boolean array indicating whether or not to include each lipid in the registration
            analysis. If the array is 1D and of shape (n_lipids), the same lipids will be used
            in the registration analysis at every frame. If the array is 2D and of shape
            (n_lipids, n_frames), the boolean value of each lipid at each frame will be
            taken into account. The default is `None`, in which case no filtering is
            performed.
        n_bins : int, optional
            The number of bins to use in each dimension for the two-dimensional density
            calculations. The default is `None`, in which case the number of bins will be given
            by the size of the system in the 'x' dimension rounded up to the nearest integer.
        gaussian_sd : float, optional
            The standard deviation of the circular Gaussian density to convolve with the
            two-dimensional densities. The spreads out the data to better represent the
            size of the lipids. The default is 15.
        """
        
        super(Registration, self).__init__(universe.trajectory)
        
        self.u = universe
        self.upper_sel = upper_sel
        self.lower_sel = lower_sel
        self.membrane = self.u.select_atoms(f"({self.upper_sel}) or ({self.lower_sel})")
        
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
        
        self.leaflets = leaflets
        
        if filter_by is not None and np.array(filter_by).ndim not in [1, 2]:
            raise ValueError("'filter_by' must either be a 1D array containing non-changing boolean"
                             "values for each lipid, or a 2D array of shape (n_residues, n_frames)"
                             " containing a boolean value for each lipid at each frame."
                             )

        elif filter_by is not None and len(filter_by) != self.membrane.n_residues:
            raise ValueError("The shape of 'filter_by' must be (n_residues,)")
        
        # determine which lipids to use in the analysis at each frame
        if filter_by is None:
            
            self.filter_by = np.full_like(
                self.leaflets,
                fill_value=True,
                dtype=bool
            )
        elif filter_by.ndim == 1:
            
            self.filter_by = np.full_like(
                self.leaflets,
                fill_value=filter_by[:, np.newaxis],
                dtype=bool
            )
        else:
            self.filter_by = filter_by.astype(bool)
                
        self.n_bins = n_bins
        self.gaussian_sd = gaussian_sd
        
        self.registration = None
    
    def _prepare(self):
        
        # Output array
        self.registration = np.full(self.n_frames, fill_value=np.NaN)
        
    def _single_frame(self):
        
        # Atoms must be inside the primary unit cell
        self.membrane.residues.atoms.wrap(inplace=True)
        
        # get the bins for the 2d histograms
        x_length = int(np.ceil(self._ts.dimensions)[0])
        if self.n_bins is not None:
            bins = np.linspace(0, x_length, self.n_bins + 1)
        else:
            bins = np.linspace(0, x_length, x_length + 1)
        
        # Upper leaflet 2d histogram
        upper_ordered_res = self.membrane.residues[np.logical_and(
            self.leaflets[:, self._frame_index] == 1,
            self.filter_by[:, self._frame_index],
        )]
        upper_ordered_atoms = upper_ordered_res.atoms.select_atoms(self.upper_sel)
        upper_hist, *_ = np.histogram2d(
            upper_ordered_atoms.positions[:, 0],
            upper_ordered_atoms.positions[:, 1],
            bins=bins
        )
        # convolve with gaussian kernel
        upper_hist = scipy.ndimage.gaussian_filter(
            upper_hist,
            sigma=self.gaussian_sd,
            mode="wrap"
        )
        
        # Find lower ordered atoms
        lower_ordered_res = self.membrane.residues[np.logical_and(
            self.leaflets[:, self._frame_index] == -1,
            self.filter_by[:, self._frame_index],
        )]
        lower_ordered_atoms = lower_ordered_res.atoms.select_atoms(self.lower_sel)
        lower_hist, *_ = np.histogram2d(
            lower_ordered_atoms.positions[:, 0],
            lower_ordered_atoms.positions[:, 1],
            bins=bins
        )
        
        # convolve with gaussian kernel
        lower_hist = scipy.ndimage.gaussian_filter(
            lower_hist,
            sigma=self.gaussian_sd,
            mode="wrap"
        )
        
        # correlation
        correlation = scipy.stats.pearsonr(
            upper_hist.flatten(),
            lower_hist.flatten()
        )
        self.registration[self._frame_index] = correlation[0]
