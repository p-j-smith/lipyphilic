# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# lipyphilic --- lipyphilic.readthedocs.io
#
# Released under the GNU Public Licence, v2 or any higher version
#

"""Plotting --- :mod:`lipyphilic.lib.plotting`
==============================================

:Author: Paul Smith
:Year: 2021
:Copyright: GNU Public License v2

This module provides methods for plotting joint probability densities and lateral
distribution maps of lipid properties projected onto the membrane plane.

The :class:`lipyphilic.lib.plotting.JointDensity` can be used, for example, to
plot a 2D PMF of cholesterol orientation and height in a lipid membrane.py. See
`Baral et al. (2020) <https://www.sciencedirect.com/science/article/pii/S0009308420300980>`_
for an example of the this 2D PMF.

The :class:`lipyphilic.lib.plotting.MembraneMap` can be used, for example, to plot
the area per lipid projected onto the membrane plane, i.e. plot the area per lipid
as a function of :math:`xy`. See `Gu et al. (2020)
<https://pubs.acs.org/doi/full/10.1021/jacs.9b11057>` for examples of these
membrane maps.

The classes and their methods
-----------------------------

.. autoclass:: JointDensity
    :members:
.. autoclass:: MembraneMap
    :members:
"""

import numpy as np
import scipy
import matplotlib.pyplot as plt
import seaborn as sns


class JointDensity():
    """Calculate and plot the joint probability density of two observables.
    
    """
    def __init__(self, ob1, ob2):
        """Set up parameters for calculating joint densities.
        
        Parameters
        ----------
        ob1: array_like
            An array containing values of the first observable.
        ob2: array_like
            An array containing values of the second observable. It **must** have the same shape as `ob1`
        """
        
        if np.array(ob1).shape != np.array(ob2).shape:
            raise ValueError("`ob1` and `ob2` must be arrays of the same shape.")
        
        self.ob1 = np.array(ob1)
        self.ob2 = np.array(ob2)
        
        self.ob1_mesh_bins = None
        self.ob2_mesh_bins = None
        self.joint_mesh_values = None
        self.temperature = None
        
    def calc_density_2D(self, bins, filter_by=None, temperature=None):
        """
        Parameters
        ----------
        bins: int or array_like or [int, int] or [array, array], optional
            The bin specification:
            * If int, the number of bins for the two dimensions (nx=ny=bins).
            * If array_like, the bin edges for the two dimensions (x_edges=y_edges=bins).
            * If [int, int], the number of bins in each dimension (nx, ny = bins).
            * If [array, array], the bin edges in each dimension (x_edges, y_edges = bins).
            * A combination [int, array] or [array, int], where int is the number of bins and array is the bin edges.
        filter_by: 2D numpy array of shape (n_residues, n_frames), optional
            ordered state of each lipid at each frame. The default is `None`, in which case
            no filtering is performed.
        temperature: float, optional
            Temperature of the system, which will be used to convert the 2D density into
            a PMF. The default is `None`, in which case the density is returned rather than
            the PMF.

        """
    
        if filter_by is not None:
            
            filter_by = np.array(filter_by)
            
            if filter_by.shape() != self.ob1.shape:
                raise ValueError("`filter_by` must be an array with the same shape as `ob1` and `ob2`.")
            
            density, ob1_bin_edges, ob2_bin_edges = np.histogram2d(
                self.ob1[filter_by].flatten(), self.ob2[filter_by].flatten(), density=True, bins=bins
            )
        
        else:
            density, ob1_bin_edges, ob2_bin_edges = np.histogram2d(
                self.ob1.flatten(), self.ob2.flatten(), density=True, bins=bins
            )

        # We need to create a grid for plotting with imshow
        self.ob1_mesh_bins, self.ob2_mesh_bins = np.meshgrid(
            ob1_bin_edges[:-1] + (ob1_bin_edges[1] - ob1_bin_edges[0]) / 2,
            ob2_bin_edges[:-1] + (ob2_bin_edges[1] - ob2_bin_edges[0]) / 2,
        )

        # determine whether we use probability density or PMF
        self.temperature = temperature
        if self.temperature is not None:

            ln_density = np.log(density)
            ln_density[~np.isfinite(ln_density)] = np.NaN
            free_energy = -(scipy.constants.Boltzmann * scipy.constants.Avogadro / 4184) * temperature * ln_density
            self.joint_mesh_values = free_energy

        else:
            self.joint_mesh_values = density
            
    def interpolate(self, method="cubic", fill_value=None, rescale=True):
        """Interpolate NaN values in the joint probability density.

        Uses scipy.interpolate.griddata to interpolate the joint density and
        optionally remove NaN values.

        Parameters
        ----------
        method: {'linear', 'nearest', 'cubic'}, optional
            Method of interpolation. One of:

            ``nearest``
              return the value at the data point closest to
              the point of interpolation. See SciPy's
              `NearestNDInterpolator` for more details.

            ``linear``
              tessellate the input point set to N-D
              simplices, and interpolate linearly on each simplex.
              See SciPy's `LinearNDInterpolator` for more details.

            ``cubic``
              return the value determined from a
              piecewise cubic, continuously differentiable (C1), and
              approximately curvature-minimizing polynomial surface. See
              SciPy's `CloughTocher2DInterpolator` for more details.

        fill_value : float, optional
            Value used to fill in for requested points outside of the
            convex hull of the input points. This option has no effect for the
            'nearest' method. If not provided, then the
            default is to use the maximum free energy value if a PMF was
            calculated, or 0 otherwise.

        rescale : bool, optional
            Rescale points to unit cube before performing interpolation.
            This is useful if some of the input dimensions have
            incommensurable units and differ by many orders of magnitude.

        """

        # this snippet is taken from: https://stackoverflow.com/a/37882746
        x, y = np.indices(self.joint_mesh_values.shape)
        
        if fill_value is None:
            if self.temperature is None:
                fill_value = 0.0
            else:
                fill_value = np.nanmax(self.joint_mesh_values)
        
        self.joint_mesh_values[np.isnan(self.joint_mesh_values)] = scipy.interpolate.griddata(
            (x[~np.isnan(self.joint_mesh_values)], y[~np.isnan(self.joint_mesh_values)]),  # points we know
            self.joint_mesh_values[~np.isnan(self.joint_mesh_values)],                     # values we know
            (x[np.isnan(self.joint_mesh_values)], y[np.isnan(self.joint_mesh_values)]),    # points to interpolate
            method=method,
            fill_value=fill_value,
            rescale=rescale
        )

    def plot_density(
        self,
        difference=None,
        ax=None,
        title=None,
        xlabel=None,
        ylabel=None,
        cmap=None,
        vmin=None, vmax=None,
        n_contours=4, contour_labels=None,
        cbar=True,
        cbar_kws={},
        imshow_kws={},
        contour_kws={},
        clabel_kws={}
    ):
        """Plot the 2D density or PMF.
        
        Use matplotlib.pyplot.imshow to plot a heatmap of the density.
        
        Optionally, add contour lines using matplotlib.pyplot.
        
        Parameters
        ----------
        difference: JointDensity, optional
            A JointDensity object for which the probability density or PMF has been calculated.
            Before ploting, the density or PMF of `difference` will be subtracted from the
            density of PMF of this object. This is useful for plotting difference in PMFs due
            to e.g a change in membrane lipid composition.
        ax: Axes, optional
            Matplotlib Axes on which to plot the 2D denstiy. The default is `None`,
            in which case a new figure and axes will be created.
        cmap : str or `~matplotlib.colors.Colormap`, optional
            The Colormap instance or registered colormap name used to map
            scalar data to colors. The default is a colormap created with
            Seaborn, seaborn.cubehelix_palette(start=.5, rot=-.75, as_cmap=True, reverse=True).
            If a difference plot is to be created, by passing a JointDensity object to `difference`,
            tha default cmap is sns.color_palette(palette="RdBu_r", n_colors=200, as_cmap=True).
        vmin, vmax : float, optional
            Define the data range that the colormap covers. By default,
            the colormap covers the complete value range of the supplied
            data.
        n_contours: float
            int or array-like, optional
            Determines the number and positions of the contour lines / regions
            plotted with matplotlib.pyplot.contour.

            If an int *n*, use `~matplotlib.ticker.MaxNLocator`, which tries
            to automatically choose no more than *n+1* "nice" contour levels
            between *vmin* and *vmax*.

            If array-like, draw contour lines at the specified levels.
            The values must be in increasing order.
            
            If 0, no contour lines are drawn.
        contour_labels : array-like, optional
            A list of contour level indices specifyig which levles should be labeled.
            The default is `None`, in which case no contours are labeled.
        cbar : bool, optional
            Whether or not to add a colorbar to the plot.
        cbar_kws : dict, optional
            A dictionary of keyword options to pass to matplotlib.pyplot.colorbar.
        imshow_kws : dict, optional
            A dictionary of keyword options to pass to matplotlib.pyplot.imshow, which
            is used to plot the 2D density map.
        contour_kws : dict, optional
            A dictionary of keyword options to pass to matplotlib.pyplot.contour, which
            is used to plot the contour lines.
        clabel_kws: dict, optional
            A dictionary of keyword options to pass to matplotlib.pyplot.contour, which
            is used to add labels to the contour lines.
        """
        
        # Determine where to plot the figure
        if ax is None:
            self.fig = plt.figure(figsize=(4, 4))
            self.ax = self.fig.add_subplot(1, 1, 1)
        else:
            self.fig = plt.gcf()
            self.ax = ax
        plt.sca(self.ax)
        
        values = self.joint_mesh_values.T - difference.joint_mesh_values.T if difference is not None else self.joint_mesh_values.T
        
        # we need vmin and vmax to be set to sensible values
        # to ensure the colorbar labels looks reasonable
        # And to clip the values
        if vmin is None and self.temperature is not None:
            vmin = min(0.0, np.floor(np.nanmin(values)))
        if vmax is None and self.temperature is not None:
            vmax = max(0.0, np.ceil(np.nanmax(values)))
        
        # make sure the colourbar is centered at zero if we're looking doing a difference plot
        if difference is not None:
            vmax = max(vmax, abs(vmin))
            vmin = -vmax
        
        values = values.clip(vmin, vmax)
        
        # Add contours if necessary
        if "colors" not in contour_kws.keys():
            contour_kws["colors"] = "xkcd:dark grey"
        contours = plt.contour(
            self.ob1_mesh_bins, self.ob2_mesh_bins, values,
            levels=n_contours,
            **contour_kws
        )
        
        # And contour labels
        if contour_labels is not None:
            levels = np.array(contours.levels)
            plt.clabel(contours, levels=levels[contour_labels], **clabel_kws)

        # Get the extent of the distributions
        # This is necessary for imshow
        dx = np.diff(self.ob1_mesh_bins[0][:2])[0]
        dy = np.diff(self.ob2_mesh_bins[0][:2])[0]
        extent = [
            self.ob1_mesh_bins[0][0] - dx / 2,
            self.ob1_mesh_bins[0][-1] + dx / 2,
            self.ob2_mesh_bins[0][0] - dy / 2,
            self.ob2_mesh_bins[-1][0] + dy / 2,
        ]
        
        # Detmine which cmap to use
        if cmap is None and difference is None:
            cmap = sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=True, reverse=True)
        elif cmap is None:
            cmap = sns.color_palette(palette="RdBu_r", n_colors=200, as_cmap=True)
        
        # Finally we can plot the density/PMF
        plt.imshow(
            values,
            origin='lower',  # this is necessary to make sure the y-axis is not inverted
            extent=extent,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            **imshow_kws
        )
        
        # And add a colourbar if necessary
        if cbar:
            
            if "label" not in cbar_kws:
                if self.temperature is not None and difference is not None:
                    cbar_kws["label"] = r"$\Delta\, \rm PMF$"
                    
                elif self.temperature is not None:
                    cbar_kws["label"] = "PMF"
                else:
                    cbar_kws["label"] = "Probability density"
            
            if "aspect" not in cbar_kws:
                cbar_kws["aspect"] = 30
                
            if "pad" not in cbar_kws:
                cbar_kws["pad"] = 0.025
                
            if "ax" not in cbar_kws:
                cbar_kws["ax"] = plt.gca()
                
            self.cbar = plt.colorbar(**cbar_kws)

            ticks = self.cbar.get_ticks()
            labels = ticks.round(2).astype(str)
            labels[-1] += "<"
            if difference is not None:
                labels[0] = "<" + labels[0]
            self.cbar.set_ticks(ticks)  # must be called before we can set the labels
            self.cbar.set_ticklabels(labels)

        plt.title(title, loc="left", weight="bold")
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.gca().set_aspect("auto")  # otherwise imshow assumes the axes have the same units
