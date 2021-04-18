# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# lipyphilic --- lipyphilic.readthedocs.io
#
# Released under the GNU Public Licence, v2 or any higher version
#

"""Flip-flop --- :mod:`lipyphilic.lib.flip_flop`
================================================

:Author: Paul Smith
:Year: 2021
:Copyright: GNU Public License v2

This module provides methods for finding flip-flop events in a lipid bilayer.

A flip-flop event occurs when a molecule - typically a sterol - moves from
one leaflet of a bilayer into the opposing leaflet.

The class :class:`lipyphilic.lib.flip_flop.FlipFlop` finds the frames at which
a flip-flop event begins and ends, as well as the direction of travel (upper-to-lower
or lower-to-upper). :class:`FlipFlop` can also determine whether each event was
successful (the molecule resides in the opposing leaflet for at least a given
length of time), or not (the molecule went to the midplane but returned to its
original leaflet).

See `Baral et al. (2020) <https://www.sciencedirect.com/science/article/pii/S0009308420300980>`_
for further discussion on flip-flop in lipid bilayers, including the affect on the flip-flop
rate of the buffer size used to assign molecules to the midplane of the bilayer.

Input
-----

Required:
  - *universe* : an MDAnalysis Universe object.
  - *lipid_sel* : atom selection for atoms to use in detecting flip-flop
  - *leaflets* : leaflet membership (-1: lower leaflet, 0: midplane, 1: upper leaflet) of each lipid in the membrane at each frame
  

Output
------

  - *resindex* : residue index of a flip-flopping molecule
  - *flip_flop_start_frame* : final frame at which the molecule was present in its original leaflet
  - *flip_flop_end_frame* : first frame at which the molecule is present in the new leaflet
  - *moves_to* : direction of travel of the molecule: equal to 1 if the upper leaflet is the new lealet, equal to -1 if the lower leaflet is the new leaflet
  
Flip-flop data area returned in a :class:`numpy.ndarray`, on a "one line, one observation" basis
and can be accessed via :attr:`FlipFlop.flip_flops`::

    flip_flops = [
        [
            <resindex (0-based)>,
            <end_frame (0-based)>,
            <start_frame (0-based)>,
            <moves_to>
        ],
        ...
    ]
    
*moves_to* is equal to 1 or -1 if the molecule flip-flops into the upper or the
lower leaflet, respectively.
    
Additionaly, the success or failure of each flip-flop event is stored in the
attribute :attr:`FlipFlop.flip_flop_success`.

Example usage of :class:`FlipFlop`
----------------------------------

An MDAnalysis Universe must first be created before using :class:`FlipFlop`::

  import MDAnalysis as mda
  from lipyphilic.lib.assign_leaflets import AssignLeaflets
  from lipyphilic.lib.flip_flop import FlipFlop

  u = mda.Universe(tpr, trajectory)

Then we need to know which leaflet each lipid is in at each frame. This may be done using
:class:`lipyphilic.lib.assign_leaflets.AssignLeaflets`::

  leaflets = AssignLeaflets(
    universe=u,
    lipid_sel="name GL1 GL2 ROH"  # assuming we are using the MARTINI forcefield
    midplane_sel="name ROH",       # only cholesterol is allowed to flip-flop
    midplane_cutoff=8.0,          # buffer size for assigning molecules to the midplane
  )
  leaflets.run()

The leaflet data are stored in the :attr:`leaflets.leaflets` attribute. We can now create our
:class:`FlipFlop` object::

  flip_flop = FlipFlop(
    universe=u,
    lipid_sel="name GL1 GL2 ROH", # this must be the same as used in AssignLeaflets
    leaflets=leaflets.leaflets
  )
  
We then select which frames of the trajectory to analyse (`None` will use every
frame)::
  
  flip_flop.run(
    start=None,
    stop=None,
    step=None
  )
  
The results are then available in the :attr:`flipflop.flip_flop` attribute as a
:class:`numpy.ndarray`. Each row corresponds to an individual flip-flop event, and
the four columns correspond, respectively, to the molecule resindex,
flip-flop start frame, flip-flop end frame, and the leaflet in which the molecule
resides after the flip-flop.

Specify minimum residence time for successful flip-flops
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can also specify the minumum number of frames a molecule must reside in its new leaflet
for the flip-flop to be considered successful. We do this using the :attr:`frame_cutoff`
parameter::

  flip_flop = FlipFlop(
    universe=u,
    lipid_sel="name GL1 GL2 ROH",
    leaflets=leaflets.leaflets,
    frame_cuotff=10,
  )

With *frame_cutoff=10*, a molecule must remain in its new leaflet for at least 10
consecutive frames for the flip-flop to be considered successful. If this condition is not met,
the flip-flop event is recorded as failing.

Calculating the flip-flop rate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The flip-flop rate can be calculatd directly from the number of successfull flip-flop evetns,
which itself can be calculated as::

  n_successful = sum(flip_flop.flip_flop_success == "Success")
  
The rate is then given by the total number of successful flip-flops divided by the total
simulations time and the number of molecules of the translocating species.

The class and its methods
-------------------------

.. autoclass:: FlipFlop
    :members:
    :exclude-members: run

"""
from tqdm.auto import tqdm
import numpy as np

from lipyphilic.lib import base


class FlipFlop(base.AnalysisBase):
    """Find flip-flop events in a lipid bilayer.
    """

    def __init__(self, universe,
                 lipid_sel,
                 leaflets,
                 frame_cutoff=1
                 ):
        """Set up parameters for finding flip-flop events.
        
        Parameters
        ----------
        universe : Universe
            MDAnalysis Universe object
        lipid_sel : str
            Selection string for atoms to use in detecting flip-flop.
        leaflets : numpy.ndarray (n_lipids,, n_frames)
            An array of leaflet membership for each lipid as each frame, in which: -1
            corresponds to the lower leaflet; 1 corresponds to the upper leaflet; and
            0 corresponds to the midplane.
        frame_cutoff : int, optional
            To be counted as a successful flip-flop, a molecule must reside in its new
            leaflet for at least 'frame_cutoff' consecutive frames. The default is `1`, in
            which case the molecule only needs to move to the opposing leaflet for a single
            frame for the flip-flop to be successful.
            
        Tip
        ---
        
        Leaflet membership can be determined using :class:`lipyphilic.lib.assign_leaflets.AssignLeaflets`.
        
        """
        super(FlipFlop, self).__init__(universe.trajectory)
        
        self.u = universe
        self.membrane = self.u.select_atoms(lipid_sel, updating=False)
        
        if (np.array(leaflets).ndim != 2) or (len(leaflets) != self.membrane.n_residues):
            raise ValueError("'leaflets' must be a 2D array of shape (n_residues, n_frames)"
                             " containing the leaflet id of each lipid at each frame."
                             )
        
        self.leaflets = np.array(leaflets)
        
        if frame_cutoff < 1:
            raise ValueError("'frame_cutoff' must be greater than or equal to 1")
           
        self.frame_cutoff = frame_cutoff
        
        self.flip_flops = None
        self.flip_flop_success = None
        
    def _setup_frames(self, trajectory, start=None, stop=None, step=None):
        """
        Pass a Reader object and define the desired iteration pattern
        through the trajectory

        Parameters
        ----------
        trajectory : MDAnalysis.Reader
            A trajectory Reader
        start : int, optional
            start frame of analysis
        stop : int, optional
            stop frame of analysis
        step : int, optional
            number of frames to skip between each analysed frame
        
        """
        self._trajectory = trajectory
        start, stop, step = trajectory.check_slice_indices(start, stop, step)
        self.start = start
        self.stop = stop
        self.step = step
        
        n_frames = len(range(start, stop, step))
        if self.leaflets.shape[1] != n_frames:
            raise ValueError("The frames to analyse must be identical to those used "
                             "in assigning lipids to leaflets."
                             )
        
        self.n_frames = n_frames
        self.frames = np.arange(start, stop, step)
          
    def _prepare(self):
        
        # Output array
        self.flip_flops = [[], [], [], []]
        self.flip_flop_success = []
        
    def _single_frame(self):
        
        # Skip if the molecule never changes leaflet
        if np.min(np.diff(self._residue_leaflets)) == np.max(np.diff(self._residue_leaflets)) == 0:
            return None
        
        # Check when the molecule leaves the upper leaflet for more than `frame_cutoff` frames
        # If it has left for at least this many frames, there is a chance it has
        # flip-flopped to the opposing leaflet
        leaflet = 1
        residue_leaflet_indices = np.nonzero(self._residue_leaflets == leaflet)[0]
        gaps = np.diff(residue_leaflet_indices) > self.frame_cutoff
        
        if gaps.size > 0:
            upper_begins = np.insert(residue_leaflet_indices[1:][gaps], 0, residue_leaflet_indices[0])
            upper_ends = np.append(residue_leaflet_indices[:-1][gaps], residue_leaflet_indices[-1])
        else:
            upper_begins = residue_leaflet_indices[:]
            upper_ends = residue_leaflet_indices[:]
        
        # Check when the molecule leaves the lower leaflet for more than `frame_cutoff` frames
        leaflet = -1
        residue_leaflet_indices = np.nonzero(self._residue_leaflets == leaflet)[0]
        gaps = np.diff(residue_leaflet_indices) > self.frame_cutoff

        if gaps.size > 0:
            lower_begins = np.insert(residue_leaflet_indices[1:][gaps], 0, residue_leaflet_indices[0])
            lower_ends = np.append(residue_leaflet_indices[:-1][gaps], residue_leaflet_indices[-1])
        else:
            lower_begins = residue_leaflet_indices[:]
            lower_ends = residue_leaflet_indices[:]

        # Combine beginnings and ending of leaflet membership
        begins = np.array(list(lower_begins) + list(upper_begins))
        ends = np.array(list(lower_ends) + list(upper_ends))
        moves_to = np.array([-1] * len(lower_begins) + [1] * len(upper_begins))

        sort = np.argsort(begins)
        begins = begins[sort]
        ends = ends[sort]
        moves_to = moves_to[sort]

        # To be considered to have flip-flopped, the molecule must remain in the leaflet for at least `frame_cutoff` frames
        keep = (ends - begins) >= (self.frame_cutoff - 1)  # if end==begin, the molecule was there for 1 frame not 0 frames
        begins = begins[keep]
        ends = ends[keep]
        moves_to = moves_to[keep]

        # Store data for when leaflet membership ends, begins, and which leaflet it has moved to
        resindex = self.membrane[self._residue_index].resindex
        
        self.flip_flops[0].extend(np.full_like(begins[1:], fill_value=resindex))
        self.flip_flops[1].extend(self.frames[ends[:-1]])
        self.flip_flops[2].extend(self.frames[begins[1:]])
        self.flip_flops[3].extend(moves_to[1:])
        
        # Check whether the flip-flop was a success or if it moved back to its
        # original leaflet
        success = np.full_like(moves_to[1:], fill_value="Success", dtype="<U10")
        success[np.diff(moves_to) == 0] = "Fail"
        self.flip_flop_success.extend(success)
        
    def _conclude(self):
    
        self.flip_flops = np.asarray(self.flip_flops).T
        self.flip_flop_success = np.asarray(self.flip_flop_success)
    
    def run(self, start=None, stop=None, step=None):
        """Perform the calculation

        Parameters
        ----------
        start : int, optional
            start frame of analysis
        stop : int, optional
            stop frame of analysis
        step : int, optional
            number of frames to skip between each analysed frame
        """
        
        self._setup_frames(self._trajectory, start, stop, step)
        self._prepare()
        
        for residue_index, residue_leaflets in tqdm(enumerate(
                self.leaflets), total=self.membrane.n_residues):
            
            self._residue_index = residue_index
            self._residue_leaflets = residue_leaflets
            self._single_frame()
            
        self._conclude()
        return self
