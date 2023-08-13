# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# lipyphilic --- lipyphilic.readthedocs.io
#
# Released under the GNU Public Licence, v2 or any higher version
#

"""Flip-flop --- :mod:`lipyphilic.analysis.flip_flop`
=====================================================

:Author: Paul Smith
:Year: 2021
:Copyright: GNU Public License v2

This module provides methods for finding flip-flop events in a lipid bilayer.

A flip-flop event occurs when a molecule - typically a sterol - moves from
one leaflet of a bilayer into the opposing leaflet.

The class :class:`lipyphilic.analysis.flip_flop.FlipFlop` finds the frames at which
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
  from lipyphilic.leaflets.assign_leaflets import AssignLeaflets
  from lipyphilic.analysis.flip_flop import FlipFlop

  u = mda.Universe(tpr, trajectory)

Then we need to know which leaflet each lipid is in at each frame. This may be done using
:class:`lipyphilic.leaflets.assign_leaflets.AssignLeaflets`::

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
    lipid_sel="name ROH",
    leaflets=leaflets.filter_leaflets("name ROH")  # pass only the relevant leaflet data
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
    lipid_sel="name ROH",
    leaflets=leaflets.filter_leaflets("name ROH")
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
from MDAnalysis.analysis.base import AnalysisBase
import numpy as np
from tqdm.auto import tqdm

from lipyphilic import _lipyferric

__all__ = [
    "FlipFlop",
]


class FlipFlop(AnalysisBase):
    """Find flip-flop events in a lipid bilayer."""

    def __init__(
        self,
        universe,
        lipid_sel,
        leaflets,
        frame_cutoff=1,
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

        Leaflet membership can be determined using :class:`lipyphilic.leaflets.assign_leaflets.AssignLeaflets`.

        """
        super().__init__(universe.trajectory)

        self.u = universe
        self.membrane = self.u.select_atoms(lipid_sel, updating=False)

        if (np.array(leaflets).ndim != 2) or (len(leaflets) != self.membrane.n_residues):
            _msg = (
                "'leaflets' must be a 2D array of shape (n_residues, n_frames)"
                " containing the leaflet id of each lipid at each frame.",
            )
            raise ValueError(_msg)

        self.leaflets = np.array(leaflets)

        if frame_cutoff < 1:
            _msg = "'frame_cutoff' must be greater than or equal to 1"
            raise ValueError(_msg)

        self.frame_cutoff = frame_cutoff

        self.results.flip_flops = None
        self.results.flip_flop_success = None

    @property
    def flip_flops(self):
        return self.results.flip_flops

    @property
    def flip_flop_success(self):
        return self.results.flip_flop_success

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
            _msg = "The frames to analyse must be identical to those used in assigning lipids to leaflets."
            raise ValueError(_msg)

        self.n_frames = n_frames
        self.frames = np.arange(start, stop, step)

    def _prepare(self):
        # Output array
        self.results.flip_flops = [[], [], [], []]
        self.results.flip_flop_success = []

    def _single_frame(self):
        # Skip if the molecule never changes leaflet
        if np.min(np.diff(self._residue_leaflets)) == np.max(np.diff(self._residue_leaflets)) == 0:
            return

        start_frames, end_frames, end_leaflets, success = _lipyferric.molecule_flip_flop(
            leaflets=self._residue_leaflets,
            frame_cutoff=self.frame_cutoff,
        )

        n_events = len(start_frames)
        resindex = self.membrane[self._residue_index].resindex
        self.results.flip_flops[0].extend([resindex for _ in range(n_events)])
        self.results.flip_flops[1].extend(start_frames)
        self.results.flip_flops[2].extend(end_frames)
        self.results.flip_flops[3].extend(end_leaflets)
        self.results.flip_flop_success.extend(success)

    def _conclude(self):
        self.results.flip_flops = np.asarray(self.results.flip_flops).T
        self.results.flip_flop_success = np.asarray(self.results.flip_flop_success)

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

        for residue_index, residue_leaflets in tqdm(enumerate(self.leaflets), total=self.membrane.n_residues):
            self._residue_index = residue_index
            self._residue_leaflets = residue_leaflets
            self._single_frame()

        self._conclude()
        return self
