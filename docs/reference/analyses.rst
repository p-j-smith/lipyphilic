 .. _Analysis-tools:

Analysis tools
==============

Here we provide a description of the analysis tools currently available in **lipyphilic**,
along with examples of how to perform each analysis.

Assign leaflets: :mod:`lipyphilic.lib.assign_leaflets`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This module provides methods for assigning lipids to leaflets in a bilayer. Leaflet
assignment is based on the distance in *z* from a lipid the midpoint of the bilayer.
Lipids may be assigned to the upper leaflet (indicated by `1`), the lower leaflet (`-1`)
or the bilayer midplane (`0`).

Below we see how to assign lipids to

.. code:: python

	import MDAnalysis as mda
	from lipyphilic.lib.assign_leaflets import AssignLeaflets

	# Load an MDAnalysis Universe
	u = mda.Universe('production.tpr','production.xtc')

	# Find which leaflet each lipid is in at each frame
	leaflets = AssignLeaflets(
	    universe=universe,
	    lipid_sel="name PO4 ROH" 
	)
	
	# Select which frames to use and perform the analysis
	leaflets.run(start=None, stop=None, step=None)  # this will use every frame in the trajectory

The results are stored as a NumPy array of shape (n_lipids, n_frames) in the
:attr:`leaflets.leaflest` attribute.

If you have used a different force field, you simply need to change the :attr:`lipid_sel` to
select the relevant headgroup atoms of your lipids. See the `MDAnalysis selection language
<https://userguide.mdanalysis.org/stable/selections.html>`__ for more info on how to select atoms.

By default, lipids are only allowed to be in the upper (`1`) or lower (`-1`) leaflet. See
:mod:`lipyphilic.lib.assign_leaflets` for more information on selecting which molecules are allowed
in the midplane.


Flip-flop: :mod:`lipyphilic.lib.flip_flop`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This module provides methods for detecting flip-flop of molecules in a lipid bilayer. A flip-flop
occurs when a molecule - typically a sterol - moves from one leaflet of a bilayer into the opposing
leaflet.

To find all flip-flop events, we first should assign lipids to leaflets as seen in the above example,
then:

.. code:: python

  flip_flops = FlipFlop(
      universe=u,
      lipid_sel="name GL1 GL2 ROH", # this must be the same as used in AssignLeaflets
      leaflets=leaflets.leaflets    # pass the NumPy array of leaflet ids
  )
    
  flip_flops.run(start=None, stop=None, step=None)

The results are stored as a NumPy array of shape (n_flip_flops, 4) in the
:attr:`flip_flops.flip_flops` attribute. Each row is a single flip-flop event, and the four columns
correspond to: the residue index of the flip-flopping molecule; the frame at which the molecule
left its original leaflet; the frame at which it entered its new leaflet; the leaflet ID to which
it moves.

See :class:`lipyphilic.lib.flip_flop` for more information on how flip-flop is detected and options such
as specifying how long a molecule must residue in the new leaflet for the flip-flop to be considered
successful.

Neighbours: :mod:`lipyphilic.lib.neighbours`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



Interlealet registration: :mod:`lipyphilic.lib.registration`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 


Area per lipid: :mod:`lipyphilic.lib.area_per_lipid`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^






