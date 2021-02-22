Usage
=====

**lipyphilic** is built on top of MDAnalysis, so to use **lipyphilic** you will also need to import MDAnalysis.

.. code:: python

   import MDAnalysis as mda
   import lipyphilic


Then you will have access to the analysis tools. The analyses are performed in the same way as the majority of those
in MDAnalysis, for example, to assign lipid to the upper or lower leaflet at each frame in a trajectory:

.. code:: python

	import MDAnalysis as mda
	from lipyphilic.lib.assign_leaflets import AssignLeaflets

	# Load an MDAnalysis Universe
	u = mda.Universe('produciton.tpr','produciotn.xtc')

	# Find which leaflet each lipid is in at each frame
	leaflets = AssignLeaflets(
	    universe=universe,
	    lipid_sel="name GL1 GL2 ROH"  # Select headgroup beads in the MARTINI forcefield
	)
	leaflets.run()


And the results will be available as a NumPy array stored in the `leaflets.leaflets` attribute.

For more details on how to use **lipyphilic**, check out all the :ref:`Analysis-tools` we
currently have available.
