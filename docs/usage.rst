Basic Usage
===========

The analysis tools in **lipyphilic** all require an `MDAnalysis Universe
<https://userguide.mdanalysis.org/stable/universe.html>`__ as input, so to use **lipyphilic** you will also
need to import MDAnalysis. The analyses are then performed in the same way as the majority of those
in MDAnalysis. For example, to assign each lipid to the upper or lower leaflet at each frame in a trajectory:

.. code:: python

	import MDAnalysis as mda
	from lipyphilic.lib.assign_leaflets import AssignLeaflets

	# Load an MDAnalysis Universe
	u = mda.Universe('production.tpr','production.xtc')

	# Find which leaflet each lipid is in at each frame
	leaflets = AssignLeaflets(
	    universe=u,
	    lipid_sel="name PO4 ROH"  # Select headgroup beads in the MARTINI forcefield
	)
	
	# Select which frames to use and perform the analysis
	leaflets.run(start=None, stop=None, step=None)  # this will use every frame in the trajectory


And the results will be available as a NumPy array stored in the `leaflets.leaflets` attribute.

For more details on how to use **lipyphilic**, check out all the :ref:`Analysis-tools` we
currently have available.
