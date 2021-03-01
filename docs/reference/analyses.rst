 .. _Analysis-tools:

Analysis tools
==============

Here we provide a brief description of the analysis tools currently available in **lipyphilic**,
along with examples of how to perform each analysis. For full information on each analysis
tool, including the APIs, see the following pages:

* Assign leaflets: :mod:`lipyphilic.lib.assign_leaflets`
* Flip-flop: :mod:`lipyphilic.lib.flip_flop`
* Interlealet registration: :mod:`lipyphilic.lib.registration`
* Neighbours: :mod:`lipyphilic.lib.neighbours`
* Area per lipid: :mod:`lipyphilic.lib.area_per_lipid`


Assign leaflets: :mod:`lipyphilic.lib.assign_leaflets`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This module provides methods for assigning lipids to leaflets in a bilayer. Leaflet
assignment is based on the distance in *z* from a lipid the midpoint of the bilayer.
Lipids may be assigned to the upper leaflet (indicated by `1`), the lower leaflet (`-1`)
or the bilayer midplane (`0`).

Below we see how to assign lipids to the upper or lower leaflet of a `MARTINI
<http://cgmartini.nl/>`__ bilayer:

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

.. note::

  Assingment of lipids to leaflets is not in itself useful, but it is required in order to calculate,
  for example, area per lipid, interleaflet correlations, and flip-flop rates.


Flip-flop: :mod:`lipyphilic.lib.flip_flop`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This module provides methods for detecting flip-flop of molecules in a lipid bilayer. A flip-flop
occurs when a molecule - typically a sterol - moves from one leaflet of a bilayer into the opposing
leaflet.

To find all flip-flop events, we first should assign lipids to leaflets as seen in the above example,
then:

.. code:: python

  from lipyphilic.lib.flip_flop import FlipFlop

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

See :mod:`lipyphilic.lib.flip_flop` for more information on how flip-flop is detected and options such
as specifying how long a molecule must residue in the new leaflet for the flip-flop to be considered
successful.


Interlealet registration: :mod:`lipyphilic.lib.registration`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 

This module provides methods for determining registration of leaflets in a bilayer. Registration is
defined by the pearson correlation coefficient of molecular densities in the two leaflets. This is
an implementation of the method described by `Thallmair et al. (2018)
<https://pubs.acs.org/doi/abs/10.1021/acs.jpclett.8b01877>`__.

To calculate the interleaflet correlation of cholesterol, we first need to calculate which leaflet each
lipid is in at each frame using :class:`lipyphilic.lib.assign_leaflets.AssignLeaflets`. Then we pass
the :class:`AssignLeaflets` object, along with atom selections for which density correlations will
be calculated, to :class:`Registration`:

.. code:: python

  from lipyphilic.lib.registration import Registration

  registration = Registration(
      leaflets=leaflets,
      upper_sel="resname CHOL and name ROH",
      lower_sel="resname CHOL and name ROH",
  )
  
  registration.run(start=None, stop=None, step=None)

The results are stored in a NumPy array of shape (n_frames), containing the pearson correlation
coefficient of cholesterol densities in the two leaflets. The data are accessible via the
:attr:`registration.registration` attribute.

As well as calcualting registration of lipid species across the two leaflets, it is also possible
to calculate the registration of arbitrary user-defined values across the two leaflets. For example,
if you have created a `Hidden Markov Model to assign lipids to the Ld or Lo phase
<https://pubs.acs.org/doi/abs/10.1021/acs.jctc.8b00828>`__, you can calculate the registration of
Lo lipids across the two leaflets. See :mod:`lipyphilic.lib.registration` for more details.


Neighbours: :mod:`lipyphilic.lib.neighbours`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This module provides methods for finding neighbouring lipids in a bilayer. Lipids are neighbours if
they are within a user-defined cutoff of one another.

Below we see how to find all neighbours in a MARTINI bilayer based on the 'GL1' and 'GL2' beads of
phospholipids and the 'ROH' bead of sterols, using a cutoff of *12* Å:

.. code:: python

	import MDAnalysis as mda
	from lipyphilic.lib.neighbours import Neighbours

	# Load an MDAnalysis Universe
	u = mda.Universe('production.tpr','production.xtc')

	# Find neighbouring lipids
	neighbours = Neighbours(
	    universe=u,
	    lipid_sel="name PO4 ROH",
		cutoff=12.0
	)
	
	# Select which frames to use and perform the analysis
	neighbours.run(start=None, stop=None, step=None)

The results are stored as a :class:`scipy.sparse.csc_matrix` in the :attr:`neighbours.neighbours`
attribute.

.. tip::

  Once the neighbour matrix has been generated, the local lipid compositions largest lipids cluster
  at each frame can be readily calculated.

See :mod:`lipyphilic.lib.neighbours` for more information on this module, including how to calculate
local lipid compositions or find the largest cluster of lipid species over time.


Area per lipid: :mod:`lipyphilic.lib.area_per_lipid`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This module provides methods for calculating the area per lipid. Areas are calculated via a 2D
Voronoi tessellation, using the `locality` module of
`Freud <https://freud.readthedocs.io/en/stable/index.html#>`_ to perform the tessellation
of atomic positions. See `Lukat et al. (2013) <https://pubs.acs.org/doi/full/10.1021/ci400172g>`_
a thorough description of calculating the area per lipid via Voronoi tessellations.

Once lipids have been assigned to leaflets, the area per lipid can be calculated as follows:

.. code:: python

  from lipyphilic.lib.area_per_lipid import AreaPerLipid

  areas = AreaPerLipid(
      universe=u,
      lipid_sel="name GL1 GL2 ROH",  # assuming we're using the MARTINI forcefield
      leaflets=leaflets.leaflets
  )
  
The above will use GL1 and GL2 beads to calculate the area of each phospholipid, and the
ROH bead to calculate the area of each sterol.

For a more complete description of calculating the area per lipid, and the API of the
analysis class, see :mod:`lipyphilic.lib.area_per_lipid`.
