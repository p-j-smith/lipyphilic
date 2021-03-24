 .. _Analysis-tools:

Analysis tools
==============

Here we provide a brief description of the analysis tools currently available in **lipyphilic**,
along with examples of how to perform each analysis. For more information on each analysis tool,
including details of all optional input parameters and further example use cases, see the :ref:`API`.


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
	  universe=u,
	  lipid_sel="name PO4 ROH" 
	)
	
	# Select which frames to use and perform the analysis
	leaflets.run(start=None, stop=None, step=None)  # this will use every frame in the trajectory


The results are stored as a NumPy array of shape (n_lipids, n_frames) in the
:attr:`leaflets.leaflets` attribute.

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

This module provides methods for detecting the flip-flop of molecules in a lipid bilayer. A flip-flop
occurs when a molecule - typically a sterol - moves from one leaflet of a bilayer into the opposing
leaflet.

To find all flip-flop events, we first should assign lipids to leaflets as seen in the above example,
then:

.. code:: python

  import MDAnalysis as mda
  from lipyphilic.lib.flip_flop import FlipFlop

  # Load an MDAnalysis Universe
  u = mda.Universe('production.tpr','production.xtc')

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

  import MDAnalysis as mda
  from lipyphilic.lib.registration import Registration

  # Load an MDAnalysis Universe
  u = mda.Universe('production.tpr','production.xtc')

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
	  lipid_sel="name GL1 GL2 ROH",
		cutoff=12.0
	)
	
	neighbours.run(start=None, stop=None, step=None)

The results are stored as a :class:`scipy.sparse.csc_matrix` in the :attr:`neighbours.neighbours`
attribute.

.. tip::

  Once the neighbour matrix has been generated, the local lipid compositions  or  the largest lipids cluster
  at each frame can be readily.

See :mod:`lipyphilic.lib.neighbours` for more information on this module, including how to calculate
local lipid compositions or the lipid enrichment/depletion index, and how to find the largest cluster of
a given lipid species over time.


Area per lipid: :mod:`lipyphilic.lib.area_per_lipid`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This module provides methods for calculating the area per lipid. Areas are calculated via a 2D
Voronoi tessellation, using the `locality` module of
`Freud <https://freud.readthedocs.io/en/stable/index.html#>`_ to perform the tessellation
of atomic positions. See `Lukat et al. (2013) <https://pubs.acs.org/doi/full/10.1021/ci400172g>`_
a thorough description of calculating the area per lipid via Voronoi tessellations.

Once lipids have been assigned to leaflets, the area per lipid can be calculated as follows:

.. code:: python

  import MDAnalysis as mda
  from lipyphilic.lib.area_per_lipid import AreaPerLipid

  # Load an MDAnalysis Universe
  u = mda.Universe('production.tpr','production.xtc')

  areas = AreaPerLipid(
    universe=u,
    lipid_sel="name GL1 GL2 ROH",  # assuming we're using the MARTINI forcefield
    leaflets=leaflets.leaflets
  )

  areas.run(start=None, stop=None, step=None)
  
The above will use GL1 and GL2 beads to calculate the area of each phospholipid, and the
ROH bead to calculate the area of each sterol.

For a more complete description of calculating the area per lipid, and the API of the
analysis class, see :mod:`lipyphilic.lib.area_per_lipid`.


Lipid order parameter --- :mod:`lipyphilic.lib.order_parameter`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This module provides methods for calculating the coarse-grained orientational order
parameter of acyl tails in a lipid bilayer. The coarse-grained order parameter, :math:`S_{CC}`,
is a measure of the degree of ordering of an acyl tail, based on the extent
to which the vector connecting two consecutive tail beads is aligned with the membrane
normal.

See `Seo et al. (2020) <https://pubs.acs.org/doi/full/10.1021/acs.jpclett.0c01317>`__ for
a definition of :math:`S_{CC}` and `Piggot et al. (2017)
<https://pubs.acs.org/doi/full/10.1021/acs.jctc.7b00643>`__ for an excellent discussion
on acyl tail order parameters in molecular dynamics simulations.

To calculate :math:`S_{CC}`, we need to provide an atom selection for the beads
in a **single** tail of lipids in the bilayer --- that is, **either** the *sn1* or *sn2*
tails, not both. If we have performed a MARTINI simulation, we can calculate the
:math:`S_{CC}` of all *sn1* tails of phospholipids as follows:

.. code:: python

  import MDAnalysis as mda
  from lipyphilic.lib.order_parameter import SCC

  # Load an MDAnalysis Universe
  u = mda.Universe('production.tpr','production.xtc')

  scc = SCC(
    universe=u,
    tail_sel="name ??A"
  )
  
The above makes use of the powerful `MDAnalysis selection language
<https://userguide.mdanalysis.org/stable/selections.html>`__. It will select beads such as
*C1A*, *C2A*, *D2A* etc. This makes it simple to quickly calculate
:math:`S_{CC}` for the *sn1* tails of all species in a bilayer.

To see how to calculate :math:`S_{CC}` using local membrane normals to define the molecular axes,
as well as the full API of the class, see :mod:`lipyphilic.lib.order_parameter`.


Lipid :math:`z` angles: :mod:`lipyphilic.lib.z_angles`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This module provides methods for calculating the angle lipids make with the
positive :math:`z` axis. If we define the orientation of MARTINI cholesterol as the
angle between the :math:`z`-axis and the vector from the the 'R5' bead to the 'ROH' bead,
we can calculate the orientation of each cholesterol molecule as follows:

.. code:: python

  import MDAnalysis as mda
  from lipyphilic.lib.z_angles import ZAngles

  # Load an MDAnalysis Universe
  u = mda.Universe('production.tpr','production.xtc')

  z_angles = ZAngles(
    universe=u,
    atom_A_sel="name R5",
    atom_B_sel="name ROH"
  )

  z_angles.run(start=None, stop=None, step=None)

The results are stored in a :class:`numpy.ndarray` of shape (n_residues, n_lipids) in the
:attr:`z_angles.z_angles` attribute.

For more information on this module, including how to return the angles in radians rather
than degrees, see :mod:`lipyphilic.lib.z_angles`.


Lipid :math:`z` positions: :mod:`lipyphilic.lib.z_positions`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This module provides methods for calculating the height in :math:`z` of lipids from the
bilayer center.

If we have used the MARTINI forcefield to study a phospholipid/cholesterol mixture,
we can calculate the height of cholesterol in the bilayer as follows:

.. code:: python

  import MDAnalysis as mda
  from lipyphilic.lib.z_positions import ZPositions

  # Load an MDAnalysis Universe
  u = mda.Universe('production.tpr','production.xtc')

  z_positions = ZPositions(
    universe=u,
    lipid_sel="name GL1 GL2 ROH",
    height_sel="name ROH",
    n_bins=10
  )

  z_positions.run(start=None, stop=None, step=None)

:attr:`lipid_sel` is an atom selection that covers all lipids in the bilayer. This
is used for calculating the membrane midpoint. :attr:`height_sel` selects which
atoms to use for caclulating the height of each lipid.

Local membrane midpoints are calculated by creating a grid of
membrane patches, with the number of grid points controlled with the :attr:`n_bins`
parameter. The distance in :math:`z` of each lipid to its local midpoint is then calculated.

Data are returned in a :class:`numpy.ndarray` of shape (n_residues, n_frames). See
:mod:`lipyphilic.lib.z_positions` for more information on this module including the
full API of the class.

Lipid :math:`z` thickness: :mod:`lipyphilic.lib.z_thickness`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This module provides methods for calculating the thickness, in :math:`z`, of lipid tails.
This is defined as the maximum distance in :math:`z` between to atoms in a tail.

If we have used the MARTINI forcefield to study a DPPC/DOPC/cholesterol mixture,
we can calculate the thickness of DPPC and DOPC *sn1* tails, as well as the thickness
of cholesterol, as follows:

.. code:: python

  import MDAnalysis as mda
  from lipyphilic.lib.z_positions import ZThickness

  # Load an MDAnalysis Universe
  u = mda.Universe('production.tpr','production.xtc')

  z_thickness = ZThickness(
    universe=u,
    lipid_sel="(name ??1 ??A) or (resname CHOL and not name ROH)"
  )

  z_thickness.run()

The above makes use of the powerful MDAnalysis atom selection language to select the DPPC
and DOPC sn1 tails along with cholesterol.

The thickness data are stored in a :class:`numpy.ndarray` of shape (n_residues, n_frames)
in the :attr:`z_thickness.z_thickness` attribute. See :mod:`lipyphilic.lib.z_thickness` for
the full API of the class.

Membrane :math:`z` thickness: :mod:`lipyphilic.lib.memb_thickness`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This module provides methods for calculating the bilayer thickness. It is defined as the
peak-to-peak distance of lipid headgroup density in :math:`z`.

Lipids must first be assigned to the upper and lower leaflets. This can be done with the
class :class:`lipyphilic.lib.assign_leaflets.AssignLeaflets`. Then, to calculate the membrane
thickness we need to define which atoms to treat as headgroup atoms and pass the leaflet
membership information to :class:`MembThickness`. If we have studied a DPPC/DOPC/cholesterol
mixture with MARTINI, we could calculate the membrane thickness as follows:

.. code:: python

  import MDAnalysis as mda
  from lipyphilic.lib.z_positions import ZThickness

  # Load an MDAnalysis Universe
  u = mda.Universe('production.tpr','production.xtc')

  memb_thickness = MembThickness(
    universe=u,
    leaflets=leaflets.filter_leaflets("resname DOPC and DPPC"),  # exclude cholesterol from thickness calculation
    lipid_sel="resname DPPC DOPC and name PO4"
  )

  memb_thickness.run()

The results are then available in the :attr:`memb_thickness.memb_thickness` attribute as a
:class:`numpy.ndarray`.

For more information on calculating membrane thickness, including options to calculating local
membrane thicknesses rather than a single global thickness, see :mod:`lipyphilic.lib.memb_thickness`.


Lateral diffusion :mod:`lipyphilic.lib.lateral_diffusion`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This module contains methods for calculating the mean squared displacement (MSD) and lateral
diffusion coefficient, :math:`D_{xy}`,of lipids in a bilayer.

The MSD of all lipids in a DPPC/DOPC/cholesterol MARTINI bilayer can be calculated using 
:class:`lipyphilic.lib.lateral_diffusion.MSD`:

.. code:: python

 import MDAnalysis as mda
 from lipyphilic.lib.lateral_diffusion import MSD

 # Load an MDAnalysis Universe
 u = mda.Universe('production.tpr','production.xtc')

 msd = MSD(
   universe=u,
   lipid_sel="name PO4 ROH"
 )

  msd.run()


The MSD of each lipid is then available in the :attr:`msd.msd` attribute as a :class:`numpy.ndarray`,
and the lagtimes are stored in the :attr:`msd.lagtimes` attribute.

For more information on this module, including how to calculate the lateral diffusion coefficient,
see :mod:`lipyphilic.lib.lateral_diffusion`.


Plotting utilities: :mod:`lipyphilic.lib.plotting`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lipyphilic** can produce joint probablity density plots (or PMFs if a temperature is provided),
as well as density maps of membrane propertes projected onto the membrane plane. The former may be
used to plot, for example, the PMF of cholesterol orientation and height in a bilayer. The latter
may be used to generate plots of, for example, the area per lipid as a function of :math:`xy` in
the membrane plane.

See :mod:`lipyphilic.lib.plotting` for the full API of :class:`lipyphilic.lib.plotting.JointDensity`
and :class:`lipyphilic.lib.plotting.ProjectionPlot`.


On-the-fly transformations :mod:`lipyphilic.transformations`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`lipyphilic` contains a module for applying on-the-fly transofrmation to atomic coordinates
while iterating over a trajectory. These are availbale in the module :mod:`lipyphilic.transformations`.

There are two transformations available in `lipyphilic`:

1. :class:`lipyphilic.transformations.nojump`, which prevents atoms from jumping across periodic boundaries. This is useful when calculating the lateral diffusion of lipids.
2. :class:`lipyphilic.transformations.center_membrane`, which can take a membrane that is split across periodic boundaries, make it whole and center it in the box.

See :mod:`lipyphilic.transformations` for full details on these transformations including how to apply
them to your trajectory.