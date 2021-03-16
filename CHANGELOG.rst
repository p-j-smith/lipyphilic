lipyphilic CHANGELOG
====================

0.5.0 (2021-03-16)
------------------

* PR#38 Add a trajectory transformation for unwrapping broken membranes (Fixes #37)
* PR#36 Add method for projecting areas onto the membrane plane (Fixes #33)
* PR#35 Added a tool for calculating membrane thickness (Fixes #34)
* PR#32 ZThickness.average() now returns a new ZThickness object rather than a NumPy array
* PR#31 SCC.weighted_average() now returns a new SCC object rather than a NumPy array
* PR#30 Add class for plotting projections of membrane properties onto the xy plane.
* PR#29 Added plotting of joint probability distributions or PMFs (Fixed #28).

0.4.0 (2021-03-05)
------------------

* PR#26 Added a tool to calculate the thickness of lipids or their tails (Fixes #25)
* PR#24 Added a tool to calculate the coarse-grained order parameter (Fixes #23)
* PR#22 Added a tool to calculate orientation of lipids in a bilayer (Fixes #20)
* PR#21 Added a tool to calculate lipid height in a bilayer (Fixes #19)
* Better description of analysis tools in the docs
* Updated installation instructions, including installing via conda-forge

0.3.2 (2021-02-27)
------------------

* Fix typo in requirements

0.3.1 (2021-02-27)
------------------

* Add support for numpy 1.20

0.3.0 (2021-02-26)
------------------

* Fix neighbour calculation for non-sequential residue indices
  Fixes #11
* Added a tool to calculate interleaflet registration

0.2.0 (2021-02-23)
------------------

* Improved documentation
* Add method to count number of each neighbour type
* Add functionality to find neighbouring lipids

0.1.0 (2021-02-17)
------------------

* Add functionality to find flip-flop events in bilayers
* Add functionality to calculate area per lipid
* Add functionality to find assign lipids to leaflets in a bilayer


0.0.0 (2021-02-08)
------------------

* First release on PyPI.
