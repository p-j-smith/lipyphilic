 .. _tutorials:

Interactive tutorials
=====================

To help you get the most out of **lipyphilic**, we have created a set of interactive tutorials
in the form of Jupyter Notebooks. There is no need to download or install anything, simply click
the link below:

.. image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/p-j-smith/lipyphilic-tutorials/main?filepath=notebooks%2F1-Introduction.ipynb

We currently have tutorials on the following topics:

1. | **Basic usage**: Illustrates basic usage of **lipyphilic**, including how to store results for later usage.
   | Also shows how to assign lipids to leaflets, which is required for many other analyses.

2. | **Flip-flop rate**: Shows how to use **lipyphilic** to calculate the rate of cholesterol flip-flop, as well as
   | identify the frames at which each flip-flop event begins and ends.

3. | **Local lipid environments**: Illustrates how to determine the local lipid environment of each lipid
   | over time, as well as the enrichment/depletion index.

4. | **Lipid domains**: Shows how to calculate the largest cluster of specific lipids over time. Examples
   | include finding the largest ganglioside cluster in a neuronal plasma membrane and identifying
   | the largest domain of |Lo| lipids in a phase separated membrane.

5. | **Interleaflet registration**: This notebook shows how to calculate the interleaflet registration over
   | time. The example shows how to calculate the registration of |Lo| lipids across leaflets.

6. | **Lateral diffusion**: Illustrates how to perform "nojump" trajectory unwrapping with LiPyphilic,
   | then use the unwrapped coordinates to calculate the mean-squared displacement and lateral
   | diffusion coefficient of lipids in a membrane.

7. | **Coarse-grained lipid order parameter**: Shows how to calculate the coarse-grained order
   | parameter, and how to create a two-dimensional projection of these values onto the membrane
   | plane.

8. | **Projection plots**: Shows how to create two-dimensional projections of arbitrary lipid properties
   | onto the membrane plane. Examples include projecting local membrane thicknesses calculated
   | using `FATSLiM <http://fatslim.github.io/>`__ onto the membrane plane, and projecting the ordered state (|Lo| and |Ld|) of lipids
   | onto the membrane plane.

9. | **Potential of mean force (PMF)**: This notebook illustrates how to use **lipyphilic** to calculate the
   | height and orientation of sterols in a membrane, and subsequently plot the two-dimensional
   | PMF of sterol height and orientation.

10. | **Hidden Markov Models (HMM)**: Learn how to use the output of **lipyphilic** to construct HMMs
    | with `HMMLearn <https://hmmlearn.readthedocs.io/en/latest/>`__. We will create a HMM based on lipid thicknesses to detect |Lo| and |Ld| lipids in a
    | phase separated membrane. The output from this is can be used as input to other
    | analyses in **lipyphilic**, such as calculating interleaflet registration or local lipid environments.

.. |Lo| replace:: L\ :sub:`o`
.. |Ld| replace:: L\ :sub:`d`