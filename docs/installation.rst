Installation
============

Conda
-----

The easiest way to install **lipyphilic** is through the `conda-forge
<https://anaconda.org/conda-forge>`__ channel of `Conda
<https://docs.conda.io/en/latest/index.html>`__::

    conda config --add channels conda-forge
    conda create -n lipyphilic -c conda-forge python=3.10 lipyphilic
    conda activate lipyphilic

This will install **lipyphilic** along with all of its dependencies into a new virtual environment.

If you do not already have Conda installed on your machine, we recommend
downloading and installing `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`__
--- a lightweight version of Conda.

PyPI
----

It's also possible to install **lipyphilic** from the `Python Package
Index <https://pypi.org/>`__. You can do this using `pip`::

    pip install lipyphilic

Alternatively, you can also install the in-development version with::

    pip install https://github.com/p-j-smith/lipyphilic/archive/main.zip

Dependencies
------------

**lipyphilic** uses `MDAnalysis <https://www.mdanalysis.org/>`__ to carry out all analysis
calculations, and `Freud <https://freud.readthedocs.io/en/stable/>`__ for performing
Voronoi tessellations.

As mentioned above, the simplest way to install these packages,
along with **lipyphilic**, is with `Conda <https://docs.conda.io/en/latest/index.html>`__.
However, it is also possible to install MDAnalysis and Freud using pip, or from source. See
the `MDAnalysis <https://userguide.mdanalysis.org/stable/installation.html>`_ and
`Freud <https://freud.readthedocs.io/en/stable/gettingstarted/installation.html>`_
installation instructions for further information.
