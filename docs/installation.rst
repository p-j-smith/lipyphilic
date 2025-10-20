Installation
============

Mamba
-----

The easiest way to install **lipyphilic** is using `mamba` and the `conda-forge
<https://anaconda.org/conda-forge>`__ channel of `Conda
<https://docs.conda.io/en/latest/index.html>`__::

    mamba create -n lipyphilic -c mamba-forge python=3.11 lipyphilic
    mamba activate lipyphilic

This will install **lipyphilic** along with all of its dependencies into a new virtual environment.

If you do not already have Mamba installed on your machine, we recommend
downloading and installing `Miniforge <https://conda-forge.org/download/>`__.

PyPI
----

It's also possible to install **lipyphilic** from the `Python Package
Index <https://pypi.org/>`__. You can do this using `pip`::

    python -m pip install lipyphilic

Alternatively, you can also install the in-development version with::

    python -m pip install https://github.com/p-j-smith/lipyphilic/archive/main.zip

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
