Installation
============

Install lipyphilic
------------------

If you already have the necessary `dependencies
<https://raw.githubusercontent.com/p-j-smith/lipyphilic/master/requirements.yml>`__ installed,
you can install **lipyphilic** using pip::

    pip install lipyphilic

You can also install the in-development version with::

    pip install https://github.com/p-j-smith/lipyphilic/archive/master.zip

Dependencies
------------

lipyphilic uses `MDAnalysis <https://www.mdanalysis.org/>`__ to carry out all analysis
calculations, and `Freud <https://freud.readthedocs.io/en/stable/>`__ for performing
Voronoi tessellations. The simplest way to install these packages,
along with lipyphilic, is through creating a new isolated environment with `Conda
<https://docs.conda.io/en/latest/index.html>`__::

    curl https://raw.githubusercontent.com/p-j-smith/lipyphilic/master/requirements.yml -o lipyphilic.yml
    conda env create -f lipyphilic.yml

You can then activate your environment::

    conda activate lipyphilic

It is also possible to install MDAnalysis and Freud using pip, or from source. See
the `MDAnalysis <https://userguide.mdanalysis.org/stable/installation.html>`_ and
`Freud <https://freud.readthedocs.io/en/stable/gettingstarted/installation.html>`_
installation instructions for further information.
