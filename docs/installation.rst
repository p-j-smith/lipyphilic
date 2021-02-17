Installation
============

You can install lipyphilic using pip::

    pip install lipyphilic

You can also install the in-development version with::

    pip install https://github.com/p-j-smith/lipyphilic/archive/master.zip

Dependencies
============

lipyphilic uses `MDAnalysis <https://www.mdanalysis.org/>`_ to carry out all analysis
calculations, and `Freud <https://freud.readthedocs.io/en/stable/>`_ for performing
Voronoi tessellations. The simplest way to install these packages,
along with lipyphilic, is through `Conda <https://docs.conda.io/en/latest/index.html>`_::

    conda create -y -n myenv -c conda-forge python=3.8 MDAnalysis>1.0 freud=2.4.1 pip
    conda activate myenv
    pip install lipyphilic

It is also possible to install MDAnalysis and Freud using pip, or from source. See
the `MDAnalysis <https://userguide.mdanalysis.org/stable/installation.html>`_ and
`Freud <https://freud.readthedocs.io/en/stable/gettingstarted/installation.html>`_
installation instructions for further information.
