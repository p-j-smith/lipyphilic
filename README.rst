==========
LiPyphilic
==========

.. start-description

**A python toolkit for the analyis of lipid membrane simulations**

.. start-badges

|conda|
|pypi|
|docs|
|travis|
|codecov|
|code-quality|
|supported-versions|
|requires|

.. |conda| image:: https://img.shields.io/conda/vn/conda-forge/lipyphilic.svg
    :alt: Conda-fogre latest release
    :target: https://anaconda.org/conda-forge/lipyphilic

.. |pypi| image:: https://img.shields.io/pypi/v/lipyphilic.svg
    :alt: PyPI Package latest release
    :target: https://pypi.org/project/lipyphilic

.. |docs| image:: https://readthedocs.org/projects/lipyphilic/badge/?style=flat
    :target: https://readthedocs.org/projects/lipyphilic
    :alt: Documentation Status

.. |travis| image:: https://api.travis-ci.com/p-j-smith/lipyphilic.svg?branch=master
    :alt: Travis-CI Build Status
    :target: https://travis-ci.com/github/p-j-smith/lipyphilic

.. |codecov| image:: https://codecov.io/gh/p-j-smith/lipyphilic/branch/master/graphs/badge.svg?branch=master
    :alt: Coverage Status
    :target: https://codecov.io/github/p-j-smith/lipyphilic

.. |code-quality| image:: https://img.shields.io/lgtm/grade/python/g/p-j-smith/lipyphilic.svg
    :alt: LGTM code quality
    :target: https://lgtm.com/projects/g/p-j-smith/lipyphilic/context:python

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/lipyphilic.svg
    :alt: Supported versions
    :target: https://pypi.org/project/lipyphilic

.. |requires| image:: https://requires.io/github/p-j-smith/lipyphilic/requirements.svg?branch=master
    :alt: Requirements Status
    :target: https://requires.io/github/p-j-smith/lipyphilic/requirements/?branch=master

.. end-badges

**lipyphilic** is free software licensed under the GNU General Public License v2 or later (GPLv2+)

Overview
========

**lipyphilic** is a set of tools for analysing MD simulations of lipid bilayers. It is an object-oriented
Python package built directly on top of `MDAnalysis <https://www.mdanalysis.org/>`__, and makes use of
`NumPy <https://numpy.org/>`__, `SciPy <https://www.scipy.org/>`__ and `pandas <https://pandas.pydata.org/>`__ for
efficient computation. The analysis classes are designed with the same interface as those of MDAnalysis -
so if you know how to `use analysis modules in
MDAnalysis <https://userguide.mdanalysis.org/stable/examples/quickstart.html#Analysis>`__ then you know how
to use **lipyphilic**!
 
Analysis tools in **lipyphilic** include: identifying sterol flip-flop events, calculating domain registration over time,
and calculating local lipid compositions. These tools position **lipyphilic** as complementary to, rather than
competing against, existing membrane analysis software such as `MemSurfer <https://github.com/LLNL/MemSurfer>`__ and
`FatSlim <http://fatslim.github.io/>`__.

Check out the `Basic Usage <https://lipyphilic.readthedocs.io/en/stable/usage.html>`__ example to see how to use
**lipyphilic**, and see the `Analysis tools <https://lipyphilic.readthedocs.io/en/stable/reference/analyses.html>`__ 
section for detailed information and exmaples on each tool.

Installation
============

The easiest way to install **lipyphilic** along with its dependencies is through `Conda
<https://docs.conda.io/en/latest/index.html>`__::

    conda config --add channels conda-forge
    conda install lipyphilic

See the `installation guide <https://lipyphilic.readthedocs.io/en/stable/installation.html>`__ for futher information.

Citing
======

If you use **lipyphilic** in your project, please cite `MDAnalysis <https://www.mdanalysis.org/pages/citations/>`__ and
if you use the Area Per Lipid tool please also cite `Freud <https://freud.readthedocs.io/en/stable/reference/citing.html>`__.

There is currently no paper describing **lipyphilic**, but we're working on it. In the meantime, if you like what we
do, please tell everyone you know to check out **lipyphilic**! And if there are things you think we could improve, features
you would like to see added, or pesky bugs that need to be fixed, please raise an issue on
`github <https://github.com/p-j-smith/lipyphilic/issues>`__.

.. end-description

Full documentation
==================

Head to `lipyphilic.readthedocs.io <https://lipyphilic.readthedocs.io/en/stable/>`__, where you will find the full
documentation of **lipyphilic**'s API as well as examples of how to use the analysis tools.

Acknowlegment
=============

The respository structure of **lipyphilic** is based on the
`PyLibrary Cookeicutter template <https://github.com/ionelmc/cookiecutter-pylibrary>`__.
