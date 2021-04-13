==========
LiPyphilic
==========

.. start-description

**A python toolkit for the analyis of lipid membrane simulations**

.. start-badges

|mdanalysis|
|conda|
|pypi|
|docs|
|travis|
|codecov|
|code-quality|
|supported-versions|
|requires|

.. |mdanalysis| image:: https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA
    :alt: Powered by MDAnalysis
    :target: https://www.mdanalysis.org

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
`NumPy <https://numpy.org/>`__ and `SciPy <https://www.scipy.org/>`__  for efficient computation.
The analysis classes are designed with the same interface as those of MDAnalysis - so if you know how to
`use analysis modules in MDAnalysis
<https://userguide.mdanalysis.org/stable/examples/quickstart.html#Analysis>`__ then learning **lipyphilic**
will be a breeze.
 
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
