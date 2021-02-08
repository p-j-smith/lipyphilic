========
Overview
========

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - tests
      - | |travis| |appveyor| |requires|
        | |codecov|
    * - package
      - | |version| |wheel| |supported-versions| |supported-implementations|
        | |commits-since|
.. |docs| image:: https://readthedocs.org/projects/lipyphilic/badge/?style=flat
    :target: https://readthedocs.org/projects/lipyphilic
    :alt: Documentation Status

.. |travis| image:: https://api.travis-ci.com/p-j-smith/lipyphilic.svg?branch=master
    :alt: Travis-CI Build Status
    :target: https://travis-ci.com/github/p-j-smith/lipyphilic

.. |appveyor| image:: https://ci.appveyor.com/api/projects/status/github/p-j-smith/lipyphilic?branch=master&svg=true
    :alt: AppVeyor Build Status
    :target: https://ci.appveyor.com/project/p-j-smith/lipyphilic

.. |requires| image:: https://requires.io/github/p-j-smith/lipyphilic/requirements.svg?branch=master
    :alt: Requirements Status
    :target: https://requires.io/github/p-j-smith/lipyphilic/requirements/?branch=master

.. |codecov| image:: https://codecov.io/gh/p-j-smith/lipyphilic/branch/master/graphs/badge.svg?branch=master
    :alt: Coverage Status
    :target: https://codecov.io/github/p-j-smith/lipyphilic

.. |version| image:: https://img.shields.io/pypi/v/lipyphilic.svg
    :alt: PyPI Package latest release
    :target: https://pypi.org/project/lipyphilic

.. |wheel| image:: https://img.shields.io/pypi/wheel/lipyphilic.svg
    :alt: PyPI Wheel
    :target: https://pypi.org/project/lipyphilic

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/lipyphilic.svg
    :alt: Supported versions
    :target: https://pypi.org/project/lipyphilic

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/lipyphilic.svg
    :alt: Supported implementations
    :target: https://pypi.org/project/lipyphilic

.. |commits-since| image:: https://img.shields.io/github/commits-since/p-j-smith/lipyphilic/v0.0.0.svg
    :alt: Commits since latest release
    :target: https://github.com/p-j-smith/lipyphilic/compare/v0.0.0...master



.. end-badges

Analyse MD simulations of lipids with python

* Free software: GNU Lesser General Public License v2.1 or later (LGPLv2+)

Installation
============

::

    pip install lipyphilic

You can also install the in-development version with::

    pip install https://github.com/p-j-smith/lipyphilic/archive/master.zip


Documentation
=============


https://lipyphilic.readthedocs.io/


Development
===========

To run all the tests run::

    tox

Note, to combine the coverage data from all the tox environments run:

.. list-table::
    :widths: 10 90
    :stub-columns: 1

    - - Windows
      - ::

            set PYTEST_ADDOPTS=--cov-append
            tox

    - - Other
      - ::

            PYTEST_ADDOPTS=--cov-append tox
