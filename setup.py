#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import io
import re
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext

from setuptools import find_packages
from setuptools import setup


def read(*names, **kwargs):
    with io.open(
        join(dirname(__file__), *names),
        encoding=kwargs.get('encoding', 'utf8')
    ) as fh:
        return fh.read()


setup(
    name='lipyphilic',
    version='0.6.2',
    license='LGPL-2.1-or-later',
    description='Analyse MD simulations of lipids with python',
    long_description='%s\n%s' % (
        re.compile('^.. start-badges.*^.. end-badges', re.M | re.S).sub('', read('README.rst')),
        re.sub(':[a-z]+:`~?(.*?)`', r'``\1``', read('CHANGELOG.rst'))
    ),
    author='Paul Smith',
    author_email='paul.smith@kcl.ac.uk',
    url='https://github.com/p-j-smith/lipyphilic',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Operating System :: Unix',
        'Operating System :: POSIX',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering'
    ],
    project_urls={
        'Documentation': 'https://lipyphilic.readthedocs.io/',
        'Changelog': 'https://lipyphilic.readthedocs.io/en/latest/changelog.html',
        'Issue Tracker': 'https://github.com/p-j-smith/lipyphilic/issues',
    },
    keywords=[
        'lipids', 'molecular dynamics'
    ],
    python_requires='>=3.6',
    install_requires=[
        'MDAnalysis>=1.0',
        'freud-analysis>=2.4.1',
        'tidynamics',
        'numpy>=1.16',
        'pandas>=1.1',  # 1.1 required by python=3.6
        'seaborn>=0.11',
    ],
    extras_require={
    },
    setup_requires=[
        'pytest-runner',
    ],
    entry_points={
    },
)
