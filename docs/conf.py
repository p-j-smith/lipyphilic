# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import os

autodoc_mock_imports = [
    "MDAnalysis",
    'MDAnalysis.analysis',
    'freud',
    'numpy',
    'scipy',
    'pandas',
    'matplotlib',
    'seaborn'
]
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx.ext.doctest',
    'sphinx.ext.extlinks',
    'sphinx.ext.ifconfig',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx'
]
source_suffix = '.rst'
master_doc = 'index'
project = 'lipyphilic'
year = '2021'
author = 'Paul Smith'
copyright = '{0}, {1}'.format(year, author)
version = release = '0.10.0'

autodoc_typehints = 'signature'
autodoc_docstring_signature = True
autoclass_content = 'both'

# pygments_style = 'trac'
pygmants_style = 'default'
templates_path = ['.']
extlinks = {
    'issue': ('https://github.com/p-j-smith/lipyphilic/issues/%s', '#'),
    'pr': ('https://github.com/p-j-smith/lipyphilic/pull/%s', 'PR #'),
    'mda': ('https://www.mdanalysis.org', 'MDAnalysis')
}
# on_rtd is whether we are on readthedocs.org
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if not on_rtd:  # only set the theme if we're building docs locally
    html_theme = 'sphinx_rtd_theme'

html_logo = 'logo/lipyphilic_logo_grey.png'
html_use_smartypants = True
html_last_updated_fmt = '%b %d, %Y'
html_split_index = False
html_sidebars = {
   '**': ['searchbox.html', 'globaltoc.html', 'sourcelink.html'],
}
html_short_title = '%s-%s' % (project, version)

napoleon_use_ivar = True
napoleon_use_rtype = False
napoleon_use_param = False
