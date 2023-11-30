import pathlib
import sys

sys.path.insert(0, pathlib.Path("../src").resolve())

import lipyphilic

project = "lipyphilic"
author = "Paul Smith"
year = "2021"
copyright = f"{year}, {author}"
version = f"v{lipyphilic.__version__}"
release = version
source_suffix = ".rst"
master_doc = "index"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.coverage",
    "sphinx.ext.doctest",
    "sphinx.ext.extlinks",
    "sphinx.ext.ifconfig",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.imgconverter",
]

autosummary_generate = True
autodoc_typehints = "signature"
autodoc_docstring_signature = True
autoclass_content = "both"

# pygments_style = 'trac'
pygmants_style = "default"
templates_path = ["."]
extlinks = {
    "issue": ("https://github.com/p-j-smith/lipyphilic/issues/%s", "#"),
    "pr": ("https://github.com/p-j-smith/lipyphilic/pull/%s", "PR #"),
    "mda": ("https://www.mdanalysis.org", "MDAnalysis"),
}

html_theme = "sphinx_rtd_theme"
html_logo = "logo/lipyphilic_logo_grey.png"
html_use_smartypants = True
html_last_updated_fmt = "%b %d, %Y"
html_split_index = False
html_sidebars = {
    "**": ["searchbox.html", "globaltoc.html", "sourcelink.html"],
}
html_short_title = f"{project}-{version}"

napoleon_use_ivar = True
napoleon_use_rtype = False
napoleon_use_param = False
