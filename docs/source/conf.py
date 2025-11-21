# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'moment_cone'
copyright = '2025, Bulois, Denis & Ressayre'
author = 'Bulois, Denis & Ressayre'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx_autodoc_typehints",
    "sphinx.ext.viewcode",
    "sphinx.ext.apidoc",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "sage": ("https://doc.sagemath.org/html/en/", None),
}

apidoc_modules = [
    {
        "path": "../../src/moment_cone",
        "destination": "moment_cone",
        'separate_modules': False,
        'include_private': True,
    },
]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
html_theme_options = {
    "github_user": "ea-icj",
    "github_repo": "moment_cone",
    "github_button": True,
    "github_banner": True,
}

autodoc_typehints = "description"

import sys
from pathlib import Path

sys.path.insert(0, str(Path('..', '..', 'src').resolve()))
print(sys.path)