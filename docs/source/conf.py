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
    "sphinx_autodoc_typehints"
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "sage": ("https://doc.sagemath.org/html/en/", None),
}

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

autodoc_typehints = "description"

import sys
from pathlib import Path

sys.path.insert(0, str(Path('..', '..', 'src').resolve()))
print(sys.path)