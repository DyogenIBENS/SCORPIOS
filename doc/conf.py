# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../'))


# -- Project information -----------------------------------------------------

project = 'SCORPiOs'
copyright = '2020, Elise Parey & Camille Berthelot'
author = 'Elise Parey & Camille Berthelot'

release = ''
version = ''

master_doc = 'index'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ["sphinx_rtd_theme", "sphinx-prompt", "sphinx.ext.autosectionlabel",
              "sphinx.ext.autodoc", 'sphinx.ext.napoleon']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']
html_static_path = ["_static"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

autodoc_mock_imports = ["matplotlib", 'seaborn', 'roman', "networkx", "svgutils",
                        "ete3", "sklearn", "numpy", "pandas", "scipy",
                        "statsmodels"]


# -- Options for HTML output -------------------------------------------------

import sphinx_rtd_theme


# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'sphinx_pdj_theme'
html_theme = "sphinx_rtd_theme"



# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
def setup(app):
    app.add_css_file('style.css')
