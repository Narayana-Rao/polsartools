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
sys.path.insert(0, os.path.abspath('../..'))
import sphinx_rtd_theme

# -- Project information -----------------------------------------------------

project = 'PolSARtools'
copyright = '2024, Narayanarao Bhogapurapu'
author = 'Narayanarao Bhogapurapu'

# The full version, including alpha/beta/rc tags
release = '0.4 (beta)'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
	'sphinx_rtd_theme',
	# 'rinoh.frontend.sphinx',
    # 'sphinx.ext.mathjax',
    # 'sphinx.ext.ifconfig',
    # 'sphinx.ext.viewcode',
    # 'sphinx.ext.autodoc',
    # 'sphinx.ext.doctest',
    # 'sphinx.ext.intersphinx',
    # 'sphinx.ext.todo',
    # 'sphinx.ext.coverage',
    # #    'nbsphinx',
    # 'sphinx.ext.napoleon',
    # 'sphinx.ext.autosectionlabel',
    #'sphinxcontrib.bibtex'
]

latex_elements = {
	'papersize':'letterpaper',
	'pointsize':'10pt',
	'preamble':'',
	'figure_align': 'htbp'
}


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# html_theme_options = {
#     'canonical_url': '',
#     'analytics_id': 'UA-XXXXXXX-1',  #  Provided by Google in your dashboard
#     'logo_only': False,
#     'display_version': True,
#     'prev_next_buttons_location': 'bottom',
#     'style_external_links': False,
#     'vcs_pageview_mode': '',
#     'style_nav_header_background': 'white',
#     # Toc options
#     'collapse_navigation': True,
#     'sticky_navigation': True,
#     'navigation_depth': 4,
#     'includehidden': True,
#     'titles_only': False
# }
bibtex_bibfiles = ['../../files/ref.bib']