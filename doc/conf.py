import os

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Interactpolate'
copyright = '2025, Lennart Klebl'
author = 'Lennart Klebl'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["hawkmoth", "sphinx.ext.autodoc", "sphinx.ext.extlinks"]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
html_theme_options = {
    'body_max_width' : '50em',
    'page_width': 'auto',
    'sidebar_width': '18em',
}

# C doc config
hawkmoth_root = os.path.abspath('../')
hawkmoth_source_uri = 'https://github.com/alcubierre-drive/interactpolate/blob/main/{source}#L{line}'

