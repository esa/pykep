# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'pykep'
copyright = '2023, 2024, Dario Izzo and Francesco Biscani'
author = 'Dario Izzo and Francesco Biscani'

import pykep as pk
release = pk.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["myst_nb", "sphinx.ext.intersphinx", "sphinx.ext.autodoc", "sphinx.ext.doctest", "sphinxcontrib.bibtex"]

bibtex_bibfiles = ['refs.bib']


intersphinx_mapping = {
    "hy": ("https://bluescarni.github.io/heyoka.py", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "python": ("https://docs.python.org/3", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None)
}

autoclass_content = 'both'

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Add 'epoch.rst' to your documentation sources
source_suffix = {
    '.rst': 'restructuredtext',
    '.ipynb': 'myst-nb',
}


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_book_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

html_logo = "_static/pykep_logo.png"

html_theme_options = {
    "repository_url": "https://github.com/esa/kep3",
    "repository_branch": "main",
    "path_to_docs": "doc",
    "use_repository_button": True,
    "use_issues_button": True,
    "launch_buttons": {
        "binderhub_url": "https://mybinder.org",
        "notebook_interface": "jupyterlab",
    "home_page_in_toc": True,
    },
}

nb_execution_mode = "force"

nb_execution_excludepatterns = [
    "udp_sf_*",
    "udp_mga*",
    "pontryagin_cartesian*",
    "approximations*",
    "sqp_solver*",
]

latex_engine = "xelatex"

myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
]

# Makes signatures more 'sane'
autodoc_preserve_defaults = True


# Make the linkcheck (and other HTTP requests) look like a real browser
linkcheck_user_agent = (
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
    "AppleWebKit/537.36 (KHTML, like Gecko) "
    "Chrome/122.0.0.0 Safari/537.36"
)

# Limit warnings as per update 6->8
suppress_warnings = [
    "RemovedInSphinx10Warning",
    "mystnb.*",
    "myst.header"
]