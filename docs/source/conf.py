# pylint: skip-file
# Copyright (c) 2024 Centre National d'Etudes Spatiales (CNES).
#
# This file is part of LIBSGM
#
#     https://github.com/CNES/Pandora_libsgm
#
# Licensed under the Apache License, Version 2.0 (the 'License');
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an 'AS IS' BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# Imports
import os
import sys
from importlib.metadata import version as get_version

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
sys.path.insert(0, os.path.abspath("."))


# -- Project information -----------------------------------------------------

project = "libSGM"
copyright = "2024, CNES"
author = "CNES"

# The full version, including alpha/beta/rc tags
try:
    version = get_version("libSGM")
    release = version
except Exception as error:
    print("WARNING: cannot find LibSGM version")
    version = "Unknown"
    release = version

# The master toctree document.
master_doc = "index"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx.ext.autodoc",
    "sphinx_rtd_theme",
    "sphinx.ext.imgmath",
    "breathe",
    "exhale",
    "autoapi.extension",
]

autoapi_dirs = ["../../src"]
autoapi_root = "api_reference"
autoapi_keep_files = True
autoapi_options = [
    "members",
    "undoc-members",
    "private-members",
    "show-inheritance",
    "show-module-summary",
    "special-members",
]

# Setup the breathe extension
breathe_projects = {"My Project": "../doxyoutput/xml"}
breathe_default_project = "My Project"

exhale_args = {
    # These arguments are required
    "containmentFolder": "./api_reference/",
    "rootFileName": "library_root.rst",
    "rootFileTitle": "Library C++ API",
    "doxygenStripFromPath": "..",
    "exhaleExecutesDoxygen": True,
    "exhaleDoxygenStdin": "INPUT =../../src/libSGM/lib/sgm.hpp",
}

# Tell sphinx what the primary language being documented is.
primary_domain = "cpp"

# Tell sphinx what the pygments highlight language should be.
highlight_language = "cpp"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named 'default.css' will overwrite the builtin 'default.css'.
html_static_path = []

latex_elements = {"papersize": "letterpaper", "pointsize": "10pt", "preamble": "", "figure_align": "htbp"}
