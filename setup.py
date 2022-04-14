# type:ignore
# pylint: disable=ungrouped-imports
#!/usr/bin/env python
# coding: utf8
#
# Copyright (c) 2020 Centre National d'Etudes Spatiales (CNES).
#
# This file is part of LIBSGM
#
#     https://github.com/CNES/Pandora_libsgm
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
"""
This module contains the required libraries and softwares allowing to execute the software,
and setup elements to configure and identify the software.
"""

import os
import shutil
from codecs import open as copen
import numpy
from setuptools import setup

try:
    from Cython.Distutils.extension import Extension
except ImportError:
    from setuptools import Extension
    from setuptools.command.build_ext import build_ext

    USING_CYTHON = False
else:
    from Cython.Distutils import build_ext

    USING_CYTHON = True

CMDCLASS = {"build_ext": build_ext}

SCR_DIR = "sources"

if USING_CYTHON:
    sources = [SCR_DIR + "/lib/sgm.cpp", SCR_DIR + "/sgm_wrapper.pyx"]
else:
    sources = [SCR_DIR + "/lib/sgm.cpp", SCR_DIR + "/sgm_wrapper.cpp"]

ext_1 = Extension(
    "libSGM.sgm_wrapper",
    sources,
    language="c++",
    library_dirs=[],
    libraries=[],
    include_dirs=[numpy.get_include(), SCR_DIR + "/lib/sgm.hpp"],
    extra_compile_args=["-O3", "-fopenmp", "-std=c++11"],
    extra_link_args=["-lgomp"],
)

extensions = [ext_1]

try:
    from sphinx.setup_command import BuildDoc

    CMDCLASS.update({"build_sphinx": BuildDoc})
except ImportError:
    print("WARNING: sphinx not available. Doc cannot be built")

os.environ["CC"] = shutil.which("gcc")
os.environ["CXX"] = shutil.which("g++")


def readme():
    with copen("README.md", "r", "utf-8") as fstream:
        return fstream.read()


setup(
    use_scm_version=True,
    long_description=readme(),
    command_options={
        "build_sphinx": {
            "build_dir": ("setup.py", "doc/build/"),
            "source_dir": ("setup.py", "doc/source/"),
            "warning_is_error": ("setup.py", True),
        }
    },
    zip_safe=False,
    ext_modules=extensions,
    cmdclass=CMDCLASS,
)
