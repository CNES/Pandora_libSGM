#!/usr/bin/env python
# coding: utf8
#
# Copyright (c) 2024 Centre National d'Etudes Spatiales (CNES).
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
import numpy
from setuptools import Extension, setup

SCR_DIR = "src/libSGM"

sources = [SCR_DIR + "/lib/sgm.cpp", SCR_DIR + "/sgm_wrapper.pyx"]

ext_1 = Extension(
    "libSGM.sgm_wrapper",
    sources,
    language="c++",
    library_dirs=[],
    libraries=[],
    include_dirs=[numpy.get_include(), SCR_DIR + "/lib/"],
    extra_compile_args=["-O3", "-fopenmp", "-std=c++11"],
    extra_link_args=["-lgomp"],
)

extensions = [ext_1]

setup(
    ext_modules=extensions,
)
