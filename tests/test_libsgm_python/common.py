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
This module contains common functions present in LibSGM's tests.
"""

import numpy as np

cv_in = np.array(
    [
        [[1, 15, 20], [14, 16, 6], [8, 19, 8]],
        [[13, 11, 3], [22, 12, 9], [16, 4, 12]],
        [[18, 2, 17], [23, 7, 1], [5, 20, 14]],
    ]
)

cv_in_nans = np.array(
    [
        [[1, 15, 20], [14, 16, 6], [8, 19, 8]],
        [[13, 11, 3], [np.nan, np.nan, np.nan], [16, 4, 12]],
        [[18, 2, 17], [23, 7, 1], [5, 20, 14]],
    ]
)
