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

import unittest
import libsgm_python.sgm_python_parall as sgm
import numpy as np


class TestSgmPythonParall(unittest.TestCase):

    """"
    Test Python version of LibSGM
    """

    ###############################################################
    # Sgm
    ###############################################################

    def test_sgm_middle_value_invalid(self):
        p1 = 8
        p2 = 32

        cv_in = np.array([[[1,15,20], [14,16,6], [8,19,8]],
                          [[13,11,3], [np.nan,np.nan,np.nan], [16,4,12]],
                          [[18,2,17], [23,7,1], [5,20,14]]])

        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))

        directions = [[0, 1], [1, 0], [1, 1], [1, -1], [0, -1], [-1, 0], [-1, -1], [-1, 1]]

        cv_out = sgm.run_sgm_parall(cv_in, p1_in, p2_in, directions, cost_paths=False, overcounting=False)

        # invalid value : nan
        self.assertTrue(np.isnan(cv_out["cv"][1,1,1]))

    def test_sgm_middle_value_invalid_overcounting(self):
        p1 = 8
        p2 = 32

        cv_in = np.array([[[1, 15, 20], [14, 16, 6], [8, 19, 8]],
                          [[13, 11, 3], [np.nan, np.nan, np.nan], [16, 4, 12]],
                          [[18, 2, 17], [23, 7, 1], [5, 20, 14]]])

        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))

        directions = [[0, 1], [1, 0], [1, 1], [1, -1], [0, -1], [-1, 0], [-1, -1], [-1, 1]]

        cv_out = sgm.run_sgm_parall(cv_in, p1_in, p2_in, directions, cost_paths=False, overcounting=True)

        # invalid value : nan
        self.assertTrue(np.isnan(cv_out["cv"][1, 1, 1]))

    def test_sgm_middle_value_lr_0_1(self):
        p1 = 8
        p2 = 32
        cv_in = np.array([[[1, 15, 20], [14, 16, 6], [8, 19, 8]],
                          [[13, 11, 3], [22, 12, 9], [16, 4, 12]],
                          [[18, 2, 17], [23, 7, 1], [5, 20, 14]]])

        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))

        directions = [[0, 1]]

        cv_out = sgm.run_sgm_parall(cv_in, p1_in, p2_in, directions, cost_paths=False, overcounting=False)

        # invalid value : nan
        self.assertEqual(cv_out["cv"][1, 1, 1], 20)

    def test_sgm_middle_value_lr_1_0(self):
        p1 = 8
        p2 = 32
        cv_in = np.array([[[1, 15, 20], [14, 16, 6], [8, 19, 8]],
                          [[13, 11, 3], [22, 12, 9], [16, 4, 12]],
                          [[18, 2, 17], [23, 7, 1], [5, 20, 14]]])

        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))

        directions = [[1, 0]]

        cv_out = sgm.run_sgm_parall(cv_in, p1_in, p2_in, directions, cost_paths=False, overcounting=False)

        # invalid value : nan
        self.assertEqual(cv_out["cv"][1, 1, 1], 20)

    def test_sgm_middle_value_lr_minus1_0(self):
        p1 = 8
        p2 = 32
        cv_in = np.array([[[1, 15, 20], [14, 16, 6], [8, 19, 8]],
                          [[13, 11, 3], [22, 12, 9], [16, 4, 12]],
                          [[18, 2, 17], [23, 7, 1], [5, 20, 14]]])

        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))

        directions = [[-1, 0]]

        cv_out = sgm.run_sgm_parall(cv_in, p1_in, p2_in, directions, cost_paths=False, overcounting=False)

        # invalid value : nan
        self.assertEqual(cv_out["cv"][1, 1, 1], 18)

    def test_sgm_middle_value_lr_0_minus1(self):
        p1 = 8
        p2 = 32
        cv_in = np.array([[[1, 15, 20], [14, 16, 6], [8, 19, 8]],
                          [[13, 11, 3], [22, 12, 9], [16, 4, 12]],
                          [[18, 2, 17], [23, 7, 1], [5, 20, 14]]])

        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))

        directions = [[0, -1]]

        cv_out = sgm.run_sgm_parall(cv_in, p1_in, p2_in, directions, cost_paths=False, overcounting=False)

        # invalid value : nan
        self.assertEqual(cv_out["cv"][1, 1, 1], 12)

    def test_sgm_middle_value_lr_1_1(self):
        p1 = 8
        p2 = 32
        cv_in = np.array([[[1, 15, 20], [14, 16, 6], [8, 19, 8]],
                          [[13, 11, 3], [22, 12, 9], [16, 4, 12]],
                          [[18, 2, 17], [23, 7, 1], [5, 20, 14]]])

        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))

        directions = [[1, 1]]

        cv_out = sgm.run_sgm_parall(cv_in, p1_in, p2_in, directions, cost_paths=False, overcounting=False)

        # invalid value : nan
        self.assertEqual(cv_out["cv"][1, 1, 1], 20)

    def test_sgm_middle_value_lr_minus1_1(self):
        p1 = 8
        p2 = 32
        cv_in = np.array([[[1, 15, 20], [14, 16, 6], [8, 19, 8]],
                          [[13, 11, 3], [22, 12, 9], [16, 4, 12]],
                          [[18, 2, 17], [23, 7, 1], [5, 20, 14]]])

        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))

        directions = [[-1, 1]]

        cv_out = sgm.run_sgm_parall(cv_in, p1_in, p2_in, directions, cost_paths=False, overcounting=False)

        # invalid value : nan
        self.assertEqual(cv_out["cv"][1, 1, 1], 12)

    def test_sgm_middle_value_lr_1_minus1(self):
        p1 = 8
        p2 = 32
        cv_in = np.array([[[1, 15, 20], [14, 16, 6], [8, 19, 8]],
                          [[13, 11, 3], [22, 12, 9], [16, 4, 12]],
                          [[18, 2, 17], [23, 7, 1], [5, 20, 14]]])

        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))

        directions = [[1, -1]]

        cv_out = sgm.run_sgm_parall(cv_in, p1_in, p2_in, directions, cost_paths=False, overcounting=False)

        # invalid value : nan
        self.assertEqual(cv_out["cv"][1, 1, 1], 20)

    def test_sgm_middle_value_lr_minus1_minus1(self):
        p1 = 8
        p2 = 32
        cv_in = np.array([[[1, 15, 20], [14, 16, 6], [8, 19, 8]],
                          [[13, 11, 3], [22, 12, 9], [16, 4, 12]],
                          [[18, 2, 17], [23, 7, 1], [5, 20, 14]]])

        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))

        directions = [[-1, -1]]

        cv_out = sgm.run_sgm_parall(cv_in, p1_in, p2_in, directions, cost_paths=False, overcounting=False)

        # invalid value : nan
        self.assertEqual(cv_out["cv"][1, 1, 1], 20)

    def test_sgm_middle_value(self):
        p1 = 8
        p2 = 32
        cv_in = np.array([[[1, 15, 20], [14, 16, 6], [8, 19, 8]],
                          [[13, 11, 3], [22, 12, 9], [16, 4, 12]],
                          [[18, 2, 17], [23, 7, 1], [5, 20, 14]]])

        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))

        directions = [[0, 1], [1, 0], [1, 1], [1, -1], [0, -1], [-1, 0], [-1, -1], [-1, 1]]

        cv_out = sgm.run_sgm_parall(cv_in, p1_in, p2_in, directions, cost_paths=False, overcounting=False)

        # invalid value : nan
        self.assertEqual(cv_out["cv"][1, 1, 1], 142)

    def test_sgm_middle_value_overcounting(self):
        p1 = 8
        p2 = 32

        cv_in = np.array([[[1, 15, 20], [14, 16, 6], [8, 19, 8]],
                          [[13, 11, 3], [22, 12, 9], [16, 4, 12]],
                          [[18, 2, 17], [23, 7, 1], [5, 20, 14]]])

        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))

        directions = [[0, 1], [1, 0], [1, 1], [1, -1], [0, -1], [-1, 0], [-1, -1], [-1, 1]]

        cv_out = sgm.run_sgm_parall(cv_in, p1_in, p2_in, directions, cost_paths=False, overcounting=True)

        # invalid value : nan
        self.assertEqual(cv_out["cv"][1, 1, 1], 58)

    def test_sgm_middle_value_min_cost(self):
        p1 = 8
        p2 = 32

        cv_in = np.array([[[1, 15, 20], [14, 16, 6], [8, 19, 8]],
                          [[13, 11, 3], [22, 12, 9], [16, 4, 12]],
                          [[18, 2, 17], [23, 7, 1], [5, 20, 14]]])

        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))

        directions = [[0, 1], [1, 0], [1, 1], [1, -1], [0, -1], [-1, 0], [-1, -1], [-1, 1]]

        cv_out = sgm.run_sgm_parall(cv_in, p1_in, p2_in, directions, cost_paths=True, overcounting=False)

        # invalid value : nan
        self.assertEqual(cv_out["cv"][1, 1, 1], 142)
        self.assertEqual(cv_out["cv_min"][1, 1, 1], 2)
