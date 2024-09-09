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
"""
This module contains tests for the sgm_python code
"""

import unittest
import warnings

import numpy as np

import libsgm_python.sgm_python as sgm
from tests.test_libsgm_python import common


class TestSgmPython(unittest.TestCase):
    """ "
    Test Python version of LibSGM
    """

    ###############################################################
    # Sgm
    ###############################################################
    @staticmethod
    def test_compute_lr():
        """ "
        Test compute lr
        """
        p1 = 3
        p2 = 4

        cv_in_2d_front = np.array(
            [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]]  # same d  # |d' - d| = 1
        )  # |d' - d| > 1

        lr_2d_previous = np.array([[1, 1, 1, 0], [1, 8, 2, 9], [9, 10, 11, 1]])

        p1_in_1d = p1 * np.ones(3)
        p2_in_1d = p2 * np.ones(3)

        # disp
        disp = 1

        # segm
        current_segm = np.ones(3)
        previous_segm = np.ones(3)

        lr_d = sgm.compute_lr(cv_in_2d_front, lr_2d_previous, disp, p1_in_1d, p2_in_1d, current_segm, previous_segm)

        real_lr_d = np.array([3, 9, 14])
        np.testing.assert_array_equal(lr_d, real_lr_d)

    @staticmethod
    def test_compute_lr_reset_history():
        """ "
        Test compute lr
        """
        p1 = 3
        p2 = 4

        cv_in_2d_front = np.array(
            [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]]  # same d  # |d' - d| = 1
        )  # |d' - d| > 1

        lr_2d_previous = np.array([[1, 1, 1, 0], [1, 8, 2, 9], [9, 10, 11, 1]])

        p1_in_1d = p1 * np.ones(3)
        p2_in_1d = p2 * np.ones(3)

        # disp
        disp = 1

        # segm
        current_segm = np.ones(3)
        previous_segm = np.zeros(3)

        lr_d = sgm.compute_lr(cv_in_2d_front, lr_2d_previous, disp, p1_in_1d, p2_in_1d, current_segm, previous_segm)

        real_lr_d = np.array([2, 6, 10])
        np.testing.assert_array_equal(lr_d, real_lr_d)

    def test_sgm_middle_value_invalid(self):
        """ "
        Test SGM middle value invalid
        """
        p1 = 8
        p2 = 32
        cv_in = common.cv_in_nans
        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))

        directions = [[0, 1], [1, 0], [1, 1], [1, -1], [0, -1], [-1, 0], [-1, -1], [-1, 1]]

        segmentation = np.ones((3, 3))

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="All-NaN slice encountered", category=RuntimeWarning)
            cv_out = sgm.run_sgm(cv_in, p1_in, p2_in, directions, segmentation, cost_paths=False, overcounting=False)

            # invalid value : nan
            self.assertTrue(np.isnan(cv_out["cv"][1, 1, 1]))

    def test_sgm_middle_value_invalid_overcounting(self):
        """ "
        Test SGM middle value invalid overcounting
        """
        p1 = 8
        p2 = 32
        cv_in = common.cv_in_nans
        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))
        directions = [[0, 1], [1, 0], [1, 1], [1, -1], [0, -1], [-1, 0], [-1, -1], [-1, 1]]

        segmentation = np.ones((3, 3))

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="All-NaN slice encountered", category=RuntimeWarning)
            cv_out = sgm.run_sgm(cv_in, p1_in, p2_in, directions, segmentation, cost_paths=False, overcounting=True)

            # invalid value : nan
            self.assertTrue(np.isnan(cv_out["cv"][1, 1, 1]))

    def test_sgm_middle_value_lr_0_1(self):
        """ "
        Test SGM middle value directions 0 1
        """
        p1 = 8
        p2 = 32
        cv_in = common.cv_in
        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))
        directions = [[0, 1]]

        segmentation = np.ones((3, 3))

        cv_out = sgm.run_sgm(cv_in, p1_in, p2_in, directions, segmentation, cost_paths=False, overcounting=False)

        # invalid value : nan
        self.assertEqual(cv_out["cv"][1, 1, 1], 20)

    def test_sgm_middle_value_lr_1_0(self):
        """ "
        Test SGM middle value directions 1 0
        """
        p1 = 8
        p2 = 32
        cv_in = common.cv_in
        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))
        directions = [[1, 0]]

        segmentation = np.ones((3, 3))

        cv_out = sgm.run_sgm(cv_in, p1_in, p2_in, directions, segmentation, cost_paths=False, overcounting=False)

        # invalid value : nan
        self.assertEqual(cv_out["cv"][1, 1, 1], 20)

    def test_sgm_middle_value_lr_minus1_0(self):
        """ "
        Test SGM middle value directions -1 1
        """
        p1 = 8
        p2 = 32
        cv_in = common.cv_in
        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))
        directions = [[-1, 0]]

        segmentation = np.ones((3, 3))

        cv_out = sgm.run_sgm(cv_in, p1_in, p2_in, directions, segmentation, cost_paths=False, overcounting=False)

        # invalid value : nan
        self.assertEqual(cv_out["cv"][1, 1, 1], 18)

    def test_sgm_middle_value_lr_0_minus1(self):
        """ "
        Test SGM middle value directions 0 -1
        """
        p1 = 8
        p2 = 32
        cv_in = common.cv_in
        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))
        directions = [[0, -1]]

        segmentation = np.ones((3, 3))

        cv_out = sgm.run_sgm(cv_in, p1_in, p2_in, directions, segmentation, cost_paths=False, overcounting=False)

        # invalid value : nan
        self.assertEqual(cv_out["cv"][1, 1, 1], 12)

    def test_sgm_middle_value_lr_1_1(self):
        """ "
        Test SGM middle value directions 1 1
        """
        p1 = 8
        p2 = 32
        cv_in = common.cv_in
        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))
        directions = [[1, 1]]

        segmentation = np.ones((3, 3))

        cv_out = sgm.run_sgm(cv_in, p1_in, p2_in, directions, segmentation, cost_paths=False, overcounting=False)

        # invalid value : nan
        self.assertEqual(cv_out["cv"][1, 1, 1], 20)

    def test_sgm_middle_value_lr_minus1_1(self):
        """ "
        Test SGM middle value directions -1 1
        """
        p1 = 8
        p2 = 32
        cv_in = common.cv_in
        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))
        directions = [[-1, 1]]

        segmentation = np.ones((3, 3))

        cv_out = sgm.run_sgm(cv_in, p1_in, p2_in, directions, segmentation, cost_paths=False, overcounting=False)

        # invalid value : nan
        self.assertEqual(cv_out["cv"][1, 1, 1], 12)

    def test_sgm_middle_value_lr_1_minus1(self):
        """ "
        Test SGM middle value directions 1 -1
        """
        p1 = 8
        p2 = 32
        cv_in = common.cv_in
        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))
        directions = [[1, -1]]

        segmentation = np.ones((3, 3))

        cv_out = sgm.run_sgm(cv_in, p1_in, p2_in, directions, segmentation, cost_paths=False, overcounting=False)

        # invalid value : nan
        self.assertEqual(cv_out["cv"][1, 1, 1], 20)

    def test_sgm_middle_value_lr_minus1_minus1(self):
        """ "
        Test SGM middle value directions -1 -1
        """
        p1 = 8
        p2 = 32
        cv_in = common.cv_in
        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))
        directions = [[-1, -1]]

        segmentation = np.ones((3, 3))

        cv_out = sgm.run_sgm(cv_in, p1_in, p2_in, directions, segmentation, cost_paths=False, overcounting=False)

        # invalid value : nan
        self.assertEqual(cv_out["cv"][1, 1, 1], 20)

    def test_sgm_middle_value(self):
        """ "
        Test SGM middle value all directions
        """
        p1 = 8
        p2 = 32
        cv_in = common.cv_in
        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))
        directions = [[0, 1], [1, 0], [1, 1], [1, -1], [0, -1], [-1, 0], [-1, -1], [-1, 1]]

        segmentation = np.ones((3, 3))

        cv_out = sgm.run_sgm(cv_in, p1_in, p2_in, directions, segmentation, cost_paths=False, overcounting=False)

        # invalid value : nan
        self.assertEqual(cv_out["cv"][1, 1, 1], 142)

    def test_sgm_middle_value_reset_middle(self):
        """ "
        Test SGM middle value all directions with a reset in middle
        """
        p1 = 8
        p2 = 32
        cv_in = common.cv_in
        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))
        directions = [[0, 1], [1, 0], [1, 1], [1, -1], [0, -1], [-1, 0], [-1, -1], [-1, 1]]

        segmentation = np.ones((3, 3))
        segmentation[1, 1] = 2

        cv_out = sgm.run_sgm(cv_in, p1_in, p2_in, directions, segmentation, cost_paths=False, overcounting=False)

        # invalid value : nan
        self.assertEqual(cv_out["cv"][1, 1, 1], 96)

    def test_sgm_middle_value_overcounting(self):
        """ "
        Test SGM middle value all directions overcounting
        """
        p1 = 8
        p2 = 32
        cv_in = common.cv_in
        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))
        directions = [[0, 1], [1, 0], [1, 1], [1, -1], [0, -1], [-1, 0], [-1, -1], [-1, 1]]

        segmentation = np.ones((3, 3))

        cv_out = sgm.run_sgm(cv_in, p1_in, p2_in, directions, segmentation, cost_paths=False, overcounting=True)

        # invalid value : nan
        self.assertEqual(cv_out["cv"][1, 1, 1], 58)

    def test_sgm_middle_value_min_cost(self):
        """ "
        Test SGM middle value minimum cost
        """
        p1 = 8
        p2 = 32
        cv_in = common.cv_in
        p1_in = p1 * np.ones((3, 3, 8))
        p2_in = p2 * np.ones((3, 3, 8))
        directions = [[0, 1], [1, 0], [1, 1], [1, -1], [0, -1], [-1, 0], [-1, -1], [-1, 1]]

        segmentation = np.ones((3, 3))

        cv_out = sgm.run_sgm(cv_in, p1_in, p2_in, directions, segmentation, cost_paths=True, overcounting=False)

        # invalid value : nan
        self.assertEqual(cv_out["cv"][1, 1, 1], 142)
        self.assertEqual(cv_out["cv_min"][1, 1, 1], 2)
