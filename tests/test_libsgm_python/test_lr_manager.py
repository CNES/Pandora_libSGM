# pylint:disable = not-callable
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
This module contains tests for the lr_manager code
"""
import unittest

import numpy as np

import libsgm_python.lr_manager as lrm


class TestSgmPythonLrManager(unittest.TestCase):
    """ "
    Test LrManager used in Python version of LibSGM
    """

    ###############################################################
    # LrManager
    ###############################################################

    def test_dir_1_0_lr_manager(self):
        """ "
        Test direction (1, 0)
        """
        shape1 = [4, 3, 3]
        direction1 = (1, 0)
        lr_manager_1 = lrm.LrManager(shape1, direction1)

        seen_pixels = np.zeros((shape1[0], shape1[1]))

        # test init
        self.assertEqual(len(lr_manager_1.planes_front), 1)
        self.assertEqual(len(lr_manager_1.planes_previous), 0)
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["i"], [0])
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["j"], [0, 1, 2])

        current_lr1 = [np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])]
        lr_manager_1.set_current_lr(current_lr1)

        curent_segm1 = [np.array([1, 2, 3])]
        lr_manager_1.set_current_segm(curent_segm1)

        # add to seen
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1

        # test next() MOVE Plan
        lr_manager_1.next()
        # add to seen
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1

        self.assertEqual(len(lr_manager_1.planes_front), 1)
        self.assertEqual(len(lr_manager_1.planes_previous), 1)
        # current
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["i"], [1])
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["j"], [0, 1, 2])
        # previous
        np.testing.assert_array_equal(lr_manager_1.planes_previous[0]["i"], [0])
        np.testing.assert_array_equal(lr_manager_1.planes_previous[0]["j"], [0, 1, 2])

        # Test lr saving + elimination
        previous_lr1 = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        np.testing.assert_array_equal(lr_manager_1.get_previous_lr(0), previous_lr1)

        # Test segm saving + elimination
        previous_segm1 = np.array([1, 2, 3])
        np.testing.assert_array_equal(lr_manager_1.get_previous_segm(0), previous_segm1)

        # test next() MOVE Plan out of dimensions AND add to seen
        lr_manager_1.next()
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1
        lr_manager_1.next()
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1
        lr_manager_1.next()

        # out of plan
        self.assertEqual(len(lr_manager_1.planes_front), 0)
        self.assertEqual(len(lr_manager_1.planes_previous), 0)

        # Test if all pixels were seen only once
        np.testing.assert_array_equal(seen_pixels, np.ones((shape1[0], shape1[1])))

    def test_dir_minus1_0_lr_manager(self):
        """ "
        Test direction (-1, 0)
        """
        shape1 = [4, 3, 3]
        direction1 = (-1, 0)
        lr_manager_1 = lrm.LrManager(shape1, direction1)

        seen_pixels = np.zeros((shape1[0], shape1[1]))

        # test init
        self.assertEqual(len(lr_manager_1.planes_front), 1)
        self.assertEqual(len(lr_manager_1.planes_previous), 0)
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["i"], [3])
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["j"], [0, 1, 2])

        current_lr1 = [np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])]
        lr_manager_1.set_current_lr(current_lr1)

        # add to seen
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1

        # test next() MOVE Plan
        lr_manager_1.next()
        # add to seen
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1

        self.assertEqual(len(lr_manager_1.planes_front), 1)
        self.assertEqual(len(lr_manager_1.planes_previous), 1)
        # current
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["i"], [2])
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["j"], [0, 1, 2])
        # previous
        np.testing.assert_array_equal(lr_manager_1.planes_previous[0]["i"], [3])
        np.testing.assert_array_equal(lr_manager_1.planes_previous[0]["j"], [0, 1, 2])

        # Test lr saving + elimination
        previous_lr1 = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        np.testing.assert_array_equal(lr_manager_1.get_previous_lr(0), previous_lr1)

        # test next() MOVE Plan out of dimensions AND add to seen
        lr_manager_1.next()
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1
        lr_manager_1.next()
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1
        lr_manager_1.next()

        # out of plan
        self.assertEqual(len(lr_manager_1.planes_front), 0)
        self.assertEqual(len(lr_manager_1.planes_previous), 0)

        # Test if all pixels were seen only once
        np.testing.assert_array_equal(seen_pixels, np.ones((shape1[0], shape1[1])))

    def test_dir_0_1_lr_manager(self):
        """ "
        Test direction (0, 1)
        """
        shape1 = [4, 3, 3]
        direction1 = (0, 1)
        lr_manager_1 = lrm.LrManager(shape1, direction1)

        seen_pixels = np.zeros((shape1[0], shape1[1]))

        # test init
        self.assertEqual(len(lr_manager_1.planes_front), 1)
        self.assertEqual(len(lr_manager_1.planes_previous), 0)
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["i"], [0, 1, 2, 3])
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["j"], [0])

        current_lr1 = [np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]])]
        lr_manager_1.set_current_lr(current_lr1)

        # add to seen
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1

        # test next() MOVE Plan
        lr_manager_1.next()
        # add to seen
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1

        self.assertEqual(len(lr_manager_1.planes_front), 1)
        self.assertEqual(len(lr_manager_1.planes_previous), 1)
        # current
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["i"], [0, 1, 2, 3])
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["j"], [1])
        # previous
        np.testing.assert_array_equal(lr_manager_1.planes_previous[0]["i"], [0, 1, 2, 3])
        np.testing.assert_array_equal(lr_manager_1.planes_previous[0]["j"], [0])

        # Test lr saving + elimination
        previous_lr1 = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]])
        np.testing.assert_array_equal(lr_manager_1.get_previous_lr(0), previous_lr1)

        # test next() MOVE Plan out of dimensions AND add to seen
        lr_manager_1.next()
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1
        lr_manager_1.next()

        # out of plan
        self.assertEqual(len(lr_manager_1.planes_front), 0)
        self.assertEqual(len(lr_manager_1.planes_previous), 0)

        # Test if all pixels were seen only once
        np.testing.assert_array_equal(seen_pixels, np.ones((shape1[0], shape1[1])))

    def test_dir_0_minus1_lr_manager(self):
        """ "
        Test direction (0, -1)
        """
        shape1 = [4, 3, 3]
        direction1 = (0, -1)
        lr_manager_1 = lrm.LrManager(shape1, direction1)

        seen_pixels = np.zeros((shape1[0], shape1[1]))

        # test init
        self.assertEqual(len(lr_manager_1.planes_front), 1)
        self.assertEqual(len(lr_manager_1.planes_previous), 0)
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["i"], [0, 1, 2, 3])
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["j"], [2])

        current_lr1 = [np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]])]
        lr_manager_1.set_current_lr(current_lr1)

        # add to seen
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1

        # test next() MOVE Plan
        lr_manager_1.next()
        # add to seen
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1

        self.assertEqual(len(lr_manager_1.planes_front), 1)
        self.assertEqual(len(lr_manager_1.planes_previous), 1)
        # current
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["i"], [0, 1, 2, 3])
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["j"], [1])
        # previous
        np.testing.assert_array_equal(lr_manager_1.planes_previous[0]["i"], [0, 1, 2, 3])
        np.testing.assert_array_equal(lr_manager_1.planes_previous[0]["j"], [2])

        # Test lr saving + elimination
        previous_lr1 = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]])
        np.testing.assert_array_equal(lr_manager_1.get_previous_lr(0), previous_lr1)

        # test next() MOVE Plan out of dimensions AND add to seen
        lr_manager_1.next()
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1
        lr_manager_1.next()

        # out of plan
        self.assertEqual(len(lr_manager_1.planes_front), 0)
        self.assertEqual(len(lr_manager_1.planes_previous), 0)

        # Test if all pixels were seen only once
        np.testing.assert_array_equal(seen_pixels, np.ones((shape1[0], shape1[1])))

    def test_dir_1_1_lr_manager(self):
        """ "
        Test direction (1, 1)
        """
        shape1 = [4, 3, 3]
        direction1 = (1, 1)
        lr_manager_1 = lrm.LrManager(shape1, direction1)

        seen_pixels = np.zeros((shape1[0], shape1[1]))

        # test init
        self.assertEqual(len(lr_manager_1.planes_front), 2)
        self.assertEqual(len(lr_manager_1.planes_previous), 0)
        # plan 1
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["i"], [0])
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["j"], [0, 1, 2])
        # plan 2
        np.testing.assert_array_equal(lr_manager_1.planes_front[1]["i"], [1, 2, 3])  # 0 is in plan 1
        np.testing.assert_array_equal(lr_manager_1.planes_front[1]["j"], [0])

        current_lr1 = [np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]), np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]) + 12]
        lr_manager_1.set_current_lr(current_lr1)

        # add to seen
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1
        seen_pixels[lr_manager_1.planes_front[1]["i"], lr_manager_1.planes_front[1]["j"]] += 1

        # test next() MOVE Plan
        lr_manager_1.next()
        # add to seen
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1
        seen_pixels[lr_manager_1.planes_front[1]["i"], lr_manager_1.planes_front[1]["j"]] += 1

        self.assertEqual(len(lr_manager_1.planes_front), 2)
        self.assertEqual(len(lr_manager_1.planes_previous), 2)
        # current
        # plan 1
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["i"], [1])
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["j"], [1, 2])
        # plan 2
        np.testing.assert_array_equal(lr_manager_1.planes_front[1]["i"], [2, 3])
        np.testing.assert_array_equal(lr_manager_1.planes_front[1]["j"], [1])
        # previous
        # plan 1
        np.testing.assert_array_equal(lr_manager_1.planes_previous[0]["i"], [0])
        np.testing.assert_array_equal(lr_manager_1.planes_previous[0]["j"], [0, 1])
        # plan 2
        np.testing.assert_array_equal(lr_manager_1.planes_previous[1]["i"], [1, 2])
        np.testing.assert_array_equal(lr_manager_1.planes_previous[1]["j"], [0])

        # Test lr saving + elimination
        previous_lr1 = np.array([[1, 2, 3], [4, 5, 6]])
        previous_lr2 = np.array([[1, 2, 3], [4, 5, 6]]) + 12

        np.testing.assert_array_equal(lr_manager_1.get_previous_lr(0), previous_lr1)
        np.testing.assert_array_equal(lr_manager_1.get_previous_lr(1), previous_lr2)

        # test next() MOVE Plan out of dimensions AND add to seen
        lr_manager_1.next()
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1
        seen_pixels[lr_manager_1.planes_front[1]["i"], lr_manager_1.planes_front[1]["j"]] += 1
        lr_manager_1.next()

        # out of plan
        self.assertEqual(len(lr_manager_1.planes_front), 0)
        self.assertEqual(len(lr_manager_1.planes_previous), 0)

        # Test if all pixels were seen only once
        np.testing.assert_array_equal(seen_pixels, np.ones((shape1[0], shape1[1])))

    def test_dir_minus1_1_lr_manager(self):
        """ "
        Test direction (-1, 1)
        """
        shape1 = [4, 3, 3]
        direction1 = (-1, 1)
        lr_manager_1 = lrm.LrManager(shape1, direction1)

        seen_pixels = np.zeros((shape1[0], shape1[1]))

        # test init
        self.assertEqual(len(lr_manager_1.planes_front), 2)
        self.assertEqual(len(lr_manager_1.planes_previous), 0)
        # plan 1
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["i"], [3])
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["j"], [0, 1, 2])
        # plan 2
        np.testing.assert_array_equal(lr_manager_1.planes_front[1]["i"], [0, 1, 2])  # 3 is in plan 1
        np.testing.assert_array_equal(lr_manager_1.planes_front[1]["j"], [0])

        current_lr1 = [np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]), np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]) + 12]
        lr_manager_1.set_current_lr(current_lr1)  #
        # add to seen
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1
        seen_pixels[lr_manager_1.planes_front[1]["i"], lr_manager_1.planes_front[1]["j"]] += 1

        # test next() MOVE Plan
        lr_manager_1.next()
        # add to seen
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1
        seen_pixels[lr_manager_1.planes_front[1]["i"], lr_manager_1.planes_front[1]["j"]] += 1

        self.assertEqual(len(lr_manager_1.planes_front), 2)
        self.assertEqual(len(lr_manager_1.planes_previous), 2)
        # current
        # plan 1
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["i"], [2])
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["j"], [1, 2])
        # plan 2
        np.testing.assert_array_equal(lr_manager_1.planes_front[1]["i"], [0, 1])
        np.testing.assert_array_equal(lr_manager_1.planes_front[1]["j"], [1])
        # previous
        # plan 1
        np.testing.assert_array_equal(lr_manager_1.planes_previous[0]["i"], [3])
        np.testing.assert_array_equal(lr_manager_1.planes_previous[0]["j"], [0, 1])
        # plan 2
        np.testing.assert_array_equal(lr_manager_1.planes_previous[1]["i"], [1, 2])
        np.testing.assert_array_equal(lr_manager_1.planes_previous[1]["j"], [0])

        # Test lr saving + elimination
        previous_lr1 = np.array([[1, 2, 3], [4, 5, 6]])
        previous_lr2 = np.array([[4, 5, 6], [7, 8, 9]]) + 12

        np.testing.assert_array_equal(lr_manager_1.get_previous_lr(0), previous_lr1)
        np.testing.assert_array_equal(lr_manager_1.get_previous_lr(1), previous_lr2)

        # test next() MOVE Plan out of dimensions AND add to seen
        lr_manager_1.next()
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1
        seen_pixels[lr_manager_1.planes_front[1]["i"], lr_manager_1.planes_front[1]["j"]] += 1
        lr_manager_1.next()

        # out of plan
        self.assertEqual(len(lr_manager_1.planes_front), 0)
        self.assertEqual(len(lr_manager_1.planes_previous), 0)

        # Test if all pixels were seen only once
        np.testing.assert_array_equal(seen_pixels, np.ones((shape1[0], shape1[1])))

    def test_dir_1_minus1_lr_manager(self):
        """ "
        Test direction (1, -1)
        """
        shape1 = [4, 3, 3]
        direction1 = (1, -1)
        lr_manager_1 = lrm.LrManager(shape1, direction1)

        seen_pixels = np.zeros((shape1[0], shape1[1]))

        # test init
        self.assertEqual(len(lr_manager_1.planes_front), 2)
        self.assertEqual(len(lr_manager_1.planes_previous), 0)
        # plan 1
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["i"], [0])
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["j"], [0, 1, 2])
        # plan 2
        np.testing.assert_array_equal(lr_manager_1.planes_front[1]["i"], [1, 2, 3])  # 0 is in plan 1
        np.testing.assert_array_equal(lr_manager_1.planes_front[1]["j"], [2])

        current_lr1 = [np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]), np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]) + 12]
        lr_manager_1.set_current_lr(current_lr1)

        # add to seen
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1
        seen_pixels[lr_manager_1.planes_front[1]["i"], lr_manager_1.planes_front[1]["j"]] += 1

        # test next() MOVE Plan
        lr_manager_1.next()
        # add to seen
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1
        seen_pixels[lr_manager_1.planes_front[1]["i"], lr_manager_1.planes_front[1]["j"]] += 1

        self.assertEqual(len(lr_manager_1.planes_front), 2)
        self.assertEqual(len(lr_manager_1.planes_previous), 2)
        # current
        # plan 1
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["i"], [1])
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["j"], [0, 1])
        # plan 2
        np.testing.assert_array_equal(lr_manager_1.planes_front[1]["i"], [2, 3])
        np.testing.assert_array_equal(lr_manager_1.planes_front[1]["j"], [1])
        # previous
        # plan 1
        np.testing.assert_array_equal(lr_manager_1.planes_previous[0]["i"], [0])
        np.testing.assert_array_equal(lr_manager_1.planes_previous[0]["j"], [1, 2])
        # plan 2
        np.testing.assert_array_equal(lr_manager_1.planes_previous[1]["i"], [1, 2])
        np.testing.assert_array_equal(lr_manager_1.planes_previous[1]["j"], [2])

        # Test lr saving + elimination
        previous_lr1 = np.array([[4, 5, 6], [7, 8, 9]])
        previous_lr2 = np.array([[1, 2, 3], [4, 5, 6]]) + 12

        np.testing.assert_array_equal(lr_manager_1.get_previous_lr(0), previous_lr1)
        np.testing.assert_array_equal(lr_manager_1.get_previous_lr(1), previous_lr2)

        # test next() MOVE Plan out of dimensions AND add to seen
        lr_manager_1.next()
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1
        seen_pixels[lr_manager_1.planes_front[1]["i"], lr_manager_1.planes_front[1]["j"]] += 1
        lr_manager_1.next()

        # out of plan
        self.assertEqual(len(lr_manager_1.planes_front), 0)
        self.assertEqual(len(lr_manager_1.planes_previous), 0)

        # Test if all pixels were seen only once
        np.testing.assert_array_equal(seen_pixels, np.ones((shape1[0], shape1[1])))

    def test_dir_minus1_minus1_lr_manager(self):
        """ "
        Test direction (-1, -1)
        """

        shape1 = [4, 3, 3]
        direction1 = (-1, -1)
        lr_manager_1 = lrm.LrManager(shape1, direction1)

        seen_pixels = np.zeros((shape1[0], shape1[1]))

        # test init
        self.assertEqual(len(lr_manager_1.planes_front), 2)
        self.assertEqual(len(lr_manager_1.planes_previous), 0)
        # plan 1
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["i"], [3])
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["j"], [0, 1, 2])
        # plan 2
        np.testing.assert_array_equal(lr_manager_1.planes_front[1]["i"], [0, 1, 2])  # 3 is in plan 1
        np.testing.assert_array_equal(lr_manager_1.planes_front[1]["j"], [2])

        current_lr1 = [np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]), np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]) + 12]
        lr_manager_1.set_current_lr(current_lr1)

        current_segm1 = [np.array([1, 2, 3]), np.array([4, 5, 6])]
        lr_manager_1.set_current_segm(current_segm1)

        # add to seen
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1
        seen_pixels[lr_manager_1.planes_front[1]["i"], lr_manager_1.planes_front[1]["j"]] += 1

        # test next() MOVE Plan
        lr_manager_1.next()
        # add to seen
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1
        seen_pixels[lr_manager_1.planes_front[1]["i"], lr_manager_1.planes_front[1]["j"]] += 1

        self.assertEqual(len(lr_manager_1.planes_front), 2)
        self.assertEqual(len(lr_manager_1.planes_previous), 2)
        # current
        # plan 1
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["i"], [2])
        np.testing.assert_array_equal(lr_manager_1.planes_front[0]["j"], [0, 1])
        # plan 2
        np.testing.assert_array_equal(lr_manager_1.planes_front[1]["i"], [0, 1])
        np.testing.assert_array_equal(lr_manager_1.planes_front[1]["j"], [1])
        # previous
        # plan 1
        np.testing.assert_array_equal(lr_manager_1.planes_previous[0]["i"], [3])
        np.testing.assert_array_equal(lr_manager_1.planes_previous[0]["j"], [1, 2])
        # plan 2
        np.testing.assert_array_equal(lr_manager_1.planes_previous[1]["i"], [1, 2])
        np.testing.assert_array_equal(lr_manager_1.planes_previous[1]["j"], [2])

        # Test lr saving + elimination
        previous_lr1 = np.array([[4, 5, 6], [7, 8, 9]])
        previous_lr2 = np.array([[4, 5, 6], [7, 8, 9]]) + 12

        np.testing.assert_array_equal(lr_manager_1.get_previous_lr(0), previous_lr1)
        np.testing.assert_array_equal(lr_manager_1.get_previous_lr(1), previous_lr2)

        # Test segm saving + elimination
        previous_segm1 = np.array([2, 3])
        previous_segm2 = np.array([5, 6])

        np.testing.assert_array_equal(lr_manager_1.get_previous_segm(0), previous_segm1)
        np.testing.assert_array_equal(lr_manager_1.get_previous_segm(1), previous_segm2)

        # test next() MOVE Plan out of dimensions AND add to seen
        lr_manager_1.next()
        seen_pixels[lr_manager_1.planes_front[0]["i"], lr_manager_1.planes_front[0]["j"]] += 1
        seen_pixels[lr_manager_1.planes_front[1]["i"], lr_manager_1.planes_front[1]["j"]] += 1
        lr_manager_1.next()

        # out of plan
        self.assertEqual(len(lr_manager_1.planes_front), 0)
        self.assertEqual(len(lr_manager_1.planes_previous), 0)

        # Test if all pixels were seen only once
        np.testing.assert_array_equal(seen_pixels, np.ones((shape1[0], shape1[1])))
