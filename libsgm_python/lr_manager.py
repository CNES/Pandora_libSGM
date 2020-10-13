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
This module contains functions to execute python SGM library.
"""

import numpy as np
import copy
from typing import List, Tuple


class LrManager:
    """
    This class manages parallel planes indices moving along given direction
    """
    # initialize
    def __init__(self, cv_shape: List, direction: Tuple):
        """
        Initialise LrManager

        :param cv_shape: shape of cost volume
        :type cv_shape:  List
        :param direction: direction used (d_i, d_j)
        :type direction:  Tuple
        """
        self.cv_shape = cv_shape
        self.d = direction

        # init working matrices
        self.planes_front = []
        self.planes_previous = []

        # store partial lr
        self.current_lr = None
        self.previous_lr = None

        # start and end indices to manage duplicated corners when two planes are used
        start = 0
        end = cv_shape[0]
        if direction[0] > 0:
            # TODO deal with directions > 1 : several planes, be careful about duplicated points
            self.planes_front.append({"i": np.array([0]), "j": np.array(list(range(0, cv_shape[1])))})
            start = 1
        if direction[0] < 0:
            self.planes_front.append({"i": np.array([cv_shape[0] - 1]), "j": np.array(list(range(0, cv_shape[1])))})
            end = end - 1
        if direction[1] > 0:
            self.planes_front.append({"i": np.array(list(range(start, end))), "j": np.array([0])})
        if direction[1] < 0:
            self.planes_front.append({"i": np.array(list(range(start, end))), "j": np.array([cv_shape[1] - 1])})

    def next(self) -> None:
        """
        Update the state of the planes, moves forward
        """
        self.planes_previous = copy.deepcopy(self.planes_front)
        self.previous_lr = self.current_lr
        delete_items_indexes = []
        for p in range(len(self.planes_front)):
            self.planes_front[p]["i"] += self.d[0]
            self.planes_front[p]["j"] += self.d[1]

            # get invalid indexes
            indexes_invalid_i = np.where(
                (self.planes_front[p]["i"] < 0) + (self.planes_front[p]["i"] >= self.cv_shape[0]))
            indexes_invalid_j = np.where(
                (self.planes_front[p]["j"] < 0) + (self.planes_front[p]["j"] >= self.cv_shape[1]))

            # remove invalid indexes
            self.planes_front[p]["i"] = np.delete(self.planes_front[p]["i"], indexes_invalid_i, axis=0)
            self.planes_previous[p]["i"] = np.delete(self.planes_previous[p]["i"], indexes_invalid_i, axis=0)

            self.planes_front[p]["j"] = np.delete(self.planes_front[p]["j"], indexes_invalid_j, axis=0)
            self.planes_previous[p]["j"] = np.delete(self.planes_previous[p]["j"], indexes_invalid_j, axis=0)

            if self.planes_front[p]["i"].shape[0] == 0 or self.planes_front[p]["j"].shape[0] == 0:
                # delete all object
                delete_items_indexes.append(p)
            else:
                # delete in stored lr, as list are of size (1) (n) or (m) (1), only one delete is effective

                if self.previous_lr is not None:
                    self.previous_lr[p] = np.delete(self.previous_lr[p], indexes_invalid_i, axis=0)
                    self.previous_lr[p] = np.delete(self.previous_lr[p], indexes_invalid_j, axis=0)

        # delete items
        for p in sorted(delete_items_indexes, reverse=True):
            self.planes_front.pop(p)
            self.planes_previous.pop(p)
            if self.previous_lr is not None:
                self.previous_lr.pop(p)

    def set_current_lr(self, lr_s: np.ndarray) -> None:
        """
        Set current lr

        :param lr_s: partial cost lr
        :type lr_s:  np.ndarray
        """
        self.current_lr = lr_s

    def get_previous_lr(self, num_plane: int) -> np.ndarray:
        """
        Get previous lr

        :param num_plane: numero of plane
        :type num_plane:  int
        :return: previous partial cost lr
        :rtype: np.ndarray
        """
        return self.previous_lr[num_plane]
