# pylint:disable=fixme
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
This module contains functions to execute python SGM library.
"""
import copy
from typing import List, Tuple

import numpy as np


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
        self.dir = direction

        # init working matrices
        self.planes_front: List = []
        self.planes_previous: List = []

        # store partial lr
        self.current_lr = None
        self.previous_lr = None

        # store partial segm
        self.current_segm = None
        self.previous_segm = None

        # start and end indices to manage duplicated corners when two planes are used
        start = 0
        end = cv_shape[0]
        if direction[0] > 0:
            # TODO deal with directions > 1 : several planes, be careful about duplicated points
            self.planes_front.append({"i": np.array([0]), "j": np.array(list(range(0, cv_shape[1])))})
            start = 1
        if direction[0] < 0:
            self.planes_front.append(
                {
                    "i": np.array([cv_shape[0] - 1]),
                    "j": np.array(list(range(0, cv_shape[1]))),
                }
            )
            end = end - 1
        if direction[1] > 0:
            self.planes_front.append({"i": np.array(list(range(start, end))), "j": np.array([0])})
        if direction[1] < 0:
            self.planes_front.append(
                {
                    "i": np.array(list(range(start, end))),
                    "j": np.array([cv_shape[1] - 1]),
                }
            )

    def next(self) -> None:
        """
        Update the state of the planes, moves forward
        """
        self.planes_previous = copy.deepcopy(self.planes_front)
        self.previous_lr = self.current_lr
        self.previous_segm = self.current_segm
        delete_items_indexes = []
        for idx, _ in enumerate(self.planes_front):
            self.planes_front[idx]["i"] += self.dir[0]
            self.planes_front[idx]["j"] += self.dir[1]

            # get invalid indexes
            indexes_invalid_i = np.where(
                (self.planes_front[idx]["i"] < 0) + (self.planes_front[idx]["i"] >= self.cv_shape[0])
            )
            indexes_invalid_j = np.where(
                (self.planes_front[idx]["j"] < 0) + (self.planes_front[idx]["j"] >= self.cv_shape[1])
            )

            # remove invalid indexes
            self.planes_front[idx]["i"] = np.delete(self.planes_front[idx]["i"], indexes_invalid_i, axis=0)
            self.planes_previous[idx]["i"] = np.delete(self.planes_previous[idx]["i"], indexes_invalid_i, axis=0)

            self.planes_front[idx]["j"] = np.delete(self.planes_front[idx]["j"], indexes_invalid_j, axis=0)
            self.planes_previous[idx]["j"] = np.delete(self.planes_previous[idx]["j"], indexes_invalid_j, axis=0)

            if self.planes_front[idx]["i"].shape[0] == 0 or self.planes_front[idx]["j"].shape[0] == 0:
                # delete all object
                delete_items_indexes.append(idx)
            else:
                # delete in stored lr, as list are of size (1) (n) or (m) (1), only one delete is effective

                if self.previous_lr is not None:
                    self.previous_lr[idx] = np.delete(self.previous_lr[idx], indexes_invalid_i, axis=0)
                    self.previous_lr[idx] = np.delete(self.previous_lr[idx], indexes_invalid_j, axis=0)

                # delete in stored segm, as list are of size (1) (n) or (m) (1), only one delete is effective

                if self.previous_segm is not None:
                    self.previous_segm[idx] = np.delete(self.previous_segm[idx], indexes_invalid_i, axis=0)
                    self.previous_segm[idx] = np.delete(self.previous_segm[idx], indexes_invalid_j, axis=0)

        # delete items
        for idx in sorted(delete_items_indexes, reverse=True):
            self.planes_front.pop(idx)
            self.planes_previous.pop(idx)
            if self.previous_lr is not None:
                self.previous_lr.pop(idx)
            if self.previous_segm is not None:
                self.previous_segm.pop(idx)

    def set_current_lr(self, lr_s: list) -> None:
        """
        Set current lr

        :param lr_s: partial cost lr
        :type lr_s:  np.ndarray
        """
        self.current_lr = lr_s  # type: ignore

    def get_previous_lr(self, num_plane: int) -> np.ndarray:
        """
        Get previous lr

        :param num_plane: numero of plane
        :type num_plane:  int
        :return: previous partial cost lr
        :rtype: np.ndarray
        """
        return self.previous_lr[num_plane]  # type:ignore

    def set_current_segm(self, segm_s: list) -> None:
        """
        Set current segmentation

        :param lr_s: partial segmentation
        :type lr_s:  np.ndarray
        """
        self.current_segm = segm_s  # type: ignore

    def get_previous_segm(self, num_plane: int) -> np.ndarray:
        """
        Get previous  segmentation

        :param num_plane: numero of plane
        :type num_plane:  int
        :return: previous partial segmentation
        :rtype: np.ndarray
        """
        return self.previous_segm[num_plane]  # type:ignore
