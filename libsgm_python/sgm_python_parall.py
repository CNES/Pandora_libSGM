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
from numba import njit, prange
from typing import List

from .lr_manager import LrManager


def run_sgm_parall(cv_in: np.ndarray, p1_in: np.ndarray, p2_in: np.ndarray, directions: np.ndarray,
                   cost_paths: bool = False, overcounting: bool = False):
    """
    Run Python LibSGM

    :param cv_in: cost volume
    :type cv_in:  np.ndarray
    :param p1_in: p1 penalties
    :type p1_in:  3D np.ndarray
    :param p2_in: p2 penalties
    :type p2_in:  3D np.ndarray
    :param directions: directions used in SGM
    :type directions: 2D np.ndarray
    :param cost_paths: True if Cost Volumes along direction are to be returned
    :type cost_paths: bool
    :param overcounting: over-counting correction option
    :type overcounting: bool
    :return: the aggregated cost volume and the minimum cost along 8 directions
    :rtype: dict of 3D np.ndarray
    """

    # Compute initial points
    # [i0, j0, di, dj, num_dir]
    starting_points = []
    for i in range(len(directions)):
        d = directions[i]
        # optimize this direction
        lr_manager = LrManager(cv_in.shape, directions[i])

        for plane in lr_manager.planes_front:
            list_indexes = [[i0, j0, d[0], d[1], i] for idx1, i0 in enumerate(plane["i"]) for idx2, j0 in
                            enumerate(plane["j"])]
            starting_points += list_indexes

    # compute costs for all
    dir_str = []
    for i in range(len(directions)):
        dir_str.append(str(i))
    dir_str = np.array(dir_str)
    cv_out = compute_costs(np.array(starting_points), cv_in, p1_in, p2_in, dir_str, cost_paths=cost_paths)

    if overcounting:
        cv_out["cv"] -= (len(directions) - 1) * cv_in

    return cv_out


@njit(parallel=True, cache=False)
def compute_costs(starting_points: List, cv_in: np.ndarray, p1_in: np.ndarray, p2_in: np.ndarray, dir_str: np.ndarray,
                  cost_paths: bool = False):
    """
    Compute cost volume, starting from given points / directions

    :param starting_points: List of starting points [i0, j0, dir_i, dir_j, num_dir]
    :type starting_points:  List
    :param cv_in: cost volume
    :type cv_in:  3D np.ndarray
    :param p1_in: p1 penalties
    :type p1_in:  3D np.ndarray
    :param p2_in: p2 penalties
    :type p2_in:  3D np.ndarray
    :param dir_str: list of direction names
    :type dir_str:  List
    :param cost_paths: True if Cost Volumes along direction are to be returned
    :type cost_paths: bool
    :return: the aggregated cost volume and the minimum cost along directions
    :rtype: dict of 3D np.ndarray
    """
    in_shape = cv_in.shape
    cv = {"cv": np.zeros(in_shape)}

    for i in range(len(dir_str)):
        name = "cv_" + dir_str[i]
        cv[name] = np.zeros(in_shape)

    cv_min = np.zeros((in_shape[0], in_shape[1], p1_in.shape[2]))

    for p in prange(starting_points.shape[0]):
        point = starting_points[p]

        d_i = point[2]
        d_j = point[3]
        i_0 = point[0]
        j_0 = point[1]
        dir = point[4]
        name = "cv_" + dir_str[dir]

        cv["cv"][i_0, j_0, :] += cv_in[i_0, j_0, :].copy()
        cv[name][i_0, j_0, :] += cv_in[i_0, j_0, :].copy()
        previous_lr = cv_in[i_0, j_0, :].copy()

        i = i_0 + d_i
        j = j_0 + d_j

        while 0 <= i < in_shape[0] and 0 <= j < in_shape[1]:

            current_lr = cv_in[i, j, :].copy()

            p1 = p1_in[i, j, dir]
            p2 = p2_in[i, j, dir]

            for d in range(in_shape[2]):
                penalties = p2 * np.ones(in_shape[2])
                penalties[max(0, d-1): min(d+1, in_shape[2]-1) + 1] = p1
                penalties[d] = 0

                val = np.nanmin(penalties + previous_lr) - np.nanmin(previous_lr)
                if not np.isnan(val):
                    current_lr[d] += val

            # save current lr
            cv["cv"][i, j, :] += current_lr
            cv[name][i, j, :] += current_lr

            # next pixel
            previous_lr = current_lr

            cv_min[i, j, dir] = np.argmin(current_lr)

            i += d_i
            j += d_j

    if cost_paths:
        cv["cv_min"] = cv_min

    return cv
