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

import numpy as np

from .lr_manager import LrManager


def run_sgm(
    cv_in: np.ndarray,
    p1_in: np.ndarray,
    p2_in: np.ndarray,
    directions: list,
    segmentation: np.ndarray,
    cost_paths: bool = False,
    overcounting: bool = False,
):
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
    :param segmentation: segmentation for piecewise optimization
    :type segmentation: 2D np.ndarray
    :param cost_paths: True if Cost Volumes along direction are to be returned
    :type cost_paths: bool
    :param overcounting: over-counting correction option
    :type overcounting: bool
    :return: the aggregated cost volume and the minimum cost along directions
    :rtype: dict of 3D np.ndarray
    """

    cv_out = {"cv": np.zeros(cv_in.shape)}

    if cost_paths:
        cv_out["cv_min"] = np.zeros((cv_in.shape[0], cv_in.shape[1], len(directions)))

    for idx_dir in range(len(directions)):  # pylint:disable=consider-using-enumerate
        # optimize this direction
        lr_manager = LrManager(cv_in.shape, directions[idx_dir])  # type: ignore

        nb_planes = len(lr_manager.planes_front)
        while nb_planes > 0:
            current_lr = []
            current_segm = []

            for p_idx in range(nb_planes):
                front_plane = lr_manager.planes_front[p_idx]
                if len(lr_manager.planes_previous) == 0:
                    # initialized, previous plane doesnt exist
                    partial_lr = cv_in[front_plane["i"], front_plane["j"], :]
                    cv_out["cv"][front_plane["i"], front_plane["j"], :] += partial_lr
                    current_lr.append(partial_lr)
                    partial_segm = segmentation[front_plane["i"], front_plane["j"]]
                    current_segm.append(partial_segm)
                else:
                    partial_lr = np.zeros(cv_in[front_plane["i"], front_plane["j"], :].shape)
                    previous_lr = lr_manager.get_previous_lr(p_idx)
                    previous_segm = lr_manager.get_previous_segm(p_idx)
                    partial_segm = segmentation[front_plane["i"], front_plane["j"]]
                    for idx in range(cv_in.shape[2]):
                        partial_lr[:, idx] = compute_lr(
                            cv_in[front_plane["i"], front_plane["j"], :],
                            previous_lr,
                            idx,
                            p1_in[front_plane["i"], front_plane["j"], idx_dir],
                            p2_in[front_plane["i"], front_plane["j"], idx_dir],
                            partial_segm,
                            previous_segm,
                        )

                    cv_out["cv"][front_plane["i"], front_plane["j"]] += partial_lr
                    current_lr.append(partial_lr)
                    current_segm.append(partial_segm)

                if cost_paths:
                    cv_out["cv_min"][front_plane["i"], front_plane["j"], idx_dir] = np.argmin(partial_lr, axis=1)

            lr_manager.set_current_lr(current_lr)
            lr_manager.set_current_segm(current_segm)

            # compute next planes
            lr_manager.next()  # pylint:disable=not-callable
            nb_planes = len(lr_manager.planes_front)

    if overcounting:
        cv_out["cv"] -= (len(directions) - 1) * cv_in

    return cv_out


def compute_lr(
    cv_in_2d_front: np.ndarray,
    lr_2d_previous: np.ndarray,
    disp: int,
    p1_in_1d: np.ndarray,
    p2_in_1d: np.ndarray,
    current_segm_1d: np.ndarray,
    previous_segm_1d: np.ndarray,
):
    """
    Compute Lr of current plane, at a given disparity

    :param cv_in_2d_front: cost volume
    :type cv_in_2d_front:  np.ndarray
    :param lr_2d_previous: previous lr computed
    :type lr_2d_previous: np.ndarray
    :param d: current disparity
    :type d:  int
    :param p1_in_1d: p1 penalties, dim=1
    :type p1_in_1d: np.ndarray
    :param p2_in_1d: p2 penalties, dim=1
    :type p2_in_1d: np.ndarray
    :param current_segm_1d: current segmentation
    :type current_segm_1d: np.ndarray
    :param previous_segm_1d:  previous segmentation
    :type previous_segm_1d: np.ndarray
    :return: partial lrs
    :rtype: np.ndarray
    """

    # Assign penalties
    penalties = np.repeat(p2_in_1d[:, np.newaxis], cv_in_2d_front.shape[1], axis=1)
    penalties[:, disp] = 0
    if (disp - 1) >= 0:
        penalties[:, disp - 1] = p1_in_1d
    if (disp + 1) < cv_in_2d_front.shape[1]:
        penalties[:, disp + 1] = p1_in_1d

    min_previous_lr_penalties = np.nanmin(lr_2d_previous + penalties, axis=1)
    min_previous_lr = np.nanmin(lr_2d_previous, axis=1)
    indexes_nan = np.where(np.isnan(min_previous_lr_penalties))
    min_previous_lr_penalties[indexes_nan] = 0
    min_previous_lr[indexes_nan] = 0

    # reset if different class
    reset = current_segm_1d == previous_segm_1d

    return cv_in_2d_front[:, disp] + reset * (min_previous_lr_penalties - min_previous_lr)
