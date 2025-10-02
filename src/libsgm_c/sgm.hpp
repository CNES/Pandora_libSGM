/*
 * Copyright (c) 2024 Centre National d'Etudes Spatiales (CNES).
 *
 * This file is part of LIBSGM
 *
 *     https://github.com/CNES/Pandora_libsgm
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef LIBSGM_SGM_HPP
#define LIBSGM_SGM_HPP

#include <stdint.h>
#include <limits>

/**
 * Structure to represent Aggregated Cost Volume and the positions of minimum costs along each direction
 */
template<typename T>
struct CostVolumes{
    T * cost_volume; /**< Aggregated Cost Volume */
    int * cost_volume_min; /**< positions of minimum costs along each direction */
};

/**
 * Structure to represent coordinates of previous point path
 */
struct Direction{
    size_t drow; /**< row coordinate */
    size_t dcol; /**< col coordinate */
};

/**
 * Structure to represent Penalty
 */
template<typename T>
struct Penalty{
    T P1; /**< Penalty 1 */
    T P2; /**< Penalty 2 */
};

// Top-left pass (first four directions)
template<typename T>
T aggregatedCostFromTopLeft0(T pixelCost, int disp, T invalid_value,
    int nb_disps, T P1, T P2, Direction direction,
    T * buff0, T min_disp0, T * pixel_0);

template<typename T>
T aggregatedCostFromTopLeft1(T pixelCost, int col, int disp, T invalid_value,
    int nb_disps, T P1, T P2, T * buff1, T min_disp1, T * pixel_1);

template<typename T>
T aggregatedCostFromTopLeft2(T pixelCost, int col, int disp, T invalid_value,
    int nb_disps, T P1, T P2, T * buff2, T * buff_disp_2, T * min_disp2, T * pixel_2);

template<typename T>
T aggregatedCostFromTopLeft3(T pixelCost, int col, int disp, T invalid_value,
    int nb_disps, T P1, T P2, Direction direction, T * buff3, T min_disp3);

// Bottom-right pass (last four directions)
template<typename T>
T aggregatedCostFromBottomRight4(T pixelCost, int row, int col, int disp, T invalid_value,
    int nb_rows, int nb_cols, int nb_disps, T P1, T P2, Direction direction,
    T * buff4, T * min_disp4, T * pixel_4, float current_class, float buff_class4, float & reset4);

template<typename T>
T aggregatedCostFromBottomRight5(T pixelCost, int row, int col, int disp, T invalid_value,
    int nb_rows, int nb_cols, int nb_disps, T P1, T P2, Direction direction,
    T * buff5, T * min_disp5, T * pixel_5, float current_class, float * buff_class5, float & reset5);

template<typename T>
T aggregatedCostFromBottomRight6(T pixelCost, int row, int col, int disp, T invalid_value,
    int nb_rows, int nb_cols, int nb_disps, T P1, T P2, Direction direction,
    T * buff6, T * buff_disp_6, T * min_disp6, T * pixel_6, float current_class, float * buff_class6, float & reset6);

template<typename T>
T aggregatedCostFromBottomRight7(T pixelCost, int row, int col, int disp, T invalid_value,
    int nb_rows, int nb_cols, int nb_disps, T P1, T P2, Direction direction, T * buff7, T * min_disp7,
    float current_class, float * buff_class7, float & reset7);

template<typename T>
void apply_penalty(Penalty<T> *penalty, T p1, T p2);

void  assignDirections(int* directions_in, Direction* dirs);

#endif // LIBSGM_SGM_HPP