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

#include <stdint.h>

/**
* Structure to represent Aggregated Cost Volume and the positions of minimum costs along each direction
*/
struct CostVolumes{
    float * cost_volume; /**< Aggregated Cost Volume */
    int * cost_volume_min; /**< positions of minimum costs along each direction */
};


/**
* Structure to represent coordinates of previous point path
*/

struct Direction{
    int drow; /**< row coordinate */
    int dcol; /**< col coordinate */
};

/**
* Structure to represent Penalty
*/
template<typename T>
struct Penalty{
    T P1; /**< Penalty 1 */
    T P2; /**< Penalty 2 */
};

/*!
 *  \brief  Compute aggregated cost volume
 *   Compute aggregated cost volume following
 *   Semi-Global algorithm by Hirschmuller
 *
 *  \param cv_in cost volume
 *  \param p1_in p1 penalty
 *  \param p2_in p2 penalty
 *  \param directions_in directions to use
 *  \param nb_rows row number of cost volume
 *  \param nb_cols column number of cost volume
 *  \param nb_disps disparity number of cost volume
 *  \param invalid_value value representing invalid cost
 *  \param segmentation segmentation map
 *  \param cost_paths True if Cost Volumes along direction are to be returned
 *  \param overcounting over-counting correction option
 *  \return cost volume aggregated, minimum cost on each direction
 */


template<typename T>
CostVolumes sgm(T * cv_in, T* p1_in, T* p2_in, int* directions_in, unsigned long int nb_rows, unsigned long int nb_cols,
 unsigned int nb_disps, T invalid_value, float* segmentation, bool cost_paths, bool overcounting);

/*!
 *  Update minimum
 *
 *  \param current_min current minimum
 *  \param value potential new minimum
 *  \param current_disp current disparity
 *  \param disp potential disparity
 *  \return new minimum, new disparity
 */

std::pair<float, int> update_minimum(float current_min, float value, int current_disp, int disp);

/*!
 *  \brief  Compute aggregated cost volume
 *   Compute aggregated cost at one point
 *   from left to right direction
 *
 *  \param pixelCost cost value 
 *  \param row row position
 *  \param col col position
 *  \param disp disp position
 *  \param invalid_value value representing invalid cost
 *  \param nb_rows row number of cost volume
 *  \param nb_cols column number of cost volume
 *  \param nb_disps disparity number of cost volume
 *  \param P1 penalty P1 from sgm equation
 *  \param P2 penalty P2 from sgm equation
 *  \param direction array of struct containing coordinates of previous point
 *  \param buff0 buffer containing previous aggregated cost (previous point)
 *  \param min_disp0 value of minimum cost at the previous pixel
 *  \param pixel_0 value of previous disp on the buffer
 *  \param current_class current class
 *  \param buff_class0 buffer containing previous class (previous point
 *  \param reset0 value of coefficient to multiply history
 *  \return cost aggregated point
 */

template<typename T>
T aggregatedCostFromTopLeft0(T pixelCost, int row, int col, int disp, T invalid_value,
    int nb_rows, int nb_cols, int nb_disps, T P1, T P2, Direction direction,
    T * buff0, T * min_disp0, T * pixel_0, float current_class, float buff_class0, float & reset0);

/*!
 *  \brief  Compute aggregated cost volume
 *   Compute aggregated cost at one point
 *   from up to down direction
 *
 *  \param pixelCost cost value 
 *  \param row row position
 *  \param col col position
 *  \param disp disp position
 *  \param invalid_value value representing invalid cost
 *  \param nb_rows row number of cost volume
 *  \param nb_cols column number of cost volume
 *  \param nb_disps disparity number of cost volume
 *  \param P1 penalty P1 from sgm equation
 *  \param P2 penalty P2 from sgm equation
 *  \param direction array of struct containing coordinates of previous point
 *  \param buff1 buffer containing previois aggregated cost (previous line)
 *  \param min_disp1 value of minimum cost at the previous pixel
 *  \param pixel_1 value of previous disp on the buffer
 *  \param current_class current class
 *  \param buff_class1 buffer containing previous class (previous point
 *  \param reset1 value of coefficient to multiply history
 *  \return cost aggregated point
 */

template<typename T> 
T aggregatedCostFromTopLeft1(T pixelCost, int row, int col, int disp, T invalid_value,
    int nb_rows, int nb_cols, int nb_disps, T P1, T P2, Direction direction,
    T * buff1, T * min_disp1, T * pixel_1, float current_class, float * buff_class1, float & reset1);
/*!
 *  \brief  Compute aggregated cost volume
 *   Compute aggregated cost at one point
 *   diagonaly from upper left
 *
 *  \param pixelCost cost value 
 *  \param row row position
 *  \param col col position
 *  \param disp disp position
 *  \param invalid_value value representing invalid cost
 *  \param nb_rows row number of cost volume
 *  \param nb_cols column number of cost volume
 *  \param nb_disps disparity number of cost volume
 *  \param P1 penalty P1 from sgm equation
 *  \param P2 penalty P2 from sgm equation
 *  \param direction array of struct containing coordinates of previous point
 *  \param buff2 buffer containing previous aggregated cost (previous line)
 *  \param buff_disp_2 buffer containing previous aggregated point
 *  \param min_disp2 value of minimum cost at the previous pixel
 *  \param pixel_2 value of previous disp on the buffer
 *  \param current_class current class
 *  \param buff_class2 buffer containing previous class (previous point
 *  \param reset2 value of coefficient to multiply history
 *  \return cost aggregated point
 */

template<typename T> 
T aggregatedCostFromTopLeft2(T pixelCost, int row, int col, int disp, T invalid_value,
    int nb_rows, int nb_cols, int nb_disps, T P1, T P2, Direction direction,
    T * buff2, T * buff_disp_2, T * min_disp2, T * pixel_2, float current_class, float * buff_class2, float & reset2);

/*!
 *  \brief  Compute aggregated cost volume
 *   Compute aggregated cost at one point
 *   diagonaly from upper right
 *
 *  \param pixelCost cost value 
 *  \param row row position
 *  \param col col position
 *  \param disp disp position
 *  \param invalid_value value representing invalid cost
 *  \param nb_rows row number of cost volume
 *  \param nb_cols column number of cost volume
 *  \param nb_disps disparity number of cost volume
 *  \param P1 penalty P1 from sgm equation
 *  \param P2 penalty P2 from sgm equation
 *  \param direction array of struct containing coordinates of previous point
 *  \param buff3 buffer containing previous aggregated cost (previous line)
 *  \param min_disp3 value of minimum cost at the previous pixel
 *  \param current_class current class
 *  \param buff_class3 buffer containing previous class (previous point
 *  \param reset3 value of coefficient to multiply history
 *  \return cost aggregated point
 */

template<typename T> 
T aggregatedCostFromTopLeft3(T pixelCost, int row, int col, int disp, T invalid_value,
    int nb_rows, int nb_cols, int nb_disps, T P1, T P2, Direction direction,
    T * buff3, T * min_disp3, float current_class, float * buff_class3, float & reset3);

/*!
 *  \brief  Compute aggregated cost volume
 *   Compute aggregated cost at one point
 *   from right to left direction
 *
 *  \param pixelCost cost value 
 *  \param row row position
 *  \param col col position
 *  \param disp disp position
 *  \param invalid_value value representing invalid cost
 *  \param nb_rows row number of cost volume
 *  \param nb_cols column number of cost volume
 *  \param nb_disps disparity number of cost volume
 *  \param P1 penalty P1 from sgm equation
 *  \param P2 penalty P2 from sgm equation
 *  \param direction array of struct containing coordinates of previous pointnt
 *  \param buff4 buffer containing previois aggregated cost (previous point)
 *  \param min_disp4 value of minimum cost at the previous pixel
 *  \param pixel_4 value of previous disp on the buffer
 *  \param current_class current class
 *  \param buff_class4 buffer containing previous class (previous point
 *  \param reset4 value of coefficient to multiply history
 *  \return cost aggregated point
 */

template<typename T> 
T aggregatedCostFromBottomRight4(T pixelCost, int row, int col, int disp, T invalid_value,
    int nb_rows, int nb_cols, int nb_disps, T P1, T P2, Direction direction,
    T * buff4, T * min_disp4, T * pixel_4, float current_class, float buff_class4, float & reset4);

/*!
 *  \brief  Compute aggregated cost volume
 *   Compute aggregated cost at one point
 *   from down to up direction
 *
 *  \param pixelCost cost value 
 *  \param row row position
 *  \param col col position
 *  \param disp disp position
 *  \param invalid_value value representing invalid cost
 *  \param nb_rows row number of cost volume
 *  \param nb_cols column number of cost volume
 *  \param nb_disps disparity number of cost volume
 *  \param P1 penalty P1 from sgm equation
 *  \param P2 penalty P2 from sgm equation
 *  \param direction array of struct containing coordinates of previous point
 *  \param buff5 buffer containing previois aggregated cost (previous line)
 *  \param min_disp5 value of minimum cost at the previous pixel
 *  \param pixel_5 value of previous disp on the buffer
 *  \param current_class current class
 *  \param buff_class5 buffer containing previous class (previous point
 *  \param reset5 value of coefficient to multiply history
 *  \return cost aggregated point
 */

template<typename T> 
T aggregatedCostFromBottomRight5(T pixelCost, int row, int col, int disp, T invalid_value,
    int nb_rows, int nb_cols, int nb_disps, T P1, T P2, Direction direction,
    T * buff5, T * min_disp5, T * pixel_5, float current_class, float * buff_class5, float & reset5);

/*!
 *  \brief  Compute aggregated cost volume
 *   Compute aggregated cost at one point
 *   diagonaly from lower left
 *
 *  \param pixelCost cost value 
 *  \param row row position
 *  \param col col position
 *  \param disp disp position
 *  \param invalid_value value representing invalid cost
 *  \param nb_rows row number of cost volume
 *  \param nb_cols column number of cost volume
 *  \param nb_disps disparity number of cost volume
 *  \param P1 penalty P1 from sgm equation
 *  \param P2 penalty P2 from sgm equation
 *  \param direction array of struct containing coordinates of previous point
 *  \param buff6 buffer containing previous aggregated cost (previous line)
 *  \param buff_disp_6 buffer containing previous aggregated point
 *  \param min_disp6 value of minimum cost at the previous pixel
 *  \param pixel_6 value of previous disp on the buffer
 *  \param current_class current class
 *  \param buff_class6 buffer containing previous class (previous point
 *  \param reset6 value of coefficient to multiply history
 *  \return cost aggregated point
 */

template<typename T>
T aggregatedCostFromBottomRight6(T pixelCost, int row, int col, int disp, T invalid_value,
    int nb_rows, int nb_cols, int nb_disps, T P1, T P2, Direction direction,
    T * buff6, T * buff_disp_6, T * min_disp6, T * pixel_6, float current_class, float * buff_class6, float & reset6);

/*!
 *  \brief  Compute aggregated cost volume
 *   Compute aggregated cost at one point
 *   diagonaly from lower right
 *
 *  \param pixelCost cost value 
 *  \param row row position
 *  \param col col position
 *  \param disp disp position
 * \param invalid_value value representing invalid cost
 *  \param nb_rows row number of cost volume
 *  \param nb_cols column number of cost volume
 *  \param nb_disps disparity number of cost volume
 *  \param P1 penalty P1 from sgm equation
 *  \param P2 penalty P2 from sgm equation
 *  \param direction array of struct containing coordinates of previous point
 *  \param buff7 buffer containing previous aggregated cost (previous line)
 *  \param min_disp7 value of minimum cost at the previous pixel
 *  \param current_class current class
 *  \param buff_class7 buffer containing previous class (previous point
 *  \param reset7 value of coefficient to multiply history
 *  \return cost aggregated point
 */

template<typename T>
T aggregatedCostFromBottomRight7(T pixelCost, int row, int col, int disp, T invalid_value,
    int nb_rows, int nb_cols, int nb_disps, T P1, T P2, Direction direction,T * buff7, T * min_disp7,
    float current_class, float * buff_class7, float & reset7);

/*!
 *  \brief  Apply penalties
 *
 *  \param penalty struct containing penalty P1,P2 for the point (row,col)
 *  \param p1 p1 value
 *  \param p2 p2 value
 */
template<typename T>
void apply_penalty(Penalty<T> *penalty, T p1, T p2);

/*!
 *  \brief  Assign Directions
 *
 *  \param directions_in directions matrix
 *  \param dirs directions
 */
void  assignDirections(int* directions_in, Direction* dirs);