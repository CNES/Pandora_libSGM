/*
 * Copyright (c) 2020 Centre National d'Etudes Spatiales (CNES).
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

#include <cstdlib>
#include <iostream>
#include <limits>
#include <algorithm>
#include <array>
#include <functional>
#include <omp.h>
#include "sgm.hpp"


template<typename T> 
CostVolumes sgm(T * cv_in, T* p1_in, T* p2_in, int* directions_in, int nb_rows, int nb_cols, int nb_disps,
    T invalid_value, bool cost_paths, bool overcounting)

{
    int nb_dir = 8;
    CostVolumes cvs;
    //Allocate final cost volume
    cvs.cost_volume = new float[nb_rows*nb_cols*nb_disps]();

    // Allocate costs
    int nb_values=1;
    if (cost_paths){
    nb_values = nb_rows*nb_cols*nb_dir;
    }
    cvs.cost_volume_min = new int[nb_values]();


    //Direction (x,y) indicating previous pixel for each path
    Direction direction[8] = {{0}};
    assignDirections(directions_in, direction);

    //Penalties
    T P1, P2;
    Penalty<T> penalty[8]={{P1, P2}, {P1, P2}, {P1, P2}, {P1, P2}, {P1, P2}, {P1, P2}, {P1, P2}, {P1, P2}};


    // Census Cost at pixel
    T pixelCost;

    /*
    Two passes: the 1st from top left , the 2nd from bottom right
    Each one aggregate 4 paths
    First pass
        0 : left-> right
        1 : up -> down
        2 : diagonal from top left
        3 : diagonal from top right
    Second pass
        4 : right->left
        5 : down->up
        6 : diagonal from lower left
        7 : diagonal from lower right
    */

    /* ---------------*/
    /* --First pass-- */
    /* ---------------*/

    //temporary buffers and variables for storing aggregated cost
    T * buff0 = new T[nb_disps]();
    T * buff1 = new T[nb_cols*nb_disps]();
    T * buff2 = new T[nb_cols*nb_disps]();
    T * buff3 = new T[nb_cols*nb_disps]();
    T pixel_0, pixel_1, pixel_2;
    T * buff_disp_2 = new T[nb_disps]();
    T  min_disp0, min_disp1, min_disp2, min_disp3;
    float costAggr;
    // temporary variables for storing outputs
    float s0, s1, s2, s3;
    // temporary variables for storing current position of minimums
    int  pos0, pos1, pos2, pos3;
    // temporary variables for storing current minimums
    float min0, min1, min2, min3;

    for (int row = 0; row < nb_rows; row++)
    {
        for (int col = 0; col < nb_cols; col++)
        {
            //left-> right
            apply_penalty(&penalty[0], p1_in[0+col*nb_dir+row*nb_dir*nb_cols], p2_in[0+col*nb_dir+row*nb_dir*nb_cols]);
            //up -> down
            apply_penalty(&penalty[1], p1_in[1+col*nb_dir+row*nb_dir*nb_cols], p2_in[1+col*nb_dir+row*nb_dir*nb_cols]);
            //diagonal from top left
            apply_penalty(&penalty[2], p1_in[2+col*nb_dir+row*nb_dir*nb_cols], p2_in[2+col*nb_dir+row*nb_dir*nb_cols]);
            //diagonal from top right
            apply_penalty(&penalty[3], p1_in[3+col*nb_dir+row*nb_dir*nb_cols], p2_in[3+col*nb_dir+row*nb_dir*nb_cols]);

            // initialize minimums
            min0 = std::numeric_limits<float>::max() , min1 = std::numeric_limits<float>::max(), min2 = std::numeric_limits<float>::max(), min3 = std::numeric_limits<float>::max();
            pos0 = 0, pos1 = 0, pos2 = 0, pos3 = 0;
            for (int disp = 0; disp < nb_disps; disp++)
            {
                costAggr = 0;
                s0= 0, s1= 0, s2= 0, s3 = 0;
                pixelCost = cv_in[disp+col*nb_disps+row*nb_disps*nb_cols];
                //left-> right
                s0 = aggregatedCostFromTopLeft0(pixelCost, row, col, disp, invalid_value, nb_rows,
                 nb_cols, nb_disps, penalty[0].P1, penalty[0].P2, direction[0], buff0, &min_disp0, &pixel_0);
                costAggr += s0;
                //up -> down
                s1 = aggregatedCostFromTopLeft1(pixelCost, row, col, disp, invalid_value, nb_rows,
                 nb_cols, nb_disps, penalty[1].P1, penalty[1].P2, direction[1], buff1, &min_disp1, &pixel_1);
                costAggr += s1;
                //diagonal from top left
                s2 = aggregatedCostFromTopLeft2(pixelCost, row, col, disp, invalid_value, nb_rows,
                 nb_cols, nb_disps, penalty[2].P1, penalty[2].P2, direction[2], buff2, buff_disp_2, &min_disp2, &pixel_2);
                costAggr += s2;
                //diagonal from top right
                s3 = aggregatedCostFromTopLeft3(pixelCost, row, col, disp, invalid_value, nb_rows,
                 nb_cols, nb_disps, penalty[3].P1, penalty[3].P2, direction[3], buff3, &min_disp3);
                costAggr += s3;

                cvs.cost_volume[disp + col*nb_disps + row*nb_disps*nb_cols] += costAggr;

                if (cost_paths){
                    std::tie(min0, pos0) = update_minimum(min0, s0, pos0, disp);
                    std::tie(min1, pos1) = update_minimum(min1, s1, pos1, disp);
                    std::tie(min2, pos2) = update_minimum(min2, s2, pos2, disp);
                    std::tie(min3, pos3) = update_minimum(min3, s3, pos3, disp);
                }
            }
            if (cost_paths){
                cvs.cost_volume_min[0 + col*nb_dir + row*nb_dir*nb_cols] = pos0;
                cvs.cost_volume_min[1 + col*nb_dir + row*nb_dir*nb_cols] = pos1;
                cvs.cost_volume_min[2 + col*nb_dir + row*nb_dir*nb_cols] = pos2;
                cvs.cost_volume_min[3 + col*nb_dir + row*nb_dir*nb_cols] = pos3;
            }
        }
    }
    //delete buffers
    delete [] buff0;
    delete [] buff1;
    delete [] buff2;
    delete [] buff3;
    delete [] buff_disp_2;
    
    /* ------------------*/
    /* -- Second pass -- */
    /* ------------------*/

    //temporary buffers and variables for storing aggregated cost
    T * buff4 = new T[nb_disps]();
    T * buff5 = new T[nb_cols*nb_disps]();
    T * buff6 = new T[nb_cols*nb_disps]();
    T * buff7 = new T[nb_cols*nb_disps]();
    T pixel_4, pixel_5, pixel_6;
    T * buff_disp_6 = new T[nb_disps]();
    T  min_disp4, min_disp5, min_disp6, min_disp7;
    // temporary variables for storing outputs
    float s4, s5, s6, s7;
    // temporary variables for storing current position of minimums
    int  pos4, pos5, pos6, pos7;
    // temporary variables for storing current minimums
    float min4, min5, min6, min7;
    int overcounting_factor;

    if (overcounting)
    {
	// Factor to correct the overcounting: number of directions [8] - 1
	overcounting_factor = 7;
    }
    else
    {
	// Factor to not correct the over-counting
	overcounting_factor = 0;
    }

    for (int row = nb_rows-1; row >= 0; row --)
    {
        for (int col = nb_cols-1; col >= 0; col --)
        {
            //right->left
            apply_penalty(&penalty[4], p1_in[4+col*nb_dir+row*nb_dir*nb_cols], p2_in[4+col*nb_dir+row*nb_dir*nb_cols]);
            //up -> down
            apply_penalty(&penalty[5], p1_in[5+col*nb_dir+row*nb_dir*nb_cols], p2_in[5+col*nb_dir+row*nb_dir*nb_cols]);
            //diagonal from lower left
            apply_penalty(&penalty[6], p1_in[6+col*nb_dir+row*nb_dir*nb_cols], p2_in[6+col*nb_dir+row*nb_dir*nb_cols]);
            //diagonal from lower right
            apply_penalty(&penalty[7], p1_in[7+col*nb_dir+row*nb_dir*nb_cols], p2_in[7+col*nb_dir+row*nb_dir*nb_cols]);

            // initialize minimums
            min4 = std::numeric_limits<float>::max(), min5 = std::numeric_limits<float>::max(), min6 = std::numeric_limits<float>::max(), min7 = std::numeric_limits<float>::max();
            pos4 = 0, pos5 = 0, pos6 = 0, pos7 = 0;

            for (int disp = 0; disp < nb_disps; disp++)
            {
                costAggr = 0;
                s4= 0, s5= 0, s6= 0, s7 = 0;

                pixelCost = cv_in[disp+col*nb_disps+row*nb_disps*nb_cols];
                //right->left
                s4 = aggregatedCostFromBottomRight4(pixelCost, row, col, disp, invalid_value, nb_rows,
                 nb_cols, nb_disps, penalty[4].P1, penalty[4].P2, direction[4], buff4, &min_disp4, &pixel_4);
                costAggr += s4;
                //up -> down
                s5 = aggregatedCostFromBottomRight5(pixelCost, row, col, disp, invalid_value, nb_rows,
                 nb_cols, nb_disps, penalty[5].P1, penalty[5].P2, direction[5], buff5, &min_disp5, &pixel_5);
                costAggr += s5;
                //diagonal from lower left
                s6 = aggregatedCostFromBottomRight6(pixelCost, row, col, disp, invalid_value, nb_rows,
                 nb_cols, nb_disps, penalty[6].P1, penalty[6].P2, direction[6], buff6, buff_disp_6, &min_disp6, &pixel_6);
                costAggr += s6;
                //diagonal from lower right
                s7 = aggregatedCostFromBottomRight7(pixelCost, row, col, disp, invalid_value, nb_rows,
                 nb_cols, nb_disps, penalty[7].P1, penalty[7].P2, direction[7], buff7, &min_disp7);
                costAggr += s7;

                cvs.cost_volume[disp + col*nb_disps + row*nb_disps*nb_cols] += costAggr;
                // Correction of the over-counting by removing (overcounting_factor * pixel cost volume) to the aggregated cost volume
                cvs.cost_volume[disp + col*nb_disps + row*nb_disps*nb_cols] -= overcounting_factor * pixelCost;

                if (cost_paths){
                    std::tie(min4, pos4) = update_minimum(min4, s4, pos4, disp);
                    std::tie(min5, pos5) = update_minimum(min5, s5, pos5, disp);
                    std::tie(min6, pos6) = update_minimum(min6, s6, pos6, disp);
                    std::tie(min7, pos7) = update_minimum(min7, s7, pos7, disp);
                }

            }
            if (cost_paths){
                cvs.cost_volume_min[4 + col*nb_dir + row*nb_dir*nb_cols] = pos4;
                cvs.cost_volume_min[5 + col*nb_dir + row*nb_dir*nb_cols] = pos5;
                cvs.cost_volume_min[6 + col*nb_dir + row*nb_dir*nb_cols] = pos6;
                cvs.cost_volume_min[7 + col*nb_dir + row*nb_dir*nb_cols] = pos7;
            }

        }
    }

    delete [] buff4;
    delete [] buff5;
    delete [] buff6;
    delete [] buff7;
    delete [] buff_disp_6;

    return cvs;
}

std::pair<float, int> update_minimum(float current_min, float value, int current_disp, int disp)
{
    if (current_min < value){
        return std::make_pair(current_min, current_disp);
    }
    else {
        return std::make_pair(value, disp);
    }
}

template<typename T> 
T aggregatedCostFromTopLeft0(T pixelCost, int row, int col, int disp, T invalid_value,
    int nb_rows, int nb_cols, int nb_disps, T P1, T P2, Direction direction,
    T * buff0, T * min_disp0, T * pixel_0)
{

    //Pixel cost at the point (row,col,disp)
    /* ------------------*/
    /* -- Left->right -- */
    /* ------------------*/

    T costAggr0 = pixelCost;

    if (! (row-direction.drow<0 || row-direction.drow > (nb_rows -1) || col-direction.dcol < 0 || col-direction.dcol > (nb_cols-1) ) )
    {   /*Minimum cost at previous point for each disparity
        Compute value just at the first disparity for a given point (row,col)
        */
        if (disp == 0)
        {
            *min_disp0 = *std::min_element(&buff0[disp], &buff0[nb_disps]);
        }
        /*If pixelCost is equal to invalid value, aggregated cost must be equal to invalid value.
        So, it's useless to compute the minimum on tmp1,tmp2,tmp3,tmp4
        */
        if(pixelCost != invalid_value)
        {
            //Previous cost
            const T tmp1 = buff0[disp];
            //Previous cost at disparity-1
            const T tmp2 = (disp>0) ? *pixel_0+P1 : std::numeric_limits<T>::max();
            //Previous cost at disparity+1
            const T tmp3 = (disp<nb_disps-1) ? buff0[disp+1]+P1 : std::numeric_limits<T>::max();

            const T tmp4 = (*min_disp0)+P2;

            //Minimum path cost
            costAggr0 += std::min({tmp1, tmp2, tmp3, tmp4}) - (*min_disp0);
        }
    }
    //sava buffer data before writing new value in buffer
    *pixel_0 = buff0[disp];
    //Save data in buffer
    buff0[disp] = costAggr0;

    return costAggr0;
}

template<typename T> 
T aggregatedCostFromTopLeft1(T pixelCost, int row, int col, int disp, T invalid_value,
    int nb_rows, int nb_cols, int nb_disps, T P1, T P2, Direction direction,
    T * buff1, T * min_disp1, T * pixel_1)
{
    /* ------------------*/
    /* --  Up->Down   -- */
    /* ------------------*/
    T costAggr1 = pixelCost;
    //Border check
    if (! (row-direction.drow<0 || row-direction.drow > (nb_rows -1) || col-direction.dcol < 0 || col-direction.dcol > (nb_cols-1) ) )
    {
        /*Minimum cost at previous point for each disparity
        Compute value just at the first disparity for a given point (row,col)
        */
        if (disp == 0)
        {
            *min_disp1 = *std::min_element(&buff1[disp+col*nb_disps], &buff1[nb_disps + col*nb_disps]);
        }
        /*If pixelCost is equal to invalid value, aggregated cost must be equal to invalid value.
        So, it's useless to compute the minimum on tmp1,tmp2,tmp3,tmp4
        */
        if(pixelCost != invalid_value)
        {
            //Previous cost
            const T tmp1 = buff1[disp+col*nb_disps];
            //Previous cost at disparity-1
            const T tmp2 = (disp>0) ? *pixel_1+P1 : std::numeric_limits<T>::max();
            //Previous cost at disparity+1
            const T tmp3 = (disp<nb_disps-1) ?
            buff1[disp+1+col*nb_disps]+P1 : std::numeric_limits<T>::max();
            //Minimum cost at previous point
            const T tmp4 = (*min_disp1)+P2;
            //Minimum path cost
            costAggr1 += std::min({tmp1, tmp2, tmp3, tmp4}) - (*min_disp1);
        }
    }
    //sava buffer data before writing new value in buffer
    *pixel_1 = buff1[disp+col*nb_disps];
    //Save data in buffer
    buff1[disp+col*nb_disps] = costAggr1;

    return costAggr1;

}

template<typename T> 
T aggregatedCostFromTopLeft2(T pixelCost, int row, int col, int disp, T invalid_value,
    int nb_rows, int nb_cols, int nb_disps, T P1, T P2, Direction direction,
    T * buff2, T * buff_disp_2, T * min_disp2, T * pixel_2)
{

    /* -----------------------------*/
    /* --Diagonal from upper left-- */
    /* -----------------------------*/
    uint8_t costAggr2 = pixelCost;
    //Border check
    if (! (row-direction.drow<0 || row-direction.drow > (nb_rows -1) || col-direction.dcol < 0 || col-direction.dcol > (nb_cols-1) ) )
    {
        /*Minimum cost at previous point for each disparity
        Compute value just at the first disparity for a given point (row,col)
        */
        if (disp == 0)
        {
            *min_disp2 = *std::min_element(&buff_disp_2[disp], &buff_disp_2[nb_disps]);
        }
        /*If pixelCost is equal to invalid value, aggregated cost must be equal to invalid value.
        So, it's useless to compute the minimum on tmp1,tmp2,tmp3,tmp4
        */
        if(pixelCost != invalid_value)
        {
            //Previous cost
            const T tmp1 = buff_disp_2[disp];
            //Previous cost at disparity-1
            const T tmp2 = (disp>0) ?*pixel_2+P1 : std::numeric_limits<T>::max();
            //Previous cost at disparity+1
            const T tmp3 = (disp<nb_disps-1) ? buff_disp_2[disp+1]+P1 : std::numeric_limits<T>::max();
            //Minimum cost at previous point
            const T tmp4 = (*min_disp2)+P2;
            //Minimum path cost
            costAggr2 += std::min({tmp1, tmp2, tmp3, tmp4}) - (*min_disp2);
        }
    }
    //sava buffer data before writing new value in buffer_disp
    *pixel_2 = buff_disp_2[disp];
    //sava buffer data before writing new value in buffer_disp
    buff_disp_2[disp] = buff2[disp+col*nb_disps];
    //Save data in buffer
    buff2[disp+col*nb_disps] = costAggr2;

    return costAggr2;
}

template<typename T> 
T aggregatedCostFromTopLeft3(T pixelCost, int row, int col, int disp, T invalid_value,
    int nb_rows, int nb_cols, int nb_disps, T P1, T P2, Direction direction,
    T * buff3, T * min_disp3)
{

    T costAggr3 = pixelCost;
    //Border check
    if (! (row-direction.drow<0 || row-direction.drow > (nb_rows -1) || col-direction.dcol < 0 || col-direction.dcol > (nb_cols-1) ) )
    {
        /*Minimum cost at previous point for each disparity
        Compute value just at the first disparity for a given point (row,col)
        */
        if (disp == 0)
        {
            *min_disp3 = *std::min_element(&buff3[disp+(col-direction.dcol)*nb_disps], &buff3[nb_disps + (col-direction.dcol)*nb_disps]);
        }
        /*If pixelCost is equal to invalid value, aggregated cost must be equal to invalid value.
        So, it's useless to compute the minimum on tmp1,tmp2,tmp3,tmp4
        */
        if(pixelCost != invalid_value)
        {
            //Previous cost
            const uint8_t tmp1 = buff3[disp+(col-direction.dcol)*nb_disps];
            //Previous cost at disparity-1
            const uint8_t tmp2 = (disp>0) ? buff3[disp-1+(col-direction.dcol)*nb_disps]+P1 : std::numeric_limits<T>::max();
            //Previous cost at disparity+1
            const uint8_t tmp3 = (disp<nb_disps-1) ? buff3[disp+1+(col-direction.dcol)*nb_disps]+P1 : std::numeric_limits<T>::max();
            //Minimum cost at previous point
            const uint8_t tmp4 = (*min_disp3)+P2;
            //Minimum path cost
            costAggr3 += std::min({tmp1, tmp2, tmp3, tmp4}) - (*min_disp3);
        }
    }
    //Save data in buffer
    buff3[disp+col*nb_disps] = costAggr3;

    return costAggr3;
}

template<typename T> 
T aggregatedCostFromBottomRight4(T pixelCost, int row, int col, int disp, T invalid_value,
    int nb_rows, int nb_cols, int nb_disps, T P1, T P2, Direction direction,
    T * buff4, T * min_disp4, T * pixel_4)
{

    /* ------------------*/
    /* -- right->left -- */
    /* ------------------*/
    T costAggr4 = pixelCost;
    //Border check
    if (! (row-direction.drow<0 || row-direction.drow > (nb_rows -1) || col-direction.dcol < 0 || col-direction.dcol > (nb_cols-1) ) )
    {
        /*Minimum cost at previous point for each disparity
        Compute value just at the first disparity for a given point (row,col)
        */
        if (disp == 0)
        {
            *min_disp4 = *std::min_element(&buff4[disp], &buff4[nb_disps]);
        }
        /*If pixelCost is equal to invalid value, aggregated cost must be equal to invalid value.
        So, it's useless to compute the minimum on tmp1,tmp2,tmp3,tmp4
        */
        if(pixelCost != invalid_value)
        {
            //Previous cost
            const uint8_t tmp1 = buff4[disp];
            //Previous cost at disparity-1
            const uint8_t tmp2 = (disp>0) ? *pixel_4+P1 : std::numeric_limits<T>::max();
            //Previous cost at disparity+1
            const uint8_t tmp3 = (disp<nb_disps-1) ? buff4[disp+1]+P1 : std::numeric_limits<T>::max();
            //Minimum cost at previous point
            const uint8_t tmp4 = (*min_disp4)+P2;
            //Minimum path cost
            costAggr4 += std::min({tmp1, tmp2, tmp3, tmp4}) - (*min_disp4);
        }
    }
    //sava buffer data before writing new value in buffer
    *pixel_4 = buff4[disp];
    //Save data in buffer
    buff4[disp] = costAggr4;

    return costAggr4;
}

template<typename T> 
T aggregatedCostFromBottomRight5(T pixelCost, int row, int col, int disp, T invalid_value,
    int nb_rows, int nb_cols, int nb_disps, T P1, T P2, Direction direction,
    T * buff5, T * min_disp5, T * pixel_5)
{

    /* ------------------*/
    /* --  Down->up   -- */
    /* ------------------*/
    T costAggr5 = pixelCost;
    //Border check
    if (! (row-direction.drow<0 || row-direction.drow > (nb_rows -1) || col-direction.dcol < 0 || col-direction.dcol > (nb_cols-1) ) )
    {
        /*Minimum cost at previous point for each disparity
        Compute value just at the first disparity for a given point (row,col)
        */
        if (disp == 0)
        {
            *min_disp5 = *std::min_element(&buff5[disp+col*nb_disps], &buff5[nb_disps + col*nb_disps]);
        }
        /*If pixelCost is equal to invalid value, aggregated cost must be equal to invalid value.
        So, it's useless to compute the minimum on tmp1,tmp2,tmp3,tmp4
        */
        if(pixelCost != invalid_value)
        {
            //Previous cost
            const T tmp1 = buff5[disp+col*nb_disps];
            //Previous cost at disparity-1
            const T tmp2 = (disp>0) ? *pixel_5+P1 : std::numeric_limits<T>::max();
            //Previous cost at disparity+1
            const T tmp3 = (disp<nb_disps-1) ? buff5[disp+1+col*nb_disps]+P1 : std::numeric_limits<T>::max();
            //Minimum cost at previous point
            const T tmp4 = (*min_disp5)+P2;
            //Minimum path cost
            costAggr5 += std::min({tmp1, tmp2, tmp3, tmp4}) - (*min_disp5);
        }

    }
    //sava buffer data before writing new value in buffer
    *pixel_5 = buff5[disp+col*nb_disps];
    //Save data in buffer
    buff5[disp+col*nb_disps] = costAggr5;

    return costAggr5;
}

template<typename T>
T aggregatedCostFromBottomRight6(T pixelCost, int row, int col, int disp, T invalid_value,
    int nb_rows, int nb_cols, int nb_disps, T P1, T P2, Direction direction,
    T * buff6, T * buff_disp_6, T * min_disp6, T * pixel_6)
{

    /* -----------------------------*/
    /* --diagonal from lower left-- */
    /* -----------------------------*/

    T costAggr6 = pixelCost;
    //Border check
    if (! (row-direction.drow<0 || row-direction.drow > (nb_rows -1) || col-direction.dcol < 0 || col-direction.dcol > (nb_cols-1) ) )
    {
        /*Minimum cost at previous point for each disparity
        Compute value just at the first disparity for a given point (row,col)
        */
        if (disp == 0)
        {
            *min_disp6 = *std::min_element(&buff_disp_6[disp], &buff_disp_6[nb_disps]);
        }
        /*If pixelCost is equal to invalid value, aggregated cost must be equal to invalid value.
        So, it's useless to compute the minimum on tmp1,tmp2,tmp3,tmp4
        */
        if(pixelCost != invalid_value)
        {
            //Previous pixels
            const T tmp1 = buff_disp_6[disp];
            //Previous cost at disparity-1
            const T tmp2 = (disp>0) ?*pixel_6+P1 : std::numeric_limits<T>::max();
            //Previous cost at disparity+1
            const T tmp3 = (disp<nb_disps-1) ? buff_disp_6[disp+1]+P1 : std::numeric_limits<T>::max();
            //Minimum cost at previous point
            const T tmp4 = (*min_disp6)+P2;
            //Minimum path cost
            costAggr6 += std::min({tmp1, tmp2, tmp3, tmp4}) - (*min_disp6);
        }
    }
    //sava buffer data before writing new value in buffer_disp
    *pixel_6 = buff_disp_6[disp];
    //sava buffer data before writing new value in buffer_disp
    buff_disp_6[disp] = buff6[disp+col*nb_disps];
    //Save data in buffer
    buff6[disp+col*nb_disps] = costAggr6;

    return costAggr6;
}

template<typename T>
T aggregatedCostFromBottomRight7(T pixelCost, int row, int col, int disp, T invalid_value,
    int nb_rows, int nb_cols, int nb_disps, T P1, T P2, Direction direction,
    T * buff7, T * min_disp7)
{

    /* ------------------------------*/
    /* --diagonal from lower right-- */
    /* ------------------------------*/
    T costAggr7 = pixelCost;
    //Border
    if (! (row-direction.drow<0 || row-direction.drow > (nb_rows -1) || col-direction.dcol < 0 || col-direction.dcol > (nb_cols-1) ) )
    {
        /*Minimum cost at previous point for each disparity
        Compute value just at the first disparity for a given point (row,col)
        */
        if (disp == 0)
        {
            *min_disp7 = *std::min_element(&buff7[disp+(col-direction.dcol)*nb_disps], &buff7[nb_disps + (col-direction.dcol)*nb_disps]);
        }
        /*If pixelCost is equal to invalid value, aggregated cost must be equal to invalid value.
        So, it's useless to compute the minimum on tmp1,tmp2,tmp3,tmp4
        */
        if(pixelCost != invalid_value)
        {
            //Previous pixels
            const T tmp1 = buff7[disp+(col-direction.dcol)*nb_disps];
            //Previous cost at disparity-1
            const T tmp2 = (disp>0) ? buff7[disp-1+(col-direction.dcol)*nb_disps]+P1 : std::numeric_limits<T>::max();
            //Previous cost at disparity+1
            const T tmp3 = (disp<nb_disps-1) ? buff7[disp+1+(col-direction.dcol)*nb_disps]+P1 : std::numeric_limits<T>::max();
            //Minimum cost at previous point
            const T tmp4 = (*min_disp7)+P2;
            //Minimum path cost
            costAggr7 += std::min({tmp1, tmp2, tmp3, tmp4}) - (*min_disp7);
        }

    }
    //sava buffer data before writing new value in buffer
    buff7[disp+col*nb_disps] = costAggr7;

    return costAggr7;
}

template<typename T>
void apply_penalty(Penalty<T> *penalty, T p1, T p2)
{
    penalty->P1 = p1 ;
    penalty->P2 = p2 ;
}

void assignDirections(int* directions_in, Direction* dirs)
{
    for (int i=0; i<8; i++)
    {
        dirs[i].drow = directions_in[2 * i];
        dirs[i].dcol = directions_in[2 * i + 1];
    }
}

/* Explicitly instantiate all the templates needed to use libSGM as an external lib */
template CostVolumes sgm(uint8_t * cv_in, uint8_t * p1_in, uint8_t * p2_in, int* directions_in,  int nb_rows,
 int nb_cols, int nb_disps, uint8_t invalid_value, bool cost_paths, bool overcounting);
template CostVolumes sgm(float * cv_in, float * p1_in, float * p2_in, int* directions_in, int nb_rows, int nb_cols,
 int nb_disps, float invalid_value, bool cost_paths, bool overcounting);