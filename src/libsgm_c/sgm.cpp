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

#include <cstdlib>
#include <iostream>
#include <limits>
#include <algorithm>
#include <array>
#include <functional>
#include "sgm.hpp"

// Forward declaration required because update_minimum is used inside template sgm
static inline void update_minimum(float& current_min, float value, int& current_disp, int disp){
  if (current_min < value){
    current_min = value;
    current_disp = disp;
  }
}


template <typename T, int colStep>
inline T costCompute(T pixelCost, int col, int disp, T invalid_value,
                             int nb_disps, Penalty<T> P, T *buff, T min_disp, T pixel)
{
  T costAggr1 = pixelCost;
  if (pixelCost != invalid_value)
  {
    // Previous cost
    const T tmp1 = buff[disp + colStep * col * nb_disps];
    // Previous cost at disparity-1
    // For disp == 0, pixel = tmp1, so tmp2 = tmp1 + P.P1 => no problem as we take the min
    const T tmp2 = pixel + P.P1;
    // Previous cost at disparity+1
    const T tmp3 = (disp < nb_disps - 1) ? buff[disp + 1 + colStep * col * nb_disps] + P.P1 : std::numeric_limits<T>::max();
    // Minimum cost at previous point
    const T tmp4 = (min_disp) + P.P2;
    // Minimum path cost
    costAggr1 += (std::min({tmp1, tmp2, tmp3, tmp4}) - (min_disp));
  }
  return costAggr1;
}

template <typename T, typename Tout>
void sgm(T *cv_in, T *p1_in, T *p2_in, int *directions_in, unsigned long int nb_rows, unsigned long int nb_cols,
                      unsigned int nb_disps, T invalid_value, float *segmentation, bool cost_paths, bool overcounting, CostVolumes<Tout> &cvs)

{
  int nb_dir = 8;

  // Direction (x,y) indicating previous pixel for each path
  Direction direction[8] = {};
  assignDirections(directions_in, direction);

  // Penalties (Fix: Initialize P1 and P2 to avoid warnings)
  Penalty<T> penalty[8] = {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}};

  // Reset cvs.cost_volume to 0
  std::fill(cvs.cost_volume, cvs.cost_volume + nb_rows * nb_cols * nb_disps, 0);

  // Reset cvs.cost_volume_min to 0
  if(cost_paths){
      std::fill(cvs.cost_volume_min, cvs.cost_volume_min + nb_rows * nb_cols * 8, 0);
  }

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

  // temporary buffers and variables for storing aggregated cost
  T *buff0 = new T[nb_disps]();
  T *buff1 = new T[nb_cols * nb_disps]();
  T *buff2 = new T[nb_cols * nb_disps]();
  T *buff3 = new T[nb_cols * nb_disps]();
  T *buff_disp_2 = new T[nb_disps]();

  // temporary variables for storing current position of minimums
  int pos0, pos1, pos2, pos3;
  // temporary variables for storing current minimums
  float min0, min1, min2, min3;
  // temporary buffers for storing current classification
  float current_class;
  // buff_class0 stores the class from the left of current pixel (previous seen pixel)
  float buff_class0 = segmentation[0];
  // tmp_buff_class1_2_3 stores the previous classes of pixels seen in the line
  float *tmp_buff_class1_2_3 = new float[nb_cols]();
  // buff_class1_2_3 stores the previous line of classes
  float *buff_class1_2_3 = new float[nb_cols]();

  for (size_t row = 0; row < nb_rows; row++)
  {
    for (size_t col = 0; col < nb_cols; col++)
    {
      // left-> right
      apply_penalty(&penalty[0], p1_in[0 + col * nb_dir + row * nb_dir * nb_cols], p2_in[0 + col * nb_dir + row * nb_dir * nb_cols]);
      // up -> down
      apply_penalty(&penalty[1], p1_in[1 + col * nb_dir + row * nb_dir * nb_cols], p2_in[1 + col * nb_dir + row * nb_dir * nb_cols]);
      // diagonal from top left
      apply_penalty(&penalty[2], p1_in[2 + col * nb_dir + row * nb_dir * nb_cols], p2_in[2 + col * nb_dir + row * nb_dir * nb_cols]);
      // diagonal from top right
      apply_penalty(&penalty[3], p1_in[3 + col * nb_dir + row * nb_dir * nb_cols], p2_in[3 + col * nb_dir + row * nb_dir * nb_cols]);

      // initialize minimums
      min0 = std::numeric_limits<float>::max(), min1 = std::numeric_limits<float>::max(), min2 = std::numeric_limits<float>::max(), min3 = std::numeric_limits<float>::max();
      pos0 = 0, pos1 = 0, pos2 = 0, pos3 = 0;
      // Get current class of pixel
      current_class = segmentation[col + row * nb_cols];

      bool isBorder0 = (row - direction[0].drow < 0 || row - direction[0].drow > (nb_rows - 1) || col - direction[0].dcol < 0 || col - direction[0].dcol > (nb_cols - 1));
      bool isBorder1 = (row - direction[1].drow < 0 || row - direction[1].drow > (nb_rows - 1) || col - direction[1].dcol < 0 || col - direction[1].dcol > (nb_cols - 1));
      bool isBorder2 = (row - direction[2].drow < 0 || row - direction[2].drow > (nb_rows - 1) || col - direction[2].dcol < 0 || col - direction[2].dcol > (nb_cols - 1));
      bool isBorder3 = (row - direction[3].drow < 0 || row - direction[3].drow > (nb_rows - 1) || col - direction[3].dcol < 0 || col - direction[3].dcol > (nb_cols - 1));
      
      bool noOp0 = isBorder0 || current_class != buff_class0;
      bool noOp1 = isBorder1 || current_class != buff_class1_2_3[col - direction[1].dcol];
      bool noOp2 = isBorder2 || current_class != buff_class1_2_3[col - direction[2].dcol];
      bool noOp3 = isBorder3 || current_class != buff_class1_2_3[col - direction[3].dcol];

      T* cv_in_base = &cv_in[col * nb_disps + row * nb_disps * nb_cols];
      Tout* cvs_out_base = &cvs.cost_volume[col * nb_disps + row * nb_disps * nb_cols];

      if(noOp0){
        for (size_t disp = 0; disp < nb_disps; disp++){
          T costAggr0 = cv_in_base[disp];

          // Save data in buffer
          buff0[disp] = costAggr0;
          cvs_out_base[disp] += costAggr0;
          if (cost_paths) update_minimum(min0, costAggr0, pos0, disp);
        }
      }else{
        T min_disp0 = *std::min_element(&buff0[0], &buff0[nb_disps]);
        T lastBuff = buff0[0];
        for (size_t disp = 0; disp < nb_disps; disp++){
          T costAggr0 = costCompute<T, 0>(cv_in_base[disp], 0, disp, invalid_value, nb_disps, penalty[0], buff0, min_disp0, lastBuff);
          // Save buffer data before writing new value in buffer
          lastBuff = buff0[disp];

          // Save data in buffer
          buff0[disp] = costAggr0;
          cvs_out_base[disp] += costAggr0;
          if (cost_paths) update_minimum(min0, costAggr0, pos0, disp);
        }
      }

      if(noOp1){
        for (size_t disp = 0; disp < nb_disps; disp++){
          T costAggr1 = cv_in_base[disp];
          // Save data in buffer
          buff1[disp + col * nb_disps] = costAggr1;
          cvs_out_base[disp] += costAggr1;
          if (cost_paths) update_minimum(min1, costAggr1, pos1, disp);
        }
      }else{
        T min_disp1 = *std::min_element(&buff1[0 + col * nb_disps], &buff1[nb_disps + col * nb_disps]);
        T lastBuff = buff1[0 + col * nb_disps];
        for (size_t disp = 0; disp < nb_disps; disp++){
          T costAggr1 = costCompute<T, 1>(cv_in_base[disp], col, disp, invalid_value, nb_disps, penalty[1], buff1, min_disp1, lastBuff);
          // sava buffer data before writing new value in buffer
          lastBuff = buff1[disp + col * nb_disps];
          // Save data in buffer
          buff1[disp + col * nb_disps] = costAggr1;
          cvs_out_base[disp] += costAggr1;
          if (cost_paths) update_minimum(min1, costAggr1, pos1, disp);
        }
      }
      
      if(noOp2){
        for (size_t disp = 0; disp < nb_disps; disp++){
          T costAggr2 = cv_in_base[disp];
          // Save data in buffer
          buff_disp_2[disp] = buff2[disp + col * nb_disps];
          buff2[disp + col * nb_disps] = costAggr2;
          cvs_out_base[disp] += costAggr2;
          if (cost_paths) update_minimum(min2, costAggr2, pos2, disp);
        }
      }else{
        T min_disp2 = *std::min_element(&buff_disp_2[0], &buff_disp_2[nb_disps]);
        T lastBuff = buff_disp_2[0];
        for (size_t disp = 0; disp < nb_disps; disp++){
          T costAggr2 = costCompute<T, 0>(cv_in_base[disp], col, disp, invalid_value, nb_disps, penalty[2], buff_disp_2, min_disp2, lastBuff);
          // Save buffer data before writing new value in buffer_disp
          lastBuff = buff_disp_2[disp];
          // Save data in buffer
          buff_disp_2[disp] = buff2[disp + col * nb_disps];
          buff2[disp + col * nb_disps] = costAggr2;
          cvs_out_base[disp] += costAggr2;
          if (cost_paths) update_minimum(min2, costAggr2, pos2, disp);
        }
      }
      
      if(noOp3){
        for (size_t disp = 0; disp < nb_disps; disp++){
          T costAggr3 = cv_in_base[disp];
          buff3[disp + col * nb_disps] = costAggr3;
          cvs_out_base[disp] += costAggr3;
          if (cost_paths) update_minimum(min3, costAggr3, pos3, disp);
        }
      }else{
        int newCol = col - direction[3].dcol;
        T min_disp3 = *std::min_element(&buff3[0 + newCol * nb_disps], &buff3[nb_disps + newCol * nb_disps]);
        T lastBuff = buff3[0 + newCol * nb_disps];
        for (size_t disp = 0; disp < nb_disps; disp++){
          T costAggr3 = costCompute<T, 1>(cv_in_base[disp], newCol, disp, invalid_value, nb_disps, penalty[3], buff3, min_disp3, lastBuff);

          // Save data in buffer
          lastBuff = buff3[disp + newCol * nb_disps];
          buff3[disp + col * nb_disps] = costAggr3;
          cvs_out_base[disp] += costAggr3;
          if (cost_paths) update_minimum(min3, costAggr3, pos3, disp);
        }
      }
      
      if (cost_paths)
      {
        cvs.cost_volume_min[0 + col * nb_dir + row * nb_dir * nb_cols] = pos0;
        cvs.cost_volume_min[1 + col * nb_dir + row * nb_dir * nb_cols] = pos1;
        cvs.cost_volume_min[2 + col * nb_dir + row * nb_dir * nb_cols] = pos2;
        cvs.cost_volume_min[3 + col * nb_dir + row * nb_dir * nb_cols] = pos3;
      }
      // Update buffer
      buff_class0 = current_class;
      tmp_buff_class1_2_3[col] = current_class;
    }
    // Update buffer
    std::copy(tmp_buff_class1_2_3, tmp_buff_class1_2_3 + nb_cols, buff_class1_2_3);
  }
  // delete buffers
  delete[] buff0;
  delete[] buff1;
  delete[] buff2;
  delete[] buff3;
  delete[] buff_disp_2;
  delete[] buff_class1_2_3;
  delete[] tmp_buff_class1_2_3;

  /* ------------------*/
  /* -- Second pass -- */
  /* ------------------*/

  // temporary buffers and variables for storing aggregated cost
  T *buff4 = new T[nb_disps]();
  T *buff5 = new T[nb_cols * nb_disps]();
  T *buff6 = new T[nb_cols * nb_disps]();
  T *buff7 = new T[nb_cols * nb_disps]();
  T *buff_disp_6 = new T[nb_disps]();
  // temporary variables for storing current position of minimums
  int pos4, pos5, pos6, pos7;
  // temporary variables for storing current minimums
  float min4, min5, min6, min7;

  // temporary buffers for storing current classification
  // buff_class4 stores the class from the right of current pixel (previous seen pixel)
  float buff_class4 = segmentation[(nb_cols - 1) + (nb_rows - 1) * nb_cols];
  // buff_class5_6_7 stores the previous line of classes
  float *buff_class5_6_7 = new float[nb_cols]();
  // tmp_buff_class5_6_7 stores the previous classes of pixels seen in the line
  float *tmp_buff_class5_6_7 = new float[nb_cols]();

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

  for (int row = nb_rows - 1; row >= 0; row--)
  {
    for (int col = nb_cols - 1; col >= 0; col--)
    {
      // right->left
      apply_penalty(&penalty[4], p1_in[4 + col * nb_dir + row * nb_dir * nb_cols], p2_in[4 + col * nb_dir + row * nb_dir * nb_cols]);
      // up -> down
      apply_penalty(&penalty[5], p1_in[5 + col * nb_dir + row * nb_dir * nb_cols], p2_in[5 + col * nb_dir + row * nb_dir * nb_cols]);
      // diagonal from lower left
      apply_penalty(&penalty[6], p1_in[6 + col * nb_dir + row * nb_dir * nb_cols], p2_in[6 + col * nb_dir + row * nb_dir * nb_cols]);
      // diagonal from lower right
      apply_penalty(&penalty[7], p1_in[7 + col * nb_dir + row * nb_dir * nb_cols], p2_in[7 + col * nb_dir + row * nb_dir * nb_cols]);

      // initialize minimums
      min4 = std::numeric_limits<float>::max(), min5 = std::numeric_limits<float>::max(), min6 = std::numeric_limits<float>::max(), min7 = std::numeric_limits<float>::max();
      pos4 = 0, pos5 = 0, pos6 = 0, pos7 = 0;
      // Get current class of pixel
      current_class = segmentation[col + row * nb_cols];

      bool isBorder4 = (row - direction[4].drow < 0 || row - direction[4].drow > (nb_rows - 1) || col - direction[4].dcol < 0 || col - direction[4].dcol > (nb_cols - 1));
      bool isBorder5 = (row - direction[5].drow < 0 || row - direction[5].drow > (nb_rows - 1) || col - direction[5].dcol < 0 || col - direction[5].dcol > (nb_cols - 1));
      bool isBorder6 = (row - direction[6].drow < 0 || row - direction[6].drow > (nb_rows - 1) || col - direction[6].dcol < 0 || col - direction[6].dcol > (nb_cols - 1));
      bool isBorder7 = (row - direction[7].drow < 0 || row - direction[7].drow > (nb_rows - 1) || col - direction[7].dcol < 0 || col - direction[7].dcol > (nb_cols - 1));

      bool noOp4 = isBorder4 || current_class != buff_class4;
      bool noOp5 = isBorder5 || current_class != buff_class5_6_7[col - direction[5].dcol];
      bool noOp6 = isBorder6 || current_class != buff_class5_6_7[col - direction[6].dcol];
      bool noOp7 = isBorder7 || current_class != buff_class5_6_7[col - direction[7].dcol];

      T* cv_in_base = &cv_in[col * nb_disps + row * nb_disps * nb_cols];
      Tout* cvs_out_base = &cvs.cost_volume[col * nb_disps + row * nb_disps * nb_cols];

      if(noOp4){
        for (size_t disp = 0; disp < nb_disps; disp++){
          T costAggr4 = cv_in_base[disp];
          // Save data in buffer
          buff4[disp] = costAggr4;
          cvs_out_base[disp] += costAggr4;
          if (cost_paths) update_minimum(min4, costAggr4, pos4, disp);
        }
      }else{
        T min_disp4 = *std::min_element(&buff4[0], &buff4[nb_disps]);
        T lastBuff = buff4[0];
        for (size_t disp = 0; disp < nb_disps; disp++){
          T costAggr4 = costCompute<T, 0>(cv_in_base[disp], row, disp, invalid_value, nb_disps, penalty[4], buff4, min_disp4, lastBuff);
          // sava buffer data before writing new value in buffer
          lastBuff = buff4[disp];
          // Save data in buffer
          buff4[disp] = costAggr4;
          cvs_out_base[disp] += costAggr4;
          if (cost_paths) update_minimum(min4, costAggr4, pos4, disp);
        }
      }

      if(noOp5){
        for (size_t disp = 0; disp < nb_disps; disp++){
          T costAggr5 = cv_in_base[disp];
          // Save data in buffer
          buff5[disp + col * nb_disps] = costAggr5;
          cvs_out_base[disp] += costAggr5;
          if (cost_paths) update_minimum(min5, costAggr5, pos5, disp);
        }
      }else{
        T min_disp5 = *std::min_element(&buff5[0 + col * nb_disps], &buff5[nb_disps + col * nb_disps]);
        T lastBuff = buff5[0 + col * nb_disps];
        for (size_t disp = 0; disp < nb_disps; disp++){
          T costAggr5 = costCompute<T, 1>(cv_in_base[disp], col, disp, invalid_value, nb_disps, penalty[5], buff5, min_disp5, lastBuff);
          // Save buffer data before writing new value in buffer
          lastBuff = buff5[disp + col * nb_disps];
          // Save data in buffer
          buff5[disp + col * nb_disps] = costAggr5;
          cvs_out_base[disp] += costAggr5;
          if (cost_paths) update_minimum(min5, costAggr5, pos5, disp);
        }
      }
      
      if(noOp6){
        for (size_t disp = 0; disp < nb_disps; disp++){
          T costAggr6 = cv_in_base[disp];
          buff_disp_6[disp] = buff6[disp + col * nb_disps];
          buff6[disp + col * nb_disps] = costAggr6;
          cvs_out_base[disp] += costAggr6;
          if (cost_paths) update_minimum(min6, costAggr6, pos6, disp);
        }
      }else{
        T min_disp6 = *std::min_element(&buff_disp_6[0], &buff_disp_6[nb_disps]);
        T lastBuff = buff_disp_6[0];
        for (size_t disp = 0; disp < nb_disps; disp++){
          T costAggr6 = costCompute<T, 0>(cv_in_base[disp], col, disp, invalid_value, nb_disps, penalty[6], buff_disp_6, min_disp6, lastBuff);
          // Save buffer data before writing new value in buffer_disp
          lastBuff = buff_disp_6[disp];
          // Save data in buffer
          buff_disp_6[disp] = buff6[disp + col * nb_disps];
          buff6[disp + col * nb_disps] = costAggr6;
          cvs_out_base[disp] += costAggr6;
          if (cost_paths) update_minimum(min6, costAggr6, pos6, disp);
        }
      }
      
      if(noOp7){
        for (size_t disp = 0; disp < nb_disps; disp++){
          T costAggr7 = cv_in_base[disp];
          buff7[disp + col * nb_disps] = costAggr7;
          cvs_out_base[disp] += costAggr7;
          cvs_out_base[disp] -= overcounting_factor*costAggr7;
          if (cost_paths) update_minimum(min7, costAggr7, pos7, disp);
        }
      }else{
        int newCol = col - direction[7].dcol;
        T min_disp7 = *std::min_element(&buff7[0 + newCol * nb_disps], &buff7[nb_disps + newCol * nb_disps]);
        T lastBuff = buff7[0 + newCol * nb_disps];
        for (size_t disp = 0; disp < nb_disps; disp++){
          T costAggr7 = costCompute<T, 1>(cv_in_base[disp], newCol, disp, invalid_value, nb_disps, penalty[7], buff7, min_disp7, lastBuff);
          lastBuff = buff7[disp + newCol * nb_disps];
          // Save buffer data before writing new value in buffer
          buff7[disp + col * nb_disps] = costAggr7;
          cvs_out_base[disp] += costAggr7;
          cvs_out_base[disp] -= overcounting_factor*cv_in_base[disp];
          if (cost_paths) update_minimum(min7, costAggr7, pos7, disp);
        }
      }
      if (cost_paths)
      {
        cvs.cost_volume_min[4 + col * nb_dir + row * nb_dir * nb_cols] = pos4;
        cvs.cost_volume_min[5 + col * nb_dir + row * nb_dir * nb_cols] = pos5;
        cvs.cost_volume_min[6 + col * nb_dir + row * nb_dir * nb_cols] = pos6;
        cvs.cost_volume_min[7 + col * nb_dir + row * nb_dir * nb_cols] = pos7;
      }

      // Update buffer
      buff_class4 = current_class;
      tmp_buff_class5_6_7[col] = current_class;
    }
    // Update buffer
    std::copy(tmp_buff_class5_6_7, tmp_buff_class5_6_7 + nb_cols, buff_class5_6_7);
  }

  delete[] buff4;
  delete[] buff5;
  delete[] buff6;
  delete[] buff7;
  delete[] buff_disp_6;
  delete[] buff_class5_6_7;
  delete[] tmp_buff_class5_6_7;
}

template <typename T>
void apply_penalty(Penalty<T> *penalty, T p1, T p2)
{
  penalty->P1 = p1;
  penalty->P2 = p2;
}

void assignDirections(int *directions_in, Direction *dirs)
{
  for (int i = 0; i < 8; i++)
  {
    dirs[i].drow = directions_in[2 * i];
    dirs[i].dcol = directions_in[2 * i + 1];
  }
}

/* Explicitly instantiate all the templates needed to use libSGM as an external lib */
template void sgm<uint8_t, uint16_t>(uint8_t *cv_in, uint8_t *p1_in, uint8_t *p2_in, int *directions_in, unsigned long int nb_rows,
                                                      unsigned long int nb_cols, unsigned int nb_disps, uint8_t invalid_value, float *segmentation, bool cost_paths, bool overcounting, CostVolumes<uint16_t> &cvs);
template void sgm<float, float>(float *cv_in, float *p1_in, float *p2_in, int *directions_in, unsigned long int nb_rows,
                                              unsigned long int nb_cols, unsigned int nb_disps, float invalid_value, float *segmentation, bool cost_paths, bool overcounting, CostVolumes<float> &cvs);