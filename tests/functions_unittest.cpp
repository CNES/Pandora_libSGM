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

#include "gtest/gtest.h"
#include "../../src/libSGM/lib/sgm.hpp"


//Global Test of sgm function: aggregation value from 8 directions on a middle point of cost volume
//with invalid middle point and without over-counting correction

TEST(sgmTest, TestMiddleValueInvalid){

    int nb_row,nb_col,nb_disp;
	nb_row = 3;
	nb_col = 3;
	nb_disp = 3;
	uint8_t P1, P2, invalid_value;
	P1 = 8;
	P2 = 32;
	invalid_value = 57;

	bool overcounting = false;
	bool cost_paths = false;

	CostVolumes cvs;
	uint8_t cv_in[27]={1,15,20,14,16,6,8,19,8,13,11,3,57,57,57,16,4,12,18,2,17,23,7,1,5,20,14};

	// method : constant
	uint8_t p1[9*8]={P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1};
	uint8_t p2[9*8]={P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2};
    int directions[2*8]={0,1, 1,0, 1,1, 1,-1, 0,-1, -1,0, -1,-1, -1,1};

    // segmentation
    float segmentation[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1}; // no piecewise optimization

	cvs = sgm(cv_in, p1, p2, directions, nb_row, nb_col, nb_disp,invalid_value, segmentation, cost_paths, overcounting);
	//middle point must stay invalid, equal to invalid value
	EXPECT_EQ(456,cvs.cost_volume[13]);
}

//with invalid middle point and with over-counting correction

TEST(sgmTest, TestMiddleValueInvalidOvercounting){

    int nb_row,nb_col,nb_disp;
	nb_row = 3;
	nb_col = 3;
	nb_disp = 3;
	uint8_t P1, P2, invalid_value;
	P1 = 8;
	P2 = 32;
	invalid_value = 57;

	bool overcounting = true;
	bool cost_paths = false;

	CostVolumes cvs;
	uint8_t cv_in[27]={1,15,20,14,16,6,8,19,8,13,11,3,57,57,57,16,4,12,18,2,17,23,7,1,5,20,14};

	// method : constant
	uint8_t p1[9*8]={P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1};
	uint8_t p2[9*8]={P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2};
    int directions[2*8]={0,1, 1,0, 1,1, 1,-1, 0,-1, -1,0, -1,-1, -1,1};

    // segmentation
    float segmentation[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1}; // no piecewise optimization

	cvs = sgm(cv_in, p1, p2, directions, nb_row, nb_col, nb_disp,invalid_value, segmentation, cost_paths, overcounting);
	//middle point must stay invalid, equal to invalid value
	EXPECT_EQ(57,cvs.cost_volume[13]);
}

//Global Test of sgm function: aggregation value from 8 directions on a middle point of cost volume with min cost

TEST(sgmTest, TestMiddleValue_cost_paths){

    int nb_row,nb_col,nb_disp;
	nb_row = 3;
	nb_col = 3;
	nb_disp = 3;
	uint8_t P1, P2, invalid_value;
	P1 = 8;
	P2 = 32;
	invalid_value = 57;

	CostVolumes cvs;
	uint8_t cv_in[27]={1,15,20,14,16,6,8,19,8,13,11,3,22,12,9,16,4,12,18,2,17,23,7,1,5,20,14};

	bool overcounting = false;
	bool cost_paths = true;

    // method : constant
	uint8_t p1[9*8]={P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1};
	uint8_t p2[9*8]={P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2};
    int directions[2*8]={0,1, 1,0, 1,1, 1,-1, 0,-1, -1,0, -1,-1, -1,1};

    // segmentation
    float segmentation[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1}; // no piecewise optimization

	cvs = sgm(cv_in, p1, p2, directions, nb_row, nb_col, nb_disp,invalid_value, segmentation, cost_paths, overcounting);
	//middle point must stay invalid, equal to invalid value
	EXPECT_EQ(142,cvs.cost_volume[13]);
	EXPECT_EQ(2,cvs.cost_volume_min[32]);
}


//Global Test of sgm function: aggregation value from 8 directions on a middle point of cost volume

TEST(sgmTest, TestMiddleValue){

    int nb_row,nb_col,nb_disp;
	nb_row = 3;
	nb_col = 3;
	nb_disp = 3;
	uint8_t P1, P2, invalid_value;
	P1 = 8;
	P2 = 32;
	invalid_value = 57;
	bool overcounting = false;
	bool cost_paths = false;

	CostVolumes cvs;
	uint8_t cv_in[27]={1,15,20,14,16,6,8,19,8,13,11,3,22,12,9,16,4,12,18,2,17,23,7,1,5,20,14};

	// method : constant
	uint8_t p1[9*8]={P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1};
	uint8_t p2[9*8]={P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2};
    int directions[2*8]={0,1, 1,0, 1,1, 1,-1, 0,-1, -1,0, -1,-1, -1,1};

    // segmentation
    float segmentation[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1}; // no piecewise optimization

	cvs = sgm(cv_in, p1, p2, directions, nb_row, nb_col, nb_disp,invalid_value, segmentation, cost_paths, overcounting);

	EXPECT_EQ(142,cvs.cost_volume[13]);
}

//Global Test of sgm function: aggregation value from 8 directions on a middle point of cost volume, and reset history in the middle

TEST(sgmTest, TestMiddleValueResetHistory){

    int nb_row,nb_col,nb_disp;
	nb_row = 3;
	nb_col = 3;
	nb_disp = 3;
	uint8_t P1, P2, invalid_value;
	P1 = 8;
	P2 = 32;
	invalid_value = 57;
	bool overcounting = false;
	bool cost_paths = false;

	CostVolumes cvs;
	uint8_t cv_in[27]={1,15,20,14,16,6,8,19,8,13,11,3,22,12,9,16,4,12,18,2,17,23,7,1,5,20,14};

	// method : constant
	uint8_t p1[9*8]={P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1};
	uint8_t p2[9*8]={P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2};
    int directions[2*8]={0,1, 1,0, 1,1, 1,-1, 0,-1, -1,0, -1,-1, -1,1};

    // segmentation
    float segmentation[9] = {1, 1, 1, 1, 2, 1, 1, 1, 1}; // piecewise optimization at the middle

	cvs = sgm(cv_in, p1, p2, directions, nb_row, nb_col, nb_disp,invalid_value, segmentation, cost_paths, overcounting);

	EXPECT_EQ(8 * cv_in[13],cvs.cost_volume[13]);
}

//with valid middle point and with over-counting correction
TEST(sgmTest, TestMiddleValueOvercounting){

    int nb_row,nb_col,nb_disp;
	nb_row = 3;
	nb_col = 3;
	nb_disp = 3;
	uint8_t P1, P2, invalid_value;
	P1 = 8;
	P2 = 32;
	invalid_value = 57;
	float alpha = 100;
	uint8_t beta = 1;
	uint8_t gamma = 2;
	bool overcounting = true;
	bool cost_paths = false;

	CostVolumes cvs;
	uint8_t cv_in[27]={1,15,20,14,16,6,8,19,8,13,11,3,22,12,9,16,4,12,18,2,17,23,7,1,5,20,14};

	// method : constant
	uint8_t p1[9*8]={P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1};
	uint8_t p2[9*8]={P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2};
    int directions[2*8]={0,1, 1,0, 1,1, 1,-1, 0,-1, -1,0, -1,-1, -1,1};

    // segmentation
    float segmentation[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1}; // no piecewise optimization

	cvs = sgm(cv_in, p1, p2, directions, nb_row, nb_col, nb_disp,invalid_value, segmentation, cost_paths, overcounting);

	EXPECT_EQ(58,cvs.cost_volume[13]);
}

// Test function aggregatedCostFromTopLeft0

// Test aggregation on cost Volume border
TEST(aggregatedCostFromTopLeft0Test, borderValue) {

	uint8_t pixel_0 ;
	uint8_t min_disp0 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	uint8_t invalid_value = 57;
	Direction direction={0,1};

	int row = 0;
	int col = 0;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff0 = new uint8_t[nb_disps]();

	// class variables
	float class0 = 1;
	float buff_class0 = 1; // same class for no piecewise optimization
	float reset0 = 1; // doesn'nt matter, disp==0

	uint8_t costAggr = aggregatedCostFromTopLeft0(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff0, &min_disp0, &pixel_0, class0, buff_class0, reset0);

	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset0);
}

// Test aggregation on cost Volume border and reset history
TEST(aggregatedCostFromTopLeft0Test, borderValueResetHistory) {

	uint8_t pixel_0 ;
	uint8_t min_disp0 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	uint8_t invalid_value = 57;
	Direction direction={0,1};

	int row = 0;
	int col = 0;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff0 = new uint8_t[nb_disps]();

	// class variables
	float class0 = 1;
	float buff_class0 = 2; // different class for piecewise optimization
	float reset0 = 1; // doesn't matter, disp==0 and not a workable pixel

	uint8_t costAggr = aggregatedCostFromTopLeft0(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff0, &min_disp0, &pixel_0, class0, buff_class0, reset0);

	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset0); // was not changed
}

// Test minimum value of the previous point
TEST(aggregatedCostFromTopLeft0Test, findMinDisp) {

	uint8_t pixel_0 ;
	uint8_t min_disp0 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={0,1};
	uint8_t invalid_value = 57;

	int row = 1;
	int col = 1;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff0 = new uint8_t[nb_disps]();
	buff0[0] = 14 ;
	buff0[1] = 28 ;
	buff0[2] = 20 ;

	// class variables
	float class0 = 1;
	float buff_class0 = 1; // same class for no piecewise optimization
	float reset0 = 1;

	uint8_t costAggr = aggregatedCostFromTopLeft0(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff0, &min_disp0, &pixel_0, class0, buff_class0, reset0);


	//pixel_0 = min(14,28,20)
	EXPECT_EQ(14,min_disp0);
	//costAggr = 15 + min(14 , 28 + 8 , 14 + 32) - 14 = 15
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset0);

}

// Test aggregation value of a non-border point
TEST(aggregatedCostFromTopLeft0Test, aggregation) {

	uint8_t pixel_0 = 14 ;
	uint8_t min_disp0 = 14 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={0,1};
	uint8_t invalid_value = 57;

	int row = 1;
	int col = 1;
	int disp = 1;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff0 = new uint8_t[nb_disps]();
	buff0[0] = 36 ;
	buff0[1] = 28 ;
	buff0[2] = 20 ;

	// class variables
	float class0 = 1;
	float buff_class0 = 1; // same class for no piecewise optimization
	float reset0 = 1;

	uint8_t costAggr = aggregatedCostFromTopLeft0(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff0, &min_disp0, &pixel_0, class0, buff_class0, reset0);

	//costAggr = 15 + min ( 14+8 , 28 , 20 + 8 , 14+32) - 14 = 15+22-14 = 23
	EXPECT_EQ(23,costAggr);
	EXPECT_EQ(1,reset0);

}


// Test aggregation value of a non-border point and reset history with reset already computed
TEST(aggregatedCostFromTopLeft0Test, aggregationResetHistory_disp1) {

	uint8_t pixel_0 = 14 ;
	uint8_t min_disp0 = 14 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={0,1};
	uint8_t invalid_value = 57;

	int row = 1;
	int col = 1;
	int disp = 1;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff0 = new uint8_t[nb_disps]();
	buff0[0] = 36 ;
	buff0[1] = 28 ;
	buff0[2] = 20 ;

	// class variables
	float class0 = 1;
	float buff_class0 = 2; // different class for  piecewise optimization and reset history
	float reset0 = 0;

	uint8_t costAggr = aggregatedCostFromTopLeft0(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff0, &min_disp0, &pixel_0, class0, buff_class0, reset0);

	//costAggr = 15 +  0 * min ( 14+8 , 28 , 20 + 8 , 14+32) - 14 = 15+22-14 = 23
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(0,reset0);

}

// Test aggregation value of a non-border point and reset history with update of reset
TEST(aggregatedCostFromTopLeft0Test, aggregationResetHistory_disp0) {

	uint8_t pixel_0 = 14 ;
	uint8_t min_disp0 = 14 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={0,1};
	uint8_t invalid_value = 57;

	int row = 1;
	int col = 1;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff0 = new uint8_t[nb_disps]();
	buff0[0] = 36 ;
	buff0[1] = 28 ;
	buff0[2] = 20 ;

	// class variables
	float class0 = 1;
	float buff_class0 = 2; // different class for  piecewise optimization and reset history
	float reset0 = 1;

	uint8_t costAggr = aggregatedCostFromTopLeft0(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff0, &min_disp0, &pixel_0, class0, buff_class0, reset0);

	//costAggr = 15 +  0 * min ( 14+8 , 28 , 20 + 8 , 14+32) - 14 = 15+22-14 = 23
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(0,reset0);

}


// Test aggregation value after an invalid point
TEST(aggregatedCostFromTopLeft0Test, aggregationAfterInvalidPoint) {

	uint8_t pixel_0 ;
	uint8_t min_disp0 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={0,1};
	uint8_t invalid_value = 57;

	int row = 1;
	int col = 1;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff0 = new uint8_t[nb_disps]();
	buff0[0] = 57 ;
	buff0[1] = 57 ;
	buff0[2] = 57 ;

	// class variables
	float class0 = 1;
	float buff_class0 = 1; // same class for no piecewise optimization
	float reset0 = 1;

	uint8_t costAggr = aggregatedCostFromTopLeft0(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff0, &min_disp0, &pixel_0, class0, buff_class0, reset0);

	//pixel_0 = min(14,28,20)
	EXPECT_EQ(57,min_disp0);
	//costAggr = 15 + min(14 , 28 + 8 , 14 + 32) - 14 = 15
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset0);

}

// Test function aggregatedCostFromTopLeft1

// Test aggregation on cost Volume border
TEST(aggregatedCostFromTopLeft1Test, borderValue) {

	uint8_t pixel_1 ;
	uint8_t min_disp1 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={1,0};
	uint8_t invalid_value = 57;

	int row = 0;
	int col = 0;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff1 = new uint8_t[nb_cols*nb_disps]();

	// class variables
	float class1 = 1;
	float buff_class1[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset1 = 1;

	uint8_t costAggr = aggregatedCostFromTopLeft1(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff1, &min_disp1, &pixel_1, class1, buff_class1, reset1) ;

	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset1);
}

// Test aggregation on cost Volume border
TEST(aggregatedCostFromTopLeft1Test, borderValueResetHistory) {

	uint8_t pixel_1 ;
	uint8_t min_disp1 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={1,0};
	uint8_t invalid_value = 57;

	int row = 0;
	int col = 0;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff1 = new uint8_t[nb_cols*nb_disps]();

	// class variables
	float class1 = 1;
	float buff_class1[5] = {0,1,1,1,1}; // different class for no piecewise optimization
	float reset1 = 1;

	uint8_t costAggr = aggregatedCostFromTopLeft1(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff1, &min_disp1, &pixel_1, class1, buff_class1, reset1) ;

	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset1);// unchanged
}

// Test minimum value of the previous point
TEST(aggregatedCostFromTopLeft1Test, findMinDisp) {

	uint8_t pixel_1 ;
	uint8_t min_disp1 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={1,0};
	uint8_t invalid_value = 57;

	int row = 1;
	int col = 1;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff1 = new uint8_t[nb_cols*nb_disps]();
	buff1[3] = 17 ;
	buff1[4] = 45 ;
	buff1[5] = 32 ;

	// class variables
	float class1 = 1;
	float buff_class1[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset1 = 1;

	uint8_t costAggr = aggregatedCostFromTopLeft1(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff1, &min_disp1, &pixel_1, class1, buff_class1, reset1) ;

	EXPECT_EQ(17,min_disp1);
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset1);

}

// Test aggregation value of a non-border point
TEST(aggregatedCostFromTopLeft1Test, aggregation) {

	uint8_t pixel_1 = 30 ;
	uint8_t min_disp1 = 28  ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8 ;
	uint8_t P2 = 32 ;
	Direction direction={1,0};
	uint8_t invalid_value = 57;

	int row = 1;
	int col = 1;
	int disp = 1;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff1 = new uint8_t[nb_disps]();
	buff1[3] = 52 ;
	buff1[4] = 28 ;
	buff1[5] = 32 ;

	// class variables
	float class1 = 1;
	float buff_class1[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset1 = 1;

	uint8_t costAggr = aggregatedCostFromTopLeft1(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff1, &min_disp1, &pixel_1, class1, buff_class1, reset1) ;

	//costAggr = 15 + min(30+8 , 28 , 32+8 , 28+32) -28
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset1);

}

// Test aggregation value of a non-border point and reset history, with reset already computed
TEST(aggregatedCostFromTopLeft1Test, aggregationResetHistory_disp1) {

	uint8_t pixel_1 = 30 ;
	uint8_t min_disp1 = 26  ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8 ;
	uint8_t P2 = 32 ;
	Direction direction={1,0};
	uint8_t invalid_value = 57;

	int row = 1;
	int col = 1;
	int disp = 1;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff1 = new uint8_t[nb_disps]();
	buff1[3] = 52 ;
	buff1[4] = 28 ;
	buff1[5] = 32 ;

	// class variables
	float class1 = 1;
	float buff_class1[5] = {1,1,1,1,1}; // unused
	float reset1 = 0;

	uint8_t costAggr = aggregatedCostFromTopLeft1(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff1, &min_disp1, &pixel_1, class1, buff_class1, reset1) ;

	//costAggr = 15 + min(30+8 , 28 , 32+8 , 28+32) -26
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(0,reset1);

}

// Test aggregation value of a non-border point and reset history, and update reset
TEST(aggregatedCostFromTopLeft1Test, aggregationResetHistory_disp0) {

	uint8_t pixel_1 = 30 ;
	uint8_t min_disp1 = 26  ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8 ;
	uint8_t P2 = 32 ;
	Direction direction={1,0};
	uint8_t invalid_value = 57;

	int row = 1;
	int col = 1;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff1 = new uint8_t[nb_disps]();
	buff1[3] = 52 ;
	buff1[4] = 28 ;
	buff1[5] = 32 ;

	// class variables
	float class1 = 1;
	float buff_class1[5] = {1,2,1,1,1}; // different class for  piecewise optimization
	float reset1 = 1;

	uint8_t costAggr = aggregatedCostFromTopLeft1(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff1, &min_disp1, &pixel_1, class1, buff_class1, reset1) ;

	//costAggr = 15 + min(30+8 , 28 , 32+8 , 28+32) -26
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(0,reset1);

}

// Test aggregation value after an invalid point
TEST(aggregatedCostFromTopLeft1Test, aggregationAfterInvalidPoint) {

	uint8_t pixel_1 ;
	uint8_t min_disp1 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={1,0};
	uint8_t invalid_value = 57;

	int row = 1;
	int col = 1;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff1 = new uint8_t[nb_cols*nb_disps]();
	buff1[3] = 57 ;
	buff1[4] = 57 ;
	buff1[5] = 57 ;

	// class variables
	float class1 = 1;
	float buff_class1[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset1 = 1;

	uint8_t costAggr = aggregatedCostFromTopLeft1(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff1, &min_disp1, &pixel_1, class1, buff_class1, reset1) ;

	EXPECT_EQ(57,min_disp1);
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset1);

}

// Test function aggregatedCostFromTopLeft2

// Test aggregation on cost Volume border
TEST(aggregatedCostFromTopLeft2Test, borderValue) {

	uint8_t pixel_2 ;
	uint8_t min_disp2 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={1,1};
	uint8_t invalid_value = 57;

	int row = 0;
	int col = 0;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff2 = new uint8_t[nb_cols*nb_disps]();
	uint8_t * buff_disp_2 = new uint8_t[nb_disps]();

	// class variables
	float class2 = 1;
	float buff_class2[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset2 = 1;

	uint8_t costAggr = aggregatedCostFromTopLeft2(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1 ,P2 ,direction, buff2, buff_disp_2, &min_disp2, &pixel_2, class2, buff_class2, reset2) ;

	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset2);
}

// Test aggregation on cost Volume border
TEST(aggregatedCostFromTopLeft2Test, borderValueResetHistory) {

	uint8_t pixel_2 ;
	uint8_t min_disp2 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={1,1};
	uint8_t invalid_value = 57;

	int row = 0;
	int col = 0;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff2 = new uint8_t[nb_cols*nb_disps]();
	uint8_t * buff_disp_2 = new uint8_t[nb_disps]();

	// class variables
	float class2 = 1;
	float buff_class2[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset2 = 0;

	uint8_t costAggr = aggregatedCostFromTopLeft2(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1 ,P2 ,direction, buff2, buff_disp_2, &min_disp2, &pixel_2, class2, buff_class2, reset2) ;

	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(0,reset2); // unchanged
}

// Test minimum value of the previous point
TEST(aggregatedCostFromTopLeft2Test, findMinDisp) {

	uint8_t pixel_2 ;
	uint8_t min_disp2;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={1,1};
	uint8_t invalid_value = 57;

	int row = 1;
	int col = 1;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff2 = new uint8_t[nb_cols*nb_disps]();
	uint8_t * buff_disp_2 = new uint8_t[nb_disps]();
	buff_disp_2[0] = 52 ;
	buff_disp_2[1] = 28 ;
	buff_disp_2[2] = 32 ;

	// class variables
	float class2 = 1;
	float buff_class2[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset2 = 1;

	uint8_t costAggr = aggregatedCostFromTopLeft2(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1 ,P2 ,direction, buff2, buff_disp_2, &min_disp2, &pixel_2, class2, buff_class2, reset2) ;

	EXPECT_EQ(28,min_disp2);
	//costAggr = 15 + min(52 , 28+8,32+28) - 28
	EXPECT_EQ(23,costAggr);
	EXPECT_EQ(1,reset2);

}

// Test aggregation value of a non-border point
TEST(aggregatedCostFromTopLeft2Test, aggregation) {

	uint8_t pixel_2 = 30 ;
	uint8_t min_disp2 = 28  ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8 ;
	uint8_t P2 = 32 ;
	Direction direction={1,1};
	uint8_t invalid_value = 57;

	int row = 1;
	int col = 1;
	int disp = 1;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff2 = new uint8_t[nb_cols*nb_disps]();
	uint8_t * buff_disp_2 = new uint8_t[nb_disps]();

	buff_disp_2[0] = 52 ;
	buff_disp_2[1] = 28 ;
	buff_disp_2[2] = 32 ;

	// class variables
	float class2 = 1;
	float buff_class2[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset2 = 1;

	uint8_t costAggr = aggregatedCostFromTopLeft2(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1 ,P2 ,direction, buff2, buff_disp_2, &min_disp2, &pixel_2, class2, buff_class2, reset2) ;

	//costAggr = 15 + min(30+8 , 28 , 32+8 , 28+32) -28
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset2);

}

// Test aggregation value of a non-border point and reset history, reset already computed
TEST(aggregatedCostFromTopLeft2Test, aggregationResetHistory_disp1) {

	uint8_t pixel_2 = 30 ;
	uint8_t min_disp2 = 26  ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8 ;
	uint8_t P2 = 32 ;
	Direction direction={1,1};
	uint8_t invalid_value = 57;

	int row = 1;
	int col = 1;
	int disp = 1;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff2 = new uint8_t[nb_cols*nb_disps]();
	uint8_t * buff_disp_2 = new uint8_t[nb_disps]();

	buff_disp_2[0] = 52 ;
	buff_disp_2[1] = 28 ;
	buff_disp_2[2] = 32 ;

	// class variables
	float class2 = 1;
	float buff_class2[5] = {1,1,1,1,1}; // unused
	float reset2 = 0;

	uint8_t costAggr = aggregatedCostFromTopLeft2(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1 ,P2 ,direction, buff2, buff_disp_2, &min_disp2, &pixel_2, class2, buff_class2, reset2) ;

	//costAggr = 15 + min(30+8 , 28 , 32+8 , 28+32) -26
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(0,reset2);

}

// Test aggregation value of a non-border point and reset history, update reset
TEST(aggregatedCostFromTopLeft2Test, aggregationResetHistory_disp0) {

	uint8_t pixel_2 = 30 ;
	uint8_t min_disp2 = 26  ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8 ;
	uint8_t P2 = 32 ;
	Direction direction={1,1};
	uint8_t invalid_value = 57;

	int row = 1;
	int col = 1;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff2 = new uint8_t[nb_cols*nb_disps]();
	uint8_t * buff_disp_2 = new uint8_t[nb_disps]();

	buff_disp_2[0] = 52 ;
	buff_disp_2[1] = 28 ;
	buff_disp_2[2] = 32 ;

	// class variables
	float class2 = 1;
	float buff_class2[5] = {2,1,1,1,1}; // different for  piecewise optimization
	float reset2 = 1;

	uint8_t costAggr = aggregatedCostFromTopLeft2(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1 ,P2 ,direction, buff2, buff_disp_2, &min_disp2, &pixel_2, class2, buff_class2, reset2) ;

	//costAggr = 15 + min(30+8 , 28 , 32+8 , 28+32) -26
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(0,reset2);

}

// Test aggregation value after an invalid point
TEST(aggregatedCostFromTopLeft2Test, aggregationAfterInvalidPoint) {

	uint8_t pixel_2 ;
	uint8_t min_disp2;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={1,1};
	uint8_t invalid_value = 57;

	int row = 1;
	int col = 1;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff2 = new uint8_t[nb_cols*nb_disps]();
	uint8_t * buff_disp_2 = new uint8_t[nb_disps]();
	buff_disp_2[0] = 57 ;
	buff_disp_2[1] = 57 ;
	buff_disp_2[2] = 57 ;

	// class variables
	float class2 = 1;
	float buff_class2[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset2 = 1;

	uint8_t costAggr = aggregatedCostFromTopLeft2(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1 ,P2 ,direction, buff2, buff_disp_2, &min_disp2, &pixel_2, class2, buff_class2, reset2) ;

	EXPECT_EQ(57,min_disp2);
	//costAggr = 15 + min(52 , 28+8,32+28) - 28
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset2);

}

// Test function aggregatedCostFromTopLeft3

// Test aggregation on cost Volume border
TEST(aggregatedCostFromTopLeft3Test, borderValue) {

	uint8_t min_disp3 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={1,-1};
	uint8_t invalid_value = 57;

	int row = 1;
	int col = 4;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff3 = new uint8_t[nb_cols*nb_disps]();

	// class variables
	float class3 = 1;
	float buff_class3[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset3 = 1;

	uint8_t costAggr = aggregatedCostFromTopLeft3(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff3, &min_disp3, class3, buff_class3, reset3) ;

	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset3);
}

// Test aggregation on cost Volume border
TEST(aggregatedCostFromTopLeft3Test, borderValueResetHistory) {

	uint8_t min_disp3 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={1,-1};
	uint8_t invalid_value = 57;

	int row = 1;
	int col = 4;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff3 = new uint8_t[nb_cols*nb_disps]();

	// class variables
	float class3 = 1;
	float buff_class3[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset3 = 0;

	uint8_t costAggr = aggregatedCostFromTopLeft3(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff3, &min_disp3, class3, buff_class3, reset3) ;

	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(0,reset3);// unchanged
}

// Test minimum value of the previous point
TEST(aggregatedCostFromTopLeft3Test, findMinDisp) {

	uint8_t min_disp3 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={1,-1};
	uint8_t invalid_value = 57;

	int row = 1;
	int col = 3;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff3 = new uint8_t[nb_cols*nb_disps]();
	buff3[12] = 14 ;
	buff3[13] = 28 ;
	buff3[14] = 20 ;

	// class variables
	float class3 = 1;
	float buff_class3[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset3 = 1;

	uint8_t costAggr = aggregatedCostFromTopLeft3(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff3, &min_disp3, class3, buff_class3, reset3) ;

	//pixel_0 = min(14,28,20)
	EXPECT_EQ(14,min_disp3);
	//costAggr = 15 + min(14 , 28 + 8 , 14 + 32) - 14 = 15
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset3);

}

// Test aggregation value of a non-border point
TEST(aggregatedCostFromTopLeft3Test, aggregation) {

	uint8_t min_disp3 = 14 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={1,-1};
	uint8_t invalid_value = 57;

	int row = 1;
	int col = 3;
	int disp = 1;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff3 = new uint8_t[nb_cols*nb_disps]();
	buff3[12] = 14 ;
	buff3[13] = 28 ;
	buff3[14] = 20 ;

	// class variables
	float class3 = 1;
	float buff_class3[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset3 = 1;

	uint8_t costAggr = aggregatedCostFromTopLeft3(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff3, &min_disp3, class3, buff_class3, reset3) ;

	//costAggr = 15 + min ( 14+8 , 28 , 20 + 8 , 14+32) - 14 = 15+22-14 = 23
	EXPECT_EQ(23,costAggr);
	EXPECT_EQ(1,reset3);

}

// Test aggregation value of a non-border point and reset history, reset already computed
TEST(aggregatedCostFromTopLeft3Test, aggregationResetHistory_disp1) {

	uint8_t min_disp3 = 14 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={1,-1};
	uint8_t invalid_value = 57;

	int row = 1;
	int col = 3;
	int disp = 1;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff3 = new uint8_t[nb_cols*nb_disps]();
	buff3[12] = 14 ;
	buff3[13] = 28 ;
	buff3[14] = 20 ;

	// class variables
	float class3 = 1;
	float buff_class3[5] = {1,1,1,1,1}; // unused
	float reset3 = 0;

	uint8_t costAggr = aggregatedCostFromTopLeft3(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff3, &min_disp3, class3, buff_class3, reset3) ;

	//costAggr = 15 + min ( 14+8 , 28 , 20 + 8 , 14+32) - 14 = 15+22-14 = 23
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(0,reset3);

}


// Test aggregation value of a non-border point and reset history, updated reset
TEST(aggregatedCostFromTopLeft3Test, aggregationResetHistory_disp0) {

	uint8_t min_disp3 = 14 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={1,-1};
	uint8_t invalid_value = 57;

	int row = 1;
	int col = 3;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff3 = new uint8_t[nb_cols*nb_disps]();
	buff3[12] = 14 ;
	buff3[13] = 28 ;
	buff3[14] = 20 ;

	// class variables
	float class3 = 1;
	float buff_class3[5] = {1,1,1,1,2}; // same class for no piecewise optimization
	float reset3 = 1;

	uint8_t costAggr = aggregatedCostFromTopLeft3(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff3, &min_disp3, class3, buff_class3, reset3) ;

	//costAggr = 15 + min ( 14+8 , 28 , 20 + 8 , 14+32) - 14 = 15+22-14 = 23
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(0,reset3);

}

// Test aggregation value after an invalid point
TEST(aggregatedCostFromTopLeft3Test, aggregationAfterInvalidPoint) {

	uint8_t min_disp3 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={1,-1};
	uint8_t invalid_value = 57;

	int row = 1;
	int col = 3;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff3 = new uint8_t[nb_cols*nb_disps]();
	buff3[12] = 57 ;
	buff3[13] = 57 ;
	buff3[14] = 57 ;

	// class variables
	float class3 = 1;
	float buff_class3[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset3 = 1;

	uint8_t costAggr = aggregatedCostFromTopLeft3(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff3, &min_disp3, class3, buff_class3, reset3) ;


	EXPECT_EQ(57,min_disp3);
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset3);

}

// Test function aggregatedCostFromBottomRight4

// Test aggregation on cost Volume border
TEST(aggregatedCostFromBottomRight4Test, borderValue) {

	uint8_t pixel_4 ;
	uint8_t min_disp4 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={0,-1};
	uint8_t invalid_value = 57;

	int row = 4;
	int col = 4;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff4 = new uint8_t[nb_disps]();

	// class variables
	float class4 = 1;
	float buff_class4 = 1; // same class for no piecewise optimization
	float reset4 = 1;

	uint8_t costAggr = aggregatedCostFromBottomRight4(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2,direction, buff4, &min_disp4, &pixel_4, class4, buff_class4, reset4) ;

	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset4);
}

// Test aggregation on cost Volume border
TEST(aggregatedCostFromBottomRight4Test, borderValueResetHistory) {

	uint8_t pixel_4 ;
	uint8_t min_disp4 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={0,-1};
	uint8_t invalid_value = 57;

	int row = 4;
	int col = 4;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff4 = new uint8_t[nb_disps]();

	// class variables
	float class4 = 1;
	float buff_class4 = 2; // different class for  piecewise optimization
	float reset4 = 1; // doesnt matter, disp ==0

	uint8_t costAggr = aggregatedCostFromBottomRight4(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2,direction, buff4, &min_disp4, &pixel_4, class4, buff_class4, reset4) ;

	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset4); // was not changed
}

// Test minimum value of the previous point
TEST(aggregatedCostFromBottomRight4Test, findMinDisp) {

	uint8_t pixel_4 ;
	uint8_t min_disp4 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={0,-1};
	uint8_t invalid_value = 57;

	int row = 3;
	int col = 3;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff4 = new uint8_t[nb_disps]();
	buff4[0] = 14 ;
	buff4[1] = 28 ;
	buff4[2] = 20 ;

	// class variables
	float class4 = 1;
	float buff_class4 = 1; // same class for no piecewise optimization
	float reset4 = 1;

	uint8_t costAggr = aggregatedCostFromBottomRight4(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2,direction, buff4, &min_disp4, &pixel_4, class4, buff_class4, reset4) ;

	//pixel_0 = min(14,28,20)
	EXPECT_EQ(14,min_disp4);
	//costAggr = 15 + min(14 , 28 + 8 , 14 + 32) - 14 = 15
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset4);

}

// Test aggregation value of a non-border point
TEST(aggregatedCostFromBottomRight4Test, aggregation) {

	uint8_t pixel_4 = 14 ;
	uint8_t min_disp4 = 14 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={0,-1};
	uint8_t invalid_value = 57;

	int row = 3;
	int col = 3;
	int disp = 1;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff4 = new uint8_t[nb_disps]();
	buff4[0] = 36 ;
	buff4[1] = 28 ;
	buff4[2] = 20 ;

	// class variables
	float class4 = 1;
	float buff_class4 = 1; // same class for no piecewise optimization
	float reset4 = 1;

	uint8_t costAggr = aggregatedCostFromBottomRight4(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2,direction, buff4, &min_disp4, &pixel_4, class4, buff_class4, reset4) ;

	//costAggr = 15 + min ( 14+8 , 28 , 20 + 8 , 14+32) - 14 = 15+22-14 = 23
	EXPECT_EQ(23,costAggr);
	EXPECT_EQ(1,reset4);

}

// Test aggregation value of a non-border point, reset history, reset already computed
TEST(aggregatedCostFromBottomRight4Test, aggregationResetHistory_disp1) {

	uint8_t pixel_4 = 14 ;
	uint8_t min_disp4 = 14 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={0,-1};
	uint8_t invalid_value = 57;

	int row = 3;
	int col = 3;
	int disp = 1;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff4 = new uint8_t[nb_disps]();
	buff4[0] = 36 ;
	buff4[1] = 28 ;
	buff4[2] = 20 ;

	// class variables
	float class4 = 1;
	float buff_class4 = 2; // different class for no piecewise optimization, disp > 0, doesnt matter
	float reset4 = 0;

	uint8_t costAggr = aggregatedCostFromBottomRight4(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2,direction, buff4, &min_disp4, &pixel_4, class4, buff_class4, reset4) ;

	//costAggr = 15 + 0 * min ( 14+8 , 28 , 20 + 8 , 14+32) - 14 = 15+22-14 = 23
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(0,reset4);

}

// Test aggregation value of a non-border point, reset history, update reset
TEST(aggregatedCostFromBottomRight4Test, aggregationResetHistory_disp0) {

	uint8_t pixel_4 = 14 ;
	uint8_t min_disp4 = 14 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={0,-1};
	uint8_t invalid_value = 57;

	int row = 3;
	int col = 3;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff4 = new uint8_t[nb_disps]();
	buff4[0] = 36 ;
	buff4[1] = 28 ;
	buff4[2] = 20 ;

	// class variables
	float class4 = 1;
	float buff_class4 = 2; // different class for no piecewise optimization
	float reset4 = 1;

	uint8_t costAggr = aggregatedCostFromBottomRight4(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2,direction, buff4, &min_disp4, &pixel_4, class4, buff_class4, reset4) ;

	//costAggr = 15 + 0 * min ( 14+8 , 28 , 20 + 8 , 14+32) - 14 = 15+22-14 = 23
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(0,reset4);

}

// Test aggregation value after an invalid point
TEST(aggregatedCostFromBottomRight4Test, aggregationAfterInvalidPoint) {

	uint8_t pixel_4 ;
	uint8_t min_disp4 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={0,-1};
	uint8_t invalid_value = 57;

	int row = 3;
	int col = 3;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff4 = new uint8_t[nb_disps]();
	buff4[0] = 57 ;
	buff4[1] = 57 ;
	buff4[2] = 57 ;

	// class variables
	float class4 = 1;
	float buff_class4 = 1; // same class for no piecewise optimization
	float reset4 = 1;

	uint8_t costAggr = aggregatedCostFromBottomRight4(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2,direction, buff4, &min_disp4, &pixel_4, class4, buff_class4, reset4) ;

	EXPECT_EQ(57,min_disp4);
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset4);

}

// Test fucntion aggregatedCostFromBottomRight5

// Test aggregation on cost Volume border
TEST(aggregatedCostFromBottomRight5Test, borderValue) {

	uint8_t pixel_5 ;
	uint8_t min_disp5 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={-1,0};
	uint8_t invalid_value = 57;

	int row = 4;
	int col = 4;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff5 = new uint8_t[nb_cols*nb_disps]();

	// class variables
	float class5 = 1;
	float buff_class5[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset5 = 1;


	uint8_t costAggr = aggregatedCostFromBottomRight5(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff5, &min_disp5, &pixel_5, class5, buff_class5, reset5) ;

	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset5);
}

// Test aggregation on cost Volume border
TEST(aggregatedCostFromBottomRight5Test, borderValueResetHistory) {

	uint8_t pixel_5 ;
	uint8_t min_disp5 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={-1,0};
	uint8_t invalid_value = 57;

	int row = 4;
	int col = 4;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff5 = new uint8_t[nb_cols*nb_disps]();

	// class variables
	float class5 = 1;
	float buff_class5[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset5 = 0;


	uint8_t costAggr = aggregatedCostFromBottomRight5(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff5, &min_disp5, &pixel_5, class5, buff_class5, reset5) ;

	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(0,reset5); // unchanged
}

// Test minimum value of the previous point
TEST(aggregatedCostFromBottomRight5Test, findMinDisp) {

	uint8_t pixel_5 ;
	uint8_t min_disp5 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={-1,0};
	uint8_t invalid_value = 57;

	int row = 3;
	int col = 3;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff5 = new uint8_t[nb_cols*nb_disps]();
	buff5[9] = 17 ;
	buff5[10] = 45 ;
	buff5[11] = 32 ;

	// class variables
	float class5 = 1;
	float buff_class5[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset5 = 1;

	uint8_t costAggr = aggregatedCostFromBottomRight5(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff5, &min_disp5, &pixel_5, class5, buff_class5, reset5) ;

	EXPECT_EQ(17,min_disp5);
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset5);

}

// Test aggregation value of a non-border point
TEST(aggregatedCostFromBottomRight5Test, aggregation) {

	uint8_t pixel_5 = 30 ;
	uint8_t min_disp5 = 28  ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8 ;
	uint8_t P2 = 32 ;
	Direction direction={-1,0};
	uint8_t invalid_value = 57;

	int row = 3;
	int col = 3;
	int disp = 1;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff5 = new uint8_t[nb_disps]();
	buff5[9] = 52 ;
	buff5[10] = 28 ;
	buff5[11] = 32 ;

	// class variables
	float class5 = 1;
	float buff_class5[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset5 = 1;

	uint8_t costAggr = aggregatedCostFromBottomRight5(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff5, &min_disp5, &pixel_5, class5, buff_class5, reset5) ;
	//costAggr = 15 + min(30+8 , 28 , 32+8 , 28+32) -28
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset5);

}

// Test aggregation value of a non-border point, reset histroy,reset already computed
TEST(aggregatedCostFromBottomRight5Test, aggregationResetHistory_disp1) {

	uint8_t pixel_5 = 30 ;
	uint8_t min_disp5 = 26  ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8 ;
	uint8_t P2 = 32 ;
	Direction direction={-1,0};
	uint8_t invalid_value = 57;

	int row = 3;
	int col = 3;
	int disp = 1;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff5 = new uint8_t[nb_disps]();
	buff5[9] = 52 ;
	buff5[10] = 28 ;
	buff5[11] = 32 ;

	// class variables
	float class5 = 1;
	float buff_class5[5] = {1,1,1,1,1}; // unused
	float reset5 = 0;

	uint8_t costAggr = aggregatedCostFromBottomRight5(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff5, &min_disp5, &pixel_5, class5, buff_class5, reset5) ;
	//costAggr = 15 + min(30+8 , 28 , 32+8 , 28+32) -26
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(0,reset5);

}

// Test aggregation value of a non-border point, reset history, update reset
TEST(aggregatedCostFromBottomRight5Test, aggregationResetHistory_disp0) {

	uint8_t pixel_5 = 30 ;
	uint8_t min_disp5 = 26  ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8 ;
	uint8_t P2 = 32 ;
	Direction direction={-1,0};
	uint8_t invalid_value = 57;

	int row = 3;
	int col = 3;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff5 = new uint8_t[nb_disps]();
	buff5[9] = 52 ;
	buff5[10] = 28 ;
	buff5[11] = 32 ;

	// class variables
	float class5 = 1;
	float buff_class5[5] = {1,1,1,2,1}; // different class
	float reset5 = 1;

	uint8_t costAggr = aggregatedCostFromBottomRight5(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff5, &min_disp5, &pixel_5, class5, buff_class5, reset5) ;
	//costAggr = 15 + min(30+8 , 28 , 32+8 , 28+32) -26
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(0,reset5);

}

// Test aggregation value after an invalid point
TEST(aggregatedCostFromBottomRight5Test, aggregationAfterInvalidPoint) {

	uint8_t pixel_5 ;
	uint8_t min_disp5 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={-1,0};
	uint8_t invalid_value = 57;

	int row = 3;
	int col = 3;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff5 = new uint8_t[nb_cols*nb_disps]();
	buff5[9] = 57 ;
	buff5[10] = 57 ;
	buff5[11] = 57 ;

	// class variables
	float class5 = 1;
	float buff_class5[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset5 = 1;

	uint8_t costAggr = aggregatedCostFromBottomRight5(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff5, &min_disp5, &pixel_5, class5, buff_class5, reset5) ;

	EXPECT_EQ(57,min_disp5);
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset5);

}

// Test function aggregatedCostFromBottomRight6

// Test aggregation on cost Volume border
TEST(aggregatedCostFromBottomRight6Test, borderValue) {

	uint8_t pixel_6 ;
	uint8_t min_disp6 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={-1,-1};
	uint8_t invalid_value = 57;

	int row = 4;
	int col = 4;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff6 = new uint8_t[nb_cols*nb_disps]();
	uint8_t * buff_disp_6 = new uint8_t[nb_disps]();

	// class variables
	float class6 = 1;
	float buff_class6[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset6 = 1;

	uint8_t costAggr = aggregatedCostFromBottomRight6(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff6, buff_disp_6, &min_disp6, &pixel_6, class6, buff_class6, reset6) ;

	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset6);
}

// Test aggregation on cost Volume border
TEST(aggregatedCostFromBottomRight6Test, borderValueResetHistory) {

	uint8_t pixel_6 ;
	uint8_t min_disp6 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={-1,-1};
	uint8_t invalid_value = 57;

	int row = 4;
	int col = 4;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff6 = new uint8_t[nb_cols*nb_disps]();
	uint8_t * buff_disp_6 = new uint8_t[nb_disps]();

	// class variables
	float class6 = 1;
	float buff_class6[5] = {1,1,1,1,1}; // unused
	float reset6 = 0;

	uint8_t costAggr = aggregatedCostFromBottomRight6(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff6, buff_disp_6, &min_disp6, &pixel_6, class6, buff_class6, reset6) ;

	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(0,reset6);// unchanged
}

// Test minimum value of the previous point
TEST(aggregatedCostFromBottomRight6Test, findMinDisp) {

	uint8_t pixel_6 ;
	uint8_t min_disp6;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={-1,-1};
	uint8_t invalid_value = 57;

	int row = 3;
	int col = 3;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff6 = new uint8_t[nb_cols*nb_disps]();
	uint8_t * buff_disp_6 = new uint8_t[nb_disps]();
	buff_disp_6[0] = 52 ;
	buff_disp_6[1] = 28 ;
	buff_disp_6[2] = 32 ;

	// class variables
	float class6 = 1;
	float buff_class6[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset6 = 1;

	uint8_t costAggr = aggregatedCostFromBottomRight6(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff6, buff_disp_6, &min_disp6, &pixel_6, class6, buff_class6, reset6) ;

	EXPECT_EQ(28,min_disp6);
	//costAggr = 15 + min(52 , 28+8,32+28) - 28
	EXPECT_EQ(23,costAggr);
	EXPECT_EQ(1,reset6);

}

// Test aggregation value of a non-border point
TEST(aggregatedCostFromBottomRight6Test, aggregation) {

	uint8_t pixel_6 = 30 ;
	uint8_t min_disp6 = 28  ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8 ;
	uint8_t P2 = 32 ;
	Direction direction={-1,-1};
	uint8_t invalid_value = 57;

	int row = 3;
	int col = 3;
	int disp = 1;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff6 = new uint8_t[nb_cols*nb_disps]();
	uint8_t * buff_disp_6 = new uint8_t[nb_disps]();

	buff_disp_6[0] = 52 ;
	buff_disp_6[1] = 28 ;
	buff_disp_6[2] = 32 ;

	// class variables
	float class6 = 1;
	float buff_class6[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset6 = 1;

	uint8_t costAggr = aggregatedCostFromBottomRight6(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff6, buff_disp_6, &min_disp6, &pixel_6, class6, buff_class6, reset6) ;

	//costAggr = 15 + min(30+8 , 28 , 32+8 , 28+32) -28
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset6);

}

// Test aggregation value of a non-border point, reset history, with reset already computed
TEST(aggregatedCostFromBottomRight6Test, aggregationResetHistory_disp1) {

	uint8_t pixel_6 = 30 ;
	uint8_t min_disp6 = 28  ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8 ;
	uint8_t P2 = 32 ;
	Direction direction={-1,-1};
	uint8_t invalid_value = 57;

	int row = 3;
	int col = 3;
	int disp = 1;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff6 = new uint8_t[nb_cols*nb_disps]();
	uint8_t * buff_disp_6 = new uint8_t[nb_disps]();

	buff_disp_6[0] = 52 ;
	buff_disp_6[1] = 28 ;
	buff_disp_6[2] = 32 ;

	// class variables
	float class6 = 1;
	float buff_class6[5] = {1,1,1,1,1}; // unused
	float reset6 = 0;

	uint8_t costAggr = aggregatedCostFromBottomRight6(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff6, buff_disp_6, &min_disp6, &pixel_6, class6, buff_class6, reset6) ;

	//costAggr = 15 + min(30+8 , 28 , 32+8 , 28+32) -28
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(0,reset6);

}

// Test aggregation value of a non-border point, reset history, update reset
TEST(aggregatedCostFromBottomRight6Test, aggregationResetHistory_disp0) {

	uint8_t pixel_6 = 30 ;
	uint8_t min_disp6 = 28  ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8 ;
	uint8_t P2 = 32 ;
	Direction direction={-1,-1};
	uint8_t invalid_value = 57;

	int row = 3;
	int col = 3;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff6 = new uint8_t[nb_cols*nb_disps]();
	uint8_t * buff_disp_6 = new uint8_t[nb_disps]();

	buff_disp_6[0] = 52 ;
	buff_disp_6[1] = 28 ;
	buff_disp_6[2] = 32 ;

	// class variables
	float class6 = 1;
	float buff_class6[5] = {1,1,1,1,2}; // same class for no piecewise optimization ( diagonal from left, but array is  reversed)
	float reset6 = 1;

	uint8_t costAggr = aggregatedCostFromBottomRight6(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff6, buff_disp_6, &min_disp6, &pixel_6, class6, buff_class6, reset6) ;

	//costAggr = 15 + min(30+8 , 28 , 32+8 , 28+32) -28
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(0,reset6);

}

// Test aggregation value after an invalid point
TEST(aggregatedCostFromBottomRight6Test, aggregationAfterInvalidPoint) {

	uint8_t pixel_6 ;
	uint8_t min_disp6;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={-1,-1};
	uint8_t invalid_value = 57;

	int row = 3;
	int col = 3;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff6 = new uint8_t[nb_cols*nb_disps]();
	uint8_t * buff_disp_6 = new uint8_t[nb_disps]();
	buff_disp_6[0] = 57 ;
	buff_disp_6[1] = 57 ;
	buff_disp_6[2] = 57 ;

	// class variables
	float class6 = 1;
	float buff_class6[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset6 = 1;

	uint8_t costAggr = aggregatedCostFromBottomRight6(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff6, buff_disp_6, &min_disp6, &pixel_6, class6, buff_class6, reset6) ;

	EXPECT_EQ(57,min_disp6);
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset6);

}

// Test function aggregatedCostFromBottomRight7

// Test aggregation on cost Volume border
TEST(aggregatedCostFromBottomRight7Test, borderValue) {

	uint8_t min_disp7 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={-1,1};
	uint8_t invalid_value = 57;

	int row = 4;
	int col = 0;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff7 = new uint8_t[nb_cols*nb_disps]();

	// class variables
	float class7 = 1;
	float buff_class7[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset7 = 1;

	uint8_t costAggr = aggregatedCostFromBottomRight7(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff7, &min_disp7, class7, buff_class7, reset7) ;

	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset7);
}


// Test aggregation on cost Volume border
TEST(aggregatedCostFromBottomRight7Test, borderValueResetHistory) {

	uint8_t min_disp7 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={-1,1};
	uint8_t invalid_value = 57;

	int row = 4;
	int col = 0;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff7 = new uint8_t[nb_cols*nb_disps]();

	// class variables
	float class7 = 1;
	float buff_class7[5] = {1,1,1,1,1}; // unused
	float reset7 = 0;

	uint8_t costAggr = aggregatedCostFromBottomRight7(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff7, &min_disp7, class7, buff_class7, reset7) ;

	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(0,reset7);// unchanged
}


// Test minimum value of the previous point
TEST(aggregatedCostFromBottomRight7Test, findMinDisp) {

	uint8_t min_disp7 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={-1,1};
	uint8_t invalid_value = 57;

	int row = 3;
	int col = 1;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff7 = new uint8_t[nb_cols*nb_disps]();
	buff7[0] = 14 ;
	buff7[1] = 28 ;
	buff7[2] = 20 ;

	// class variables
	float class7 = 1;
	float buff_class7[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset7 = 1;

	uint8_t costAggr = aggregatedCostFromBottomRight7(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff7, &min_disp7, class7, buff_class7, reset7) ;
	//pixel_0 = min(14,28,20)
	EXPECT_EQ(14,min_disp7);
	//costAggr = 15 + min(14 , 28 + 8 , 14 + 32) - 14 = 15
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset7);

}

// Test aggregation value of a non-border point
TEST(aggregatedCostFromBottomRight7Test, aggregation) {

	uint8_t min_disp7 = 14 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={-1,1};
	uint8_t invalid_value = 57;

	int row = 3;
	int col = 1;
	int disp = 1;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff7 = new uint8_t[nb_cols*nb_disps]();
	buff7[0] = 14 ;
	buff7[1] = 28 ;
	buff7[2] = 20 ;

	// class variables
	float class7 = 1;
	float buff_class7[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset7 = 1;

	uint8_t costAggr = aggregatedCostFromBottomRight7(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff7, &min_disp7, class7, buff_class7, reset7) ;

	//costAggr = 15 + min ( 14+8 , 28 , 20 + 8 , 14+32) - 14 = 15+22-14 = 23
	EXPECT_EQ(23,costAggr);
	EXPECT_EQ(1,reset7);

}

// Test aggregation value of a non-border point, reset history, reset already computed
TEST(aggregatedCostFromBottomRight7Test, aggregationResetHistory_disp1) {

	uint8_t min_disp7 = 14 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={-1,1};
	uint8_t invalid_value = 57;

	int row = 3;
	int col = 1;
	int disp = 1;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff7 = new uint8_t[nb_cols*nb_disps]();
	buff7[0] = 14 ;
	buff7[1] = 28 ;
	buff7[2] = 20 ;

	// class variables
	float class7 = 1;
	float buff_class7[5] = {1,1,1,1,1}; // unused
	float reset7 = 0;

	uint8_t costAggr = aggregatedCostFromBottomRight7(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff7, &min_disp7, class7, buff_class7, reset7) ;

	//costAggr = 15 +  0 * min ( 14+8 , 28 , 20 + 8 , 14+32) - 14 = 15
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(0,reset7);

}

// Test aggregation value of a non-border point, reset history, update reset
TEST(aggregatedCostFromBottomRight7Test, aggregationResetHistory_disp0) {

	uint8_t min_disp7 = 14 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={-1,1};
	uint8_t invalid_value = 57;

	int row = 3;
	int col = 1;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff7 = new uint8_t[nb_cols*nb_disps]();
	buff7[0] = 14 ;
	buff7[1] = 28 ;
	buff7[2] = 20 ;

	// class variables
	float class7 = 1;
	float buff_class7[5] = {2,1,1,1,1}; // same class for no piecewise optimization
	float reset7 = 1;

	uint8_t costAggr = aggregatedCostFromBottomRight7(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff7, &min_disp7, class7, buff_class7, reset7) ;

	//costAggr = 15 + 0 *  min ( 14+8 , 28 , 20 + 8 , 14+32) - 14 = 15
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(0,reset7);

}

// Test aggregation value after an invalid point
TEST(aggregatedCostFromBottomRight7Test, aggregationAfterInvalidPoint) {

	uint8_t min_disp7 ;

	uint8_t pixelCost = 15 ;
	uint8_t P1 = 8;
	uint8_t P2 = 32;
	Direction direction={-1,1};
	uint8_t invalid_value = 57;

	int row = 3;
	int col = 1;
	int disp = 0;

	int nb_rows = 5;
	int nb_cols = 5;
	int nb_disps = 3;
	uint8_t * buff7 = new uint8_t[nb_cols*nb_disps]();
	buff7[0] = 57 ;
	buff7[1] = 57 ;
	buff7[2] = 57 ;

	// class variables
	float class7 = 1;
	float buff_class7[5] = {1,1,1,1,1}; // same class for no piecewise optimization
	float reset7 = 1;

	uint8_t costAggr = aggregatedCostFromBottomRight7(pixelCost, row, col, disp, invalid_value, nb_rows, nb_cols, nb_disps,
		P1, P2, direction, buff7, &min_disp7, class7, buff_class7, reset7) ;

	EXPECT_EQ(57,min_disp7);
	EXPECT_EQ(15,costAggr);
	EXPECT_EQ(1,reset7);


}

//Global Test of sgm function, input float cost volume: aggregation value from 8 directions on a middle point of cost volume
//with invalid middle point and without over-counting correction

TEST(sgmTestFloat, TestMiddleValueInvalid){

    int nb_row,nb_col,nb_disp;
	nb_row = 3;
	nb_col = 3;
	nb_disp = 3;
	float P1, P2, invalid_value;
	P1 = 1.0;
	P2 = 3.0;
	invalid_value = 5.0;

	bool overcounting = false;
	bool cost_paths = false;

	CostVolumes cvs;
	float cv_in[27]={1,0.5,0.3,-0.2,0.4,0.1,-0.8,0.6,0.2,-0.3,-0.7,0.4,5.0,5.0,5.0,0.4,0.4,-0.1,-0.8,0.7,-0.6,0.1,-0.1,0.1,0.5,0.2,0.4};

	// method : constant
	float p1[9*8]={P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1};
	float p2[9*8]={P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2};
    int direction[2*8]={0,1, 1,0, 1,1, 1,-1, 0,-1, -1,0, -1,-1, -1,1};

    float segmentation[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1}; // no piecewise optimization

	cvs = sgm(cv_in, p1, p2, direction, nb_row, nb_col, nb_disp,invalid_value, segmentation, cost_paths, overcounting);
	//middle point must stay invalid, equal to invalid value
	EXPECT_FLOAT_EQ(40.0,cvs.cost_volume[13]);
}

//Global Test of sgm function, input float cost volume: aggregation value from 8 directions on a middle point of cost volume
//with invalid middle point and with over-counting correction

TEST(sgmTestFloat, TestMiddleValueInvalidOvercounting){

    int nb_row,nb_col,nb_disp;
	nb_row = 3;
	nb_col = 3;
	nb_disp = 3;
	float P1, P2, invalid_value;
	P1 = 1.0;
	P2 = 3.0;
	invalid_value = 5.0;

	bool overcounting = true;
	bool cost_paths = false;

	CostVolumes cvs;
	float cv_in[27]={1,0.5,0.3,-0.2,0.4,0.1,-0.8,0.6,0.2,-0.3,-0.7,0.4,5.0,5.0,5.0,0.4,0.4,-0.1,-0.8,0.7,-0.6,0.1,-0.1,0.1,0.5,0.2,0.4};

	// method : constant
	float p1[9*8]={P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1};
	float p2[9*8]={P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2};
    int direction[2*8]={0,1, 1,0, 1,1, 1,-1, 0,-1, -1,0, -1,-1, -1,1};

    // segmentation
    float segmentation[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1}; // no piecewise optimization

	cvs = sgm(cv_in, p1, p2, direction, nb_row, nb_col, nb_disp,invalid_value, segmentation, cost_paths, overcounting);
	//middle point must stay invalid, equal to invalid value
	EXPECT_FLOAT_EQ(5.0,cvs.cost_volume[13]);
}

//Global Test of sgm function: aggregation value from 8 directions on a middle point of cost

TEST(sgmTestFloat, TestMiddleValue){

    int nb_row,nb_col,nb_disp;
	nb_row = 3;
	nb_col = 3;
	nb_disp = 3;
	float P1, P2, invalid_value;
	P1 = 1.0;
	P2 = 3.0;
	invalid_value = 5.0;

	CostVolumes cvs;
	float cv_in[27]={1,0.5,0.3,-0.2,0.4,0.1,-0.8,0.6,0.2,-0.3,-0.7,0.4,0.8,-0.7,0.2,0.4,0.4,-0.1,-0.8,0.7,-0.6,0.1,-0.1,0.1,0.5,0.2,0.4};

	bool overcounting = false;
	bool cost_paths = false;

    // method : constant
	float p1[9*8]={P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1};
	float p2[9*8]={P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2};
    int direction[2*8]={0,1, 1,0, 1,1, 1,-1, 0,-1, -1,0, -1,-1, -1,1};

    // segmentation
    float segmentation[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1}; // no piecewise optimization

	cvs = sgm(cv_in, p1, p2, direction, nb_row, nb_col, nb_disp,invalid_value, segmentation, cost_paths, overcounting);
	//middle point must stay invalid, equal to invalid value
	// L0 = -0.7 + min(-0.7 ; -0.3+1; 0,4+1) - (-0.7) = - 0.7
	// L1 = -0.7 + min(0.4 ; -0.2+1; 0.1+1) - (-0.2) = - 0.1
	// L2 = -0.7 + min(0.5 ; 1+1 ; 0.3 +1) - 0.3 = -0.5
	// L3 = -0.7 + min(0.6 ; 0.2+1 ; -0.8 + 1) - (-0.8) = 0.3
	// L4 = -0.7 + min(0.4 ; 0.4+1 ; -0.1+1) - (-0.1) = - 0.2
	// L5 = -0.7 + min(-0.1 ; 0.1+1 ; 0.1+1) - (-0.1) = -0.7
	// L6 = -0.7 + min(0.2 ; 0.4 +1 ; 0.5 +1) - 0.2 = -0.7
	// L7 = -0.7 + min(0.7 ; -0.6+1 ; -0.8+1) - (-0.8) = 0.3

	EXPECT_FLOAT_EQ(-2.3 ,cvs.cost_volume[13]);

}

//Global Test of sgm function: aggregation value from 8 directions on a middle point of cost with overcounting

TEST(sgmTestFloat, TestMiddleValueOvercounting){

    int nb_row,nb_col,nb_disp;
	nb_row = 3;
	nb_col = 3;
	nb_disp = 3;
	float P1, P2, invalid_value;
	P1 = 1.0;
	P2 = 3.0;
	invalid_value = 5.0;

	CostVolumes cvs;
	float cv_in[27]={1,0.5,0.3,-0.2,0.4,0.1,-0.8,0.6,0.2,-0.3,-0.7,0.4,0.8,-0.7,0.2,0.4,0.4,-0.1,-0.8,0.7,-0.6,0.1,-0.1,0.1,0.5,0.2,0.4};

	bool overcounting = true;
	bool cost_paths = false;

    // method : constant
	float p1[9*8]={P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1};
	float p2[9*8]={P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2};
    int penalties[2*8]={0,1, 1,0, 1,1, 1,-1, 0,-1, -1,0, -1,-1, -1,1};

    // segmentation
    float segmentation[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1}; // no piecewise optimization

	cvs = sgm(cv_in, p1, p2, penalties, nb_row, nb_col, nb_disp,invalid_value, segmentation, cost_paths, overcounting);
	//middle point must stay invalid, equal to invalid value
	// L0 = -0.7 + min(-0.7 ; -0.3+1; 0,4+1) - (-0.7) = - 0.7
	// L1 = -0.7 + min(0.4 ; -0.2+1; 0.1+1) - (-0.2) = - 0.1
	// L2 = -0.7 + min(0.5 ; 1+1 ; 0.3 +1) - 0.3 = -0.5
	// L3 = -0.7 + min(0.6 ; 0.2+1 ; -0.8 + 1) - (-0.8) = 0.3
	// L4 = -0.7 + min(0.4 ; 0.4+1 ; -0.1+1) - (-0.1) = - 0.2
	// L5 = -0.7 + min(-0.1 ; 0.1+1 ; 0.1+1) - (-0.1) = -0.7
	// L6 = -0.7 + min(0.2 ; 0.4 +1 ; 0.5 +1) - 0.2 = -0.7
	// L7 = -0.7 + min(0.7 ; -0.6+1 ; -0.8+1) - (-0.8) = 0.3
	// Overcounting : -2.3 + (-0.7 x 7 ) = 2.6
	EXPECT_FLOAT_EQ(2.6 ,cvs.cost_volume[13]);

}

//Global Test of sgm function: aggregation value from 8 directions on a middle point of cost with cost path option

TEST(sgmTestFloat, TestMiddleValueCostPaths){

    int nb_row,nb_col,nb_disp;
	nb_row = 3;
	nb_col = 3;
	nb_disp = 3;
	float P1, P2, invalid_value;
	P1 = 1.0;
	P2 = 3.0;
	invalid_value = 5.0;

	CostVolumes cvs;
	float cv_in[27]={1,0.5,0.3,-0.2,0.4,0.1,-0.8,0.6,0.2,-0.3,-0.7,0.4,0.8,-0.7,0.2,0.4,0.4,-0.1,-0.8,0.7,-0.6,0.1,-0.1,0.1,0.5,0.2,0.4};

	bool overcounting = false;
	bool cost_paths = true;

    // method : constant
	float p1[9*8]={P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1};
	float p2[9*8]={P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2,P2};
    int direction[2*8]={0,1, 1,0, 1,1, 1,-1, 0,-1, -1,0, -1,-1, -1,1};

    // segmentation
    float segmentation[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1}; // no piecewise optimization

	cvs = sgm(cv_in, p1, p2, direction, nb_row, nb_col, nb_disp,invalid_value, segmentation, cost_paths, overcounting);
	//middle point must stay invalid, equal to invalid value
	// L0 = -0.7 + min(-0.7 ; -0.3+1; 0,4+1) - (-0.7) = - 0.7
	// L1 = -0.7 + min(0.4 ; -0.2+1; 0.1+1) - (-0.2) = - 0.1
	// L2 = -0.7 + min(0.5 ; 1+1 ; 0.3 +1) - 0.3 = -0.5
	// L3 = -0.7 + min(0.6 ; 0.2+1 ; -0.8 + 1) - (-0.8) = 0.3
	// L4 = -0.7 + min(0.4 ; 0.4+1 ; -0.1+1) - (-0.1) = - 0.2
	// L5 = -0.7 + min(-0.1 ; 0.1+1 ; 0.1+1) - (-0.1) = -0.7
	// L6 = -0.7 + min(0.2 ; 0.4 +1 ; 0.5 +1) - 0.2 = -0.7
	// L7 = -0.7 + min(0.7 ; -0.6+1 ; -0.8+1) - (-0.8) = 0.3
	// Overcounting : -2.3 + (-0.7 x 7 ) = 2.6
	EXPECT_EQ(1,cvs.cost_volume_min[32]);

}


int main(int argc, char **argv) {

	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}