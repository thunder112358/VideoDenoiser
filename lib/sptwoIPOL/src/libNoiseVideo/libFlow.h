/*----------------------------------------------------------------------------

Copyright (c) 2016-2018 A. Buades

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

----------------------------------------------------------------------------*/

/*****
REMARK:

THIS CODE INCLUDES THE FUNCTIONS FROM THE FOLLOWING FILES 
tvl1flow_lib.c, mask.c, zoom.c and bicubic_interpolation.c, 
by Javier Sanchez,
available at http://www.ipol.im/pub/art/2013/26/

****/

#ifndef DUAL_TVL1_OPTIC_FLOW_H
#define DUAL_TVL1_OPTIC_FLOW_H

#include "../library/libBasic.h"


namespace libUSTGFLOW
{

#define MAX_ITERATIONS 300
#define PRESMOOTHING_SIGMA 0.8
#define GRAD_IS_ZERO 1E-10


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
    
    
    
struct flow_params {
    float tau;     // time step
    float lambda;  // weight parameter for the data term
    float theta;   // weight parameter for (u - v)²
    int   nscales; // number of scales
    float zfactor; // factor for building the image piramid
    int   warps;   // number of warpings per scale
    float epsilon; // tolerance for numerical convergence
    bool  verbose;  // enable/disable the verbose mode
    int iflagMedian;
};


void check_flow_reciprocity(float *uflow, float *vflow, float *uflow2, float *vflow2, float *omask,float fThreshold, int width, int height);




void divergence(
    const float *v1, // x component of the vector field
    const float *v2, // y component of the vector field
    float *div,      // output divergence
    const int nx,    // image width
    const int ny     // image height
);


float bicubic_interpolation_at(
    const float *input, //image to be interpolated
    const float  uu,    //x component of the vector field
    const float  vv,    //y component of the vector field
    const int    nx,    //image width
    const int    ny,    //image height
    bool         border_out //if true, return zero outside the region
);



/**
 *
 * Function to compute the optical flow using multiple scales
 *
 **/
void Dual_TVL1_optic_flow_multiscale(
    float *I0,           // source image
    float *I1,           // target image
    float *u1,           // x component of the optical flow
    float *u2,           // y component of the optical flow
    const int   nxx,     // image width
    const int   nyy,     // image height
    struct flow_params &params);




}

#endif//DUAL_TVL1_OPTIC_FLOW_H
