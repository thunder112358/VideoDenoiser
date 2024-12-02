// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.


#ifndef ZOOM_H
#define ZOOM_H

#include <cmath>
#include <cstring>
#include <iostream>

#include "gaussian.h"
#include "bicubic_interpolation.h"

#define ZOOM_SIGMA_ZERO 0.6

/**
  *
  * Compute the size of a zoomed image from the zoom factor
  *
**/
void zoom_size(
    int nx,             // width of the original image
    int ny,             // height of the original image
    int *nxx,          // width of the zoomed image
    int *nyy,          // height of the zoomed image
    float factor = 0.5  // zoom factor between 0 and 1
)
{
    *nxx = (int)((float) nx * factor + 0.5);
    *nyy = (int)((float) ny * factor + 0.5);
}

/**
  *
  * Function to downsample the image
  *
**/
void zoom_out(
    const float *I,          // input image
    float *Iout,            // output image
    const int nx,           // image width
    const int ny,           // image height
    const float factor = 0.5 // zoom factor between 0 and 1
)
{
    std::cout << "zoom_out: Input dimensions: " << nx << "x" << ny << std::endl;
    
    int nxx, nyy;
    zoom_size(nx, ny, &nxx, &nyy, factor);
    std::cout << "zoom_out: Output dimensions: " << nxx << "x" << nyy << std::endl;

    try {
        // Allocate temporary buffer
        float *Is = new float[nx * ny];
        memcpy(Is, I, nx * ny * sizeof(float));

        // Compute the Gaussian sigma for smoothing
        const float sigma = ZOOM_SIGMA_ZERO * sqrt(1.0f/(factor*factor) - 1.0f);
        std::cout << "zoom_out: Smoothing with sigma=" << sigma << std::endl;

        // Pre-smooth the image
        gaussian(Is, nx, ny, sigma);

        // Re-sample the image
        std::cout << "zoom_out: Resampling image..." << std::endl;
        for (int i1 = 0; i1 < nyy; i1++) {
            for (int j1 = 0; j1 < nxx; j1++) {
                const float i2 = (float) i1 / factor;
                const float j2 = (float) j1 / factor;
                Iout[i1 * nxx + j1] = bicubic_interpolation_at(Is, j2, i2, nx, ny);
            }
        }

        delete[] Is;
        std::cout << "zoom_out: Complete" << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Exception in zoom_out: " << e.what() << std::endl;
        throw;
    }
}


/**
  *
  * Function to upsample the image
  *
**/
void zoom_in(
    const float *I, //input image
    float *Iout,    //output image
    int nx,         //width of the original image
    int ny,         //height of the original image
    int nxx,        //width of the zoomed image
    int nyy         //height of the zoomed image
)
{
    // compute the zoom factor
    const float factorx = ((float)nxx / nx);
    const float factory = ((float)nyy / ny);

    // re-sample the image using bicubic interpolation
	for (int i1 = 0; i1 < nyy; i1++)
	    for (int j1 = 0; j1 < nxx; j1++)
	    {
		float i2 =  (float) i1 / factory;
		float j2 =  (float) j1 / factorx;

		Iout[i1 * nxx + j1] = bicubic_interpolation_at(I, j2, i2, nx, ny); 
	    }
}

#endif
