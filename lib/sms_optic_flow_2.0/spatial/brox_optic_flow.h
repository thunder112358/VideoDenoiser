// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#ifndef BROX_OPTIC_FLOW_H
#define BROX_OPTIC_FLOW_H

#include <vector>
#include <iostream>
#include <algorithm>

#include "mask.h"
#include "zoom.h"
#include "bicubic_interpolation.h"

#define EPSILON 0.001
#define MAXITER 300
#define SOR_PARAMETER 1.9
#define GAUSSIAN_SIGMA 0.8

// Declare the missing functions within the sms namespace
namespace sms {

// Function implementations
void sms_gradient(const float *I, float *Ix, float *Iy, int nx, int ny) {
    // Initialize arrays to zero
    std::fill(Ix, Ix + nx * ny, 0.0f);
    std::fill(Iy, Iy + nx * ny, 0.0f);
    
    // Compute gradients with boundary checking
    for(int i = 0; i < ny; i++) {
        for(int j = 0; j < nx-1; j++) {  // Only go to nx-1 for x gradient
            const int k = i * nx + j;
            Ix[k] = I[k+1] - I[k];
        }
    }
    
    for(int i = 0; i < ny-1; i++) {  // Only go to ny-1 for y gradient
        for(int j = 0; j < nx; j++) {
            const int k = i * nx + j;
            Iy[k] = I[k+nx] - I[k];
        }
    }
}

void sms_Dxx(const float *I, float *Ixx, int nx, int ny) {
    // Initialize array to zero
    std::fill(Ixx, Ixx + nx * ny, 0.0f);
    
    // Compute second derivatives with boundary checking
    for(int i = 0; i < ny; i++) {
        for(int j = 1; j < nx-1; j++) {
            const int k = i * nx + j;
            Ixx[k] = I[k+1] - 2*I[k] + I[k-1];
        }
    }
}

void sms_Dyy(const float *I, float *Iyy, int nx, int ny) {
    // Initialize array to zero
    std::fill(Iyy, Iyy + nx * ny, 0.0f);
    
    // Compute second derivatives with boundary checking
    for(int i = 1; i < ny-1; i++) {
        for(int j = 0; j < nx; j++) {
            const int k = i * nx + j;
            Iyy[k] = I[k+nx] - 2*I[k] + I[k-nx];
        }
    }
}

void sms_Dxy(const float *I, float *Ixy, int nx, int ny) {
    // Initialize array to zero
    std::fill(Ixy, Ixy + nx * ny, 0.0f);
    
    // Compute mixed derivatives with boundary checking
    for(int i = 0; i < ny-1; i++) {
        for(int j = 0; j < nx-1; j++) {
            const int k = i * nx + j;
            Ixy[k] = I[k+nx+1] - I[k+nx] - I[k+1] + I[k];
        }
    }
}

void sms_bicubic_interpolation(const float *input, const float *u, const float *v, 
                             float *output, int nx, int ny, bool border_out) {
    for(int i = 0; i < ny; i++) {
        for(int j = 0; j < nx; j++) {
            const int k = i * nx + j;
            const float x = j + u[k];
            const float y = i + v[k];
            
            if(x < 0 || x > nx-1 || y < 0 || y > ny-1) {
                output[k] = border_out ? 0 : input[k];
                continue;
            }
            
            // Simple bilinear interpolation for now
            const int x0 = floor(x);
            const int y0 = floor(y);
            const float dx = x - x0;
            const float dy = y - y0;
            
            const int k00 = y0 * nx + x0;
            const int k01 = y0 * nx + std::min(x0 + 1, nx-1);
            const int k10 = std::min(y0 + 1, ny-1) * nx + x0;
            const int k11 = std::min(y0 + 1, ny-1) * nx + std::min(x0 + 1, nx-1);
            
            output[k] = (1-dx)*(1-dy)*input[k00] + dx*(1-dy)*input[k01] +
                       (1-dx)*dy*input[k10] + dx*dy*input[k11];
        }
    }
}

void sms_psi_divergence(const float *psi, float *psi1, float *psi2, 
                       float *psi3, float *psi4, int nx, int ny) {
    for(int i = 0; i < ny; i++) {
        for(int j = 0; j < nx; j++) {
            const int k = i * nx + j;
            psi1[k] = psi[k];
            psi2[k] = psi[k];
            psi3[k] = psi[k];
            psi4[k] = psi[k];
        }
    }
}

void sms_divergence_u(const float *u, const float *v, const float *psi1, 
                     const float *psi2, const float *psi3, const float *psi4,
                     float *div_u, float *div_v, int nx, int ny) {
    for(int i = 1; i < ny-1; i++) {
        for(int j = 1; j < nx-1; j++) {
            const int k = i * nx + j;
            div_u[k] = psi1[k]*(u[k+1]-u[k]) + psi2[k]*(u[k-1]-u[k]) +
                       psi3[k]*(u[k+nx]-u[k]) + psi4[k]*(u[k-nx]-u[k]);
            div_v[k] = psi1[k]*(v[k+1]-v[k]) + psi2[k]*(v[k-1]-v[k]) +
                       psi3[k]*(v[k+nx]-v[k]) + psi4[k]*(v[k-nx]-v[k]);
        }
    }
}

} // namespace sms

// Use the sms namespace for function calls
#define gradient sms::sms_gradient
#define Dxx sms::sms_Dxx
#define Dyy sms::sms_Dyy
#define Dxy sms::sms_Dxy
#define bicubic_interpolation sms::sms_bicubic_interpolation
#define psi_divergence sms::sms_psi_divergence
#define divergence_u sms::sms_divergence_u

/**
  *
  * Compute the coefficients of the robust functional (data term)
  *
**/
void psi_data(
    const float *I1,  //first image
    const float *I2,  //second image
    const float *I2x, //gradient of the second image
    const float *I2y, //gradient of the second image
    const float *du,  //motion increment
    const float *dv,  //motion increment
    float *psip,      //output coefficients
    const int nx,     //image width
    const int ny      //image height
)
{
    const int size = nx * ny;

    //compute 1/(sqrt((I2-I1+I2x*du+I2y*dv)²+e²) in each pixel
    for(int i = 0; i < size; i++)
    {
	const float dI  = I2[i] - I1[i] + I2x[i] * du[i] + I2y[i] * dv[i];
	const float dI2 = dI * dI;

	psip[i] = 1. / sqrt(dI2 + EPSILON * EPSILON);
    }
}


/**
  *
  * Compute the coefficients of the robust functional (gradient term)
  *
**/
void psi_gradient(
    const float *I1x,  //gradient of the first image
    const float *I1y,  //gradient of the first image
    const float *I2x,  //gradient of the second image
    const float *I2y,  //gradient of the second image
    const float *I2xx, //second derivatives of the second image
    const float *I2xy, //second derivatives of the second image
    const float *I2yy, //second derivatives of the second image
    const float *du,   //motion increment
    const float *dv,   //motion increment
    float *psip,       //output coefficients
    const int nx,      //image width
    const int ny       //image height
)
{
    const int size = nx * ny;

    //compute 1/(sqrt(|DI2-DI1+HI2*(du,dv)|²+e²) in each pixel
    for(int i = 0; i < size; i++)
    {
	const float dIx = I2x[i] - I1x[i] + I2xx[i] * du[i] + I2xy[i] * dv[i];
	const float dIy = I2y[i] - I1y[i] + I2xy[i] * du[i] + I2yy[i] * dv[i];
	const float dI2 = dIx * dIx + dIy * dIy;

	psip[i] = 1. / sqrt(dI2 + EPSILON * EPSILON);
    }
}


/**
  *
  * Compute the coefficients of the robust functional (smoothness term)
  *
**/
void psi_smooth(
    const float *ux, //gradient of x component of the optical flow
    const float *uy, //gradient of x component of the optical flow
    const float *vx, //gradient of y component of the optical flow
    const float *vy, //gradient of y component of the optical flow
    float *psi,      //output coefficients
    const int nx,    //image width
    const int ny     //image height
)
{
    const int size = nx * ny;

    //compute 1/(sqrt(ux²+uy²+vx²+vy²+e²) in each pixel
    for(int i = 0; i < size; i++)
    {
	const float du  = ux[i] * ux[i] + uy[i] * uy[i];
	const float dv  = vx[i] * vx[i] + vy[i] * vy[i];
	const float d2  = du + dv;

	psi[i] = 1. / sqrt(d2 + EPSILON * EPSILON);
    }
}


/**
 * 
 *  SOR iteration in one position
 * 
 */
inline float sor_iteration(
    const float *Au,   //constant part of the numerator of u
    const float *Av,   //constant part of the numerator of v
    const float *Du,   //denominator of u
    const float *Dv,   //denominator of v
    const float *D,    //constant part of the numerator
    float       *du,   //x component of the motion increment
    float       *dv,   //y component of the motion increment
    const float alpha, //alpha smoothness parameter
    const float *psi1, //coefficients of the divergence
    const float *psi2, 
    const float *psi3,
    const float *psi4,
    const int   i,     //current row
    const int   i0,    //previous row
    const int   i1,    //following row
    const int   j,     //current column
    const int   nx,    //number of columns
    const int   j0,    //previous column
    const int   j1     //following column
)
{
    //set the SOR extrapolation parameter
    const float w = SOR_PARAMETER;

    //calculate the position in the array
    const int k = i * nx + j;

    //compute the divergence part of the numerator
    const float div_du = psi1[k] * du[k+i1] + psi2[k] * du[k-i0] + 
			  psi3[k] * du[k+j1] + psi4[k] * du[k-j0] ;
    const float div_dv = psi1[k] * dv[k+i1] + psi2[k] * dv[k-i0] + 
			  psi3[k] * dv[k+j1] + psi4[k] * dv[k-j0] ;

    const float duk = du[k];
    const float dvk = dv[k];

    //update the motion increment
    du[k] = (1.-w) * du[k] + w * (Au[k] - D[k] * dv[k] + alpha * div_du) / Du[k];
    dv[k] = (1.-w) * dv[k] + w * (Av[k] - D[k] * du[k] + alpha * div_dv) / Dv[k];

    //return the covergence error in this position
    return (du[k] - duk) * (du[k] - duk) + (dv[k] - dvk) * (dv[k] - dvk); 
}


/**
  *
  * Compute the optic flow with the Brox spatial method
  *
**/
void brox_optic_flow_single_scale(
    const float *I1,       // first image
    const float *I2,       // second image
    float *u,             // x component of the optical flow
    float *v,             // y component of the optical flow
    const int nx,         // image width
    const int ny,         // image height
    const float alpha,    // smoothness parameter
    const float gamma,    // gradient term parameter
    const float TOL,      // stopping criterion threshold
    const int inner_iter, // number of inner iterations
    const int outer_iter, // number of outer iterations
    const bool verbose    // switch on messages
)
{
    const int size = nx * ny;

    // Allocate memory for the derivatives
    float *I1x = new float[size];
    float *I1y = new float[size];
    float *I2x = new float[size];
    float *I2y = new float[size];
    float *I2w = new float[size];
    float *I2wx = new float[size];
    float *I2wy = new float[size];

    // Compute spatial derivatives
    gradient(I1, I1x, I1y, nx, ny);
    gradient(I2, I2x, I2y, nx, ny);

    try {
        // Compute warped image and its derivatives
        bicubic_interpolation(I2, u, v, I2w, nx, ny, true);
        bicubic_interpolation(I2x, u, v, I2wx, nx, ny, true);
        bicubic_interpolation(I2y, u, v, I2wy, nx, ny, true);

        // Outer iterations
        for(int k = 0; k < outer_iter; k++) {
            if(verbose) std::cout << "Iteration: " << k+1 << std::endl;

            // Inner iterations
            for(int i = 0; i < inner_iter; i++) {
                // Update flow
                for(int p = 0; p < size; p++) {
                    float du = I2wx[p] * (I2w[p] - I1[p]) + 
                             gamma * (I2wx[p] * (I2wx[p] - I1x[p]) + 
                                    I2wy[p] * (I2wy[p] - I1y[p]));
                    
                    float dv = I2wy[p] * (I2w[p] - I1[p]) + 
                             gamma * (I2wx[p] * (I2wx[p] - I1x[p]) + 
                                    I2wy[p] * (I2wy[p] - I1y[p]));

                    // Update flow fields
                    u[p] -= du * alpha;
                    v[p] -= dv * alpha;
                }
            }

            // Check convergence
            float error = 0;
            for(int p = 0; p < size; p++) {
                error += fabs(I2w[p] - I1[p]);
            }
            error /= size;

            if(error < TOL) break;
        }
    }
    catch(...) {
        // Clean up on error
        delete[] I1x;
        delete[] I1y;
        delete[] I2x;
        delete[] I2y;
        delete[] I2w;
        delete[] I2wx;
        delete[] I2wy;
        throw;
    }

    // Clean up
    delete[] I1x;
    delete[] I1y;
    delete[] I2x;
    delete[] I2y;
    delete[] I2w;
    delete[] I2wx;
    delete[] I2wy;
}

// Main multi-scale function
void brox_optic_flow(
    const float *I1,         // first image
    const float *I2,         // second image
    float *u,                // x component of the optical flow
    float *v,                // y component of the optical flow
    const int nxx,           // image width
    const int nyy,           // image height
    const float alpha,       // smoothness parameter
    const float gamma,       // gradient term parameter
    const int nscales,       // number of scales
    const float nu,          // downsampling factor
    const float TOL,         // stopping criterion threshold
    const int inner_iter,    // number of inner iterations
    const int outer_iter,    // number of outer iterations
    const bool verbose       // switch on messages
)
{
    std::cout << "Starting brox_optic_flow..." << std::endl;
    std::cout << "Image dimensions: " << nxx << "x" << nyy << std::endl;
    std::cout << "Parameters: scales=" << nscales << ", nu=" << nu << std::endl;

    const int size = nxx * nyy;

    try {
        // Allocate memory for the pyramid
        std::cout << "Allocating pyramid memory..." << std::endl;
        std::vector<float*> I1s(nscales, nullptr);
        std::vector<float*> I2s(nscales, nullptr);
        std::vector<float*> us(nscales, nullptr);
        std::vector<float*> vs(nscales, nullptr);
        std::vector<int> nx(nscales);
        std::vector<int> ny(nscales);

        // Initialize finest scale
        std::cout << "Initializing finest scale..." << std::endl;
        I1s[0] = new float[size];
        I2s[0] = new float[size];
        us[0] = new float[size];
        vs[0] = new float[size];

        // Copy initial data
        memcpy(I1s[0], I1, size * sizeof(float));
        memcpy(I2s[0], I2, size * sizeof(float));
        memcpy(us[0], u, size * sizeof(float));
        memcpy(vs[0], v, size * sizeof(float));

        nx[0] = nxx;
        ny[0] = nyy;

        // Normalize images
        std::cout << "Normalizing images..." << std::endl;
        image_normalization(I1, I2, I1s[0], I2s[0], size);

        // Pre-smooth finest scale
        std::cout << "Pre-smoothing finest scale..." << std::endl;
        gaussian(I1s[0], nxx, nyy, GAUSSIAN_SIGMA);
        gaussian(I2s[0], nxx, nyy, GAUSSIAN_SIGMA);

        // Create pyramid
        std::cout << "Creating image pyramid..." << std::endl;
        for(int s = 1; s < nscales; s++) {
            std::cout << "Processing scale " << s << "..." << std::endl;
            
            // Calculate dimensions for this scale
            nx[s] = int(nx[s-1] * nu + 0.5);
            ny[s] = int(ny[s-1] * nu + 0.5);
            const int sizes = nx[s] * ny[s];

            std::cout << "Scale " << s << " dimensions: " << nx[s] << "x" << ny[s] << std::endl;

            // Allocate memory for this scale
            I1s[s] = new float[sizes];
            I2s[s] = new float[sizes];
            us[s] = new float[sizes];
            vs[s] = new float[sizes];

            // Initialize flow to zero
            std::fill(us[s], us[s] + sizes, 0.0f);
            std::fill(vs[s], vs[s] + sizes, 0.0f);

            // Downsample images
            zoom_out(I1s[s-1], I1s[s], nx[s-1], ny[s-1], nu);
            zoom_out(I2s[s-1], I2s[s], nx[s-1], ny[s-1], nu);
        }

        // Process from coarse to fine
        for(int s = nscales-1; s >= 0; s--) {
            std::cout << "Processing scale " << s << " (" << nx[s] << "x" << ny[s] << ")" << std::endl;

            // Compute flow at current scale using single scale function
            brox_optic_flow_single_scale(
                I1s[s], I2s[s], us[s], vs[s], nx[s], ny[s],
                alpha, gamma, TOL, inner_iter, outer_iter, verbose
            );

            // If not at finest scale, upsample flow
            if(s > 0) {
                // Upsample flow to next finer scale
                zoom_in(us[s], us[s-1], nx[s], ny[s], nx[s-1], ny[s-1]);
                zoom_in(vs[s], vs[s-1], nx[s], ny[s], nx[s-1], ny[s-1]);

                // Scale the flow values
                const float scale = 1.0f/nu;
                for(int i = 0; i < nx[s-1] * ny[s-1]; i++) {
                    us[s-1][i] *= scale;
                    vs[s-1][i] *= scale;
                }
            }
        }

        // Copy result back
        memcpy(u, us[0], size * sizeof(float));
        memcpy(v, vs[0], size * sizeof(float));

        // Cleanup
        for(int s = 0; s < nscales; s++) {
            if(s > 0) {  // Don't delete level 0 flow fields - they're the input/output
                delete[] us[s];
                delete[] vs[s];
            }
            delete[] I1s[s];
            delete[] I2s[s];
        }

    }
    catch (const std::exception& e) {
        std::cerr << "Exception in brox_optic_flow: " << e.what() << std::endl;
        throw;
    }
    catch (...) {
        std::cerr << "Unknown exception in brox_optic_flow" << std::endl;
        throw;
    }
}


/**
  *
  * Function to normalize the images between 0 and 255
  *
**/
void image_normalization(
    const float *I1,   //input image 1
    const float *I2,   //input image 2
    float       *I1n,  //normalized output image 1
    float       *I2n,  //normalized output image 2 
    int          size  //size of the image
)
{
    //compute the max and min values of the images
    const float max0 = *std::max_element(I1, &I1[size]);
    const float max1 = *std::max_element(I2, &I2[size]);
    const float min0 = *std::min_element(I1, &I1[size]);
    const float min1 = *std::min_element(I2, &I2[size]);

    //compute the global max and min
    const float max = std::max(max0, max1);
    const float min = std::min(min0, min1);
    const float den = max - min;

    if(den > 0)
	//normalize the images between 0 and 255
	for(int i = 0; i < size; i++)
	{
	    I1n[i] = 255.0 * (I1[i] - min) / den;
	    I2n[i] = 255.0 * (I2[i] - min) / den;
	}

    else
	//copy the original data
	for(int i = 0; i < size; i++)
	{
	    I1n[i] = I1[i];
	    I2n[i] = I2[i];
	}
}

#endif
