// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Jorge Jara <jjara@dcc.uchile.cl>, Jose Delpiano
// <jdelpian@uandes.cl> and Mauricio Cerda <mauriciocerda@med.uchile.cl>.
// All rights reserved.

#ifndef CLG_OF_C
#define CLG_OF_C

/**
 * @file clg_of.c
 * @brief Implementation of CLG optical flow for 2D grayscale images.
 * @author Jorge Jara <jjara@dcc.uchile.cl>, Jose Delpiano <jdelpian@uandes.cl>
 */
#include "clg_of.h"
#include <math.h>   // Used for "exp".
#include <stdio.h>  // Used for "printf".
#include <stdlib.h> // Used for "free".
#include <string.h>  // For memcpy and memset

// Add these function prototypes before calcMSCLG_OF
static void clg_zoom_out(const double *I, double *Iout, const int nx, const int ny, const double factor);
int clg_gaussian(double* image, int width, int height, double sigma);

/**
 * @brief boundaryCondition
 *
 * Neumann boundary conditions (derivatives are set to zero).
 *
 * Parameters:
 *
 * u               Pointer to the optical flow horizontal component
 *                 matrix/array.
 * v               Pointer to the optical flow vertical component matrix/array.
 * nRows           Number of rows of the optical flow arrrays.
 * nCols           Number of columns of the optical flow arays.
 *
 */
double boundaryCondition(double **u, double **v, int nRows, int nCols) {

    int i, j;
    double error = 0.0;

    // First and last rows.
    for (j=1; j<nCols-1; j++) {

        i = 0;
        error += (u[i+1][j]-u[i][j]) * (u[i+1][j]-u[i][j]);
        error += (v[i+1][j]-v[i][j]) * (v[i+1][j]-v[i][j]);

        u[i][j] = u[i+1][j];
        v[i][j] = v[i+1][j];

        i = nRows-1;
        error += (u[i-1][j]-u[i][j]) * (u[i-1][j]-u[i][j]);
        error += (v[i-1][j]-v[i][j]) * (v[i-1][j]-v[i][j]);

        u[i][j] = u[i-1][j];
        v[i][j] = v[i-1][j];
    }

    // First and last columns.
    for (i=1; i<nRows-1; i++) {

        j = 0;
        error += (u[i][j+1]-u[i][j]) * (u[i][j+1]-u[i][j]);
        error += (v[i][j+1]-v[i][j]) * (v[i][j+1]-v[i][j]);

        u[i][j] = u[i][j+1];
        v[i][j] = v[i][j+1];

        j = nCols-1;
        error += (u[i][j-1]-u[i][j]) * (u[i][j-1]-u[i][j]);
        error += (v[i][j-1]-v[i][j]) * (v[i][j-1]-v[i][j]);

        u[i][j] = u[i][j-1];
        v[i][j] = v[i][j-1];
    }

    // Corners.
    i=0; j=0;
    error += (u[i+1][j+1]-u[i][j]) * (u[i+1][j+1]-u[i][j]);
    error += (v[i+1][j+1]-v[i][j]) * (v[i+1][j+1]-v[i][j]);

    u[i][j] = u[i+1][j+1];
    v[i][j] = v[i+1][j+1];

    j = nCols-1;
    error += (u[i+1][j-1]-u[i][j]) * (u[i+1][j-1]-u[i][j]);
    error += (v[i+1][j-1]-v[i][j]) * (v[i+1][j-1]-v[i][j]);

    u[i][j] = u[i+1][j-1];
    v[i][j] = v[i+1][j-1];

    i=nRows-1; j=0;
    error += (u[i-1][j+1]-u[i][j]) * (u[i-1][j+1]-u[i][j]);
    error += (v[i-1][j+1]-v[i][j]) * (v[i-1][j+1]-v[i][j]);

    u[i][j] = u[i-1][j+1];
    v[i][j] = v[i-1][j+1];

    j = nCols-1;
    error += (u[i-1][j-1]-u[i][j]) * (u[i-1][j-1]-u[i][j]);
    error += (v[i-1][j-1]-v[i][j]) * (v[i-1][j-1]-v[i][j]);

    u[i][j] = u[i-1][j-1];
    v[i][j] = v[i-1][j-1];

    return error;
}


/**
 * @brief SOR_at
 *
 * SOR iteration at location (i,j).
 * This method return the new values.
 *
 * Parameters:
 *
 * u               Pointer to the optical flow horizontal component
 *                 matrix/array.
 * v               Pointer to the optical flow vertical component matrix/array.
 * J[JROWS][JCOLS] Pointer to the array that stores the computed (and possibly
 *                 smoothed) derivatives for each image/optical flow pixel.
 * i               Row of the location to compute.
 * j               Column of the location to compute.
 * alpha           Optical flow global smoothing coefficient.
 * wFactor         Relaxation parameter (if w=1.0 then the function becomes the
 *                 Gauss-Seidel method).
 *
 */
double SOR_at(double **u,
              double **v,
              double **J[JROWS][JCOLS],
              int i,
              int j,
              double alpha,
              double wFactor) {

    double h2a, numerator, denominator, ulocal, vlocal, error;

    // h, could be less than 1.0.
    h2a = 1.0 / alpha;

    // SOR formula.
    numerator   = u[i][j-1] + u[i-1][j] + u[i][j+1] + u[i+1][j] -
                  h2a * (J[0][1][i][j] * v[i][j]+J[0][2][i][j]);

    denominator = 4.0 + h2a * J[0][0][i][j];

    ulocal      = (1.0-wFactor) * u[i][j] + wFactor * numerator / denominator;

    numerator   = v[i][j-1] + v[i-1][j] + v[i][j+1] + v[i+1][j] -
                  h2a * (J[0][1][i][j] * ulocal+J[1][2][i][j]);

    denominator = 4.0 + h2a * J[1][1][i][j];

    vlocal      = (1.0-wFactor) * v[i][j] + wFactor * numerator / denominator;

    error       = (ulocal - u[i][j]) * (ulocal - u[i][j]);
    error      += (vlocal - v[i][j]) * (vlocal - v[i][j]);

    u[i][j]     = ulocal;
    v[i][j]     = vlocal;

    return error;
}


/**
 * @brief relaxPointwiseCoupledGaussSeidel
 *
 * Pointwise coupled Gauss-Seidel relaxation iteration for CLG-OF equations.
 * Each call to this function updates the current value of the solution,
 * u[1..m][1..n], v[1..m][1..n], using the motion tensor
 * J[JROWS][JCOLS][1..m][1..n].
 * Neumann boundary conditions are used (derivatives are set to zero).
 * The return value is the total error of the current iteration.
 *
 * Parameters:
 *
 * u               Pointer to the optical flow horizontal component
 *                 matrix/array.
 * v               Pointer to the optical flow vertical component matrix/array.
 * J[JROWS][JCOLS] Pointer to the array that stores the computed (and possibly
 *                 smoothed) derivatives for each image/optical flow pixel.
 * nRows           Number of rows of the optical flow arrrays.
 * nCols           Number of columns of the optical flow arays.
 * alpha           Optical flow global smoothing coefficient.
 *
 */
double relaxPointwiseCoupledGaussSeidel(double **u,
                                        double **v,
                                        double **J[JROWS][JCOLS],
                                        int nRows,
                                        int nCols,
                                        double alpha) {
    int i, j;
    double hx, hy, hx2, hy2, discr, detU, detV, error, ulocal, vlocal;

    // h, could be less than 1.0.
    hx    = 1.0;
    hy    = hx;
    hx2   = hx * hx;
    hy2   = hy * hy;
    error = 0.0;

    for (i=1; i<nRows-1; i++) {
        for (j=1; j<nCols-1; j++) {
            // Point-wise Coupled Gauss-Seidel formula.

            // Coupled system for two variables solved by Cramer's rule.
            // This is the discriminant of that system
            discr = lin2by2det(alpha/hx2*2.0 + alpha/hy2*2.0 + J[0][0][i][j],
                               J[0][1][i][j],
                               J[0][1][i][j],
                               alpha/hx2*2.0 + alpha/hy2*2.0 + J[1][1][i][j]);

            if (abs(discr) > EPS) {
                // Determinant replacing first column to solve
                // the 2x2 system by using Cramer's rule.
                detU = lin2by2det(alpha / hx2 * (u[i-1][j] + u[i+1][j])
                                  + alpha / hy2 * (u[i][j-1] + u[i][j+1])
                                  - J[0][2][i][j],
                                    J[0][1][i][j],
                                    alpha / hx2 * (v[i-1][j] + v[i+1][j])
                                  + alpha / hy2 * (v[i][j-1] + v[i][j+1])
                                  - J[1][2][i][j],
                                    alpha / hx2 * 2.0 + alpha / hy2 * 2.0
                                  + J[1][1][i][j]);

                // Determinant replacing second column to solve
                // the 2x2 system by using Cramer's rule.
                detV = lin2by2det(alpha / hx2 * 2.0 + alpha / hy2 * 2.0
                                  + J[0][0][i][j],
                                    alpha / hx2 * (u[i-1][j] + u[i+1][j])
                                  + alpha / hy2 * (u[i][j-1] + u[i][j+1])
                                  - J[0][2][i][j],
                                    J[0][1][i][j],
                                    alpha / hx2 * (v[i-1][j] + v[i+1][j])
                                  + alpha / hy2 * (v[i][j-1] + v[i][j+1])
                                  - J[1][2][i][j]);

                // Division of two discriminants (Cramer's rule).
                ulocal = detU / discr;
                vlocal = detV / discr;

                error += (ulocal-u[i][j])*(ulocal-u[i][j]);
                error += (vlocal-v[i][j])*(vlocal-v[i][j]);

                u[i][j]=ulocal;
                v[i][j]=vlocal;

                // Normal Gauss-Seidel iteration.
            } else {

                // w_Factor=1.0, thus a Gauss-Seidel iteration.
                error += SOR_at(u, v, J, i, j, alpha, 1.0);

            } // End if.
        } // End columns.
    } // End rows.

    error += boundaryCondition(u, v, nRows, nCols);

    return sqrt(error / (nRows*nCols));;
}


/**
 * @brief relaxSOR
 *
 * SOR relaxation iteration for CLG-OF equations.
 * Each call to this function updates the current value of the solution,
 * u[1..m][1..n], v[1..m][1..n], using the motion tensor
 * J[JROWS][JCOLS][1..m][1..n].
 * Neumann boundary conditions are used (derivatives are set to zero).
 *
 * Parameters:
 *
 * u               Pointer to the optical flow horizontal component
 *                 matrix/array.
 * v               Pointer to the optical flow vertical component matrix/array.
 * J[JROWS][JCOLS] Pointer to the array that stores the computed (and possibly
 *                 smoothed) derivatives for each image/optical flow pixel.
 * nRows           Number of rows of the optical flow arrrays.
 * nCols           Number of columns of the optical flow arays.
 * alpha           Optical flow global smoothing coefficient.
 * wFactor         Relaxation parameter (if w=1.0 then the function is the
 *                 Gauss-Seidel method).
 *
 */
double relaxSOR(double **u,
                double **v,
                double **J[JROWS][JCOLS],
                int nRows,
                int nCols,
                double alpha,
                double wFactor) {
    int i, j;
    double error = 0.0;

    for (i=1; i<nRows-1; i++)
        for (j=1; j<nCols-1; j++)
            error += SOR_at(u, v, J, i, j, alpha, wFactor);

    error += boundaryCondition(u, v, nRows, nCols);

    return sqrt(error / (nRows*nCols));
}


/**
 * @brief calcCLG_OF
 *
 * Main CLG-optical flow (CLG-OF) computation function.
 *
 * Parameters:
 *
 * image1          Pointer to the first image (the "previous" time frame).
 * image2          Pointer to the second image (the "current" time frame).
 * uOut            Pointer to the horizontal component of the CLG-OF solution.
 * vOut            Pointer to the vertical component of the CLG-OF solution.
 * nRows           Number of image rows (same for the CLG-OF vector field).
 * nCols           Number of image columns (same for the CLG-OF vector field).
 * iterations      Number of iterations for iterative solution.
 * alpha           Global smoothing coefficient of the CLG-OF.
 * rho             Local spatio-temporal smoothing coefficient of the CLG-OF.
 * wFactor         SOR relaxation factor, between 0 and 2.
 * verbose         Display/hide messages to the stdout.
 * coupledMode     Iteration type. 1->Pointwise-Coupled Gauss-Seidel, 0->SOR.
 *
 */
int calcCLG_OF(double* image1,
               double* image2,
               double* uOut,
               double* vOut,
               int nRows,
               int nCols,
               int iterations,
               double alpha,
               double rho,
               double wFactor,
               int verbose,
               int coupledMode) {

    if (verbose) {
        printf("calc_clg\n");
        printf("  setting up variables\n");
    }

    int i=0, j=0;

    // h, could be less than 1.0.
    double h = 1.0;

    // Matrix to vector.
    double **prevFrame, **currFrame;
    prevFrame = pMatrix(nRows, nCols);
    currFrame = pMatrix(nRows, nCols);

    for (i=0; i<nRows; i++) {
        for (j=0; j<nCols; j++) {
            prevFrame[i][j] = image1[j + i*nCols];
            currFrame[i][j] = image2[j + i*nCols];
        }
    }

    if (verbose)
        printf("  allocating memory for arrays\n");

    double **u, **v;
    double **dfdx, **dfdxw;
    double **dfdy, **dfdyw;
    double **dfdt, **dfdtw;

    u = pMatrix(nRows, nCols);
    v = pMatrix(nRows, nCols);

    // Derivatives and their warped versions.
    dfdx = pMatrix(nRows, nCols);
    dfdy = pMatrix(nRows, nCols);
    dfdt = pMatrix(nRows, nCols);

    for (i=0; i<nRows; i++) {
        for (j=0; j<nCols; j++) {
            u[i][j] = 0.0;
            v[i][j] = 0.0;
            dfdt[i][j] = currFrame[i][j] - prevFrame[i][j];
        }
    }

    if (verbose)
        printf("  allocating memory for derivatives matrices\n");

    double **J[JROWS][JCOLS];

    // Because of symmetry, only the upper part is allocated.
    int k, l;
    for (k=0; k<JROWS; k++)
        for (l=k; l<JCOLS; l++)
            J[k][l] = pMatrix(nRows, nCols);

    // Spatial derivatives obtention.
    computeDerivatives(prevFrame, currFrame, dfdx, dfdy, nRows, nCols, verbose);
    // Compute J tensor.
    computeJTensor(dfdx, dfdy, dfdt, J, nRows, nCols);

    if (verbose)
        printf("  local spatio temporal smoothing\n");

    if (rho > 0) {

        k=0, l=0;
        for (k=0; k<JROWS; k++)
            for (l=k; l<JCOLS; l++)
                matrixSmooth(J[k][l], nRows, nCols, rho);

    }

    if (iterations == 0)
        iterations = (int) (nRows * nCols / 8.0);

    if (verbose)
        printf("  performing %i relax iterations\n", iterations);

    int count = 0;
    double error = 1000000;
    double convergenceError = 0.0;

    for (count=0; count<iterations && convergenceError*1.01 < error; count++) {

        if (count > 0)
            error = convergenceError;

        if (coupledMode == 1) {

            if (verbose && count % 50 == 0 && count > 0)
                printf("  iteration %d/%d (P-C Gauss-Seidel), error=%f\n",
                       count, iterations, error);

            convergenceError = relaxPointwiseCoupledGaussSeidel(u, v, J,
                                                                nRows, nCols,
                                                                alpha);
        } else {
            if (verbose && count % 50 == 0 && count > 0)
                printf("  iteration %d/%d (SOR), error=%f\n",
                       count, iterations, error);

            convergenceError = relaxSOR(u, v, J, nRows, nCols, alpha, wFactor);
        }
    }

    // Show debug information.
    if (verbose)
        printf("  filling output after %d iterations, error=%f\n", count, error);

    // Fill output variables.
    for (i=0; i<nRows; i++) {
        for (j=0; j<nCols; j++) {
            uOut[j + i*nCols] += u[i][j];
            vOut[j + i*nCols] += v[i][j];
        }
    }

    // Free memory.
    if (verbose)
        printf("  freeing memory\n");

    freePmatrix(u, nRows);
    freePmatrix(v, nRows);

    for (k=0; k<JROWS; k++)
        for (l=k; l<JCOLS; l++)
            freePmatrix(J[k][l], nRows);

    freePmatrix(prevFrame, nRows);
    freePmatrix(currFrame, nRows);
    freePmatrix(dfdx,  nRows);
    freePmatrix(dfdy,  nRows);
    freePmatrix(dfdt,  nRows);

    if (verbose)
        printf("calc_clg: done\n");
    return 1;
}


/**
 * @brief calcMSCLG_OF
 *
 * Main Multiscale CLG-optical flow (CLG-OF) computation function.
 *
 * Parameters:
 *
 * image1          Pointer to the first image (the "previous" time frame).
 * image2          Pointer to the second image (the "current" time frame).
 * uOut            Pointer to the horizontal component of the CLG-OF solution.
 * vOut            Pointer to the vertical component of the CLG-OF solution.
 * nRows           Number of image rows (same for the CLG-OF vector field).
 * nCols           Number of image columns (same for the CLG-OF vector field).
 * iterations      Number of iterations for iterative solution.
 * alpha           Global smoothing coefficient of the CLG-OF.
 * rho             Local spatio-temporal smoothing coefficient of the CLG-OF.
 * sigma           Standard deviation for an optional gaussian smoothing kernel,
 *                 applied to the input images prior to the CLG-OF computation.
 * wFactor         SOR relaxation factor, between 0 and 2.
 * nScales         Total number of scales (at least 1).
 * scaleFactor     Downsampling factor, between 0.5 and 1.0.
 * coupledMode     Iteration type. 1->Pointwise-Coupled Gauss-Seidel, 0->SOR.
 * verbose         Display/hide messages to the stdout.
 *
 */
int calcMSCLG_OF(double* image1,
                 double* image2,
                 double* uOut,
                 double* vOut,
                 int nRows,
                 int nCols,
                 int iterations,
                 double alpha,
                 double rho,
                 double sigma,
                 double wFactor,
                 int nScales,
                 const double scaleFactor,
                 int coupledMode,
                 int verbose) {

    // Input validation
    if (!image1 || !image2 || !uOut || !vOut) {
        printf("Error: Null input pointers\n");
        return 0;
    }

    if (nRows <= 0 || nCols <= 0 || nScales <= 0) {
        printf("Error: Invalid dimensions or scales (rows=%d, cols=%d, scales=%d)\n", 
               nRows, nCols, nScales);
        return 0;
    }

    if (verbose) {
        printf("calcMSCLG_OF\n");
        printf("  parameters: nScales=%d, scaleFactor=%f\n", nScales, scaleFactor);
    }

    // Allocate memory for scale parameters
    double **I1s = (double**)malloc(nScales * sizeof(double*));
    double **I2s = (double**)malloc(nScales * sizeof(double*));
    double **us = (double**)malloc(nScales * sizeof(double*));
    double **vs = (double**)malloc(nScales * sizeof(double*));
    int *nxx = (int*)malloc(nScales * sizeof(int));
    int *nyy = (int*)malloc(nScales * sizeof(int));

    if (!I1s || !I2s || !us || !vs || !nxx || !nyy) {
        printf("Error: Failed to allocate scale arrays\n");
        free(I1s);
        free(I2s);
        free(us);
        free(vs);
        free(nxx);
        free(nyy);
        return 0;
    }

    // Initialize first scale
    const int size = nCols * nRows;
    I1s[0] = (double*)malloc(size * sizeof(double));
    I2s[0] = (double*)malloc(size * sizeof(double));
    
    if (!I1s[0] || !I2s[0]) {
        printf("Error: Failed to allocate first scale images\n");
        free(I1s[0]);
        free(I2s[0]);
        free(I1s);
        free(I2s);
        free(us);
        free(vs);
        free(nxx);
        free(nyy);
        return 0;
    }

    // Copy input data
    memcpy(I1s[0], image1, size * sizeof(double));
    memcpy(I2s[0], image2, size * sizeof(double));
    us[0] = uOut;
    vs[0] = vOut;
    nxx[0] = nCols;
    nyy[0] = nRows;

    // Pre-smoothing at finest scale
    if (verbose) printf("  apply gaussian smoothing for each frame\n");
    
    if (!clg_gaussian(I1s[0], nCols, nRows, sigma) ||
        !clg_gaussian(I2s[0], nCols, nRows, sigma)) {
        printf("Error: Gaussian smoothing failed\n");
        free(I1s[0]);
        free(I2s[0]);
        free(I1s);
        free(I2s);
        free(us);
        free(vs);
        free(nxx);
        free(nyy);
        return 0;
    }

    // Create the scales
    for (int s = 1; s < nScales; s++) {
        if (verbose) printf("  scales %d init\n", s);

        // Calculate dimensions for this scale
        nxx[s] = (int)(nxx[s-1] * scaleFactor);
        nyy[s] = (int)(nyy[s-1] * scaleFactor);

        // Ensure minimum dimensions
        if (nxx[s] < 32 || nyy[s] < 32) {
            printf("Warning: Scale %d too small (%dx%d), stopping at scale %d\n", 
                   s, nxx[s], nyy[s], s-1);
            nScales = s;
            break;
        }

        const int sizes = nxx[s] * nyy[s];
        
        // Allocate memory for this scale
        I1s[s] = (double*)malloc(sizes * sizeof(double));
        I2s[s] = (double*)malloc(sizes * sizeof(double));
        us[s] = (double*)malloc(sizes * sizeof(double));
        vs[s] = (double*)malloc(sizes * sizeof(double));

        if (!I1s[s] || !I2s[s] || !us[s] || !vs[s]) {
            printf("Error: Memory allocation failed for scale %d\n", s);
            // Clean up all allocated memory
            for (int j = 0; j <= s; j++) {
                free(I1s[j]);
                free(I2s[j]);
                if (j > 0) {
                    free(us[j]);
                    free(vs[j]);
                }
            }
            free(I1s);
            free(I2s);
            free(us);
            free(vs);
            free(nxx);
            free(nyy);
            return 0;
        }

        // Initialize flow fields to zero
        memset(us[s], 0, sizes * sizeof(double));
        memset(vs[s], 0, sizes * sizeof(double));

        // Compute the zoom from the previous finer scale
        if (verbose) {
            printf("  computing zoom for scale %d\n", s);
            printf("  source dims: %dx%d, target dims: %dx%d\n", 
                   nxx[s-1], nyy[s-1], nxx[s], nyy[s]);
        }

        // Verify memory and dimensions before zoom
        if (!I1s[s-1] || !I1s[s] || !I2s[s-1] || !I2s[s]) {
            printf("Error: Invalid pointers for zoom operation at scale %d\n", s);
            // ... cleanup code ...
            return 0;
        }

        clg_zoom_out(I1s[s-1], I1s[s], nxx[s-1], nyy[s-1], scaleFactor);
        clg_zoom_out(I2s[s-1], I2s[s], nxx[s-1], nyy[s-1], scaleFactor);
    }

    // Process scales from coarse to fine
    // ... rest of the processing code ...

    // Cleanup
    for (int s = 0; s < nScales; s++) {
        free(I1s[s]);
        free(I2s[s]);
        if (s > 0) {
            free(us[s]);
            free(vs[s]);
        }
    }
    free(I1s);
    free(I2s);
    free(us);
    free(vs);
    free(nxx);
    free(nyy);

    return 1;
}

// Change the return type of the gaussian function from void to int
int clg_gaussian(double* image, int width, int height, double sigma) {
    if (!image || width <= 0 || height <= 0 || sigma <= 0) {
        printf("Error: Invalid input parameters for Gaussian smoothing\n");
        return 0;
    }

    // Calculate kernel size (make sure it's odd)
    int ksize = (int)(2.0 * ceil(3.0 * sigma) + 1.0);
    ksize = (ksize % 2 == 0) ? ksize + 1 : ksize;
    
    if (ksize > width || ksize > height) {
        printf("Error: Gaussian kernel size (%d) too large for image dimensions (%dx%d)\n", 
               ksize, width, height);
        return 0;
    }

    // Allocate temporary buffer
    double* temp = (double*)malloc(width * height * sizeof(double));
    if (!temp) {
        printf("Error: Failed to allocate memory for Gaussian smoothing\n");
        return 0;
    }

    // Create Gaussian kernel
    double* kernel = (double*)malloc(ksize * sizeof(double));
    if (!kernel) {
        printf("Error: Failed to allocate memory for Gaussian kernel\n");
        free(temp);
        return 0;
    }

    // Fill kernel with Gaussian values
    int half = ksize / 2;
    double sum = 0.0;
    for (int x = 0; x < ksize; x++) {
        double diff = (x - half);
        kernel[x] = exp(-(diff * diff) / (2.0 * sigma * sigma));
        sum += kernel[x];
    }

    // Normalize kernel
    for (int x = 0; x < ksize; x++) {
        kernel[x] /= sum;
    }

    // Horizontal pass
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            double sum = 0.0;
            for (int k = -half; k <= half; k++) {
                int xx = x + k;
                // Mirror boundary conditions
                if (xx < 0) xx = -xx;
                if (xx >= width) xx = 2 * (width - 1) - xx;
                sum += image[y * width + xx] * kernel[k + half];
            }
            temp[y * width + x] = sum;
        }
    }

    // Vertical pass
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            double sum = 0.0;
            for (int k = -half; k <= half; k++) {
                int yy = y + k;
                // Mirror boundary conditions
                if (yy < 0) yy = -yy;
                if (yy >= height) yy = 2 * (height - 1) - yy;
                sum += temp[yy * width + x] * kernel[k + half];
            }
            image[y * width + x] = sum;
        }
    }

    // Clean up
    free(kernel);
    free(temp);
    return 1;
}

static void clg_zoom_out(const double *I, double *Iout, const int nx, const int ny, const double factor) {
    if (!I || !Iout || nx <= 0 || ny <= 0 || factor <= 0.0 || factor >= 1.0) {
        printf("Error: Invalid parameters in zoom_out\n");
        return;
    }

    // Calculate new dimensions
    const int nxx = (int)(nx * factor);
    const int nyy = (int)(ny * factor);

    if (nxx <= 0 || nyy <= 0) {
        printf("Error: Invalid output dimensions in zoom_out (%dx%d)\n", nxx, nyy);
        return;
    }

    printf("  zoom_out: input dims=%dx%d, output dims=%dx%d, factor=%f\n", 
           nx, ny, nxx, nyy, factor);

    // Compute the zoom factor
    const double factorx = ((double)nx-1.0) / ((double)nxx-1.0);
    const double factory = ((double)ny-1.0) / ((double)nyy-1.0);

    // For each pixel in the output image
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < nyy; i++) {
        for (int j = 0; j < nxx; j++) {
            // Get the real position in the input image
            const double i_in = i * factory;
            const double j_in = j * factorx;

            // Get the four neighboring pixels
            const int i1 = (int)floor(i_in);
            const int i2 = (i1 < ny-1) ? i1+1 : i1;
            const int j1 = (int)floor(j_in);
            const int j2 = (j1 < nx-1) ? j1+1 : j1;

            // Get the interpolation weights
            const double di = i_in - i1;
            const double dj = j_in - j1;

            // Compute the interpolated value using bilinear interpolation
            const double value = 
                (1.0-di) * (1.0-dj) * I[i1*nx + j1] +
                (1.0-di) * dj * I[i1*nx + j2] +
                di * (1.0-dj) * I[i2*nx + j1] +
                di * dj * I[i2*nx + j2];

            // Store the result
            Iout[i*nxx + j] = value;
        }
    }
}

#endif //CLG_OF_C
