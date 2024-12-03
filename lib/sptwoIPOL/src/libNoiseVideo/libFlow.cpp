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


#include "libFlow.h"

namespace libUSTGFLOW
{


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>



// this function is like "malloc", but it returns always a valid pointer
void *xmalloc(size_t size)
{
    void *p = malloc(size);
    if (!p)
        exit(fprintf(stderr, "out of memory\n"));
    return p;
}



//!
//! Rejects points for which the corresponding point has not been matches or the matching does not correspond to the reference point
//!
void check_flow_reciprocity(float *uflow, float *vflow, float *uflow2, float *vflow2, float *omask,float fThreshold, int width, int height)
{

    for (int ii=0; ii < width*height; ii++)  omask[ii]=1.0f;

    for(int ipy = 0 ;  ipy < height ; ipy++)
        for(int ipx = 0 ;  ipx < width ; ipx++) {

            int l = ipy * width + ipx;

            // pixel in input2 for ipx
            int ipxRev = ipx +  (int) rintf(uflow[l]);
            int ipyRev = ipy +  (int) rintf(vflow[l]);

            if (ipxRev >= 0 && ipxRev < width && ipyRev >= 0 && ipyRev < height) {

                float aux1 = fabsf( uflow[l] + uflow2[ipyRev * width + ipxRev]);
                float aux2 = fabsf( vflow[l] + vflow2[ipyRev * width + ipxRev]);
                float aux = MAX(aux1,aux2);


                if (aux > fThreshold)   omask[l] = 0.0;

            } else
                omask[l] = 0.0;

        }


}






/*****

THE FOLLOWING FUNCTIONS COME FROM bicubic_interpolation.c, by Javier Sanchez,
available at http://www.ipol.im/pub/art/2013/26/

****/


#define BOUNDARY_CONDITION 0
//0 Neumann
//1 Periodic
//2 Symmetric

/**
 *
 * Neumann boundary condition test
 *
 **/
static int neumann_bc(int x, int nx, bool *out)
{
    if(x < 0) {
        x = 0;
        *out = true;
    } else if (x >= nx) {
        x = nx - 1;
        *out = true;
    }

    return x;
}

/**
 *
 * Periodic boundary condition test
 *
 **/
static int periodic_bc(int x, int nx, bool *out)
{
    if(x < 0) {
        const int n   = 1 - (int)(x/(nx+1));
        const int ixx = x + n * nx;

        x =   ixx% nx;
        *out = true;
    } else if(x >= nx) {
        x = x % nx;
        *out = true;
    }

    return x;
}


/**
 *
 * Symmetric boundary condition test
 *
 **/
static int symmetric_bc(int x, int nx, bool *out)
{
    if(x < 0) {
        const int border = nx - 1;
        const int xx = -x;
        const int n  = (int)(xx/border) % 2;

        if ( n ) x = border - ( xx % border );
        else x = xx % border;
        *out = true;
    } else if ( x >= nx ) {
        const int border = nx - 1;
        const int n = (int)(x/border) % 2;

        if ( n ) x = border - ( x % border );
        else x = x % border;
        *out = true;
    }

    return x;
}


/**
 *
 * Cubic interpolation in one dimension
 *
 **/
static double cubic_interpolation_cell (
    double v[4],  //interpolation points
    double x      //point to be interpolated
)
{
    return  v[1] + 0.5 * x * (v[2] - v[0] +
                              x * (2.0 *  v[0] - 5.0 * v[1] + 4.0 * v[2] - v[3] +
                                   x * (3.0 * (v[1] - v[2]) + v[3] - v[0])));
}


/**
 *
 * Bicubic interpolation in two dimensions
 *
 **/
static double bicubic_interpolation_cell (
    double p[4][4], //array containing the interpolation points
    double x,       //x position to be interpolated
    double y        //y position to be interpolated
)
{
    double v[4];
    v[0] = cubic_interpolation_cell(p[0], y);
    v[1] = cubic_interpolation_cell(p[1], y);
    v[2] = cubic_interpolation_cell(p[2], y);
    v[3] = cubic_interpolation_cell(p[3], y);
    return cubic_interpolation_cell(v, x);
}

/**
 *
 * Compute the bicubic interpolation of a point in an image.
 * Detect if the point goes outside the image domain.
 *
 **/
float bicubic_interpolation_at(
    const float *input, //image to be interpolated
    const float  uu,    //x component of the vector field
    const float  vv,    //y component of the vector field
    const int    nx,    //image width
    const int    ny,    //image height
    bool         border_out //if true, return zero outside the region
)
{
    const int sx = (uu < 0)? -1: 1;
    const int sy = (vv < 0)? -1: 1;

    int x, y, mx, my, dx, dy, ddx, ddy;
    bool out[1] = {false};

    //apply the corresponding boundary conditions
    switch(BOUNDARY_CONDITION) {

    case 0:
        x   = neumann_bc((int) uu, nx, out);
        y   = neumann_bc((int) vv, ny, out);
        mx  = neumann_bc((int) uu - sx, nx, out);
        my  = neumann_bc((int) vv - sx, ny, out);
        dx  = neumann_bc((int) uu + sx, nx, out);
        dy  = neumann_bc((int) vv + sy, ny, out);
        ddx = neumann_bc((int) uu + 2*sx, nx, out);
        ddy = neumann_bc((int) vv + 2*sy, ny, out);
        break;

    case 1:
        x   = periodic_bc((int) uu, nx, out);
        y   = periodic_bc((int) vv, ny, out);
        mx  = periodic_bc((int) uu - sx, nx, out);
        my  = periodic_bc((int) vv - sx, ny, out);
        dx  = periodic_bc((int) uu + sx, nx, out);
        dy  = periodic_bc((int) vv + sy, ny, out);
        ddx = periodic_bc((int) uu + 2*sx, nx, out);
        ddy = periodic_bc((int) vv + 2*sy, ny, out);
        break;

    case 2:
        x   = symmetric_bc((int) uu, nx, out);
        y   = symmetric_bc((int) vv, ny, out);
        mx  = symmetric_bc((int) uu - sx, nx, out);
        my  = symmetric_bc((int) vv - sx, ny, out);
        dx  = symmetric_bc((int) uu + sx, nx, out);
        dy  = symmetric_bc((int) vv + sy, ny, out);
        ddx = symmetric_bc((int) uu + 2*sx, nx, out);
        ddy = symmetric_bc((int) vv + 2*sy, ny, out);
        break;

    default:
        x   = neumann_bc((int) uu, nx, out);
        y   = neumann_bc((int) vv, ny, out);
        mx  = neumann_bc((int) uu - sx, nx, out);
        my  = neumann_bc((int) vv - sx, ny, out);
        dx  = neumann_bc((int) uu + sx, nx, out);
        dy  = neumann_bc((int) vv + sy, ny, out);
        ddx = neumann_bc((int) uu + 2*sx, nx, out);
        ddy = neumann_bc((int) vv + 2*sy, ny, out);
        break;
    }

    if(*out && border_out)
        return 0.0;

    else {
        //obtain the interpolation points of the image
        const float p11 = input[mx  + nx * my];
        const float p12 = input[x   + nx * my];
        const float p13 = input[dx  + nx * my];
        const float p14 = input[ddx + nx * my];

        const float p21 = input[mx  + nx * y];
        const float p22 = input[x   + nx * y];
        const float p23 = input[dx  + nx * y];
        const float p24 = input[ddx + nx * y];

        const float p31 = input[mx  + nx * dy];
        const float p32 = input[x   + nx * dy];
        const float p33 = input[dx  + nx * dy];
        const float p34 = input[ddx + nx * dy];

        const float p41 = input[mx  + nx * ddy];
        const float p42 = input[x   + nx * ddy];
        const float p43 = input[dx  + nx * ddy];
        const float p44 = input[ddx + nx * ddy];

        //create array
        double pol[4][4] = {
            {p11, p21, p31, p41},
            {p12, p22, p32, p42},
            {p13, p23, p33, p43},
            {p14, p24, p34, p44}
        };

        //return interpolation
        return bicubic_interpolation_cell(pol, uu-x, vv-y);
    }
}






/**
 *
 * Compute the bicubic interpolation of an image.
 *
 **/
void bicubic_interpolation_warp(
    const float *input,     // image to be warped
    const float *u,         // x component of the vector field
    const float *v,         // y component of the vector field
    float       *output,    // image warped with bicubic interpolation
    const int    nx,        // image width
    const int    ny,        // image height
    bool         border_out // if true, put zeros outside the region
)
{


    #pragma omp parallel shared(input,u,v,border_out)
    {
        #pragma omp for //schedule(dynamic) nowait
        for(int i = 0; i < ny; i++)
        {
            int   p;
            float uu;
            float vv;


            for(int j = 0; j < nx; j++) {
                p  = i * nx + j;
                uu = (float) (j + u[p]);
                vv = (float) (i + v[p]);

                // obtain the bicubic interpolation at position (uu, vv)
                output[p] = bicubic_interpolation_at(input,
                uu, vv, nx, ny, border_out);
            }
        }
    }
}


/*****
END OF FUNCTIONS FROM bicubic_interpolation.c, by Javier Sanchez,
available at http://www.ipol.im/pub/art/2013/26/

****/


/*****

THE FOLLOWING FUNCTIONS COME FROM mask.c, by Javier Sanchez,
available at http://www.ipol.im/pub/art/2013/26/

****/

#define BOUNDARY_CONDITION_DIRICHLET 0
#define BOUNDARY_CONDITION_REFLECTING 1
#define BOUNDARY_CONDITION_PERIODIC 2

#define DEFAULT_GAUSSIAN_WINDOW_SIZE 5
#define DEFAULT_BOUNDARY_CONDITION BOUNDARY_CONDITION_REFLECTING


/**
 *
 * Details on how to compute the divergence and the grad(u) can be found in:
 * [2] A. Chambolle, "An Algorithm for Total Variation Minimization and
 * Applications", Journal of Mathematical Imaging and Vision, 20: 89-97, 2004
 *
 **/


/**
 *
 * Function to compute the divergence with backward differences
 * (see [2] for details)
 *
 **/
void divergence(
    const float *v1, // x component of the vector field
    const float *v2, // y component of the vector field
    float *div,      // output divergence
    const int nx,    // image width
    const int ny     // image height
)
{




    #pragma omp parallel shared(v1,v2,div)
    {
        #pragma omp for //schedule(dynamic) nowait
        for (int i = 1; i < ny-1; i++)
        {

            int p, p1, p2;
            double v1x, v2y;

            for(int j = 1; j < nx-1; j++) {
                p  = i * nx + j;
                p1 = p - 1;
                p2 = p - nx;

                v1x = v1[p] - v1[p1];
                v2y = v2[p] - v2[p2];

                div[p] = v1x + v2y;
            }
        }
    }


    // compute the divergence on the first and last rows
    for (int j = 1; j < nx-1; j++) {
        const int p = (ny-1) * nx + j;

        div[j] = v1[j] - v1[j-1] + v2[j];
        div[p] = v1[p] - v1[p-1] - v2[p-nx];
    }

    // compute the divergence on the first and last columns
    for (int i = 1; i < ny-1; i++) {
        const int p1 = i * nx;
        const int p2 = (i+1) * nx - 1;

        div[p1] =  v1[p1]   + v2[p1] - v2[p1 - nx];
        div[p2] = -v1[p2-1] + v2[p2] - v2[p2 - nx];

    }

    div[0]         =  v1[0] + v2[0];
    div[nx-1]      = -v1[nx - 2] + v2[nx - 1];
    div[(ny-1)*nx] =  v1[(ny-1)*nx] - v2[(ny-2)*nx];
    div[ny*nx-1]   = -v1[ny*nx - 2] - v2[(ny-1)*nx - 1];
}




/**
 *
 * Function to compute the gradient with forward differences
 * (see [2] for details)
 *
 **/
void forward_gradient(
    const float *f, //input image
    float *fx,      //computed x derivative
    float *fy,      //computed y derivative
    const int nx,   //image width
    const int ny    //image height
)
{


    #pragma omp parallel shared(f,fx,fy)
    {
        #pragma omp for //schedule(dynamic) nowait
        for (int i = 0; i < ny-1; i++)
        {

            int p,p1,p2;
            for(int j = 0; j < nx-1; j++) {
                p  = i * nx + j;
                p1 = p + 1;
                p2 = p + nx;

                fx[p] = f[p1] - f[p];
                fy[p] = f[p2] - f[p];
            }
        }
    }



    // compute the gradient on the last row
    for (int j = 0; j < nx-1; j++) {
        const int p = (ny-1) * nx + j;

        fx[p] = f[p+1] - f[p];
        fy[p] = 0;
    }

    // compute the gradient on the last column
    for (int i = 1; i < ny; i++) {
        const int p = i * nx-1;

        fx[p] = 0;
        fy[p] = f[p+nx] - f[p];
    }

    fx[ny * nx - 1] = 0;
    fy[ny * nx - 1] = 0;
}







/**
 *
 * Function to compute the gradient with centered differences
 *
 **/
void centered_gradient(
    const float *input,  //input image
    float *dx,           //computed x derivative
    float *dy,           //computed y derivative
    const int nx,        //image width
    const int ny         //image height
)
{
    // compute the gradient on the center body of the image
    //printf("agma4\n");

    for (int i = 1; i < ny-1; i++) {
        for(int j = 1; j < nx-1; j++) {
            const int k = i * nx + j;
            dx[k] = 0.5*(input[k+1] - input[k-1]);
            dy[k] = 0.5*(input[k+nx] - input[k-nx]);
        }
    }

    // compute the gradient on the first and last rows
    for (int j = 1; j < nx-1; j++) {
        dx[j] = 0.5*(input[j+1] - input[j-1]);
        dy[j] = 0.5*(input[j+nx] - input[j]);

        const int k = (ny - 1) * nx + j;

        dx[k] = 0.5*(input[k+1] - input[k-1]);
        dy[k] = 0.5*(input[k] - input[k-nx]);
    }

    // compute the gradient on the first and last columns
    for(int i = 1; i < ny-1; i++) {
        const int p = i * nx;
        dx[p] = 0.5*(input[p+1] - input[p]);
        dy[p] = 0.5*(input[p+nx] - input[p-nx]);

        const int k = (i+1) * nx - 1;

        dx[k] = 0.5*(input[k] - input[k-1]);
        dy[k] = 0.5*(input[k+nx] - input[k-nx]);
    }

    // compute the gradient at the four corners
    dx[0] = 0.5*(input[1] - input[0]);
    dy[0] = 0.5*(input[nx] - input[0]);

    dx[nx-1] = 0.5*(input[nx-1] - input[nx-2]);
    dy[nx-1] = 0.5*(input[2*nx-1] - input[nx-1]);

    dx[(ny-1)*nx] = 0.5*(input[(ny-1)*nx + 1] - input[(ny-1)*nx]);
    dy[(ny-1)*nx] = 0.5*(input[(ny-1)*nx] - input[(ny-2)*nx]);

    dx[ny*nx-1] = 0.5*(input[ny*nx-1] - input[ny*nx-1-1]);
    dy[ny*nx-1] = 0.5*(input[ny*nx-1] - input[(ny-1)*nx-1]);
}



/**
 *
 * In-place Gaussian smoothing of an image
 *
 */
void gaussian(
    float *I,             // input/output image
    const int xdim,       // image width
    const int ydim,       // image height
    const double sigma    // Gaussian sigma
)
{
    const int boundary_condition = DEFAULT_BOUNDARY_CONDITION;
    const int window_size = DEFAULT_GAUSSIAN_WINDOW_SIZE;

    const double den  = 2*sigma*sigma;
    const int   size = (int) (window_size * sigma) + 1 ;
    const int   bdx  = xdim + size;
    const int   bdy  = ydim + size;

    if (boundary_condition && size > xdim) {
        fprintf(stderr, "GaussianSmooth: sigma too large\n");
        abort();
    }

    // compute the coefficients of the 1D convolution kernel
    double *B = new double[size];
    for(int i = 0; i < size; i++)
        B[i] = 1 / (sigma * sqrt(2.0 * 3.1415926)) * exp(-i * i / den);

    // normalize the 1D convolution kernel
    double norm = 0;
    for(int i = 0; i < size; i++)
        norm += B[i];
    norm *= 2;
    norm -= B[0];
    for(int i = 0; i < size; i++)
        B[i] /= norm;

    // convolution of each line of the input image

    #pragma omp parallel shared(I, B)

    {
        #pragma omp for
        for (int k = 0; k < ydim; k++) {
            double *R = (double *) xmalloc((size + xdim + size)*sizeof*R);

            int i, j;
            for (i = size; i < bdx; i++)
                R[i] = I[k * xdim + i - size];

            switch (boundary_condition) {
            case BOUNDARY_CONDITION_DIRICHLET:
                for(i = 0, j = bdx; i < size; i++, j++)
                    R[i] = R[j] = 0;
                break;

            case BOUNDARY_CONDITION_REFLECTING:
                for(i = 0, j = bdx; i < size; i++, j++) {
                    R[i] = I[k * xdim + size-i];
                    R[j] = I[k * xdim + xdim-i-1];
                }
                break;

            case BOUNDARY_CONDITION_PERIODIC:
                for(i = 0, j = bdx; i < size; i++, j++) {
                    R[i] = I[k * xdim + xdim-size+i];
                    R[j] = I[k * xdim + i];
                }
                break;
            }

            for (i = size; i < bdx; i++) {
                double sum = B[0] * R[i];
                for (j = 1; j < size; j++ )
                    sum += B[j] * ( R[i-j] + R[i+j] );
                I[k * xdim + i - size] = sum;
            }

            free(R);
        }

    }



    // convolution of each column of the input image
    #pragma omp parallel shared(I, B)

    {
        #pragma omp for
        for (int k = 0; k < xdim; k++) {
            double *T = (double *)xmalloc((size + ydim + size)*sizeof*T);
            int i, j;
            for (i = size; i < bdy; i++)
                T[i] = I[(i - size) * xdim + k];

            switch (boundary_condition) {
            case BOUNDARY_CONDITION_DIRICHLET:
                for (i = 0, j = bdy; i < size; i++, j++)
                    T[i] = T[j] = 0;
                break;

            case BOUNDARY_CONDITION_REFLECTING:
                for (i = 0, j = bdy; i < size; i++, j++) {
                    T[i] = I[(size-i) * xdim + k];
                    T[j] = I[(ydim-i-1) * xdim + k];
                }
                break;

            case BOUNDARY_CONDITION_PERIODIC:
                for( i = 0, j = bdx; i < size; i++, j++) {
                    T[i] = I[(ydim-size+i) * xdim + k];
                    T[j] = I[i * xdim + k];
                }
                break;
            }

            for (i = size; i < bdy; i++) {
                double sum = B[0] * T[i];
                for (j = 1; j < size; j++ )
                    sum += B[j] * (T[i-j] + T[i+j]);
                I[(i - size) * xdim + k] = sum;
            }

            free(T);
        }
    }





    delete[] B;
}


/*****
END OF FUNCTIONS FROM mask.c, by Javier Sanchez,
available at http://www.ipol.im/pub/art/2013/26/

****/


/*****

THE FOLLOWING FUNCTIONS COME FROM zoom.c, by Javier Sanchez,
available at http://www.ipol.im/pub/art/2013/26/

****/

#define ZOOM_SIGMA_ZERO 0.6

/**
 *
 * Compute the size of a zoomed image from the zoom factor
 *
 **/
void zoom_size(
    int nx,      // width of the orignal image
    int ny,      // height of the orignal image
    int *nxx,    // width of the zoomed image
    int *nyy,    // height of the zoomed image
    float factor // zoom factor between 0 and 1
)
{
    //compute the new size corresponding to factor
    //we add 0.5 for rounding off to the closest number
    *nxx = (int)((float) nx * factor + 0.5);
    *nyy = (int)((float) ny * factor + 0.5);
}

/**
 *
 * Downsample an image
 *
 **/
void zoom_out(
    const float *I,    // input image
    float *Iout,       // output image
    const int nx,      // image width
    const int ny,      // image height
    const float factor // zoom factor between 0 and 1
)
{
    // temporary working image
    float *Is = (float *) xmalloc(nx * ny * sizeof*Is);
    for(int i = 0; i < nx * ny; i++)
        Is[i] = I[i];

    // compute the size of the zoomed image
    int nxx, nyy;
    zoom_size(nx, ny, &nxx, &nyy, factor);

    // compute the Gaussian sigma for smoothing
    const float sigma = ZOOM_SIGMA_ZERO * sqrt(1.0/(factor*factor) - 1.0);

    // pre-smooth the image
    gaussian(Is, nx, ny, sigma);

    // re-sample the image using bicubic interpolation
    //#pragma omp parallel for
    #pragma omp parallel shared(Is, nxx, nyy, Iout)

    {
        #pragma omp for //schedule(dynamic) nowait
        for (int i1 = 0; i1 < nyy; i1++) {
            float i2, j2;

            for (int j1 = 0; j1 < nxx; j1++) {
                i2  = (float) i1 / factor;
                j2  = (float) j1 / factor;

                float g = bicubic_interpolation_at(Is, j2, i2, nx, ny, false);
                Iout[i1 * nxx + j1] = g;
            }
        }
    }


    free(Is);
}



/**
 *
 * Function to upsample the image
 *
 **/
void zoom_in(
    const float *I, // input image
    float *Iout,    // output image
    int nx,         // width of the original image
    int ny,         // height of the original image
    int nxx,        // width of the zoomed image
    int nyy         // height of the zoomed image
)
{
    // compute the zoom factor
    const float factorx = ((float)nxx / nx);
    const float factory = ((float)nyy / ny);

    // re-sample the image using bicubic interpolation
//#pragma omp parallel shared(I, Iout)
    {
//#pragma omp for //schedule(dynamic) nowait
        for (int i1 = 0; i1 < nyy; i1++)
            for (int j1 = 0; j1 < nxx; j1++) {
                float i2 =  (float) i1 / factory;
                float j2 =  (float) j1 / factorx;

                float g = bicubic_interpolation_at(I, j2, i2, nx, ny, false);
                Iout[i1 * nxx + j1] = g;
            }
    }

}


/*****
END OF FUNCTIONS FROM zoom.c, by Javier Sanchez,
available at http://www.ipol.im/pub/art/2013/26/

****/





/*****

THE FOLLOWING FUNCTIONS COME FROM tvl1flow_lib.c, by Javier Sanchez,
available at http://www.ipol.im/pub/art/2013/26/

****/


#define MAX_ITERATIONS 300
#define PRESMOOTHING_SIGMA 0.8
#define GRAD_IS_ZERO 1E-10


/**
 * Implementation of the Zach, Pock and Bischof dual TV-L1 optic flow method
 *
 * see reference:
 *  [1] C. Zach, T. Pock and H. Bischof, "A Duality Based Approach for Realtime
 *      TV-L1 Optical Flow", In Proceedings of Pattern Recognition (DAGM),
 *      Heidelberg, Germany, pp. 214-223, 2007
 *
 *
 * Details on the total variation minimization scheme can be found in:
 *  [2] A. Chambolle, "An Algorithm for Total Variation Minimization and
 *      Applications", Journal of Mathematical Imaging and Vision, 20: 89-97, 2004
 **/


/**
 *
 * Function to compute the optical flow in one scale
 *
 **/
void Dual_TVL1_optic_flow(
    float *I0,           // source image
    float *I1,           // target image
    float *u1,           // x component of the optical flow
    float *u2,           // y component of the optical flow
    const int   nx,      // image width
    const int   ny,      // image height
    const float tau,     // time step
    const float lambda,  // weight parameter for the data term
    const float theta,   // weight parameter for (u - v)²
    const int   warps,   // number of warpings per scale
    const float epsilon, // tolerance for numerical convergence
    const bool  verbose,  // enable/disable the verbose mode
    int iflagMedian)
{
    const int   size = nx * ny;
    const float l_t = lambda * theta;

    size_t sf = sizeof(float);
    float *I1x    = (float *) malloc(size*sf);
    float *I1y    = (float *) xmalloc(size*sf);
    float *I1w    = (float *) xmalloc(size*sf);
    float *I1wx   = (float *) xmalloc(size*sf);
    float *I1wy   = (float *) xmalloc(size*sf);
    float *rho_c  = (float *) xmalloc(size*sf);
    float *v1     = (float *) xmalloc(size*sf);
    float *v2     = (float *) xmalloc(size*sf);
    float *p11    = (float *) xmalloc(size*sf);
    float *p12    = (float *) xmalloc(size*sf);
    float *p21    = (float *) xmalloc(size*sf);
    float *p22    = (float *) xmalloc(size*sf);
    float *div    = (float *) xmalloc(size*sf);
    float *grad   = (float *) xmalloc(size*sf);
    float *div_p1 = (float *) xmalloc(size*sf);
    float *div_p2 = (float *) xmalloc(size*sf);
    float *u1x    = (float *) xmalloc(size*sf);
    float *u1y    = (float *) xmalloc(size*sf);
    float *u2x    = (float *) xmalloc(size*sf);
    float *u2y    = (float *) xmalloc(size*sf);


    //! Compute gradient of second image
    centered_gradient(I1, I1x, I1y, nx, ny);


    // initialization of p
    for (int i = 0; i < size; i++) {
        p11[i] = p12[i] = 0.0;
        p21[i] = p22[i] = 0.0;
    }

    for (int warpings = 0; warpings < warps; warpings++) {

        //! compute the warping of the target image and its derivatives
        //#pragma omp parallel sections
        {
            //#pragma omp section
            {
                bicubic_interpolation_warp(I1,  u1, u2, I1w,  nx, ny, true);
            }

            //#pragma omp section
            {
                bicubic_interpolation_warp(I1x, u1, u2, I1wx, nx, ny, true);
            }

            //#pragma omp section
            {
                bicubic_interpolation_warp(I1y, u1, u2, I1wy, nx, ny, true);
            }
        }


        for (int i = 0; i < size; i++) {
            const float Ix2 = I1wx[i] * I1wx[i];
            const float Iy2 = I1wy[i] * I1wy[i];

            // store the |Grad(I1)|^2
            grad[i] = (Ix2 + Iy2);

            // compute the constant part of the rho function
            rho_c[i] = (I1w[i] - I1wx[i] * u1[i]
                        - I1wy[i] * u2[i] - I0[i]);
        }

        int n = 0;
        float error = INFINITY;
        while (error > epsilon * epsilon && n < MAX_ITERATIONS) {
            n++;
            // estimate the values of the variable (v1, v2)
            // (thresholding opterator TH)

            for (int i = 0; i < size; i++) {
                const float rho = rho_c[i]
                                  + (I1wx[i] * u1[i] + I1wy[i] * u2[i]);

                float d1, d2;

                if (rho < - l_t * grad[i]) {
                    d1 = l_t * I1wx[i];
                    d2 = l_t * I1wy[i];
                } else {
                    if (rho > l_t * grad[i]) {
                        d1 = -l_t * I1wx[i];
                        d2 = -l_t * I1wy[i];
                    } else {
                        if (grad[i] < GRAD_IS_ZERO)
                            d1 = d2 = 0;
                        else {
                            float fi = -rho/grad[i];
                            d1 = fi * I1wx[i];
                            d2 = fi * I1wy[i];
                        }
                    }
                }

                v1[i] = u1[i] + d1;
                v2[i] = u2[i] + d2;
            }

            // compute the divergence of the dual variable (p1, p2)
            //#pragma omp parallel sections
            {
                //#pragma omp section
                {
                    divergence(p11, p12, div_p1, nx ,ny);
                }

                //#pragma omp section
                {

                    divergence(p21, p22, div_p2, nx ,ny);
                }

            }


            // estimate the values of the optical flow (u1, u2)
            error = 0.0;


            for (int i = 0; i < size; i++) {
                const float u1k = u1[i];
                const float u2k = u2[i];

                u1[i] = v1[i] + theta * div_p1[i];
                u2[i] = v2[i] + theta * div_p2[i];

                error += (u1[i] - u1k) * (u1[i] - u1k) +
                         (u2[i] - u2k) * (u2[i] - u2k);
            }
            error /= size;

            // compute the gradient of the optical flow (Du1, Du2)
            //#pragma omp parallel sections
            {
                //#pragma omp section
                {
                    forward_gradient(u1, u1x, u1y, nx ,ny);
                }

                //#pragma omp section
                {
                    forward_gradient(u2, u2x, u2y, nx ,ny);

                }
            }


            for (int i = 0; i < size; i++) {
                const float taut = tau / theta;
                const float g1   = hypot(u1x[i], u1y[i]);
                const float g2   = hypot(u2x[i], u2y[i]);
                const float ng1  = 1.0 + taut * g1;
                const float ng2  = 1.0 + taut * g2;

                p11[i] = (p11[i] + taut * u1x[i]) / ng1;
                p12[i] = (p12[i] + taut * u1y[i]) / ng1;
                p21[i] = (p21[i] + taut * u2x[i]) / ng2;
                p22[i] = (p22[i] + taut * u2y[i]) / ng2;
            }
        }

        if (verbose)
            fprintf(stderr, "Warping: %d, "
                    "Iterations: %d, "
                    "Error: %f\n", warpings, n, error);



        //! If selected, filter flow (u1, u2) by median filtering
        if (iflagMedian) {
            float fRadius=2.50f;
            float *tmp = new float[nx*ny];

            libUSTG::fiPatchMedian(u1, tmp,fRadius, nx, ny);
            libUSTG::fpCopy(tmp, u1, nx*ny);

            libUSTG::fiPatchMedian(u2, tmp,fRadius, nx, ny);
            libUSTG::fpCopy(tmp, u2, nx*ny);

        }


    }

    // delete allocated memory
    free(I1x);
    free(I1y);
    free(I1w);
    free(I1wx);
    free(I1wy);
    free(rho_c);
    free(v1);
    free(v2);
    free(p11);
    free(p12);
    free(p21);
    free(p22);
    free(div);
    free(grad);
    free(div_p1);
    free(div_p2);
    free(u1x);
    free(u1y);
    free(u2x);
    free(u2y);
}

/**
 *
 * Compute the max and min of an array
 *
 **/
static void getminmax(
    float *min,     // output min
    float *max,     // output max
    const float *x, // input array
    int n           // array size
)
{
    *min = *max = x[0];
    for (int i = 1; i < n; i++) {
        if (x[i] < *min)
            *min = x[i];
        if (x[i] > *max)
            *max = x[i];
    }
}

/**
 *
 * Function to normalize the images between 0 and 255
 *
 **/
void image_normalization(
    const float *I0,  // input image0
    const float *I1,  // input image1
    float *I0n,       // normalized output image0
    float *I1n,       // normalized output image1
    int size          // size of the image
)
{
    float max0, max1, min0, min1;

    // obtain the max and min of each image
    getminmax(&min0, &max0, I0, size);
    getminmax(&min1, &max1, I1, size);

    // obtain the max and min of both images
    const float max = (max0 > max1)? max0 : max1;
    const float min = (min0 < min1)? min0 : min1;
    const float den = max - min;

    if (den > 0)
        // normalize both images
        for (int i = 0; i < size; i++) {
            I0n[i] = 255.0 * (I0[i] - min) / den;
            I1n[i] = 255.0 * (I1[i] - min) / den;
        }

    else
        // copy the original images
        for (int i = 0; i < size; i++) {
            I0n[i] = I0[i];
            I1n[i] = I1[i];
        }
}





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
    struct flow_params &params)
{
    
    const float tau = params.tau;     // time step
    const float lambda = params.lambda;  // weight parameter for the data term
    const float theta = params.theta;   // weight parameter for (u - v)²
    const int   nscales = params.nscales; // number of scales
    const float zfactor = params.zfactor; // factor for building the image piramid
    const int   warps = params.warps;   // number of warpings per scale
    const float epsilon = params.epsilon; // tolerance for numerical convergence
    const bool  verbose = params.verbose;  // enable/disable the verbose mode
    int iflagMedian = params.iflagMedian;
    
    int size = nxx * nyy;

    // allocate memory for the pyramid structure
    float **I0s = (float **) xmalloc(nscales * sizeof(float*));
    float **I1s = (float **) xmalloc(nscales * sizeof(float*));
    float **u1s = (float **) xmalloc(nscales * sizeof(float*));
    float **u2s = (float **) xmalloc(nscales * sizeof(float*));
    int    *nx  = (int *) xmalloc(nscales * sizeof(int));
    int    *ny  = (int *) xmalloc(nscales * sizeof(int));

    I0s[0] = (float *)xmalloc(size*sizeof(float));
    I1s[0] = (float *)xmalloc(size*sizeof(float));

    u1s[0] = u1;
    u2s[0] = u2;
    nx [0] = nxx;
    ny [0] = nyy;

    // normalize the images between 0 and 255
    image_normalization(I0, I1, I0s[0], I1s[0], size);

    // pre-smooth the original images
    gaussian(I0s[0], nx[0], ny[0], PRESMOOTHING_SIGMA);
    gaussian(I1s[0], nx[0], ny[0], PRESMOOTHING_SIGMA);

    // create the scales
    for (int s = 1; s < nscales; s++) {
        zoom_size(nx[s-1], ny[s-1], &nx[s], &ny[s], zfactor);
        const int sizes = nx[s] * ny[s];

        // allocate memory
        I0s[s] = (float *)xmalloc(sizes*sizeof(float));
        I1s[s] = (float *)xmalloc(sizes*sizeof(float));
        u1s[s] = (float *)xmalloc(sizes*sizeof(float));
        u2s[s] = (float *)xmalloc(sizes*sizeof(float));

        // zoom in the images to create the pyramidal structure
        zoom_out(I0s[s-1], I0s[s], nx[s-1], ny[s-1], zfactor);
        zoom_out(I1s[s-1], I1s[s], nx[s-1], ny[s-1], zfactor);
    }


    //! initialize the flow at the coarsest scale
    for (int i = 0; i < nx[nscales-1] * ny[nscales-1]; i++)
        u1s[nscales-1][i] = u2s[nscales-1][i] = 0.0;

    //! pyramidal structure for computing the optical flow
    for (int s = nscales-1; s >= 0; s--) {
        if (verbose)
            fprintf(stderr, "Scale %d: %dx%d\n", s, nx[s], ny[s]);

        // compute the optical flow at the current scale
        Dual_TVL1_optic_flow(
            I0s[s], I1s[s], u1s[s], u2s[s], nx[s], ny[s],
            tau, lambda, theta, warps, epsilon, verbose, iflagMedian
        );

        // if this was the last scale, finish now
        if (!s) break;

        // otherwise, upsample the optical flow

        // zoom the optical flow for the next finer scale
        //#pragma omp parallel sections
        {

            //#pragma omp section
            {

                zoom_in(u1s[s], u1s[s-1], nx[s], ny[s], nx[s-1], ny[s-1]);
            }

            //#pragma omp section
            {
                zoom_in(u2s[s], u2s[s-1], nx[s], ny[s], nx[s-1], ny[s-1]);
            }

        }


        // scale the optical flow with the appropriate zoom factor
        for (int i = 0; i < nx[s-1] * ny[s-1]; i++) {
            u1s[s-1][i] *= (float) 1.0 / zfactor;
            u2s[s-1][i] *= (float) 1.0 / zfactor;
        }
    }

    // delete allocated memory
    for (int i = 1; i < nscales; i++) {
        free(I0s[i]);
        free(I1s[i]);
        free(u1s[i]);
        free(u2s[i]);
    }
    free(I0s[0]);
    free(I1s[0]);

    free(I0s);
    free(I1s);
    free(u1s);
    free(u2s);
    free(nx);
    free(ny);
}


/*****
END OF FUNCTIONS FROM tvl1flow_lib.c, by Javier Sanchez,
available at http://www.ipol.im/pub/art/2013/26/

****/



}


