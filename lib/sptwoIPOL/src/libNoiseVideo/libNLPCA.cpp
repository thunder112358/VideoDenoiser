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

#include "libNLPCA.h"
#include <math.h>


namespace libUSTG
{




//put patch data in vector format
void  denoising_block_to_vector_nsq(cflimage &fpI, float *fpV, int i, int j,
                                    int ixWin, int iyWin)
{
    int kk0 = 0;
    for (int ii=0; ii < fpI.c(); ii++) {	
        int kk=kk0;	
        for(int r=0; r < iyWin; r++) {
            for(int s= 0; s < ixWin; s++) {
                int m = (j + r) * fpI.w() + i + s;
                fpV[kk] = fpI[ii * fpI.wh() + m];
                kk+=fpI.c();
            }
        }
        kk0++;
    }
}


//accumulate data from vectors to patch format
void denoising_add_vector_to_block_nsq(cflimage &fpO, float *fpV, float *fpW,  int ci, int cj, int ixWin, int iyWin)
{


    int kk = 0;
    for(int r=0; r < iyWin; r++)
        for(int s=0; s < ixWin; s++) {
            int m = (cj + r) * fpO.w() + ci + s;

            for (int ic = 0; ic < fpO.c(); ic++, kk++) fpO[ic * fpO.wh() + m]   +=  fpV[kk];

            if (fpW) fpW[m] += 1.0f;
        }
}





    
//!
//! Center X data, n row vectors of length p to its barycenter
//! That is, substract to each column its average
//!
//! X              matrix of data n x p
//! barycenter     output the barycenter, average of each column
void  lmmse_center_data(laMatrix  &X,float *barycenter)
{

    int n = X.nrows();
    int p = X.ncols();

    for(int j=0; j<p; j++) {

        float sum=0.0;
        for(int i=0; i<n; i++) sum += X[i][j];

        sum/=(float) n;
        barycenter[j]=sum;

        for(int i=0; i<n; i++) X[i][j] -= barycenter[j];
    }

}


void hard_principal_values(float *S, float **M, float rsigma2, int p, int npts)
{
    for(int jj=0; jj < p; jj++)
        for(int ii=0; ii < npts; ii++)
            if (S[jj] < rsigma2)
                M[ii][jj] = 0.0f;
            else
                M[ii][jj] = 1.0f;
}


void wiener_oracle_coefficients(float **U, float **M, float rsigma2, int p, int npts)
{

    float fAux;
    for(int ii=0; ii < npts; ii++) {
        for(int jj=0; jj < p; jj++) {
            fAux =  U[ii][jj] * U[ii][jj];
            M[ii][jj] = fAux / (fAux + rsigma2);

        }
    }

}



float scalar_product(float *u, float *v, int n)
{

    double aux = 0.0f;
    for(int i=0; i < n; i++)
        aux += (double) u[i] * (double) v[i];

    return (float) aux;

}




//Algorithm 9
void  NLPCAdenoising(laMatrix &X2res, laMatrix &X02res, laMatrix &XR2,
                     int flagX02, float fSigma, float fRMult)
{


    int npts = X2res.nrows();
    int p = X2res.ncols();

    float fSigma2 = fSigma * fSigma;

    //! Center X data, n row vectors of length p to its barycenter (output)
    laVector barycenter(p);
    laVector barycenter0(p);

    //! Center X02 data
    lmmse_center_data(X2res, barycenter.v());
    if (flagX02)
        lmmse_center_data(X02res, barycenter0.v());

    // NLPCA

    //! Compute PCA by SVD decomposition
    laVector S(p);          //! Variance or eigenvalues
    laMatrix V(p, p);       //! PCs or eigenvectors
    laMatrix A(npts, p);    //! New coefficients (coefficients in eigenvectors base)

    laMatrix A02;           //! Coefficients of auxiliary image if used



    if (flagX02) {          //! If oracle

        A02.create(npts, p);

        //! Compute PCA on auxiliary image
        compute_pca_svd(X02res, S, V, A02);

        for(int ii=0; ii < npts; ii++) {            //! Compute U coefficients with PCA basis V
            for(int jj=0; jj < p; jj++) {

                float aux = 0.0;
                for (int kk=0; kk <p; kk++)
                    aux += X2res[ii][kk] * V[kk][jj];

                A[ii][jj] =  aux;

            }

        }


    } else {                                //! If  non secondary image, just compute PCA
        compute_pca_svd(X2res, S, V, A);
    }


    //! Compute coefficient modifiers F
    laMatrix F(npts, p);  
    float fRSigma2 = fRMult * fRMult * fSigma2;

    if (flagX02) {          //! If oracle
        wiener_oracle_coefficients(A02.v(),F.v(),fRSigma2,p,npts);
    } else {
        hard_principal_values(S.v(),F.v(),fRSigma2,p,npts);
    }


    //! Modify coefficients
    for(int jj=0; jj < p; jj++)
        for(int ii=0; ii < npts; ii++)
            A[ii][jj] *= F[ii][jj];



    /// Reconstruct denoised patches and save in XR2
    for(int ii=0; ii < npts; ii++) {
        for (int kk=0; kk < p; kk++) {
            XR2[ii][kk] = barycenter[kk] + (float) scalar_product(A[ii], V[kk],p);

        }
    }

}



}


