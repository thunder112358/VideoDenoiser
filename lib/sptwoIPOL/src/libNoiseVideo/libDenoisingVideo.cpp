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

#include "libDenoisingVideo.h"
#include <math.h>
#include <time.h>


namespace libUSTG
{

struct diffpatches {
    int i, j;
    double diff;
};

int ixWin, iyWin, iTemp, iBloc;
int iWidth, iHeight, iChannels, aChannels, iFrames, ixDWin, iyDWin, iWin2, p;
int iDTemp, iDBloc, bsize, ifK, tframe;
unsigned char useOracle;
float fFixedThrDist;
    
float *mask, *mask2;


void SetGlobalParameters(cflimage *fpI, cflimage *fpA,
                         struct video_denoising_params &params)
{
    iWidth = fpI[0].w();
    iHeight = fpI[0].h();
    iChannels = fpI[0].c();    
    
    ixWin=params.iWin;
    iyWin=params.iWin;
    iTemp=params.iTemp;
    iBloc=params.iBloc;
    tframe=params.tframe;
    ifK=params.ifK;
    iFrames=params.iFrames;
    useOracle=params.useOracle;
    fFixedThrDist=params.fFixedThrDist;

    if (useOracle) aChannels = fpA[0].c();

    // Half side of patch for filtering (Averaging and PCA)
    ixDWin = (ixWin-1) / 2;
    // Half side of patch for filtering (Averaging and PCA)
    iyDWin = (iyWin-1) / 2;

    // Number of pixels per patch for filtering (Averaging and PCA)
    iWin2 =  ixWin * iyWin;
    // Length of patch vector for filtering (Averaging and PCA)
    p     =  ixWin * iyWin * iChannels;


    // Half side of temporal block
    iDTemp = (iTemp-1)/2;
    // Half side of research block
    iDBloc = (iBloc-1) / 2;
    // Number of pixel in research block
    bsize = iBloc * iBloc * iTemp;
}


void AllocateMemory()
{
    //! Mask: sum of weights in average per pixel
    mask = new float[iWidth*iHeight];
    fpClear(mask,0.0f,iWidth*iHeight);


    //! Mask2: treated pixels being central pixels of the window
    //!        or copied by copy paste strategy being center of the window
    mask2 = new float[iWidth*iHeight];
    fpClear(mask2,0.0f,iWidth*iHeight);
}

void FreeMemory()
{
    delete[] mask;
    delete[] mask2;
}

void  InitializeAuxiliaryVariables(laMatrix &X2, laMatrix &X02)
{
    // Memory: data matrix with as maximum bsize samples and length p
    //(sample per row)
    X2.create(bsize,p);
    
    //  Memory: data matrix for a secondary image,
    //(computing distances in second image and/or computing PCA)
    if (useOracle) X02.create(bsize, p);
}

void UpdateAuxiliaryVariables(laMatrix &X2, laMatrix &X02, laMatrix &X2res,
                              laMatrix &X02res, laMatrix &XR2, int npts)
{
    // Resize X2
    X2res = X2.copyBlock(0, 0, npts, p);
    
    // Resize secondary matrix
    if (useOracle) X02res = X02.copyBlock(0, 0, npts, p);
    
    // Memory for matrix of denoised vectors
    XR2.create(npts, p);
    XR2 = 0.0;
}



//auxiliary function to sort 3D blocks differences in increasing order
int order_increasing_diffs(const void *pVoid1, const void *pVoid2)
{
    struct diffpatches *node1, *node2;

    node1=(struct diffpatches *) pVoid1;
    node2=(struct diffpatches *) pVoid2;

    if (node1->diff < node2->diff) return -1;
    if (node1->diff > node2->diff) return 1;
    return 0;
}

//check if a patch contains occluded pixels
unsigned char check_zerovalues(float *u0,int i0,int j0,int xradius,
                               int yradius, int width0)
{

    //if some value in the compared patches is 0, then discard the match
    unsigned char valid=1;
    for (int s=-yradius; (s<= yradius) && valid; s++) {

        int l = (j0+s)*width0 + (i0-xradius);
        float *ptr0 = &u0[l];

        for(int r=-xradius; (r<=xradius) && valid; r++,ptr0++) {
            if (*ptr0 == 0) valid=0;
        }

    }

    return valid;
}

//assign label 1 to frames not containing occluded pixels
unsigned char *get_valid_frames(flimage *fpM, int tmin, int tmax, int x, int y,
                                int &nvalidt)
{
    unsigned char *validt=new unsigned char[tmax-tmin+1];
    nvalidt=0;
    for (int t=tmin; t <= tmax; t++) {
        //check for occluded pixels in patch
        validt[t-tmin]=check_zerovalues(fpM[t].v(), x, y, ixDWin, iyDWin, iWidth);
        if (validt[t-tmin]) nvalidt++;
    }
    return validt;
}


//Lines 8 and 9 of Algorithm 5: compute differences between 3D blocks (Algorithm 6)
//and sort the results
struct diffpatches *get_block_differences(cflimage *fpI, cflimage *fpA,
                                          unsigned char *validt, int x, int y,
                                          int tmin, int tmax, int imin, int imax,
                                          int jmin, int jmax)
{
    struct diffpatches *pdiff=new struct diffpatches[(imax-imin+1)*(jmax-jmin+1)];
    int ikt = 0;
    for(int j=jmin; j <= jmax; j++)
        for(int i=imin; i <= imax; i++,ikt++) {
            float fDif = 0.0f;
            for(int t=tmin; t <= tmax; t++) {
                if (validt[t-tmin]) {
                    if (useOracle)
                        for (int k=0; k < aChannels; k++)
                            fDif += fiL2FloatDist(fpA[t].v(k),fpA[t].v(k),x,y,
                                                  i,j,ixDWin,iyDWin,
                                                  iWidth,iWidth);
                    else
                        for (int k=0; k < iChannels; k++)
                            fDif += fiL2FloatDist(fpI[t].v(k),fpI[t].v(k),x,y,
                                                  i,j,ixDWin,iyDWin,
                                                  iWidth,iWidth);
                }
            }
            pdiff[ikt].diff=fDif;
            pdiff[ikt].i=i;
            pdiff[ikt].j=j;
        }
    
    //sort differences
    qsort(pdiff, (imax-imin+1)*(jmax-jmin+1), sizeof(struct diffpatches), order_increasing_diffs);
    
    return pdiff;
}

//Algorithm 5
void Get3DBlock(cflimage *fpI, cflimage *fpA, flimage *fpM, int x,int y,
                               laMatrix &X2, laMatrix &X02,
                               int *indexes, int &scurrent, int &npts)
{

    // Learning zone
    int tmin = MAX(tframe-iDTemp,0);
    int tmax = MIN(tframe+iDTemp,iFrames-1);

    int imin=MAX(x-iDBloc,ixDWin);
    int jmin=MAX(y-iDBloc,iyDWin);
    int imax=MIN(x+iDBloc,iWidth-1-ixDWin);
    int jmax=MIN(y+iDBloc,iHeight-1-iyDWin);

    // Spatiotemporal Block comparison with occlusion handling

    //check for which of the frames the reference block does not contain
    //occlusions
    int nvalidt;
    unsigned char *validt=get_valid_frames(fpM, tmin, tmax, x, y, nvalidt);

    //compute sorted block differences (Algorithm 6)
    struct diffpatches *pdiff=get_block_differences(fpI, fpA, validt, x, y,
                                                    tmin, tmax, imin, imax,
                                                    jmin, jmax);

    unsigned char addpatches=1;
    unsigned char validpatch;
    npts=0;
    scurrent=0;

    int iAmp = nvalidt;
    float fThresh = fFixedThrDist * fFixedThrDist * (float) p * (float) iAmp;

    for (int ikt=0; (ikt < (imax-imin+1)*(jmax-jmin+1)) && addpatches; ikt++) {
        int i=pdiff[ikt].i;
        int j=pdiff[ikt].j;
        float adif = pdiff[ikt].diff;

        for(int t=tmin; t<=tmax; t++) {
            if (validt[t-tmin]) {
                validpatch=check_zerovalues(fpM[t].v(), i,j,
                                            ixDWin,iyDWin,iWidth);
                if (validpatch) {
                    // Fill X2 with noisy patch
                    denoising_block_to_vector_nsq(fpI[t], X2[npts],
                                                  i-ixDWin, j-iyDWin,
                                                  ixWin, iyWin);
                    if (useOracle)
                        denoising_block_to_vector_nsq(fpA[t], X02[npts],
                                                      i-ixDWin, j-iyDWin,
                                                      ixWin, iyWin);

                    // Position in image of current patch
                    indexes[npts] = t * iWidth*iHeight + j*iWidth+i;

                    // Index of pixel being denoised in data matrix X2
                    if (t == tframe && y == j && x == i) scurrent = npts;

                    npts++;
                }

            }
        }
        if (npts >= ifK && adif > fThresh) addpatches=0;
    }

    delete[] validt;
    delete[] pdiff;
}

//Algorithm 8
//compute average and standard deviation of 3D block, per channel
//return average per channel,
//and average standard deviation of the three channels
void AveragePatches(laMatrix &XR2, const laMatrix & X2, int npts,
                   int p, int iChannels, float *vmean, float &stdAvg)
{
    int iSize = p / iChannels;      //! Number of pixel per channel
    double fValueAux = 0.0f;         //! Auxiliary
    
    //compute average of 3D block, per channel
    //! For each channel
    stdAvg=0.0f;
    for (int cc = 0; cc < iChannels; cc++) {
        double mean=0.0f, std=0.0f;
        //! For all patches
        for (int ii=0; ii < npts; ii++)
        for (int kk=0; kk < iSize ; kk++) {
            //! Average only values corresponding to cc channel
            fValueAux = (double) X2[ii][iChannels * kk + cc];
            mean += fValueAux;
            std += (fValueAux*fValueAux);
        }
        
        mean /= (double) (npts * iSize);
        std /= (double) (npts * iSize);
        std -= (mean*mean);
        std = sqrt(std);
        
        vmean[cc]=mean;
        stdAvg += (float) std;
    }
    
    stdAvg /= iChannels;
    
}
    
//set 3D block to average value
void SetMeanValue(laMatrix &XR2, const laMatrix & X2, int npts,
                      int p, int iChannels, float *vmean)
{
    
    for(int ii=0; ii < npts; ii++)
        for (int cc=0; cc < iChannels ; cc++)
            for (int kk=0; kk < iWin2; kk++)
                XR2[ii][ iChannels * kk + cc] = vmean[cc];

}

    

    
//Algorithm 7
void Denoise3DBlock(laMatrix &X2, laMatrix &X02, laMatrix &X2res,
                    laMatrix &X02res, laMatrix &XR2,
                    int npts, struct video_denoising_params &params)
{
    UpdateAuxiliaryVariables(X2, X02, X2res, X02res, XR2, npts);

    float *vmean = new float[iChannels];
    float stdAvg;
    AveragePatches(XR2, X2res, npts, p, iChannels, vmean, stdAvg); //Algorithm 8
    
    //Check if flat zone
    if (stdAvg < params.useFlatPar * params.fSigma ) { // is flat
        SetMeanValue(XR2, X2res, npts, p, iChannels, vmean);
    } else {
        int flagX02 = (params.useOracle)?(1):(0);
        //Algorithm 9
        NLPCAdenoising(X2res, X02res, XR2, flagX02, params.fSigma,
                       params.fRMult);
    }

    delete[] vmean;

}
    

void Aggregation(laMatrix &XR2, int *indexes, int npts, cflimage &fpO)
{
    // Put back denoised values into image
    for(int ii=0; ii < npts; ii++) {
        int index = indexes[ii];
        int ct = index / (iWidth*iHeight);
        index = index % (iWidth*iHeight);
        if (ct == tframe) {
            int ci = index % iWidth;
            int cj = index / iWidth;
            #pragma omp flush(fpO, mask)
            #pragma omp critical
            denoising_add_vector_to_block_nsq(fpO, XR2[ii], mask,
                                              ci-ixDWin, cj-iyDWin,
                                              ixWin, iyWin);

            // Current pixel denoised as center of a block
            // Mark center of patch so the patch is no longer used for
            // denoising the current frame
            #pragma omp flush(mask2)
            #pragma omp atomic
            mask2[cj * iWidth + ci] += 1.0;
        }
    }

}

void NormalizeOutput(cflimage *fpI, cflimage &fpO)
{
    // Normalize values
    for(int i=0; i < iWidth*iHeight; i++)
        if (mask[i] > 0.0) {
            for (int ic=0; ic < iChannels; ic++)
                fpO[ic * fpO.wh() + i] /= mask[i];

        } else {
            for (int ic=0; ic < iChannels; ic++)
                fpO[ic * fpO.wh() + i] = fpI[tframe][ic * fpO.wh() + i];

        }
}


void DenoiseFrame(cflimage *fpI, cflimage *fpA, flimage *fpM, cflimage &fpO,
                                     struct video_denoising_params &params)
    
{
    SetGlobalParameters(fpI, fpA, params);

    AllocateMemory();

    fpO=0.0f; // clear output

    #pragma omp parallel for schedule(dynamic)
    //for each image row
    for (int y = ixDWin; y < iHeight - ixDWin - 1; y++) {
        //Auxiliary variables
        laMatrix X2, X02, X2res, X02res, XR2;
        unsigned char mask2iszero;
        // scurrent: position of reference pixel in the X2 table
        int npts, scurrent;
        // indexes[k]: image coordinate (j*width+i) of block in k-entry of X2
        int *indexes = new int[bsize];
        InitializeAuxiliaryVariables(X2, X02);
        
        //for each pixel in row
        for (int x = iyDWin; x < iWidth - iyDWin - 1; x++) {
            #pragma omp flush(mask2)
            #pragma omp critical
            mask2iszero=(mask2[y*iWidth+x] == 0.0f);
            
            if (mask2iszero) { //not processed pixel
                
                //Algorithm 5:
                Get3DBlock(fpI,fpA,fpM,x,y, X2, X02, indexes, scurrent, npts);
                //Algorithm 7:
                Denoise3DBlock(X2, X02, X2res, X02res, XR2, npts, params);
                Aggregation(XR2, indexes, npts, fpO);
            }
            
        } // END for pixel x

        delete[] indexes;
    }  // END for pixel y

    NormalizeOutput(fpI, fpO);
    FreeMemory();
}




}
