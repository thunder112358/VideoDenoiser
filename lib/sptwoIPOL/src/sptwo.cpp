/*----------------------------------------------------------------------------

Copyright (c) 2016-2018 A. Buades and J.L. Lisani

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <stdbool.h>
#include <time.h>
#include <math.h>

#include "library/libImage.h"
#include "library/libBasic.h"
#include "libNoiseVideo/libFlow.h"
#include "libNoiseVideo/libDenoisingVideo.h"

//Default parameters for optical flow estimation
#define PAR_DEFAULT_NPROC   0
#define PAR_DEFAULT_TAU     0.25
#define PAR_DEFAULT_LAMBDA  0.15
#define PAR_DEFAULT_THETA   0.3
#define PAR_DEFAULT_NSCALES 5 //100
#define PAR_DEFAULT_ZFACTOR 0.5
#define PAR_DEFAULT_NWARPS  5 //5
#define PAR_DEFAULT_EPSILON 0.01 //0.01

using namespace std;

//Global variables
int width, height, nchannels, nframes;
unsigned char useOracle;
float lrdist, hBin, h02;


libUSTG::flimage *uImages, *vImages, *mwImages;
libUSTG::cflimage *iImages, *oImages, *aImages, *iwImages, *awImages;


//Read input data and check if they are color or gray level images
unsigned char ReadInputData(char *inputname)
{
    libUSTG::cflmovie input(inputname);

    nframes=input.n();

    //load first frame
    libUSTG::cflimage auximage = input.getframe(0);
    
    //get info and check color
    width = auximage.w();
    height = auximage.h();
    nchannels = auximage.c();
    
    unsigned char iscolor=(nchannels >= 3)?(1):(0);
    if (iscolor) {
        //check that the 3 channels are different
        iscolor=0;
        for (int n=0; (n < width*height) && !iscolor; n++) {
            if (auximage.v(0)[n] != auximage.v(1)[n]) iscolor=1;
            if (auximage.v(0)[n] != auximage.v(2)[n]) iscolor=1;
            if (auximage.v(1)[n] != auximage.v(2)[n]) iscolor=1;
        }
    }

    nchannels=(iscolor)?(3):(1);
    
    iImages = new libUSTG::cflimage[nframes];
    for (int ii=0; ii < nframes; ii++) {
        if (iscolor) iImages[ii] = input.getframe(ii);
        else {
            libUSTG::cflimage auximage = input.getframe(ii);
            iImages[ii].create(width, height, nchannels); //nchannels=1
            //copy first channel
            memcpy(iImages[ii].v(0), auximage.v(0), width*height*sizeof(float));
        }
    }

    return iscolor;
}

//write result in PNG format (8-bits)
void WriteOutput(char *outputname, int iFrame)
{
    for (int ii=0; ii < nframes; ii++) {
        if ((iFrame == -1) || (ii == iFrame)) { 
            string name(outputname);
            name += libUSTG::int2string(ii);
            name += ".png";
            oImages[ii].save(name.c_str());
        }
    }

}


void MemoryAllocation()
{

    //! Output
    oImages = new libUSTG::cflimage[nframes];

    aImages = new libUSTG::cflimage[nframes];

    for (int ii=0; ii < nframes; ii++) {
        oImages[ii] = iImages[ii];
        aImages[ii] = iImages[ii];
    }


    //! Flow memory
    uImages = new libUSTG::flimage[nframes];
    vImages = new libUSTG::flimage[nframes];
    for (int ii=0; ii < nframes; ii++)     uImages[ii].create(width, height);
    for (int ii=0; ii < nframes; ii++)     vImages[ii].create(width, height);

    //! Memory for warped stuff
    iwImages = new libUSTG::cflimage[nframes];
    awImages = new libUSTG::cflimage[nframes];

    mwImages = new libUSTG::flimage[nframes];

    for (int ii=0; ii < nframes; ii++)     iwImages[ii].create(width, height,
                                                               nchannels);

    for (int ii=0; ii < nframes; ii++)     awImages[ii].create(width, height,
                                                               nchannels);

    for (int ii=0; ii < nframes; ii++)     mwImages[ii].create(width, height);

}

void FreeMemory()
{
    delete[] uImages;
    delete[] vImages;
    delete[] iwImages;
    delete[] awImages;
    delete[] mwImages;

    delete[] iImages;
    delete[] aImages;
    delete[] oImages;
}


void SetGlobalParameters(libUSTGFLOW::flow_params &fparams,
                         libUSTG::video_denoising_params &dparams,
                         float lrdistTh,
                         float hBinTh, float factorocc,
                         int iBloc, int iWin, int iTemp, float factorflat,
                         int iKnn)
{

    //Optical flow parameters
    fparams.tau     = PAR_DEFAULT_TAU;
    fparams.theta   = PAR_DEFAULT_THETA;
    fparams.nscales = PAR_DEFAULT_NSCALES;
    fparams.zfactor = PAR_DEFAULT_ZFACTOR;
    fparams.warps  = PAR_DEFAULT_NWARPS;
    fparams.epsilon = PAR_DEFAULT_EPSILON;
    fparams.verbose = 0;
    fparams.iflagMedian = 0;
    const float N = 1 + log(hypot(width, height)/16.0) / log(1/fparams.zfactor);
    if (N < fparams.nscales) fparams.nscales = N;
    
    //Occlussion mask parameters
    lrdist = lrdistTh;
    hBin = hBinTh;
    h02 = factorocc * dparams.fSigma;
    h02*=h02;

    //Denoising parameters
    dparams.iFrames=nframes;
    dparams.iBloc = iBloc;
    dparams.iWin  = iWin;
    dparams.iTemp = iTemp;
    dparams.useFlatPar=factorflat;
    dparams.ifK = iKnn;
}


void SetParameters_1st_iteration(libUSTGFLOW::flow_params &fparams,
                                 libUSTG::video_denoising_params &dparams,
                                 float thdist, float lambda, float factorPCA)
{
    //strcpy(dparams.algoType,"pcaThr");
    dparams.useOracle = 0;

    dparams.fRMult =  factorPCA;
    dparams.fFixedThrDist = thdist;
    fparams.lambda = lambda;

    useOracle=dparams.useOracle;
}

void SetParameters_2nd_iteration(libUSTGFLOW::flow_params &fparams,
                                 libUSTG::video_denoising_params &dparams,
                                 float thdist, float lambda, float factorPCA)
{
    //strcpy(dparams.algoType,"pcaCWie");
    dparams.useOracle = 1;
    
    dparams.fRMult =  factorPCA;
    dparams.fFixedThrDist = thdist;
    fparams.lambda = lambda;

    useOracle=dparams.useOracle;

    //set oracle sequence = output of 1st iteration
    for (int ii=0; ii < nframes; ii++) aImages[ii]=oImages[ii];
}

//Compute flow from frame ii to frame jj
void ComputeFlow(int ii, int jj, libUSTGFLOW::flow_params &fparams)
{

    // reference frame
    libUSTG::flimage gray1;
    if (useOracle)
        gray1 = aImages[ii].getGray();
    else {
        gray1 = iImages[ii].getGray();
    }

    //target frame
    libUSTG::flimage gray2;
    if (useOracle)
        gray2 = aImages[jj].getGray();
    else {
        gray2 = iImages[jj].getGray();
    }

    //compute flow
    // compute optical flow in u[jj], v[jj]
    libUSTGFLOW::Dual_TVL1_optic_flow_multiscale(gray1.v(), gray2.v(),
                                                 uImages[jj].v(), vImages[jj].v(),
                                                 width, height, fparams);

}





//Warp image
void WarpImage(libUSTG::cflimage &in, libUSTG::flimage &xflow,
               libUSTG::flimage &yflow, libUSTG::cflimage *out)
{
    for (int ii=0; ii < nchannels; ii++) {
        float *pin = in.v(ii);
        float *pout = out->v(ii);
        for(int j=0; j < height; j++)
            for(int i=0; i < width; i++) {
                float p[2] = {i + xflow[j*width + i], j + yflow[j * width + i]};
                float value = libUSTGFLOW::bicubic_interpolation_at(pin, p[0],
                                                                    p[1], width,
                                                                    height,false);
                pout[j * width + i] = value;
            }
    }
}

void WarpFrame(int jj)
{
    // apply warping
    WarpImage(iImages[jj], uImages[jj], vImages[jj], &(iwImages[jj]));
    if (useOracle)
        WarpImage(aImages[jj], uImages[jj], vImages[jj], &(awImages[jj]));
}

//Check Left-Right coherence of the flow
void CheckLeftRightFlow(int ii, int jj, float *mask,
                        libUSTGFLOW::flow_params &fparams)
{
    // reference frame
    libUSTG::flimage gray1;
    if (useOracle)
        gray1 = aImages[ii].getGray();
    else {
        gray1 = iImages[ii].getGray();
    }

    //target frame
    libUSTG::flimage gray2;
    if (useOracle)
        gray2 = aImages[jj].getGray();
    else {
        gray2 = iImages[jj].getGray();
    }

    libUSTG::flimage ur(width, height);
    libUSTG::flimage vr(width, height);

    libUSTGFLOW::Dual_TVL1_optic_flow_multiscale( gray2.v(), gray1.v(),
                                                 ur.v(), vr.v(), width, height,
                                                 fparams);

    //get mask
    libUSTGFLOW::check_flow_reciprocity(uImages[jj].v(), vImages[jj].v(),
                                        ur.v(), vr.v(), mask, lrdist,
                                        width,  height);
}

//Algorithm 4:
//Occlusions mask, depending on left-right coherence and divergence of the flow
void GetOcclusionsMask(int ii, int jj, libUSTGFLOW::flow_params &fparams)
{
    //compute divergence
    libUSTG::flimage divergence(width,height);
    divergence=0.0f;
    libUSTGFLOW::divergence(uImages[jj].v(), vImages[jj].v(), divergence.v(),
                            width, height);

    // get mask
    for (int iky=0; iky < height; iky++) 
        for (int ikx=0; ikx < width; ikx++) {
            int iind = iky*width+ikx;

            //gray level term
            float fDif = distanceL2(iImages[ii], iwImages[jj],ikx,iky,ikx,iky);
            fDif /= (float) nchannels;
            mwImages[jj][iind] = expf(-fDif / h02);

            //divergence term
            float divV = divergence[iind];
            if (divV < 0) mwImages[jj][iind] *= exp(-divV*divV);

            //Mask binarization
            if (mwImages[jj][iind] > hBin) mwImages[jj][iind] = 1.0f;
            else mwImages[jj][iind] = 0.0f;
        }


    //check left-right coherence of the flow
    libUSTG::flimage lrmask(width, height);
    CheckLeftRightFlow(ii, jj, lrmask.v(), fparams);

    //intersection of masks
    for (int n=0; n < width*height; n++) mwImages[jj][n] *= lrmask[n];

}


//Algorithm 2: align frames and denoise
void denoise_function(libUSTGFLOW::flow_params &fparams,
                      libUSTG::video_denoising_params &dparams,
                      int iFrame)
{
    int iDTemp = (dparams.iTemp-1)/2;    // Half side of temporal block

    //do, for each frame
    for (int ii=0; ii < nframes; ii++) {
        if ((iFrame == -1) || (ii == iFrame)) { 
            // Define temporal neighborhood
            int jjmin = MAX(ii-iDTemp,0);
            int jjmax = MIN(ii+iDTemp,nframes-1);

            //Align frames (Algorithm 3)
            #pragma omp parallel for schedule(dynamic)
            for (int jj=jjmin; jj <= jjmax; jj++) {
                if (jj != ii) {
                    ComputeFlow(ii, jj, fparams);
                    WarpFrame(jj);
                    GetOcclusionsMask(ii, jj, fparams); //Algorithm 4
                } else {
                    //flow = 0 and no occlusions for reference frame
                    iwImages[ii] = iImages[ii];
                    if (useOracle)
                        awImages[ii] = aImages[ii];
                    mwImages[ii]=1.0f;
                }
            }

            // Denoise frame
            dparams.tframe=ii;
            libUSTG::DenoiseFrame(iwImages, awImages, mwImages, oImages[ii],
                                  dparams);
        }

    }


}




int main(int argc, char **argv)
{

    std::vector <OptStruct *> options;
    OptStruct osigma = 	{"s:", 0, "5.0", NULL, "noise standard deviation"};
    options.push_back(&osigma);
    OptStruct oiFrame = {"i:", 0, "-1", NULL, "frame to denoise (-1: denoise all frames)"};
    options.push_back(&oiFrame);
    
    //Parameters
    OptStruct oibloc = 	{"b:", 0, "12", NULL, "radius of search region"};
    options.push_back(&oibloc);
    OptStruct oiwin = 	{"w:", 0, "2", NULL, "radius of patch"};
    options.push_back(&oiwin);
    OptStruct oitemp = 	{"t:", 0, "7", NULL, "radius of temporal neighborhood"};
    options.push_back(&oitemp);
    OptStruct oiknn = 	{"k:", 0, NULL, NULL, "minimum number of patches (recommended: 55 gray images, 95 color images)"};
    options.push_back(&oiknn);
    OptStruct oflat = 	{"f:", 0, "0.85f", NULL, "flat parameter"};
    options.push_back(&oflat);
    OptStruct olrdist = {"c:", 0, "1.0f", NULL, "threshold for left-right coherence in occlusions mask"};
    options.push_back(&olrdist);
    OptStruct ohbin = 	{"h:", 0, "0.5f", NULL, "occlusion binarization threshold"};
    options.push_back(&ohbin);
    OptStruct ofocc = 	{"o:", 0, "5.5f", NULL, "occlusion factor"};
    options.push_back(&ofocc);
    OptStruct ofpca1 = 	{"p:", 0, "1.8f", NULL, "PCA factor 1st step"};
    options.push_back(&ofpca1);
    OptStruct ofpca2 = 	{"q:", 0, "1.45f", NULL, "PCA factor 2nd step"};
    options.push_back(&ofpca2);
    OptStruct odist1 = 	{"d:", 0, "0.0f", NULL, "3D blocks distances 1st step"};
    options.push_back(&odist1);
    OptStruct odist2 = 	{"e:", 0, "2.0f", NULL, "3D blocks distances 2nd step"};
    options.push_back(&odist2);
    OptStruct olambda1 = 	{"l:", 0, "0.075f", NULL, "optical flow lambda 1st step"};
    options.push_back(&olambda1);
    OptStruct olambda2 = 	{"m:", 0, "0.15f", NULL, "optical flow lambda 2nd step"};
    options.push_back(&olambda2);
   

    std::vector<ParStruct *> parameters;
    ParStruct pinput = {"input", NULL, "input file"};
    parameters.push_back(&pinput);
    ParStruct pout = {"out", NULL, "output file"};
    parameters.push_back(&pout);

    if (!parsecmdline("sptwo","Denoise video sequence with optical flow estimation",
                      argc, argv, options, parameters))
        return EXIT_FAILURE;

    //Set parameters
    
    libUSTGFLOW::flow_params fparams;
    libUSTG::video_denoising_params dparams;
    dparams.fSigma = atof(osigma.value);
    int iFrame = atoi(oiFrame.value);

    unsigned char iscolor=ReadInputData(pinput.value);
    
    int iBloc, iWin, iTemp, iKnn;
    float hBinTh, factorocc, factorflat, lrdistTh;
    iBloc = 2*atoi(oibloc.value)+1;
    iWin  = 2*atoi(oiwin.value)+1;
    iTemp = 2*atoi(oitemp.value)+1;
    
    if (oiknn.flag) iKnn = atoi(oiknn.value);
    else {
        iKnn=(iscolor)?(95):(55);
    }

    factorflat= atof(oflat.value);
    
    lrdistTh=atof(olrdist.value);
    hBinTh=atof(ohbin.value);
    factorocc=atof(ofocc.value);
    
    float thdist_1st, lambda_1st, factorPCA_1st;
    factorPCA_1st = atof(ofpca1.value);
    thdist_1st = atof(odist1.value);
    lambda_1st = atof(olambda1.value);
    
    float thdist_2nd, lambda_2nd, factorPCA_2nd;
    factorPCA_2nd = atof(ofpca2.value);
    thdist_2nd = atof(odist2.value);
    lambda_2nd = atof(olambda2.value);
    

    //Global Parameters
    MemoryAllocation();
    SetGlobalParameters(fparams, dparams, lrdistTh, hBinTh, factorocc, iBloc, iWin,
                        iTemp, factorflat, iKnn);

    //First Step
    SetParameters_1st_iteration(fparams, dparams, thdist_1st, lambda_1st,
                                factorPCA_1st);
    denoise_function(fparams, dparams, -1); //Algorithm 2

    //Second Step
    SetParameters_2nd_iteration(fparams, dparams, thdist_2nd, lambda_2nd,
                                factorPCA_2nd);
    denoise_function(fparams, dparams, iFrame); //Algorithm 2

    //Save Result
    WriteOutput(pout.value, iFrame);

    FreeMemory();

    return EXIT_SUCCESS;
}








