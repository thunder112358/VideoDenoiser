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

#ifndef _LIBDENOISINGVIDEO_H_
#define _LIBDENOISINGVIDEO_H_

#include "../library/libBasic.h"
#include "../library/libImage.h"
#include "libNLPCA.h"

#define MNPCA 2
#define VERBOSE 0


namespace libUSTG
{

struct video_denoising_params {
    int iTemp;
    int iWin;
    int iBloc;
    float fSigma;
    int ifK;
    unsigned char useOracle;
    float fRMult;
    float fFixedThrDist;
    float useFlatPar;
    int iFrames;
    int tframe;
};
    

void DenoiseFrame(cflimage *fpI, cflimage *fpA, flimage *fpM, cflimage &fpO,
                                     struct video_denoising_params &params);

}


#endif
