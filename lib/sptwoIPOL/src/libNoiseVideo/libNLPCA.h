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

#ifndef _LIBDENOISING_H_
#define _LIBDENOISING_H_

#include "../library/libBasic.h"
#include "../library/libImage.h"

#define MNPCA 2
#define VERBOSE 0



namespace libUSTG
{

void  denoising_block_to_vector_nsq(cflimage &fpI, float *fpV, int i, int j, int ixWin, int iyWin);

void denoising_add_vector_to_block_nsq(cflimage &fpO, float *fpV, float *fpW,  int ci, int cj, int ixWin, int iyWin);

void  NLPCAdenoising(laMatrix &X2res, laMatrix &X02res, laMatrix &XR2, int flagX02, float fSigma, float fRMult);



}


#endif
