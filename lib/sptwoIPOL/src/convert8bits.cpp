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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "library/libImage.h"
#include "library/libBasic.h"

using namespace std;

int main(int argc, char **argv)
{
	
	
  vector <OptStruct *> options;
  OptStruct oC = {"c:", 0, "0", NULL, "clip option"}; options.push_back(&oC);
	vector<ParStruct *> parameters;
	ParStruct pinput = {"input", NULL, "float image"}; parameters.push_back(&pinput);
	ParStruct pout = {"out.png", NULL, "output PNG"}; parameters.push_back(&pout);
	
	
	if (!parsecmdline("convert8bits","Convert float image to 0-255 PNG image", argc, argv, options, parameters))
		return EXIT_FAILURE;

  unsigned char clip=(unsigned char) atoi(oC.value);
  
  //! Input
  libUSTG::cflimage input;
  input.load(pinput.value);
	
  libUSTG::cflimage output = input;
  
  //convert
  int nchannels=input.c();
  int size=input.wh();
  float *in, *out;
  in=input.v(0);
  float min=in[0];
  float max=in[0];
  for (int c = 0; c < nchannels; c++) {
    in=input.v(c);
    for (int n = 0; n < size; n++) {
      if (in[n] < min) min=in[n];
      if (in[n] > max) max=in[n];
    }
  }
  //printf("min=%2.2f  max=%2.2f\n", min, max);
  float v;
  for (int c = 0; c < nchannels; c++) {
    in=input.v(c);
    out=output.v(c);
    for (int n = 0; n < size; n++) {
      if (clip) {
        //clip values
        v=in[n];
        if (v < 0) v=0.0;
        if (v > 255) v=255.0;
      } else {
        //for each channel: linear transform from [min, max] to [0, 255]
        v=(in[n]-min)*255.0f/(max-min);
      }
      //png image value format: RRR..GGG...BBB...
      out[n]=(int) (v+0.5f);
    }
  }

	
	//! Save result
	output.save( pout.value);
	return EXIT_SUCCESS;
	
}
