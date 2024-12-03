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
	OptStruct oa = {"a:", 0,  "1.0", NULL,"a value"};  options.push_back(&oa);	
	
	vector<ParStruct *> parameters;
	ParStruct pinput = {"image", NULL, "image"}; parameters.push_back(&pinput);
	ParStruct pout = {"out", NULL, "output file"}; parameters.push_back(&pout);
	
	
	if (!parsecmdline("addGaussianNoise","\n## Adds Gaussian noise of std a", argc, argv, options, parameters))
		return EXIT_FAILURE;

		
	//! Input
	libUSTG::cflimage auximage;
	auximage.load(pinput.value);
	

    //get info and check color
    int width = auximage.w();
    int height = auximage.h();
    int nchannels = auximage.c();
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


    libUSTG::cflimage input;
    if (!iscolor && (nchannels > 1)) {
        //gray image with 3 channels -> convert it to 1 channel image
        input.create(width, height, 1); //nchannels=1
        //copy first channel
        memcpy(input.v(0), auximage.v(0), width*height*sizeof(float));
    } else {
        input = auximage;
    }
    
	
	//! Parameters
	float a = atof(oa.value);
	
	input.addGaussianNoise(a);
	
	
	//! Save result
	input.save( pout.value); 
	return EXIT_SUCCESS;
	
}
