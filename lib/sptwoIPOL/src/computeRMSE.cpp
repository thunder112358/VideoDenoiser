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
	OptStruct oP = {"p:", 0,  "2", NULL,"L^p distance (p=0,1,2)"};  options.push_back(&oP);
	OptStruct oM = {"m:", 0,  NULL, NULL, "mask>0 points used in distance"};  options.push_back(&oM);
	OptStruct oV = {"v", 0,  NULL, NULL, "verbose mode"};  options.push_back(&oV);
	OptStruct oB = {"b:", 0,  "0", NULL, "boundary elimination "};  options.push_back(&oB);
	OptStruct oC = {"c:", 0,  "0", NULL, "distance for specific channels"};  options.push_back(&oC);
	
    
	
	vector<ParStruct *> parameters;
	ParStruct pinput = {"image", NULL, "image"}; parameters.push_back(&pinput);
	ParStruct pinput2 = {"image", NULL, "image"}; parameters.push_back(&pinput2);
	
	if (!parsecmdline("computeRMSE","\n## Distance l^p (default: RMSE)", argc, argv, options, parameters))
		return EXIT_FAILURE;
	
	
    
	// Parameters
	int p = atoi(oP.value);
    
    
	// Input
    libUSTG::cflimage input;
	input.load(pinput.value);
	
    
	
	// Input
	libUSTG::cflimage input2;
	input2.load(pinput2.value);
	
    
	
	// Mask
	libUSTG::flimage imask;
	if (oM.flag)
		imask.load(oM.value);
	
    
    
	// process
	assert(input.c() == input2.c() && input.w() == input2.w() && input.h() == input2.h());
	
	
	// remove boundary
	int boundary = atoi(oB.value);
    
	libUSTG::cflimage cinput = input.copy(boundary/2, boundary/2, input.w() - boundary, input.h()-boundary);
	libUSTG::cflimage cinput2 = input2.copy(boundary/2, boundary/2, input2.w() - boundary, input2.h()-boundary);
	
    
    
    
	libUSTG::flimage cimask;
	if (oM.flag)
		cimask = imask.copy(boundary/2, boundary/2, imask.w() - boundary, imask.h()-boundary);
    
	
    if (!oC.flag)
    {
        
        float fDist = 0.0f;
        for (int ii=0; ii < input.c(); ii++)
        {
            
            float fDif = 0.0f;
            
            if (oM.flag)
                fDif = libUSTG::fpDistLp(cinput.v(ii), cinput2.v(ii), cimask.v(), p, cinput.wh());
            
            else
                fDif = libUSTG::fpDistLp(cinput.v(ii), cinput2.v(ii), p, cinput.wh());
            
            fDist += fDif;
            
            
            if (oV.flag)
                printf("channel %d: %3.3f\n", ii, fDif);
        }
        
        
        fDist /= (float) cinput.c();
        
        
        if (oV.flag)
            printf("mean error: %3.3f\n",fDist);
        else
            printf("%3.3f\n", fDist);
    } else
    {
        
        int ii = atoi(oC.value);
        assert(ii>=0 && ii<input.c());
        
        float fDif = 0.0f;
        if (oM.flag)
			fDif = libUSTG::fpDistLp(cinput.v(ii), cinput2.v(ii), cimask.v(), p, cinput.wh());
		
		else
			fDif = libUSTG::fpDistLp(cinput.v(ii), cinput2.v(ii), p, cinput.wh());
		
        printf("channel %d: %3.3f\n", ii, fDif);
        
    }
    
	
	return EXIT_SUCCESS;
	
	
}
