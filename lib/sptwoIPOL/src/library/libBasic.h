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

#ifndef __library__
#define __library__


#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cmath>
#include <cassert>
#include <vector>
#include <getopt.h>

#include <fftw3.h>

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include <sys/types.h>
#include <unistd.h>


namespace libUSTG
{


#ifndef PI
#define PI 3.14159265358979323846264338327
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

#define dTiny 1e-10
#define fTiny 0.0000000001f
#define fLarge 10000000000.0f
#define dLarge 1e+10

#define MAX(i,j) ( (i)<(j) ? (j):(i) )
#define MIN(i,j) ( (i)<(j) ? (i):(j) )


//! Forward declarations
class laMatrix;
class laVector;

/*

void src_exit(const char* message);

*/

std::string int2string(int number);



////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//! Float Value operations
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////



void fpClear(float *fpI,float fValue, int iLength);
void fpCopy(float *fpI,float *fpO, int iLength);

float fpVar(float *u,int size);
float fpMean(float *u,int size);

float fpDistLp(float *u, float *v, int p, int size);
float fpDistLp(float *u, float *v, float *m, int p, int size);



////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//! Float pointer ordering
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////


void fpQuickSort(float *fpI, int iLength, int inverse = 0 );


////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//! Patch Statistics
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////


void fiPatchStatistics(float *fpIn, float *fpMinV, float *fpMaxV, float *fpMeanV, float *fpVarV, float *fpMedianV, float fRadius, int iWidth, int iHeight);


void fiPatchMedian(float *fpIn, float *fpMedianV, float fRadius, int iWidth, int iHeight);




////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//! Add Noise
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////


void fpAddNoiseGaussian(float *u, float *v, float std, long int randinit, int size);



////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//! Patch Distances
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////



//! Distances are not normalized by number of pixels or channels
float fiL2FloatDist(float *u0,float *u1,int i0,int j0,int i1,int j1,int xradius, int yradius, int width0, int width1);





////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//! Numerical Analysis Classes
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////


class laVector
{
protected:
    int d_n;	// size of array. upper index is nn-1
    float *d_v;
public:

    laVector();
    explicit laVector(int n);								// Zero-based array
    laVector(const float &a, int n);						// initialize to constant value
    laVector(const float *a, int n);						// Initialize to array

    laVector(const laVector &rhs);							// Copy constructor
    laVector & operator=(const laVector &rhs);				// assignment
    laVector & operator=(const float &a);					// assign a to every element

    float & operator[](const int i);				// i'th element
    const float & operator[](const int i) const;

    void create(int n);
    void erase();

    laVector copyBlock(int i0, int length);

    void sort(int decreasing_order);

    float * v();

    int size() const;
    ~laVector();

    friend class laMatrix;

};



class laMatrix
{

protected:
    int d_n;			// nrows
    int d_m;			// ncols
    float **d_v;		// pointer

public:


    //! Construction / Destruction

    laMatrix();
    laMatrix(int n, int m);
    laMatrix(const float &a, int n, int m);
    laMatrix(const float *a, int n, int m);
    laMatrix(const laMatrix &rhs);

    laMatrix & operator=(const laMatrix &rhs);
    laMatrix & operator=(const float &a);

    ~laMatrix();



    //! Basic operators

    float * operator[](const int i);	//subscripting: pointer to row i
    const float * operator[](const int i) const;

    float ** v();

    int nrows() const;
    int ncols() const;


    void create(int n, int m);



    //! Non member Arithmetic Operations

    friend laMatrix operator*  (float a, const laMatrix& rhs);                               // scalar matrix product
    friend laMatrix operator/  (const laMatrix& lhs, float a);                               // matrix scalar division
    friend laMatrix operator+  (const laMatrix& lhs, const laMatrix& rhs);                   // matrix sum
    friend laMatrix operator-  (const laMatrix& lhs, const laMatrix& rhs);                   // matrix subtraction
    friend laMatrix operator*  (const laMatrix& lhs, const laMatrix& rhs);                   // matrix product

    friend laVector operator*  (const laMatrix& lhs, const laVector & rhs);                  // matrix vector product



    //! Other
    laMatrix transposed();

    laMatrix copyBlock(int i0, int j0, int rowb, int colb);

    friend class laVector;

};

void compute_svd_double(laMatrix &Ain, laMatrix &m_Uin, laMatrix &m_Vin, laVector &m_Win);
void compute_pca_svd(laMatrix &X, laVector &S, laMatrix &V, laMatrix &U);




}



//
//! Parser
//



//! structure for parameters and options which are
//! optional in the code or they already have a default value
typedef struct optstruct {
    char *gp;           //! string of two letters  "a:"  as necessary for using getopt afterwards
    //! the ":" indicates that the activation of the option requires a value for the parameter
    //! and "a" that this option is activated by "-a " in the command

    int flag;           //! flag indicating that the option has been activated

    char *defvalue;     //! default value for this parameter if not modified by console

    char *value;        //! value of the associated parameter to current option

    char *comment;      //! comment that appears by console

} OptStruct;



//! structure for necessary parameters of the method
typedef struct parstruct {
    char * name;
    char * value;       //! value of the parameter
    char * comment;     //! comment that appears by console

} ParStruct;



int parsecmdline(char *pname,
                 char *function,
                 int argc, char **argv,
                 std::vector <OptStruct*> & opt,
                 std::vector <ParStruct*> & par);



#endif
