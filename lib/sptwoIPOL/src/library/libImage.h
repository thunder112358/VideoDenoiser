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

#ifndef __libImage__
#define __libImage__



extern "C" {
#include "libImageFormats.h"
}


#define USE_FFTW
#ifdef USE_FFTW
#include <fftw3.h>
#endif


#include "libBasic.h"
#include "libImageFormatPM.h"


#define iipHorizontal 0
#define iipVertical 1



namespace libUSTG
{


//
//! Class definitions
//

class cflimage;
class flimage;


////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//! Image Classes
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////


class cflimage
{


protected:



public:

    int d_c, d_w, d_h, d_wh, d_whc;		// channels, width, height, width*height, width*height*channels
    float *d_v;							// pointer


    char name[128];
    float visuMin, visuMax;				// useful only for visualization with mxview



    //
    //! Construction / Memory
    //

    cflimage();
    cflimage(int w, int h, int c);
    cflimage(int w, int h, float *igray);
    cflimage(int w, int h, float *ired, float *igreen, float *iblue);

    cflimage(int w, int h, unsigned char *igray);
    cflimage(int w, int h, unsigned char *ired, unsigned char *igreen, unsigned char *iblue);

    cflimage(const flimage &red, const flimage &green, const flimage &blue);


    void create(int w, int h, int c);

    cflimage(const cflimage& im);

    void erase();
    virtual ~cflimage();



    //
    //! Operators
    //

    cflimage& operator= (const cflimage& im);



    //
    //! Load / Save
    //
    void load(const char *filename);
    void save(const char *filename);



    //
    //! Get Basic Data
    //

    int c() const
    {
        return d_c;
    }
    int w() const
    {
        return d_w;
    }
    int h() const
    {
        return d_h;
    }
    int wh() const
    {
        return d_wh;
    }
    int whc() const
    {
        return d_whc;
    }


    float * v() const
    {
        return d_v;
    }
    float * v(int i) const
    {
        return &d_v[i * d_wh];
    }

    inline float operator[](int i) const
    {
        return d_v[i];
    }
    inline float& operator[](int i)
    {
        return d_v[i];
    }



    flimage getChannel(int i);
    void getChannel(int i, flimage *out);

    operator  flimage();



    int isSameSize(const  cflimage &inIm);

    char * getName()
    {
        return name;
    }
    void  setName(const char *iname)
    {
        strcpy(name, iname);
    }


    float getVisuMin()
    {
        return visuMin;
    }
    float getVisuMax()
    {
        return visuMax;
    }

    void setVisuMin(float a)
    {
        visuMin = a;
    }
    void setVisuMax(float a)
    {
        visuMax = a;
    }


    //
    //! Math
    //

    cflimage& operator= (float a);
    void operator-=(float a);
    void operator+=(float a);
    void operator*=(float a);

    float  min();
    float  max();

    float  min_channel (int i);
    float  max_channel (int i);

    void  max (float M);
    void  min (float m);

    void normalizeL1();

    void rint();
    void abs();

    void thre(float m, float M);



    //
    //! Color Conversion
    //

    flimage getGray();
    flimage getGray(float wr, float wg, float wb);

    void getGray(flimage * out);
    void getGray(float wr, float wg, float wb, flimage *out);

    int isGray();

    void Rgb2Yuv(int iflagOrto);
    void Yuv2Rgb(int iflagOrto);

    cflimage binarize(float value, int inverse);
    void binarize(float value, int inverse, cflimage *out);


    //
    //! Block operations
    //

    cflimage copy(int ipx, int ipy, int iw, int ih);
    void paste(const cflimage &im, int x, int y);

    cflimage padding(int w, int h, float fValue);
    cflimage append(const cflimage &imIn, int extension);



    //
    //! Value operations
    //


    void addGaussianNoise(float std);


    //
    //! Patch Processing
    //


    friend float distanceL2(const cflimage &input1,const  cflimage &input2, int ipx, int ipy, int iqx, int iqy);





};



class flimage : public cflimage
{

public:


    flimage();
    flimage(int w, int h);
    flimage(int w, int h, float *ptr);
    flimage(int w, int h, unsigned char *ptr);

    flimage(const flimage& im);

    void load(const char* filename);


    void create(int w, int h);


    using cflimage::operator=;
    flimage& operator=(const flimage& im);


    float operator()(int i, int j) const
    {
        return this->d_v[j * this->d_w + i];
    }
    float& operator()(int i, int j)
    {
        return this->d_v[j * this->d_w + i];
    }



    virtual ~flimage() {};

};









////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//! Movie Classes
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////

class cflmovie
{

public:


    // Constructor && Destructor
    cflmovie();
    cflmovie(const char * ifilename);							/// Reads a cflmovie file
    cflmovie(const char * ifilename,  int inframes);   			/// Opens a cflmovie file in writable format


    cflmovie(const cflmovie & imovie);

    ~cflmovie();


    cflmovie& operator= (const cflmovie& im);


    // Load
    void load(const char * ifilename);											/// Reads a cflmovie file


    // Main variables
    int n()
    {
        return d_n;
    }



    // Getting image at position fpos
    cflimage getframe(int fpos);


    // Write image into film
    void write(cflimage &frame);


private:

    int d_n;

    char* filename;
    char* fileadress;				/// Used for reading/writting image names in format .mov

    bool  writable;				/// Is current cflmovie writable
    char ** strings;   			/// In case we read each time the image from the file

    int pos;        				/// Current position for writting
    std::ofstream outfile;     		/// Output mov file for writing

};






}





#endif


