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

#include "libImage.h"


namespace libUSTG
{



//////////////////////////////////////////////////////////////////
//! Begin Cflimage
//////////////////////////////////////////////////////////////////

//! Constructors
cflimage::cflimage() : d_c(0), d_w(0), d_h(0), d_wh(0), d_whc(0), d_v(0),  visuMin(0.0f), visuMax(255.f)
{
}



cflimage::cflimage(int w, int h, int c) : d_c(c), d_w(w), d_h(h),d_wh(w*h), d_whc(c*w*h), d_v(new float[c*w*h]),  visuMin(0.0f), visuMax(255.f)
{
    for (int ii=0; ii < d_whc; ii++) d_v[ii] = 0.0f;

}



cflimage::cflimage(int w, int h, float *igray) : d_c(1), d_w(w), d_h(h),d_wh(w*h), d_whc(w*h), d_v(new float[w*h]), visuMin(0.0f), visuMax(255.f)
{
    memcpy(d_v, igray, w * h * sizeof(float));
}


cflimage::cflimage(int w, int h, unsigned char *igray) : d_c(1), d_w(w), d_h(h),d_wh(w*h), d_whc(w*h), d_v(new float[w*h]), visuMin(0.0f), visuMax(255.f)
{
    for (int ii=0; ii < d_wh; ii++) d_v[ii] = (float) igray[ii];
}


cflimage::cflimage(int w, int h, float *ired, float *igreen, float *iblue) : d_c(3), d_w(w), d_h(h),d_wh(w*h), d_whc(3*w*h), d_v(new float[3*w*h]),  visuMin(0.0f), visuMax(255.f)
{
    memcpy(d_v, ired, w * h * sizeof(float));
    memcpy(d_v + w*h, igreen, w * h * sizeof(float));
    memcpy(d_v + 2*w*h, iblue, w * h * sizeof(float));
}


cflimage::cflimage(int w, int h, unsigned char *ired, unsigned char *igreen, unsigned char *iblue) : d_c(3), d_w(w), d_h(h),d_wh(w*h), d_whc(3*w*h), d_v(new float[3*w*h]),  visuMin(0.0f), visuMax(255.f)
{
    for (int ii=0; ii < d_wh; ii++) {
        d_v[ii] = (float) ired[ii];
        d_v[d_wh + ii] = (float) igreen[ii];
        d_v[2*d_wh + ii] = (float) iblue[ii];
    }
}




cflimage::cflimage(const flimage &red, const  flimage &green,const  flimage &blue) : d_c(0), d_w(0), d_h(0), d_wh(0), d_whc(0), d_v(0),  visuMin(0.0f), visuMax(255.f)
{

    assert(red.d_w == green.d_w && green.d_w == blue.d_w);
    assert(red.d_h == green.d_h && green.d_h == blue.d_h);

    d_w = red.d_w;
    d_h = red.d_h;
    d_c = 3;

    d_wh = d_w * d_h;
    d_whc = d_wh * d_c;

    d_v = new float[3 * d_wh];
    memcpy(d_v, red.d_v, d_wh * sizeof(float));
    memcpy(d_v + d_wh, green.d_v, d_wh * sizeof(float));
    memcpy(d_v + 2 * d_wh, blue.d_v, d_wh * sizeof(float));


}


cflimage::cflimage(const cflimage& im) : d_c(im.d_c), d_w(im.d_w), d_h(im.d_h),d_wh(im.d_wh), d_whc(im.d_whc), d_v(0),  visuMin(0.0f), visuMax(255.f)
{

    if (d_whc > 0) {
        d_v = new float[d_whc];
        memcpy((void *) d_v, (const void *) im.d_v, d_whc * sizeof(float));

    }

}



void cflimage::create(int w, int h, int c)
{
    erase();
    d_c = c;
    d_w = w;
    d_h=h;
    d_wh = w*h;
    d_whc = c*w*h;
    d_v = new float[d_whc];
    visuMin=0.0f;
    visuMax=255.0f;

    for (int ii=0; ii < d_whc; ii++) d_v[ii] = 0.0f;
}



void cflimage::erase()
{
    d_w = d_h = d_wh = d_whc = 0;
    if (d_v) delete[] d_v;
    d_v=0;
}




cflimage::~cflimage()
{
    erase();
}






//////////////////////////////////////////////////////////////////
//! Begin Operators
//////////////////////////////////////////////////////////////////




cflimage&  cflimage::operator= (const cflimage& im)
{
    if (&im == this) {
        return *this;
    }


    if (d_c != im.d_c || d_w != im.d_w || d_h != im.d_h) {
        erase();
        d_c = im.d_c;
        d_w = im.d_w;
        d_h=im.d_h;
        d_wh = d_w * d_h;
        d_whc=d_c * d_w * d_h;
        d_v = new float[d_whc];
    }


    memcpy((void *) d_v, (const void *) im.d_v, d_whc * sizeof(float));

    visuMin = im.visuMin;
    visuMax = im.visuMax;

    return *this;

}








//////////////////////////////////////////////////////////////////
//! Begin Load/Save
//////////////////////////////////////////////////////////////////



void cflimage::load(const char *filename)
{



    // erase current image
    erase();



    // look if it exists
    std::ifstream fileImage(filename);
    if (!fileImage.good()) {
        std::cout << "... failed to read image " << filename << std::endl;
        exit(-1);

    } else fileImage.close();


    // try pmf format
    d_v = _load_pm(filename, d_c, d_w, d_h);
    if (d_v) {

        d_wh = d_w * d_h;
        d_whc = d_c * d_w * d_h;
        return;
    }



    // the sure bet
    d_v = iio_read_image_float_split(filename, &d_w, &d_h, &d_c);

    if (d_v) {
        if (d_c == 2) d_c = 1;
        //if (d_c > 3)  d_c = 3;

        d_wh = d_w * d_h;
        d_whc = d_c * d_w * d_h;



        return;
    }




    std::cout << "... failed to read image " << filename << std::endl;
    exit(-1);


}



void cflimage::save(const char *filename)
{

    if (!d_v) {
        std::cout << "... failed to save image " << filename << std::endl;
        exit(-1);
    }

    std::string strname(filename);
    size_t pos = strname.find_last_of ('.');
    std::string extension = strname.substr(pos+1);

    if ( extension == "png" || extension == "PNG" || extension == "tif" || extension == "TIF" || extension == "tiff" || extension == "TIFF" ) {

        iio_save_image_float_split((char*)filename, d_v, d_w, d_h, d_c);

    } else {

        if (_create_pm(filename,d_c,d_w,d_h,d_v) == 0) {
            std::cout << "... failed to save pmf image " << filename << std::endl;
            exit(-1);

        } else return;

    }

}






//////////////////////////////////////////////////////////////////
//! Begin Get Basic Data
//////////////////////////////////////////////////////////////////




cflimage::operator  flimage()
{
    return getGray();
}





int  cflimage::isSameSize(const  cflimage &inIm)
{

    if (d_c != inIm.d_c || d_w != inIm.d_w || d_h != inIm.d_h) return 0;
    else return 1;

}






flimage cflimage::getChannel(int i)
{

    assert(i < d_c);

    flimage image(d_w,d_h);

    for (int jj=0; jj < d_wh; jj++) image.d_v[jj] = d_v[ i * d_wh + jj];

    return image;
}



void cflimage::getChannel(int i, flimage *out)
{

    assert(i < d_c);
    assert(d_v != NULL);

    if (out->d_w != d_w || out->d_h != d_h) {
        out->erase();
        out->create(d_w, d_h);
    }


    for (int jj=0; jj < d_wh; jj++) out->d_v[jj] = d_v[ i * d_wh + jj];


}




//////////////////////////////////////////////////////////////////
//! Begin Math
//////////////////////////////////////////////////////////////////



cflimage& cflimage::operator= (float a)
{
    if (d_v)  for (int j=0; j<d_whc ; j++) d_v[j] = a;

    return *this;
}


void cflimage::operator-= (float a)
{
    if (d_v)  for (int j=0; j<d_whc ; j++) d_v[j] -= a;

}


void cflimage::operator+= (float a)
{
    if (d_v)  for (int j=0; j<d_whc ; j++) d_v[j] += a;

}


void cflimage::operator*= (float a)
{
    if (d_v)  for (int j=0; j<d_whc ; j++) d_v[j] *= a;

}




float  cflimage::max ()
{

    assert(d_v != NULL);

    float fmax = d_v[0];
    for (int j=0; j < d_whc ; j++)  if (d_v[j] > fmax)  fmax = d_v[j];
    return fmax;


}



float  cflimage::min ()
{
    assert(d_v != NULL);

    float fmin = d_v[0];
    for (int j=0; j < d_whc ; j++)  if (d_v[j] < fmin)  fmin = d_v[j];
    return fmin;
}



float  cflimage::min_channel (int i)
{
    assert(d_v != NULL);
    assert(i>= 0 && i < d_c);

    float *ptr = &d_v[i * d_wh];
    float fmin = *ptr;

    for (int j=0; j < d_wh ; j++, ptr++)  if (*ptr < fmin)  fmin = *ptr;
    return fmin;
}



float  cflimage::max_channel (int i)
{
    assert(d_v != NULL);
    assert(i>= 0 && i < d_c);

    float *ptr = &d_v[i * d_wh];
    float fmax = *ptr;

    for (int j=0; j < d_wh ; j++, ptr++)  if (*ptr > fmax)  fmax = *ptr;
    return fmax;
}



void  cflimage::min (float m)
{
    assert(d_v != NULL);
    for (int j=0; j < d_whc ; j++)  if (d_v[j] > m)  d_v[j] = m;
}



void  cflimage::max (float M)
{
    assert(d_v != NULL);
    for (int j=0; j < d_whc ; j++)  if (d_v[j] < M)  d_v[j] = M;
}




void cflimage::thre(float m, float M)
{

    assert(d_v != NULL);
    for (int ii=0; ii < d_whc; ii++) {
        if (d_v[ii] >= M) 	d_v[ii]= M;
        else if (d_v[ii] <= m)  d_v[ii]= m;
    }

}



void cflimage::normalizeL1()
{

    float fSum = 0.0f;
    for (int j=0; j < d_whc ; j++)    fSum += d_v[j];

    assert(fSum != 0.0f);
    float dfSum = 1.0 / fSum;
    for (int j=0; j < d_whc ; j++)    d_v[j] *= dfSum;

}



void cflimage::rint()
{

    assert(d_v != NULL);
    for (int ii=0; ii < d_whc; ii++)		d_v[ii]= rintf(d_v[ii]);

}


void cflimage::abs()
{

    assert(d_v != NULL);
    for (int ii=0; ii < d_whc; ii++)		d_v[ii]= fabsf(d_v[ii]);

}





///////////////////////////////////////////////
//! Begin Color Conversion
///////////////////////////////////////////////


flimage cflimage::getGray()
{



    flimage image(d_w, d_h);
    image=0.0f;

    getGray(&image);

    return image;
}





void cflimage::getGray(flimage * out)
{




    assert(d_v != NULL);
    if (out->d_w != d_w || out->d_h != d_h) {
        out->erase();
        out->create(d_w, d_h);
    }


    for (int i=0; i < d_whc; i++) {
        out->d_v[i % d_wh] += d_v[i];
    }

    for (int i=0; i < d_wh; i++)
        out->d_v[i] /= (float) d_c;

}



flimage cflimage::getGray(float wr, float wg, float wb)
{

    //assert(d_c == 1  || d_c == 3);

    flimage image(d_w, d_h);
    image=0.0f;

    getGray(wr,wg,wb,&image);

    return image;
}





void cflimage::getGray(float wr, float wg, float wb, flimage *out)
{
    //assert(d_c == 1  || d_c == 3);
    assert(d_v != NULL);

    if (out->d_w != d_w || out->d_h != d_h) {
        out->erase();
        out->create(d_w, d_h);
    }


    if (d_c == 1)  *out = (flimage) (*this);
    else {
        for (int i=0; i < d_wh; i++) out->d_v[i] = wr * d_v[i] + wg * d_v[d_wh + i] + wb * d_v[2*d_wh+i];
    }


}



int cflimage::isGray()
{
    if (d_c == 1) return 1;

    for (int i=0; i < d_wh; i++) {
        float value = d_v[i];
        for (int j=1; j < d_c; j++)
            if ( d_v[ j * d_wh + i] != value  ) return 0;
    }

    return 1;

}

cflimage cflimage::binarize(float value, int inverse)
{
    assert(d_v != NULL);
    cflimage binary(d_w,d_h,d_c);

    binarize(value, inverse, &binary);

    return binary;
}




void cflimage::binarize(float value, int inverse, cflimage *out)
{

    assert(d_v != NULL);
    if (!isSameSize(*out)) {
        out->erase();
        out->create(d_w, d_h, d_c);
    }


    for (int ii=0; ii < d_whc; ii++) {
        if (d_v[ii] >= value && !inverse) 	out->d_v[ii]= 1.0;
        else if (d_v[ii] < value && inverse)  out->d_v[ii]= 1.0;
        else out->d_v[ii]= 0.0;

    }

}









void cflimage::Rgb2Yuv(int iflagOrto)
{

    assert(d_c==3);
    float *r,*g,*b;
    r=d_v;
    g=&d_v[d_wh];
    b=&d_v[2*d_wh];

    float vr, vg,vb;

    if (iflagOrto) {
        for(int i=0; i<d_wh; i++) {
            vr = r[i];
            vg = g[i];
            vb = b[i];

            r[i] =    0.577350 * vr + 0.577350 * vg + 0.577350  * vb;
            g[i] =    0.707106 * vr					- 0.707106  * vb;
            b[i] =    0.408248 * vr - 0.816496 * vg	+ 0.408248  * vb;
        }
    } else {
        for(int i=0; i<d_wh; i++) {
            vr = r[i];
            vg = g[i];
            vb = b[i];

            r[i] = ( 0.299 *  vr + 0.587 * vg + 0.114 * vb);
            g[i] =  ( vr - r[i]);
            b[i] =  ( vb - r[i]);

        }
    }

}






void  cflimage::Yuv2Rgb(int iflagOrto)
{

    assert(d_c==3);
    assert(d_c==3);
    float *r,*g,*b;
    r=d_v;
    g=&d_v[d_wh];
    b=&d_v[2*d_wh];

    float vr, vg,vb;

    if (iflagOrto) {
        for(int i=0; i<d_wh; i++) {
            vr = r[i];
            vg = g[i];
            vb = b[i];

            r[i] =    0.577350 * vr + 0.707106 * vg + 0.408248  * vb;
            g[i] =    0.577350 * vr					- 0.816496  * vb;
            b[i] =    0.577350 * vr - 0.707106 * vg	+ 0.408248  * vb;
        }
    } else {
        for(int i=0; i<d_wh; i++) {
            vr = r[i];
            vg = g[i];
            vb = b[i];

            g[i] =  ( vr - 0.299 * (vg + vr) - 0.114 * (vb +  vr) ) / 0.587;
            r[i] =  ( vg + vr);
            b[i] =  ( vb + vr);

        }
    }



}









//////////////////////////////////////////////////////////////////
//! Begin Block Operations
//////////////////////////////////////////////////////////////////



cflimage cflimage::padding(int w, int h, float fValue)
{

    assert(w >= d_w  && h >= d_h);

    cflimage image(w,h,d_c);
    image=fValue;

    for (int ii=0; ii < d_c; ii++) {

        for(int j=0; j < d_h; j++)
            for(int i=0; i < d_w; i++)
                image.d_v[ii * image.d_wh + j * image.d_w + i] = d_v[ ii * d_wh + j * d_w + i];

    }

    return image;
}






cflimage cflimage::copy(int ipx, int ipy, int iw, int ih)
{

    assert(iw>0 && ih>0);
    assert(ipx>=0 && ipy>=0);
    assert(ipx + iw - 1 < d_w && ipy + ih - 1 < d_h);

    cflimage image(iw, ih, d_c);

    int nn=0;
    for (int ii=0; ii < d_c; ii++) {

        int l = ii * d_wh +  ipy * d_w + ipx;

        for (int jj = 0; jj < ih; jj++) {

            for (int kk = 0; kk < iw; kk++,nn++,l++)
                image[nn] = d_v[l];

            l += d_w - iw;
        }


    }

    return image;
}




void cflimage::paste(const cflimage &im, int ipx, int ipy)
{

    assert(ipx>=0 && ipy>=0);
    assert(ipx + im.d_w - 1 < d_w && ipy + im.d_h - 1 < d_h);
    assert(d_c == im.d_c);


    for (int ii=0; ii < d_c; ii++) {


        int ll = ii * im.d_wh;
        int nn = ii * d_wh +  ipy * d_w + ipx;


        for (int jj = 0; jj < im.d_h; jj++) {

            for (int kk = 0; kk < im.d_w; ll++, nn++, kk++)
                d_v[nn] = im.d_v[ll];

            nn += d_w - im.d_w;
        }


    }

}






cflimage cflimage::append(const  cflimage &imIn, int extension)
{

    assert(d_c == imIn.d_c);

    //! create image
    cflimage image;
    if (extension == iipHorizontal)
        image.create(d_w + imIn.d_w, MAX(d_h, imIn.d_h), d_c);
    else
        image.create(MAX(d_w, imIn.d_w), d_h + imIn.d_h, d_c);


    image = 0.0f;
    image.paste(*this, 0, 0);
    if (extension == iipHorizontal)
        image.paste(imIn, d_w, 0);
    else
        image.paste(imIn, 0, d_h);

    return image;
}





//////////////////////////////////////////////////////////////////
//! Begin value operations
//////////////////////////////////////////////////////////////////



void  cflimage::addGaussianNoise(float std)
{
    assert(d_v != NULL);
    fpAddNoiseGaussian(d_v, d_v, std, 0, d_whc);
}


//////////////////////////////////////////////////////////////////
//! Begin block processing
//////////////////////////////////////////////////////////////////


float distanceL2(const  cflimage &input1, const  cflimage &input2, int ipx, int ipy, int iqx, int iqy)
{

    assert(ipx >=0 && ipy>=0 && iqx>=0 && iqy>=0);
    assert(ipx < input1.w() && ipy  < input1.h() && iqx  < input2.w() && iqy  < input2.h() );

    float fDif = 0.0f;
    float fDist = 0.0f;
    for (int ii=0; ii < input1.d_c; ii++) {

        float *ptr1 = &input1.d_v[ii * input1.d_wh +  ipy  * input1.d_w + ipx];
        float *ptr2 = &input2.d_v[ii * input2.d_wh +  iqy  * input2.d_w + iqx];

        fDif = *ptr1 - *ptr2;
        fDist += fDif * fDif;

    }

    return fDist;
}






/// Class flimage

flimage::flimage() : cflimage()
{
}




flimage::flimage(int w, int h) : cflimage(w, h, 1)
{
}


flimage::flimage(int w, int h, float *ptr) : cflimage(w, h, 1)
{
    memcpy(this->d_v, ptr, w * h * sizeof(float));
}

flimage::flimage(int w, int h, unsigned char *ptr) : cflimage(w, h, ptr)
{
}




flimage::flimage(const flimage& im)
    : cflimage(im)
{}




flimage& flimage::operator=(const flimage & im)
{
    cflimage::operator=(im);
    return *this;
}




void flimage::create(int w, int h)
{
    cflimage::create(w,h,1);
}


void flimage::load(const char* filename)
{
    cflimage image;
    image.load(filename);
    image.getGray(this);
}





cflmovie::cflmovie():
    d_n(0), filename(NULL), fileadress(NULL), writable(0), strings(NULL), pos(-1), outfile(0)
{}




cflmovie::cflmovie(const char* ifilename,  int inframes)
    :  d_n(inframes), filename(NULL), fileadress(NULL), writable(0), strings(NULL), pos(-1), outfile(0)
{

    filename = new char[128];
    strcpy(filename,ifilename);


    outfile.open(filename);
    if (!outfile.is_open()) {
        printf("cflmovie file not writtable\n");
        exit(-1);
    }

    outfile << "cflmovie" << std::endl;
    outfile << d_n << std::endl;


    writable = true;
    pos = -1;

}



cflmovie::cflmovie(const cflmovie & imovie)
    :  d_n(imovie.d_n), filename(NULL), fileadress(NULL), writable(imovie.writable), strings(NULL), pos(imovie.pos), outfile(0)
{

    strcpy(filename, imovie.filename);
    strcpy(fileadress, imovie.fileadress);

    if (imovie.strings) {
        strings = new char*[d_n];
        for (int i=0; i < d_n; i++) {
            strcpy(strings[i], imovie.strings[i]);
        }

    }


}



cflmovie::cflmovie(const char * ifilename)
    :  d_n(0), filename(NULL), fileadress(NULL), writable(0), strings(NULL), pos(-1), outfile(0)
{

    /// Reading
    writable = false;
    pos = -1;
    fileadress = new char[128];
    filename = new char[128];

    strcpy(filename,ifilename);



    std::ifstream file(ifilename);
    if (!file.is_open()) {
        printf("cflmovie %s file not found or impossible to open\n", ifilename);
        exit(-1);
    }

    /// Reading format
    char  iinfo[128];
    file >> iinfo;

    if( strcmp(iinfo, "cflmovie") == 0) {

        file >> d_n;


        strings = new char*[d_n];
        for(int i=0; i < d_n; i++) strings[i] = new char[128];

        for(int i = 0; i < d_n ; i++ ) {
            file >>iinfo;
            strcpy(strings[i],iinfo);

        }



    } else {

        printf("Error: not a correct cflmovie list file\n");
        exit(-1);

    }

}


void cflmovie::load(const char * ifilename)
{

    /// Reading
    writable = false;
    pos = -1;
    fileadress = new char[128];
    filename = new char[128];

    strcpy(filename,ifilename);



    std::ifstream file(ifilename);
    if (!file.is_open()) {
        printf("cflmovie %s file not found or impossible to open\n", ifilename);
        exit(-1);
    }

    /// Reading format
    char  iinfo[128];
    file >> iinfo;

    if( strcmp(iinfo, "cflmovie") == 0) {

        file >> d_n;

        strings = new char*[d_n];
        for(int i=0; i < d_n; i++) strings[i] = new char[128];

        for(int i = 0; i < d_n ; i++ ) {
            file >>iinfo;
            strcpy(strings[i],iinfo);

        }



    } else {

        printf("Error: not a correct cflmovie list file\n");
        exit(-1);

    }

}





cflmovie& cflmovie::operator= (const cflmovie& im)
{
    printf("warning :: using cflmovie operator = which is not defined properly\n");

    if (&im == this) {
        return *this;
    }

    return *this;
}




cflmovie::~cflmovie()
{

    if (outfile.is_open()) outfile.close();

}





cflimage cflmovie::getframe(int fpos)
{


    cflimage image;


    image.load(strings[fpos]);
    return image;

 
}





void cflmovie::write(cflimage &frame)
{

    if (writable) {

        pos++;

        char* imfilename = new char[128];
        strcpy(imfilename,filename);

        char buf[128];
        sprintf(buf, "%d", pos);

        strcat(imfilename, "_");
        strcat(imfilename, buf);
        //strcat(imfilename, ".png");


        frame.save(imfilename);

        outfile <<  imfilename << std::endl;

    }

}





}
