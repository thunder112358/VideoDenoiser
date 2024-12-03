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

#include "libBasic.h"


namespace libUSTG
{


std::string int2string(int number)
{
    std::stringstream s;
    s << number;
    return s.str();
}




//
//! Value operations
//

void fpClear(float *fpI,float fValue, int iLength)
{
    assert(iLength > 0);
    for (int ii=0; ii < iLength; ii++) fpI[ii] = fValue;
}

void fpCopy(float *fpI,float *fpO, int iLength)
{
    assert(iLength > 0);
    if (fpI != fpO)  memcpy((void *) fpO, (const void *) fpI, iLength * sizeof(float));
}


float  fpMean(float *u,int size)
{
    assert(size > 0);
    float *ptru=&u[0];
    float mean=0.0;
    for(int i=0; i<size; i++,ptru++)  mean +=  *ptru;
    mean/=(float) size;
    return mean;
}



float  fpVar(float *u,int size)
{

    assert(size > 0);
    float *ptru=&u[0];
    float mean=0.0;
    float mean2 = 0.0;
    for(int i=0; i<size; i++,ptru++) {
        mean +=  *ptru;
        mean2 +=  *ptru *  (*ptru);
    }

    mean/=(float) size;
    mean2/=(float) size;
    float var = mean2- mean*mean;

    var = fabsf(var);
    return var;
}








float fpDistLp(float *u, float *v, int p, int size)
{

    float fDist = 0.0f;

    for (int ii=0; ii < size; ii++) {

        float dif = fabsf(u[ii] - v[ii]);


        if (p == 0) fDist = MAX(fDist, dif);
        else if (p == 1)  fDist += dif;
        else if (p == 2)
            fDist += dif*dif;
        else {
            fDist += powf(dif, (float) p);

        }

    }

    fDist /= (float) size;

    if (p>0)
        fDist = powf(fDist, 1.0 / (float) p);


    return fDist;

}



float fpDistLp(float *u, float *v, float *m, int p, int size)
{

    float fDist = 0.0f;
    int iCount = 0;

    for (int ii=0; ii < size; ii++)
        if (m[ii] > 0.0f) {

            float dif = fabsf(u[ii] - v[ii]);


            if (p == 0) fDist = MAX(fDist, dif);
            else if (p == 1)  fDist += dif;
            else if (p == 2)
                fDist += dif*dif;
            else {
                fDist += powf(dif, (float) p);

            }

            iCount++;
        }

    fDist /= (float) iCount;

    if (p>0)
        fDist = powf(fDist, 1.0 / (float) p);


    return fDist;

}



//
//! Float pointer ordering
//



int order_float_increasing(const void *a, const void *b)
{
    if ( *(float*)a  > *(float*)b ) return 1;
    else if ( *(float*)a  < *(float*)b ) return -1;

    return 0;
}




int order_float_decreasing(const void *a, const void *b)
{
    if ( *(float*)a  > *(float*)b ) return -1;
    else if ( *(float*)a  < *(float*)b ) return 1;

    return 0;
}





void fpQuickSort(float *fpI, int iLength, int inverse)
{

    if (inverse)
        qsort(fpI, iLength, sizeof(float), order_float_decreasing);
    else
        qsort(fpI, iLength, sizeof(float), order_float_increasing);

}






////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//! Patch Statistics
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////


void fiPatchStatistics(float *fpIn, float *fpMinV, float *fpMaxV, float *fpMeanV, float *fpVarV, float *fpMedianV, float fRadius, int iWidth, int iHeight)
{

    //! Parameters
    int iRadius = (int)(fRadius+1.0);
    int iNeigSize = (2*iRadius+1)*(2*iRadius+1);
    float fRadiusSqr = fRadius * fRadius;


    float * vector = new float[iNeigSize];



    //! For each pixel
    for(int x=0; x < iWidth; x++)
        for(int y=0; y< iHeight; y++) {

            int iCount=0;
            float fMin = fLarge;
            float fMax = -fLarge;


            for(int i=-iRadius; i<=iRadius; i++)
                for(int j=-iRadius; j<=iRadius; j++)
                    if ((float) (i*i + j*j) <= fRadiusSqr) {

                        int x0=x+i;
                        int y0=y+j;

                        if (x0 >= 0 && y0 >= 0 && x0 < iWidth && y0 < iHeight) {
                            float fValue= fpIn[y0*iWidth+x0];

                            if (fValue > fMax) fMax = fValue;
                            if (fValue < fMin) fMin = fValue;

                            vector[iCount] = fValue;
                            iCount++;

                        }

                    }

            int l = y*iWidth+x;
            if (fpMinV) fpMinV[l] = fMin;
            if (fpMaxV) fpMaxV[l] = fMax;

            if (fpMeanV) fpMeanV[l] = fpMean(vector, iCount);
            if (fpVarV)  fpVarV[l] = fpVar(vector, iCount);


            if (fpMedianV) {
                fpQuickSort(vector, iCount, 0);
                fpMedianV[l] = vector[iCount / 2];
            }



        }

    delete[] vector;
}



void fiPatchMedian(float *fpIn, float *fpMedianV, float fRadius, int iWidth, int iHeight)
{
    fiPatchStatistics( fpIn, NULL, NULL,  NULL,  NULL, fpMedianV,  fRadius,  iWidth, iHeight);
}




//
//! Noise
//



void fpAddNoiseGaussian(float *u, float *v, float std, long int randinit, int size)
{

    srand48( (long int) time (NULL) + (long int) getpid()  + (long int) randinit);
    //srand48( (long int) time (NULL)  + (long int) randinit);

    for (int i=0; i< size; i++) {

        float a = drand48();
        float b = drand48();
        float z = (float)(std)*sqrt(-2.0*log(a))*cos(2.0*M_PI*b);

        v[i] =  u[i] + (float) z;

    }

}




//
//! Patch distances
//


float fiL2FloatDist(float *u0,float *u1,int i0,int j0,int i1,int j1,int xradius, int yradius, int width0, int width1)
{



    float dist=0.0;
    for (int s=-yradius; s<= yradius; s++) {

        int l = (j0+s)*width0 + (i0-xradius);
        float *ptr0 = &u0[l];

        l = (j1+s)*width1 + (i1-xradius);
        float *ptr1 = &u1[l];

        for(int r=-xradius; r<=xradius; r++,ptr0++,ptr1++) {
            float dif = (*ptr0 - *ptr1);
            dist += (dif*dif);
        }

    }

    return dist;
}











////////////////////////////////////////////////////////////////////////////////
//! BEGIN Numerics
////////////////////////////////////////////////////////////////////////////////





laVector::laVector() : d_n(0), d_v(0) {}



laVector::laVector(int n) : d_n(n), d_v(new float[n]) {}


laVector::laVector(const float& a, int n) : d_n(n), d_v(new float[n])
{
    for(int i=0; i<n; i++)
        d_v[i] = a;
}


laVector::laVector(const float *a, int n) : d_n(n), d_v(new float[n])
{
    for(int i=0; i<n; i++)
        d_v[i] = *a++;
}


laVector::laVector(const laVector &rhs) : d_n(rhs.d_n), d_v(new float[d_n])
{
    for(int i=0; i<d_n; i++)
        d_v[i] = rhs[i];
}


laVector & laVector::operator=(const laVector &rhs)
{
    if (this != &rhs) {
        if (d_n != rhs.d_n) {
            if (d_v != 0) delete [] d_v;
            d_n =rhs.d_n;
            d_v = new float[d_n];
        }

        for (int i=0; i<d_n; i++)
            d_v[i]=rhs[i];
    }
    return *this;
}


laVector & laVector::operator=(const float &a)	//assign a to every element
{
    for (int i=0; i<d_n; i++)
        d_v[i]=a;
    return *this;
}


float & laVector::operator[](const int i)	//subscripting
{
    return d_v[i];
}


const float & laVector::operator[](const int i) const	//subscripting
{
    return d_v[i];
}


float * laVector::v()
{
    return d_v;
}



int laVector::size() const
{
    return d_n;
}


laVector::~laVector()
{
    if (d_v != 0)
        delete[] d_v;
}


void laVector::erase()
{
    d_n = 0;
    if (d_v) delete[] d_v;
    d_v=0;
}


void laVector::create(int n)
{
    erase();
    d_n = n;
    d_v = new float[d_n];
    for (int ii=0; ii < d_n; ii++) d_v[ii] = 0.0f;
}



laVector laVector::copyBlock(int i0, int length)
{

    laVector block(length);
    for (int ii=0; ii < length; ii++)  block.d_v[ii] = d_v[i0+ii];

    return block;
}


void laVector::sort(int decreasing_order)
{

    if (decreasing_order)
        qsort(d_v, d_n, sizeof(float), order_float_decreasing);
    else
        qsort(d_v, d_n, sizeof(float), order_float_increasing);

}




laMatrix::laMatrix() : d_n(0), d_m(0), d_v(0) {}


laMatrix::laMatrix(int n, int m) : d_n(n), d_m(m), d_v(new float*[n])
{
    d_v[0] = new float[m*n];
    for (int i=1; i< n; i++)
        d_v[i] = d_v[i-1] + m;
}



void laMatrix::create(int n, int m)
{

    if (d_v != 0) {
        delete[] d_v[0];
        delete[] d_v;
    }

    d_n = n;
    d_m = m;

    d_v = new float*[n];
    d_v[0] = new float[m*n];
    for (int i=1; i< n; i++)
        d_v[i] = d_v[i-1] + m;
}





laMatrix::laMatrix(const float &a, int n, int m) : d_n(n), d_m(m), d_v(new float*[n])
{
    int i,j;
    d_v[0] = new float[m*n];
    for (i=1; i< n; i++)
        d_v[i] = d_v[i-1] + m;
    for (i=0; i< n; i++)
        for (j=0; j<m; j++)
            d_v[i][j] = a;
}


laMatrix::laMatrix(const float *a, int n, int m) : d_n(n), d_m(m), d_v(new float*[n])
{
    int i,j;
    d_v[0] = new float[m*n];
    for (i=1; i< n; i++)
        d_v[i] = d_v[i-1] + m;
    for (i=0; i< n; i++)
        for (j=0; j<m; j++)
            d_v[i][j] = *a++;
}


laMatrix::laMatrix(const laMatrix &rhs) : d_n(rhs.d_n), d_m(rhs.d_m), d_v(new float*[rhs.d_n])
{
    assert(d_m>0 && d_n>0);

    d_v[0] = new float[d_m * d_n];

    for (int i=1; i< d_n; i++)
        d_v[i] = d_v[i-1] + d_m;

    for (int i=0; i< d_n; i++)
        for (int j=0; j<d_m; j++)
            d_v[i][j] = rhs[i][j];
}


laMatrix & laMatrix::operator=(const laMatrix &rhs)
{
    if (this != &rhs) {
        int i,j;
        if (d_n != rhs.d_n || d_m != rhs.d_m) {
            if (d_v != 0) {
                delete[] (d_v[0]);
                delete[] (d_v);
            }

            d_n=rhs.d_n;
            d_m=rhs.d_m;
            d_v = new float*[d_n];
            d_v[0] = new float[d_m*d_n];


        }


        for (i=1; i< d_n; i++)
            d_v[i] = d_v[i-1] + d_m;
        for (i=0; i< d_n; i++)
            for (j=0; j<d_m; j++)
                d_v[i][j] = rhs[i][j];
    }
    return *this;
}


laMatrix & laMatrix::operator=(const float &a)	//assign a to every element
{
    for (int i=0; i< d_n; i++)
        for (int j=0; j<d_m; j++)
            d_v[i][j] = a;
    return *this;
}



float* laMatrix::operator[](const int i)	//subscripting: pointer to row i
{
    return d_v[i];
}


const float* laMatrix::operator[](const int i) const
{
    return d_v[i];
}


float ** laMatrix::v()
{
    return d_v;
}





int laMatrix::nrows() const
{
    return d_n;
}


int laMatrix::ncols() const
{
    return d_m;
}


laMatrix::~laMatrix()
{
    if (d_v != 0) {
        delete[] (d_v[0]);
        delete[] (d_v);
    }
}







laMatrix operator*(float a, const laMatrix& rhs)
{

    laMatrix r(rhs.d_n, rhs.d_m);

    for (int k=r.d_n * r.d_m - 1; k>=0; k--)
        r.d_v[0][k] = a * rhs.d_v[0][k];

    return r;
}



laMatrix operator/(const laMatrix& lhs, float a)
{
    laMatrix r(lhs.d_n, lhs.d_m);

    for (int k=r.d_n * r.d_m - 1; k>=0; k--)
        r.d_v[0][k] = lhs.d_v[0][k] / a;
    return r;

}



laMatrix operator+(const laMatrix& lhs, const laMatrix& rhs)
{
    laMatrix r(lhs.d_n, lhs.d_m);
    for (int k=r.d_m * r.d_n -1; k>=0; k--)
        r.d_v[0][k] = lhs.d_v[0][k] + rhs.d_v[0][k];
    return r;
}



laMatrix operator-(const laMatrix& lhs, const laMatrix& rhs)
{
    laMatrix r(lhs.d_m, lhs.d_n);
    for (int k=r.d_m * r.d_n - 1; k>=0; k--)
        r.d_v[0][k] = lhs.d_v[0][k] - rhs.d_v[0][k];
    return r;
}



laMatrix operator*(const laMatrix& lhs, const laMatrix& rhs)
{
    float aux;

    laMatrix r(lhs.d_n, rhs.d_m);
    for (int i=0; i< r.d_n; i++)
        for (int j=0; j< r.d_m; j++) {
            aux = 0.0;
            for (int k=0; k< lhs.d_m; k++)
                aux += lhs.d_v[i][k] * rhs.d_v[k][j];

            r.d_v[i][j] = aux;
        }

    return r;
}



laVector operator*(const laMatrix& lhs, const laVector& rhs)
{


    laVector r(lhs.d_n);
    for (int i=0; i < r.size(); i++) {

        r[i] = 0;
        for (int k=0; k< rhs.size(); k++) {
            r[i] += lhs.d_v[i][k] * rhs[k];
        }
    }

    return r;
}



laMatrix laMatrix::transposed()
{
    laMatrix r(d_m, d_n);

    for (int ii=0; ii < d_m; ii++)
        for (int jj=0; jj < d_n; jj++) {
            r.d_v[ii][jj] = d_v[jj][ii];
        }

    return r;
}






laMatrix laMatrix::copyBlock(int i0, int j0, int rowb, int colb)
{

    laMatrix block(rowb, colb);

    for (int i=0; i < rowb; i++)
        for (int j=0; j < colb; j++)
            block.d_v[i][j] = d_v[i0+i][j0+j];

    return block;
}


float withSignOf(float a, float b)
{
    return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}


float svdhypot(float a, float b)
{
    a = fabsf(a);
    b = fabsf(b);
    if(a > b) {
        b /= a;
        return a*sqrt(1.0 + b*b);
    } else if(b) {
        a /= b;
        return b*sqrt(1.0 + a*a);
    }
    return 0.0;
}



void svdrotate_double(double& a, double& b, double c, double s)
{
    double d = a;
    a = +d*c +b*s;
    b = -d*s +b*c;
}





double **laMatrix2double(laMatrix &Ain)
{
    int rows = Ain.nrows();
    int cols = Ain.ncols();

    double **A=new double *[rows];
    A[0] = new double[rows*cols];
    for (int i=1; i< rows; i++) A[i]=A[i-1]+cols;

    for (int i=0; i< rows; i++)
        for (int j=0; j< cols; j++) A[i][j]=(double) Ain[i][j];

    return A;
}

double *laVector2double(laVector &Vin)
{
    int size = Vin.size();

    double *V=new double[size];
    for (int i=0; i< size; i++) V[i]=(double) Vin[i];

    return V;
}

//includes delete[]
void double2laMatrix(double **A, laMatrix &Ain)
{
    int rows = Ain.nrows();
    int cols = Ain.ncols();

    for (int i=0; i< rows; i++)
        for (int j=0; j< cols; j++) Ain[i][j]=(float) A[i][j];

    delete[] A[0];
    delete[] A;
}

//includes delete[]
void double2laVector(double *V, laVector &Vin)
{
    int size = Vin.size();

    for (int i=0; i< size; i++) Vin[i]=(float) V[i];

    delete[] V;
}

//internal computations in double precision
//input and output: float
void compute_svd_double(laMatrix &Ain, laMatrix &m_Uin, laMatrix &m_Vin, laVector &m_Win)
{

    int rows = Ain.nrows();
    int cols = Ain.ncols();

    double **A=laMatrix2double(Ain);
    double **m_U=laMatrix2double(m_Uin);
    double **m_V=laMatrix2double(m_Vin);
    double *m_W=laVector2double(m_Win);

    const double	EPSILON = 0.00001;
    const int SVD_MAX_ITS = 100;

    double g, scale, anorm;
    double * RV1 = new double[cols];


    for(int i=0; i < rows; i++)
        for(int j=0; j < cols; j++)
            m_U[i][j] = A[i][j];

    // Householder reduction to bidiagonal form:
    anorm = g = scale = 0.0;
    for (int i=0; i< cols; i++) {
        int l = i + 1;
        RV1[i] = scale*g;
        g = scale = 0.0;
        if(i < rows) {
            for (int k=i; k< rows; k++)
                scale += fabsf(m_U[k][i]);
            if (scale != 0.0) {
                double invScale=1.0/scale, s=0.0;
                for (int k=i; k< rows; k++) {
                    m_U[k][i] *= invScale;
                    s += m_U[k][i] * m_U[k][i];
                }
                double f = m_U[i][i];
                g = - withSignOf(sqrt(s),f);
                double h = 1.0 / (f*g - s);
                m_U[i][i] = f - g;
                for (int j=l; j< cols; j++) {
                    s = 0.0;
                    for (int k=i; k< rows; k++)
                        s += m_U[k][i] * m_U[k][j];
                    f = s * h;
                    for (int k=i; k< rows; k++)
                        m_U[k][j] += f * m_U[k][i];
                }
                for (int k=i; k< rows; k++)
                    m_U[k][i] *= scale;
            }
        }


        m_W[i] = scale * g;
        g = scale = 0.0;
        if ( i< rows && i< cols-1 ) {
            for (int k=l; k< cols; k++)
                scale += fabsf(m_U[i][k]);
            if (scale != 0.0) {
                double invScale=1.0/scale, s=0.0;
                for (int k=l; k< cols; k++) {
                    m_U[i][k] *= invScale;
                    s += m_U[i][k] * m_U[i][k];
                }
                double f = m_U[i][l];
                g = - withSignOf(sqrt(s),f);
                double h = 1.0 / (f*g - s);
                m_U[i][l] = f - g;
                for (int k=l; k< cols; k++)
                    RV1[k] = m_U[i][k] * h;
                for (int j=l; j< rows; j++) {
                    s = 0.0;
                    for (int k=l; k< cols; k++)
                        s += m_U[j][k] * m_U[i][k];
                    for (int k=l; k< cols; k++)
                        m_U[j][k] += s * RV1[k];
                }
                for (int k=l; k< cols; k++)
                    m_U[i][k] *= scale;
            }
        }
        anorm = MAX(anorm, fabsf(m_W[i]) + fabsf(RV1[i]) );
    }

    // Accumulation of right-hand transformations:
    m_V[cols-1][cols-1] = 1.0;
    for (int i= cols-2; i>=0; i--) {
        m_V[i][i] = 1.0;
        int l = i+1;
        g = RV1[l];
        if (g != 0.0) {
            double invgUil = 1.0 / (m_U[i][l]*g);
            for (int j=l; j< cols; j++)
                m_V[j][i] = m_U[i][j] * invgUil;
            for (int j=l; j< cols; j++) {
                double s = 0.0;
                for (int k=l; k< cols; k++)
                    s += m_U[i][k] * m_V[k][j];
                for (int k=l; k< cols; k++)
                    m_V[k][j] += s * m_V[k][i];
            }
        }
        for (int j=l; j< cols; j++)
            m_V[i][j] = m_V[j][i] = 0.0;
    }

    // Accumulation of left-hand transformations:
    for (int i=MIN(rows,cols)-1; i>=0; i--) {
        int l = i+1;
        g = m_W[i];
        for (int j=l; j< cols; j++)
            m_U[i][j] = 0.0;
        if (g != 0.0) {
            g = 1.0 / g;
            double invUii = 1.0 / m_U[i][i];
            for (int j=l; j< cols; j++) {
                double s = 0.0;
                for (int k=l; k< rows; k++)
                    s += m_U[k][i] * m_U[k][j];
                double f = (s * invUii) * g;
                for (int k=i; k< rows; k++)
                    m_U[k][j] += f * m_U[k][i];
            }
            for (int j=i; j< rows; j++)
                m_U[j][i] *= g;
        } else
            for (int j=i; j< rows; j++)
                m_U[j][i] = 0.0;
        m_U[i][i] = m_U[i][i] + 1.0;
    }

    // Diagonalization of the bidiagonal form:
    for (int k=cols-1; k>=0; k--) { // Loop over singular values
        for (int its=1; its<=SVD_MAX_ITS; its++) {
            bool flag = false;
            int l  = k;
            int nm = k-1;
            while(l>0 && fabsf(RV1[l]) > EPSILON*anorm) { // Test for splitting
                if(fabsf(m_W[nm]) <= EPSILON*anorm) {
                    flag = true;
                    break;
                }
                l--;
                nm--;
            }
            if (flag) {	// Cancellation of RV1[l], if l > 0
                double c=0.0, s=1.0;
                for (int i=l; i< k+1; i++) {
                    double f = s * RV1[i];
                    RV1[i] = c * RV1[i];
                    if (fabsf(f)<=EPSILON*anorm)
                        break;
                    g = m_W[i];
                    double h = svdhypot(f,g);
                    m_W[i] = h;
                    h = 1.0 / h;
                    c = g * h;
                    s = - f * h;
                    for (int j=0; j< rows; j++)
                        svdrotate_double(m_U[j][nm],m_U[j][i], c,s);
                }
            }
            double z = m_W[k];
            if (l==k) {		// Convergence of the singular value
                if (z< 0.0) {	// Singular value is made nonnegative
                    m_W[k] = -z;
                    for (int j=0; j< cols; j++)
                        m_V[j][k] = - m_V[j][k];
                }
                break;
            }

            // Exception if convergence to the singular value not reached:
            if(its==SVD_MAX_ITS) {
                printf("svd::convergence_error\n");
                delete[] RV1;
                exit(-1);
            }
            double x = m_W[l]; // Get QR shift value from bottom 2x2 minor
            nm = k-1;
            double y = m_W[nm];
            g = RV1[nm];
            double h = RV1[k];
            double f = ( (y-z)*(y+z) + (g-h)*(g+h) ) / ( 2.0*h*y );
            g = svdhypot(f,1.0);
            f = ( (x-z)*(x+z) + h*(y/(f+withSignOf(g,f)) - h) ) / x;
            // Next QR transformation (through Givens reflections)
            double c=1.0, s=1.0;
            for (int j=l; j<=nm; j++) {
                int i = j+1;
                g = RV1[i];
                y = m_W[i];
                h = s * g;
                g = c * g;
                z = svdhypot(f,h);
                RV1[j] = z;
                z = 1.0 / z;
                c = f * z;
                s = h * z;
                f = x*c + g*s;
                g = g*c - x*s;
                h = y * s;
                y *= c;
                for(int jj=0; jj < cols; jj++)
                    svdrotate_double(m_V[jj][j],m_V[jj][i], c,s);
                z = svdhypot(f,h);
                m_W[j] = z;
                if (z!=0.0) { // Rotation can be arbitrary if z = 0.0
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = c*g + s*y;
                x = c*y - s*g;
                for(int jj=0; jj < rows; jj++)
                    svdrotate_double(m_U[jj][j],m_U[jj][i], c,s);
            }
            RV1[l] = 0.0;
            RV1[k] = f;
            m_W[k] = x;
        }
    }


    double2laMatrix(A, Ain); //include delete double array
    double2laMatrix(m_U, m_Uin); //include delete double array
    double2laMatrix(m_V, m_Vin); //include delete double array
    double2laVector(m_W, m_Win); //include delete double array
    delete[] RV1;

}



void compute_pca_svd(laMatrix &X, laVector &S, laMatrix &V, laMatrix &U)
{

    int n = X.nrows();
    int p = X.ncols();

    compute_svd_double(X,U,V,S);


    // U contain new coefficients, which must be normalized by eigenvalues
    for(int i=0; i < n; i++)
        for(int j=0; j < p; j++)
            U[i][j] *= S[j];

    // Normalize eigenvalues
    float norm = (float) (n-1);
    for(int i=0; i < p; i++)
        S[i] = S[i] * S[i] / norm;

    // If n < p, principal component should be zero from n to p-1
    // Coefficients of these principal components should be zero
    if (n < p) {
        for(int i=n-1; i < p; i++) S[i] = 0.0f;

        for(int j=0; j < n; j++)
            for(int i=n-1; i < p; i++)
                U[j][i] = 0.0f;
    }


}





}



void  printusage(char *pname,
                 char *gp,
                 std::vector<OptStruct*>  &opt,
                 std::vector<ParStruct*>  &par)
{


    int npar = par.size();

    //! print function prototype
    //! "usage:" function_name [-r r] [-g g] [-b b]  par1 par2 ... parn
    printf("\nusage: %s ", pname);
    for(int i=0; i < (int) strlen(gp); i++)
        if (gp[i] != ':') {
            printf("[-%c",gp[i]);

            if (i+1 < (int) strlen(gp) && gp[i+1] ==  ':') printf(" %c] ", gp[i]);
            else printf("] ");

        }

    for(int i=0; i < npar; i++)
        printf(" %s ", par[i]->name);

    printf("\n");


    //! print options with associated descriptions and defaulted values
    int j=0;
    for(int i=0; i < (int) strlen(gp); i++)
        if (gp[i] != ':') {
            printf("\t-%c",gp[i]);

            if (i+1 < (int) strlen(gp) && gp[i+1] ==  ':') {

                printf("  %c\t %s ", gp[i], opt[j]->comment);
                if (opt[j]->defvalue != NULL) printf("(Default: %s)",opt[j]->defvalue);

                printf("\n");

            } else printf("\t %s \n", opt[j]->comment);

            j++;
        }

    //! print mandatory parameters with associated descriptions
    for(int i=0; i < npar; i++) {
        printf("\t%s",par[i]->name);
        printf("\t %s\n", par[i]->comment);
    }



}


void print_call(int argc, char **argv)
{

    printf("\nCall:  ");
    for (int ii=0; ii < argc; ii++)  printf("%s ", argv[ii]);
    printf("\n");

}


int parsecmdline(char *pname,
                 char *function,
                 int argc, char **argv,
                 std::vector <OptStruct*> & opt,
                 std::vector <ParStruct*> & par)
{

    //! number of options and obligatory parameters
    int nopt = opt.size();
    int npar = par.size();


    //! size of auxiliar parameter needed by getopt
    //! the length is at maximum 2*nopt since each option contains one caracter
    //! or two if the option requires a values
    char *gp = new char[2*nopt+1];
    gp[0]='\0';


    //! clear options and concatenate option identifiers into gp vector
    for(int i=0; i < nopt; i++) {
        opt[i]->flag = 0;
        opt[i]->value=NULL;
        strcat(gp, opt[i]->gp);
    }

    //! clear necessary parameter values
    for(int i=0; i < npar; i++) {
        par[i]->value = NULL;
    }


    //! in this way getopt doesn't print any information and errors are
    //! treated and printed by this program.
    opterr = 0;


    //! getopt return the following option of the concatenated vector gp
    //! or -1 if all options have already been treated
    //! it also fills "optarg" with the associated value to this option as it is given in argv
    int c;
    while ((c = getopt (argc, argv, gp)) != -1) {

        //! set current option given by getopt
        //! if getopt finds an option which was not in the option of the program (vector gp)
        //! returns ?
        int j=0;
        for(unsigned int i=0; i < strlen(gp); i++)
            if (c == gp[i]) {
                opt[j]->flag = 1;
                opt[j]->value = optarg;  // getopt fills "optarg" with the associated value to this option as it is given in argv

                break;

            } else if (gp[i] != ':') j++;


        //! option found in console command was not one of the parameters of our program
        //! or should have a mandatory values which is not provided
        if (c == '?') {


            //! when getopt encounters an unknown option character or an option with a missing required argument
            //! it stores that option character in optopt variable
            unsigned int i = 0;
            for(i=0; i < strlen(gp); i++)
                if (optopt == gp[i]) {
                    printf("\n%s: %s\n", pname, function);
                    printf("\nerror: option -%c requires an argument.\n", optopt);
                    break;
                }

            if (i == strlen(gp)) {
                printf("\n%s: %s\n", pname, function);
                printf ("\nerror: unknown option `-%c'.\n", optopt);
            }

            print_call(argc,argv);
            printusage(pname, gp,  opt,  par);
            delete[] gp;
            return 0;

        }



    }


    //! Setting default values for non selected options
    for(int j=0; j < nopt; j++)
        if (opt[j]->flag == 0 && opt[j]->defvalue != NULL) opt[j]->value =  opt[j]->defvalue;


    //! Check remaining words in command after reading option
    if (argc - optind != npar) {
        //printf("\n%s: %s\n", pname, function);
        printf("%s\n", function);
        print_call(argc,argv);
        fprintf (stderr, "\nerror: incorrect number of parameters\n");
        printusage(pname, gp,  opt,par);
        delete[] gp;
        return 0;
    }

    //! Read mandatory parameter values
    int i=0;
    for (int index = optind; index < argc ; index++, i++) {
        par[i]->value = argv[index];
    }


    delete[] gp;
    return 1;


}






