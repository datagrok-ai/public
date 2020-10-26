#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <time.h>
clock_t start, end;
#ifndef PI
#define PI 3.14159265358979323846264338327950288
#endif

//Old Discrete Fourier Transform
void dft(float *inreal, float *inimag, float *outreal, float *outimag, int length)
{
    for (int k = 0; k < length; k++)
    { // For each output element
        float sumreal = 0;
        float sumimag = 0;
        for (int t = 0; t < length; t++)
        {                                                     // For each input element
            float angle = 2 * M_PI * t * k / ((float)length); // 2*PI*t*k / length;
            sumreal += inreal[t] * cos(angle) + inimag[t] * sin(angle);
            sumimag += -inreal[t] * sin(angle) + inimag[t] * cos(angle);
        }
        outreal[k] = sumreal;
        outimag[k] = sumimag;
    }
}

//Old Inverse Discrete Fourier Transform
void idft(float *inreal, float *inimag, float *outreal, float *outimag, int length)
{
    for (int k = 0; k < length; k++)
    { // For each output element
        float sumreal = 0;
        float sumimag = 0;
        for (int t = 0; t < length; t++)
        {                                                     // For each input element
            float angle = 2 * M_PI * t * k / ((float)length); // 2*PI*t*k / length;
            sumreal += inreal[t] * cos(angle) - inimag[t] * sin(angle);
            sumimag += inreal[t] * sin(angle) + inimag[t] * cos(angle);
        }
        outreal[k] = sumreal / ((float)length);
        outimag[k] = sumimag / ((float)length);
    }
}

// New Fast Discrete Fourier Transform
//len == power of 2 is expected
void fdft(float *ar, float *ai, bool invert, int len)
{
    int n = len;
    if (n != 1)
    {
        float a0r[n / 2], a0i[n / 2], a1r[n / 2], a1i[n / 2];
        for (int i = 0; 2 * i < n; i++)
        {
            a0r[i] = ar[2 * i];
            a0i[i] = ai[2 * i];
            a1r[i] = ar[2 * i + 1];
            a1i[i] = ai[2 * i + 1];
        }
        float *at1 = a0r;
        float *at2 = a0i;
        float *at3 = a1r;
        float *at4 = a1i;

        fdft(at1, at2, invert, n / 2);
        fdft(at3, at4, invert, n / 2);

        float angle = 2 * PI / n * (invert ? -1 : 1);
        float wr = 1;
        float wi = 0;
        float wnr = cos(angle);
        float wni = sin(angle);
        for (int i = 0; 2 * i < n; i++)
        {
            ar[i] = a0r[i] + (wr * a1r[i] - wi * a1i[i]);
            ai[i] = a0i[i] + (wi * a1r[i] + wr * a1i[i]);
            // a[i] = a0[i] + w * a1[i];
            ar[i + n / 2] = a0r[i] - (wr * a1r[i] - wi * a1i[i]);
            ai[i + n / 2] = a0i[i] - (wi * a1r[i] + wr * a1i[i]);
            //a[i + n/2] = a0[i] - w * a1[i];
            if (invert)
            {
                ar[i] /= 2;
                ai[i] /= 2;
                ar[i + n / 2] /= 2;
                ai[i + n / 2] /= 2;
            }
            float tmp = wr;
            wr = wr * wnr - wi * wni;
            wi = tmp * wni + wi * wnr;
            //w *= wn;
        }
    }
}

//Fast Fourier Transform
void fft(float *inreal, float *inimag, float *outreal, float *outimag, int length, bool invert)
{
    //Length Equal to power of 2 is expected
    // use {invert = true} for inverse transform
    float aRe[length];
    float aIm[length];
    for (int i = 0; i < length; i++)
    {
        aRe[i] = inreal[i];
        aIm[i] = inimag[i];
    }
    float *ar = aRe;
    float *ai = aIm;
    fdft(ar, ai, invert, length);

    for (int i = 0; i < length; i++)
    {
        outreal[i] = ar[i];
        outimag[i] = ai[i];
    }
}