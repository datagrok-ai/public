#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "dft.h"
#include <emscripten.h> // note we added the emscripten header

void EMSCRIPTEN_KEEPALIVE sma(float *buffer, int bufSize, int window_size)
{
    float sum = 0;
    float res[bufSize];
    for (int i = 0; i < window_size - 1; i++)
    {
        sum = sum + buffer[i];
        res[i] = 0;
    }
    for (int i = window_size - 1; i < bufSize; i++)
    {
        sum = sum + buffer[i];
        res[i] = sum / window_size;
        sum = sum - buffer[i - window_size + 1];
    }
    for (int i = 0; i < bufSize; i++)
    {
        buffer[i] = res[i];
    }
}

void EMSCRIPTEN_KEEPALIVE exps(float *buffer, int bufSize, float alpha)
{
    float res[bufSize];
    res[0] = buffer[0];
    for (int i = 1; i < bufSize; i++)
    {
        res[i] = alpha * buffer[i] + (1 - alpha) * res[i - 1];
    }
    for (int i = 0; i < bufSize; i++)
    {
        buffer[i] = res[i];
    }
}

void EMSCRIPTEN_KEEPALIVE kalman(float *buffer, int bufSize, float q, float r, float p)
{
    float res[bufSize];
    res[0] = buffer[0];
    float estimate = buffer[0];
    float p_prev = 0;
    float k = p / (p + r);
    for (int i = 1; i < bufSize; i++)
    {
        estimate = res[i - 1];
        p = p_prev + q;
        k = p / (p + r);
        res[i] = estimate + k * (buffer[i] - estimate);
        p_prev = (1 - k) * p;
    }
    for (int i = 0; i < bufSize; i++)
    {
        buffer[i] = res[i];
    }
}

void EMSCRIPTEN_KEEPALIVE minmax(float *buffer, int bufSize)
{
    float min = buffer[0];
    float max = buffer[0];
    for (int i = 1; i < bufSize; i++)
    {
        if (buffer[i] > max)
        {
            max = buffer[i];
        }
        if (buffer[i] < min)
        {
            min = buffer[i];
        }
    }
    for (int i = 0; i < bufSize; i++)
    {
        buffer[i] = (buffer[i] - min) / (max - min);
    }
}

void EMSCRIPTEN_KEEPALIVE zscore(float *buffer, int bufSize)
{
    float mean = 0;
    float deviation = 0;
    for (int i = 0; i < bufSize; i++)
    {
        mean = mean + buffer[i];
    }
    mean = mean / bufSize;
    for (int i = 0; i < bufSize; i++)
    {
        deviation = deviation + pow((buffer[i] - mean), 2);
    }
    deviation = sqrt(deviation / bufSize);
    for (int i = 0; i < bufSize; i++)
    {
        buffer[i] = (buffer[i] - mean) / deviation;
    }
}

void EMSCRIPTEN_KEEPALIVE boxcox(float *buffer, int bufSize, float lambda, float ofset)
{
    if (lambda == 0)
    {
        for (int i = 0; i < bufSize; i++)
        {
            buffer[i] = log(buffer[i] + ofset);
        }
    }
    else
    {
        for (int i = 0; i < bufSize; i++)
        {
            buffer[i] = (pow((buffer[i] + ofset), lambda) - 1) / lambda;
        }
    }
}

void EMSCRIPTEN_KEEPALIVE ffilterold(float *buffer, int bufSize, float lowcuthz, float hicuthz)
{
    //void EMSCRIPTEN_KEEPALIVE ffilter(float *buffer, int bufSize, float lowcut, float hicut){
    float imbuffer[bufSize];
    float rout[bufSize];
    float imout[bufSize];
    for (int i = 0; i < bufSize; i++)
    {
        imbuffer[i] = 0;
    }
    float *imin = imbuffer;
    float *routp = rout;
    float *imoutp = imout;
    dft(buffer, imbuffer, routp, imoutp, bufSize);
    int low_frq_number = (int)lowcuthz;
    int hi_frq_number = (int)hicuthz;
    //int low_frq_number = (int) (bufSize*(1-lowcut))/2;
    //int hi_frq_number = (int) (bufSize*hicut)/2;
    int center = (int)bufSize / 2.0;
    for (int i = low_frq_number; i < center; i++)
    {
        rout[i] = 0;
        imout[i] = 0;
        rout[bufSize - i - 1] = 0;
        imout[bufSize - i - 1] = 0;
    }
    for (int i = 0; i < hi_frq_number; i++)
    {
        rout[i] = 0;
        imout[i] = 0;
        rout[bufSize - i - 1] = 0;
        imout[bufSize - i - 1] = 0;
    }
    idft(routp, imoutp, buffer, imbuffer, bufSize);
}

void EMSCRIPTEN_KEEPALIVE sdensityold(float *buffer, int bufSize)
{
    float imbuffer[bufSize];
    float rout[bufSize];
    float imout[bufSize];
    for (int i = 0; i < bufSize; i++)
    {
        imbuffer[i] = 0;
    }
    float *imin = imbuffer;
    float *routp = rout;
    float *imoutp = imout;
    dft(buffer, imbuffer, routp, imoutp, bufSize);
    for (int i = 0; i < bufSize; i++)
    {
        buffer[i] = pow(rout[i], 2) + pow(imout[i], 2);
    }
}

void EMSCRIPTEN_KEEPALIVE samplitude(float *buffer, int bufSize)
{
    float imbuffer[bufSize];
    float rout[bufSize];
    float imout[bufSize];
    for (int i = 0; i < bufSize; i++)
    {
        imbuffer[i] = 0;
    }
    float *imin = imbuffer;
    float *routp = rout;
    float *imoutp = imout;
    dft(buffer, imbuffer, routp, imoutp, bufSize);
    for (int i = 0; i < bufSize; i++)
    {
        buffer[i] = rout[i];
    }
}

void EMSCRIPTEN_KEEPALIVE gettrend(float *buffer, int bufSize)
{
    float sumXY = 0;
    float sumXs = 0;
    float sumX = bufSize;
    float sumY = 0;
    for (int i = 0; i < bufSize; i++)
    {
        sumXY = sumXY + i * buffer[i];
        sumXs = sumXs + pow(i, 2);
        sumY = sumY + buffer[i];
    }
    float beta = (bufSize * sumXY - sumX * sumY) / (bufSize * sumXs - pow(sumX, 2));
    float alpha = (sumY - beta * sumX) / bufSize;
    for (int i = 0; i < bufSize; i++)
    {
        buffer[i] = alpha + beta * i;
    }
}

void EMSCRIPTEN_KEEPALIVE removetrend(float *buffer, int bufSize)
{
    float sumXY = 0;
    float sumXs = 0;
    float sumX = bufSize;
    float sumY = 0;
    for (int i = 0; i < bufSize; i++)
    {
        sumXY = sumXY + i * buffer[i];
        sumXs = sumXs + pow(i, 2);
        sumY = sumY + buffer[i];
    }
    float beta = (bufSize * sumXY - sumX * sumY) / (bufSize * sumXs - pow(sumX, 2));
    float alpha = (sumY - beta * sumX) / bufSize;
    for (int i = 0; i < bufSize; i++)
    {
        buffer[i] = buffer[i] - (alpha + beta * i);
    }
}

void EMSCRIPTEN_KEEPALIVE ffilter(float *buffer, int bufSize, float lowcuthz, float hicuthz)
{
    int powerOf2 = pow(2, ceil(log(bufSize) / log(2)));
    float reInput[powerOf2];
    float imInput[powerOf2];
    float reOutput[powerOf2];
    float imOutput[powerOf2];
    for (int i = 0; i < bufSize; i++)
    {
        imInput[i] = 0;
        reInput[i] = buffer[i];
    }
    for (int i = bufSize; i < powerOf2; i++)
    {
        imInput[i] = 0;
        reInput[i] = 0;
    }
    float *rInPtr = reInput;
    float *imInPtr = imInput;
    float *rOutPtr = reOutput;
    float *imOutPtr = imOutput;
    fft(rInPtr, imInPtr, rOutPtr, imOutPtr, powerOf2, false);
    int low_frq_number = (int)lowcuthz;
    int hi_frq_number = (int)hicuthz;
    int center = (int)powerOf2 / 2.0;
    for (int i = low_frq_number; i < center; i++)
    {
        reOutput[i] = 0;
        imOutput[i] = 0;
        reOutput[bufSize - i - 1] = 0;
        imOutput[bufSize - i - 1] = 0;
    }
    for (int i = 0; i < hi_frq_number; i++)
    {
        reOutput[i] = 0;
        imOutput[i] = 0;
        reOutput[bufSize - i - 1] = 0;
        imOutput[bufSize - i - 1] = 0;
    }
    fft(rOutPtr, imOutPtr, rInPtr, imInPtr, powerOf2, true);
    for (int i = 0; i < bufSize; i++)
    {
        buffer[i] = reInput[i];
    }
}

void EMSCRIPTEN_KEEPALIVE sdensity(float *buffer, int bufSize)
//HIGHER FREQUENCIES MAY BE CUT OUT
{
    int powerOf2 = pow(2, ceil(log(bufSize) / log(2)));
    float reInput[powerOf2];
    float imInput[powerOf2];
    float reOutput[powerOf2];
    float imOutput[powerOf2];
    for (int i = 0; i < bufSize; i++)
    {
        imInput[i] = 0;
        reInput[i] = buffer[i];
    }
    for (int i = bufSize; i < powerOf2; i++)
    {
        imInput[i] = 0;
        reInput[i] = 0;
    }
    float *rInPtr = reInput;
    float *imInPtr = imInput;
    float *rOutPtr = reOutput;
    float *imOutPtr = imOutput;
    fft(rInPtr, imInPtr, rOutPtr, imOutPtr, powerOf2, false);
    //dft(buffer, imbuffer, routp, imoutp, bufSize);
    for (int i = 0; i < bufSize; i++)
    {
        buffer[i] = pow(reOutput[i], 2) + pow(imOutput[i], 2);
    }
}