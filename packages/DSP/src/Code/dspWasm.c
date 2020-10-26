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

void EMSCRIPTEN_KEEPALIVE ffilter(float *buffer, int bufSize, int lowcuthz, int hicuthz)
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
    int low_frq_number = (int)lowcuthz/2.0;
    int hi_frq_number = (int) hicuthz/2.0;
    int center = (int)powerOf2 / 2.0;
    for (int i = low_frq_number; i < center; i++)
    {
        reOutput[i] = 0;
        imOutput[i] = 0;
        reOutput[powerOf2 - i - 1] = 0;
        imOutput[powerOf2 - i - 1] = 0;
    }
    for (int i = 0; i < hi_frq_number; i++)
    {
        reOutput[i] = 0;
        imOutput[i] = 0;
        reOutput[powerOf2 - i - 1] = 0;
        imOutput[powerOf2 - i - 1] = 0;
    }
    fft(rOutPtr, imOutPtr, rInPtr, imInPtr, powerOf2, true);
    for (int i = 0; i < bufSize; i++)
    {
        buffer[i] = reInput[i];
    }
}

int EMSCRIPTEN_KEEPALIVE sdensity(float *buffer, int bufSize)
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
    return ((int)powerOf2/2);
}

int EMSCRIPTEN_KEEPALIVE asamp(float *buffer, int bufSize, int windowSize, int ofset)
{
    // Simple Moving average resampling
    // windowSize - resampling window size;
    // ofset - additional parameter, allowing downsampling to start not from zero
    int newSize = floor((bufSize - ofset) / ((float)windowSize));
    float out[newSize];
    for (int i = 0; i < newSize; i++)
    {
        float currentValue = 0;
        for (int j = 0; j < windowSize; j++)
        {
            currentValue = currentValue + buffer[windowSize * i + ofset + j];
        }
        out[i] = currentValue / windowSize;
    }
    for (int i = 0; i < newSize; i++)
    {
        buffer[i] = out[i];
    }
    return newSize;
}

void EMSCRIPTEN_KEEPALIVE subsamp(float *buffer, int bufSize, int nbufSize, int ofset)
{
    // Simple subsample resampling
    // nbufSize - desired new sample size;
    // ofset - additional parameter, allowing downsampling to start not from zero
    float out[nbufSize];
    int windowSize = floor((bufSize - ofset) / nbufSize);
    for (int i = 0; i < nbufSize; i++)
    {
        out[i] = buffer[windowSize * i + ofset];
    }
    for (int i = 0; i < nbufSize; i++)
    {
        buffer[i] = out[i];
    }
}