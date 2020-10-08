#include <stdio.h>
#include <math.h>
#include "dft.h"
#include <emscripten.h> // note we added the emscripten header

void EMSCRIPTEN_KEEPALIVE sma(float *buffer, int bufSize, int window_size) {
        float sum=0;
        float res[bufSize];
        for(int i =0; i<window_size-1;i++){
            sum = sum+ buffer[i];
            res[i]=0;
        }
        for(int i=window_size-1;i<bufSize;i++){
            sum = sum + buffer[i];
            res[i]=sum/window_size;
            sum = sum - buffer[i-window_size+1];
        }
        for (int i = 0; i<bufSize;i++){
            buffer[i]=res[i];
        }
    }

void EMSCRIPTEN_KEEPALIVE exps(float *buffer, int bufSize, float alpha) {
        float res[bufSize];
        res[0]=buffer[0];
        for (int i =1; i<bufSize;i++){
            res[i]=alpha*buffer[i]+(1-alpha)*res[i-1];
        }
        for(int i =0;i<bufSize;i++){
            buffer[i]=res[i];
        }
    }


void EMSCRIPTEN_KEEPALIVE kalman(float *buffer, int bufSize, float q, float r, float p) {
    float res[bufSize];
    res[0] = buffer[0];
    float estimate = buffer[0];
    float p_prev = 0;
    float k = p/(p+r);
    for(int i =1; i<bufSize;i++){
        estimate = res[i-1];
        p = p_prev+q;
        k = p/(p+r);
        res[i]=estimate+k*(buffer[i]-estimate);
        p_prev = (1-k)*p;
    }
    for(int i =0;i<bufSize;i++){
        buffer[i]=res[i];
    }
}

void EMSCRIPTEN_KEEPALIVE minmax(float *buffer, int bufSize) {
//void minmax(float *buffer, int bufSize) {
        float min=buffer[0];
        float max=buffer[0];
        for(int i =1; i<bufSize;i++){
            if(buffer[i]>max){max = buffer[i];}
            if(buffer[i]<min){min = buffer[i];}
        }
        for (int i=0;i<bufSize;i++){
            buffer[i] = (buffer[i] - min)/(max-min);
        }
    }


void EMSCRIPTEN_KEEPALIVE zscore(float *buffer, int bufSize) {
  //void zscore(float *buffer, int bufSize){
        float mean=0;
        float deviation=0;
        for(int i =0; i<bufSize;i++){
            mean = mean + buffer[i];
        }
        mean = mean/bufSize;
        for(int i =0; i<bufSize;i++){
            deviation = deviation + pow((buffer[i]-mean),2);
        }
        deviation = sqrt(deviation/bufSize);
        for(int i = 0;i<bufSize;i++){
            buffer[i]=(buffer[i]-mean)/deviation;
        }
    }

void EMSCRIPTEN_KEEPALIVE boxcox(float *buffer, int bufSize, float lambda, float ofset) {
    //void boxcox(float *buffer, int bufSize, float lambda, float ofset){
    if (lambda == 0){
        for(int i = 0;i<bufSize;i++){
            buffer[i]=log(buffer[i]+ofset);
        }
    }
    else{
        for(int i = 0;i<bufSize;i++){
            buffer[i]=(pow((buffer[i]+ofset),lambda)-1)/lambda;
        }
    }
    }

void EMSCRIPTEN_KEEPALIVE ffilter(float *buffer, int bufSize, float lowcut, float hicut){
//void ffilter(float *buffer, int bufSize, float lowcut, float hicut){
            float imbuffer[bufSize];
            float rout[bufSize];
            float imout[bufSize];
            for(int i =0;i<bufSize;i++){
                imbuffer[i]=0;
            }
            float * imin = imbuffer;
            float * routp = rout;
            float * imoutp = imout;
            dft(buffer,imbuffer,routp,imoutp,bufSize);
            int low_frq_number = (int) (bufSize*(1-lowcut))/2;
            int hi_frq_number = (int) (bufSize*hicut)/2;
            int center = (int) bufSize/2.0;
            for (int i =low_frq_number; i<center;i++){
                rout[i]=0;
                imout[i]=0;
                rout[bufSize-i-1]=0;
                imout[bufSize-i-1]=0;
            }
            for (int i =0; i <hi_frq_number;i++){
                rout[i]=0;
                imout[i]=0;
                rout[bufSize-i-1]=0;
                imout[bufSize-i-1]=0;
            }
            idft(routp,imoutp,buffer,imbuffer,bufSize);
    }