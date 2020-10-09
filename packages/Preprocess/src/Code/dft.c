#include <stdio.h>
#include <math.h>

void dft( float* inreal , float* inimag,float* outreal, float* outimag, int length) {
    for (int k = 0; k < length; k++) {  // For each output element
        float sumreal = 0;
        float sumimag = 0;
        for (int t = 0; t < length; t++) {  // For each input element
            float angle =2*M_PI*t*k/((float)length);// 2*PI*t*k / length;
            sumreal +=  inreal[t] * cos(angle) + inimag[t] * sin(angle);
            sumimag += -inreal[t] * sin(angle) + inimag[t] * cos(angle);
        }
        outreal[k] = sumreal;
        outimag[k] = sumimag;
    }
}

void idft(float* inreal , float* inimag,float* outreal, float* outimag, int length){
    for (int k = 0; k < length; k++) {  // For each output element
        float sumreal = 0;
        float sumimag = 0;
        for (int t = 0; t < length; t++) {  // For each input element
            float angle =2*M_PI*t*k/((float)length);// 2*PI*t*k / length;
            sumreal +=  inreal[t] * cos(angle) - inimag[t] * sin(angle);
            sumimag +=  inreal[t] * sin(angle) + inimag[t] * cos(angle);
        }
        outreal[k] = sumreal/((float)length);
        outimag[k] = sumimag/((float)length);
    }


}