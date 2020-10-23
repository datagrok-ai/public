//#ifndef FUNCTIONS_H_INCLUDED
//#define FUNCTIONS_H_INCLUDED
//Forward Fourier Transform
void dft( float* inreal , float* inimag,float* outreal, float* outimag, int length);

//Inverse Fourier Transform
void idft(float* inreal , float* inimag,float* outreal, float* outimag, int length);

//Fast Fourier Transform
void fft(float *inreal, float *inimag, float *outreal, float *outimag, int length, bool invert);

//#endif