// lib2.c
// Another simple lib with several test functions

#include <emscripten.h>

EMSCRIPTEN_KEEPALIVE
int myProd(int n, int m){
    return n * m;
} // myProd

EMSCRIPTEN_KEEPALIVE
float myCube(float n) {
    return n * n * n;
} // myAdd

int functionNotToInclude(){
    return 2023;
} // functionNotToInclude

EMSCRIPTEN_KEEPALIVE
int myDif(int a, int b) {
    return a - b;
} // myDif

