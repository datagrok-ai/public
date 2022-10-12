// lib1.c
// Simple lib with several test functions

#include <emscripten.h>

EMSCRIPTEN_KEEPALIVE
int mySqr(int n){
    return n * n;
} // mySqr

EMSCRIPTEN_KEEPALIVE
float myAdd(float n, float m) {
    return n + m;
} // myAdd

EMSCRIPTEN_KEEPALIVE
int tripleProduct(int a, int b, int c) {
    return a * b * c;
} // tripleProduct

int someFunction(){
    return 2022;
} // someFunction