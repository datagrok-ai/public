// simple.cpp

// Simple library example

#include <emscripten.h>  // <-- this is required for Emscripten Tool

extern "C" {  // <-- this is required in order to get an access to the function srt further
    int sqr(int num);
}

EMSCRIPTEN_KEEPALIVE // <-- this is required for Emscripten Tool
int sqr(int num)
{
    return num * num;
} // sqr