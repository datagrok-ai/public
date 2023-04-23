#include <emscripten.h>
#include <cstdlib>
extern "C"
{
    void add(int *a, int *b, int len, int *c);
};

EMSCRIPTEN_KEEPALIVE
void add(int *a, int *b, int len, int *c)
{
    for (int i = 0; i < len; i++)
    {
        c[i] = a[i] + b[i];
    }
}
