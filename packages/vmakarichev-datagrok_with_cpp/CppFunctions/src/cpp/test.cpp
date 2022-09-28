// test.cpp

// Test of cpp-functions

#include <iostream>
#include "simple.cpp"  // <-- file with simple library example

using namespace std;

// test of the simple cpp-library
void testSimpleLib();

int main()
{

    testSimpleLib();

    return 0;

} // main


// test of the simple cpp-library
void testSimpleLib() {
    for(int i = 0; i <= 10; i++)
       cout << i << "^2 = " << sqr(i) << endl;
} // testSimpleLib