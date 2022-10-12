// test.cpp

// Simple console app for lib1.c, lib2.c testing

#include <iostream>
#include "lib1.c"
#include "lib2.c"

using namespace std;

int main()
{
    cout << "func1: " << mySqr(11) << endl
         << "func2: " << myAdd(1.3, 1.7) << endl
         << "func3: " << tripleProduct(2, 3, 4) << endl
         << "func4: " << someFunction() << endl
         << "func5: " << myProd(12, 13) << endl
         << "func6: " << myCube(9) << endl
         << "func7: " << functionNotToInclude() << endl
         << "func8: " << myDif(90, 100) << endl;

    return 0;
}