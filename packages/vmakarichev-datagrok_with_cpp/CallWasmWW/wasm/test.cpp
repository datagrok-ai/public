// Simple test of functions from lib.cpp

#include <iostream>
using namespace std;

#include "lib.cpp"

template<typename T>
void printArr(T * arr, int arrLength) {
    for(int i = 0; i < arrLength; i++)
        cout << "  " << arr[i];
    cout << "\n";
}

int main()
{
    cout << "Test of lib.cpp\n\n";

    cout << "1. Test of sum:\n  10 + 20 = " << sum(10, 20) << "\n\n";

    const int rowCount = 3;
    const int columnCount = 5;

    float floatCol[rowCount] = {1.2, 2.2, 3.3};
    cout << "2. Test of maxFloatCol:\n arr:";
    printArr(floatCol, rowCount);
    cout << " maxFloatCol: " << maxFloatCol(floatCol, rowCount) << "\n\n";

    int intCol[rowCount] = {-3, 5, 0};
    cout << "3. Test of maxIntCol:\n arr:";
    printArr(intCol, rowCount);
    cout << " maxIntCol: " << maxIntCol(intCol, rowCount) << "\n\n";

    float otherFloatCol[rowCount] = {4.2, 9.2, 1.3};
    cout << "3. Test of addFloatCols:\n col1:";
    printArr(floatCol, rowCount);
    cout << " col2:";
    printArr(otherFloatCol, rowCount);
    float sumOfFloatCols[rowCount] = {0.0};
    addFloatCols(floatCol, rowCount, otherFloatCol, rowCount, sumOfFloatCols, rowCount);
    cout << " sumOfFloatCols: ";
    printArr(sumOfFloatCols, rowCount);
    cout << "\n\n";

    int otherIntCol[rowCount] = {4, 9, 13};
    cout << "4. Test of addIntCols:\n col1:";
    printArr(intCol, rowCount);
    cout << " col2:";
    printArr(otherIntCol, rowCount);
    int sumOfIntCols[rowCount] = {0};
    addIntCols(intCol, rowCount, otherIntCol, rowCount, sumOfIntCols, rowCount);
    cout << " sumOfIntFloatCols: ";
    printArr(sumOfIntCols, rowCount);
    cout << "\n\n";

    int intCols[rowCount * columnCount] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    int doubledIntCols[rowCount * columnCount] = { 0 };
    cout << "5. Test of doubledInts:\n columns:";
    printArr(intCols, rowCount * columnCount);
    doubledInts(intCols, rowCount, columnCount, doubledIntCols, rowCount, columnCount);
    cout << " doubled columns: ";
    printArr(doubledIntCols, rowCount * columnCount);
    cout << "\n\n";

    float floatCols[rowCount * columnCount] 
        = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 
           0.8, 0.9, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15};
    float doubledFloatCols[rowCount * columnCount] = { 0 };
    cout << "6. Test of doubledFloats:\n columns:";
    printArr(floatCols, rowCount * columnCount);
    doubledFloats(floatCols, rowCount, columnCount, doubledFloatCols, rowCount, columnCount);
    cout << " doubled columns: ";
    printArr(doubledFloatCols, rowCount * columnCount);
    cout << "\n\n";

    return 0;
}