#include <iostream>
using namespace std;

#include "sum.c"

int main()
{
    cout << "Hi" << endl;
    
    int arr[] = {1,2,3,4,5};
    int len = 5;

    cout << "sum is " << sum(arr, len) << endl;
    
    return 0;
}