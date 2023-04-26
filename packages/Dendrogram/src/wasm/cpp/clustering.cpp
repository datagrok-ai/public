#include <math.h>
#include <string>
#include <iostream>
#include "fastcluster.h"
#include <iomanip>
#include <sstream>
// #include <emscripten.h>

using namespace std;

extern "C"{
    void getDendrogram(float *distmat, int npoints, int optMethod, int *merge, float *heights);
}

// string toNewick(int *merge, double *heights, const int n) {

//     string *clusters = new string[n - 1]{};
//     double *clusterHeights = new double[n - 1]{};
//     string a = "(";
//     int left,right;
//     for (int i = 0; i < n - 1; i++) {
//         double height = 0;
//         a = "(";
//         left = merge[i];
//         right = merge[i + n - 1];
//         stringstream stream;
//         stream << fixed << setprecision(2) << heights[i];
//         string curHeight = stream.str();
        
//         if(left<0){
//             // its in R format so indexing starts with 1..
//             height = curHeight;
//             a += to_string(-left - 1) + ":" + height + ",";
            
//         }
//         else{
//             height = curHeight + heights[left - 1];
//             a += clusters[left - 1] + ":" + height + ",";
            
//         }

//         if(right<0){
//             height = curHeight;
//             a += to_string(-right - 1) + ":" + height + ")";
//         }
//         else{
//             height = curHeight + heights[right - 1];
//             a += clusters[right - 1] + ":" + curHeight + ")";
//         }
//         clusters[i] = a;
//     }
//     string out = clusters[n - 2] + ";";
//     delete[] clusters;
//     return out;
// }

// EMSCRIPTEN_KEEPALIVE
void getDendrogram(float *distmat, int npoints, int optMethod, int *merge, float *heights){

    hclust_fast(npoints, distmat, optMethod, merge, heights);
    delete[] distmat;
    // don't forget to delete merges and heights from the client;
    return;
}
