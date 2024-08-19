#include <stdio.h>
#include <iostream>
#include "dbscan.h"
#include <emscripten.h>

#define MINIMUM_POINTS 4     // minimum number of cluster
#define EPSILON (0.75*0.75)  // distance for clustering

extern "C"{
    void dbscan(float* embedX, float* embedY, int numPoints, double epsilon, int min_points, int* clusterIds);
}

EMSCRIPTEN_KEEPALIVE
void dbscan (float* embedX, float* embedY, int numPoints, double epsilon, int min_points, int* clusterIds) {
    vector<Point> points;

    for (int i = 0; i < numPoints; i++) {
        Point p = Point();
        p.x = embedX[i];
        p.y = embedY[i];
        p.clusterID = UNCLASSIFIED;
        points.push_back(p);
    }
    DBSCAN ds(min_points, epsilon, points);
    ds.run();
    for (int i = 0; i < numPoints; i++) {
        clusterIds[i] = ds.m_points[i].clusterID;
    }
    return;
}