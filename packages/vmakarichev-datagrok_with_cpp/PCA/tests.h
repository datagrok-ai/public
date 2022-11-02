// test.h
// Headers of test-functions

#ifndef TESTS_H
#define TESTS_H

// Manipulating matrices of different types: PROBLEM when different types!
void differentTypeValuesManipulation();

// Test of PCA that is performed using correlation matrix
void pcaUsingCorrelationMatrixTest();

// Performance comparison of PCAs implemented using correlation matrix
void comparePerformanceOfPCAwithCorMatr();

// Performance investiagtion of PCA using correlation matrix: data is given by void **
void investigatePCAwithCorMatrVoidDataRepresentation(int numOfLaunches = 10);

// Performance investiagtion of PCA using correlation matrix: data is given by float *
void investigatePCAwithCorMatrFloatDataRepresentation(int numOfLaunches = 10);

#endif 

