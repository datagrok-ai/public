// test.h
// Declarations of functions for testing Solver

#ifndef TEST_H
#define TEST_H

// try solver function
void trySolver();

// test of solver wrapper
void testSolverWrapper();

// test computation of k-s using fixedPoint
void testFixedPoint();

// test computation of k-s using fixedPoint, another example
void testFixedPointAnother();

// test of solver of the equation x = g(x) using Newton's method: one-dimensional case
void testNewtonSolver1D();

// test of solver of the equation F(x) = 0 using Newton's method: multi-dimensional case
void testNewtonSolverND();

#endif // !TEST_H

