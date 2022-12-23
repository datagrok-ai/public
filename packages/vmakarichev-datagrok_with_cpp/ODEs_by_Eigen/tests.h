// tests.h

// Declaration of functions for testing methods

#ifndef TESTS_H
#define TESTS_H

// test of one-step method: one-dimensional case
void testOneStepMethodOneDim();

// test of one-step method: multi-dimensional case
void testOneStepMethodMultiDim();

// analysis of solving stiff problem: one-dimensional case
void analyzeStiffProblemOneDim();

// analysis of solving stiff problem: multi-dimensional case
void analyzeStiffProblemMultiDim();

// test ode23s solver
void testODE23sSolver();

#endif // !TESTS_H

