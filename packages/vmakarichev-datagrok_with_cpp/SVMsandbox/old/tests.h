// test.h

// Headers of funtions for testing methods

#ifndef TESTS_H
#define TESTS_H

void testCreateDataset();

void testCreateNormalizedDataset();

void testTrainModelSimpleLinear();

void testGeneratorLinearSeparable();

void testGeneratorLinearNonSeparable();

void testTrainModelComplexLinear();

void testTrainModelNormalizedDataLinear();

void testTrainModelNormalizedDataLinearHighDim(bool);

void testTrainModelNormalizedDataLinearHighDimDouble();

#endif // TESTS_H

