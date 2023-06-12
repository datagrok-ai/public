#include <iostream>
#include "../../../../../Eigen/Eigen/Dense"
#include <emscripten.h>

using namespace Eigen;
using namespace std;

/* One-way Analysis of Variance (ANOVA) is a statistical method for testing for
   differences in the means of three or more groups.
   https://www.cuemath.com/anova-formula
   columns - arrays of observations or measurements taken at different time points or conditions;
   rowCount - number of levels (groups) of the factor;
   columnCount - number of observations in each group. */

extern "C" {
  double oneWayAnova(float * columns, int rowCount, int columnCount);
}

//name: oneWayAnova
//input: dataframe table
//input: column_list columns
//output: double ftest [_callResult]
EMSCRIPTEN_KEEPALIVE
double oneWayAnova(float * columns, int rowCount, int columnCount) {
  // Input data
  Map<Matrix<float, Dynamic, Dynamic, ColMajor> > data(columns, rowCount, columnCount);

  // Calculate overall mean
  double overallMean = data.mean();

  // Calculate group means
  VectorXf groupMeans = data.rowwise().mean();

  // Calculate group sample sizes
  VectorXf groupSampleSizes = VectorXf::Zero(rowCount);
  groupSampleSizes.fill(columnCount);

  // Calculate sum of squares between groups
  VectorXf ssb = groupSampleSizes.array() * (groupMeans.array() - overallMean).square();
  double ssbTotal = ssb.sum();

  // Calculate sum of squares within groups
  VectorXf ssr = VectorXf::Zero(rowCount);
  for (int i = 0; i < rowCount; ++i) {
      VectorXf group_data = data.row(i);
      ssr(i) = (group_data.array() - groupMeans(i)).square().sum();
  }
  double ssrTotal = ssr.sum();

  // Calculate degrees of freedom
  int dfBetween = rowCount - 1;
  int dfWithin = rowCount * (columnCount - 1);

  // Calculate mean squares
  double msb = ssbTotal / dfBetween;
  double msr = ssrTotal / dfWithin;

  // Calculate F-statistic
  double fStatistic = msb / msr;

  return fStatistic;
}
