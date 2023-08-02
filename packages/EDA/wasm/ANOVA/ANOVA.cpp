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
  double oneWayAnova(float * columns, int rowCount, int columnCount, float * predictionColumn, int resLength);
}

//name: oneWayAnova
//input: dataframe table
//input: column_list columns
//output: column res [new(columns.rowCount)]
EMSCRIPTEN_KEEPALIVE
double oneWayAnova(float * columns, int rowCount, int columnCount, float * res, int resRowCount) {
  // Input data
  Map<Matrix<float, Dynamic, Dynamic, ColMajor>> data(columns, rowCount, columnCount);
  // Calculate overall mean
  double overallMean = data.mean();

  // Calculate group means
  VectorXf groupMeans = data.colwise().mean();

  // Calculate group sample sizes
  VectorXf groupSampleSizes = VectorXf::Zero(columnCount);
  groupSampleSizes.fill(rowCount);

  // Calculate sum of squares between groups
  VectorXf ssb = groupSampleSizes.array() * (groupMeans.array() - overallMean).square();
  double ssbTotal = ssb.sum();

  // Calculate sum of squares within groups
  VectorXf ssr = VectorXf::Zero(columnCount);
  for (int i = 0; i < columnCount; ++i) {
      VectorXf group_data = data.col(i);
      ssr(i) = (group_data.array() - groupMeans(i)).square().sum();
  }
  double ssrTotal = ssr.sum();

  // Calculate degrees of freedom
  int dfBetween = columnCount - 1;
  int dfWithin = columnCount * (rowCount - 1);

  // Calculate mean squares
  double msb = ssbTotal / dfBetween;
  double msr = ssrTotal / dfWithin;

  // Calculate F-statistic
  double fStatistic = msb / msr;
  Map<Vector<float, Dynamic>> column(res, resRowCount);
  column[0] = ssbTotal;
  column[1] = ssrTotal;
  column[2] = dfBetween;
  column[3] = dfWithin;
  column[4] = msb;
  column[5] = msr;
  column[6] = fStatistic;

  return fStatistic;
}
