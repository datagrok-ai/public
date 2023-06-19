/* inputs-generation.ts

   Generator of inputs for variance-based sensitivity analysis.

   Computation of Sobol indeces requires the set of matrices {A, B and AB[i]} (see [1, 2] for more details).
   Generation of these matrices is provided below.

   REMARK. We use notations from the following sources. It provides simplicity of the further
           computations implementation.

   References
     [1] A. Saltelli, et. al., "Variance based sensitivity analysis of model output. Design and estimator
         for the total sensitivity index", DOI: https://doi.org/10.1016/j.cpc.2009.09.018

     [2] Variance-based sensitivity analysis,
         LINK: https://en.wikipedia.org/wiki/Variance-based_sensitivity_analysis
*/

import {checkSize} from './utils';

// Inputs generator
class InputsForVarianceBasedSensitivityAnalysis {
  private N: number; // samplesCount
  private d: number; // dimension, i.e. model inputs count

  private A: Array<Float32Array>; // columns of the matrix A
  private B: Array<Float32Array>; // columns of the matrix B
  private AB: Array<Array<Float32Array>>; // each element is an array of columns of the matrices AB

  constructor(samplesCount: number = 1, dimension: number = 1) {
    // check sizes
    checkSize(samplesCount);
    checkSize(dimension);

    this.N = samplesCount;
    this.d = dimension;

    this.A = this.generateMatrix(samplesCount, dimension);
    this.B = this.generateMatrix(samplesCount, dimension);
    this.AB = this.createMixMatrax(this.A, this.B);
  }

  // Returns array of columns of randomly generated elements with uniform distribution on the segemnt [0, 1]
  private generateMatrix(rowCount: number, colCount: number): Array<Float32Array> {
    const M = new Array<Float32Array>(colCount);

    for (let j = 0; j < colCount; ++j) {
      M[j] = new Float32Array(rowCount);

      for (let i = 0; i < rowCount; ++i)
        M[j][i] = Math.random();
    }

    return M;
  } // generateMatrix

  // Returns a set of matrices AB (see [1, 2] for more details).
  private createMixMatrax(A: Array<Float32Array>, B: Array<Float32Array>): Array<Array<Float32Array>> {
    const d = this.d;
    const AB = new Array<Array<Float32Array>>(d);

    /* For any i, the matrix AB[i] has the same columns as the matrix A except i-th column
       that is taken from the matrix B. */
    for (let i = 0; i < d; ++i) {
      AB[i] = new Array<Float32Array>(0);

      for (let j = 0; j < d; ++j) {
        if (i !== j)
          AB[i].push(A[j]); // copy all non i-th columns of A
        else
          AB[i].push(B[j]);
      } // i-th column is taken from B
    }

    return AB;
  } // createMixMatrax

  // Returns the specified row of the matrix, which is given as an array of columns.
  private getRow(M: Array<Float32Array>, index: number): Float32Array {
    const row = new Float32Array(this.d);

    for (let j = 0; j < this.d; ++j)
      row[j] = M[j][index];

    return row;
  }

  // Returns the specified row of A
  getRowOfMatrixA(index: number): Float32Array {
    return this.getRow(this.A, index);
  }

  // Returns the specified row of B
  getRowOfMatrixB(index: number): Float32Array {
    return this.getRow(this.B, index);
  }

  // Returns the specified row of the specified matrix AB
  getRowOfMatrixAB(matrixIndex: number, rowIndex: number): Float32Array {
    return this.getRow(this.AB[matrixIndex], rowIndex);
  }

  // Returns concatenated column of A,B and AB's
  getConcatenatedColumns(): Array<Float32Array> {
    const columns = new Array<Float32Array>(0);

    const N = this.N;
    const d = this.d;

    for (let i = 0; i < d; ++i) {
      const col = new Float32Array(N * (d + 2));

      col.set(this.A[i]);

      col.set(this.B[i], N);

      for (let j = 0; j < d; ++j)
        col.set(this.AB[j][i], 2 * N + N * j);

      columns.push(col);
    }

    return columns;
  }
} // InputsForVarianceBasedSensitivityAnalysis

// Returns input columns data for variance-based sensitivity analysis
export function getGeneratedColumnsData(samplesCount: number, dimension: number): Array<Float32Array> {
  const inputs = new InputsForVarianceBasedSensitivityAnalysis(samplesCount, dimension);

  return inputs.getConcatenatedColumns();
}

