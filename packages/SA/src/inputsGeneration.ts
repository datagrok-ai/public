

export class InputsForVarianceBasedSensitivityAnalysis {
  private N: number; // samplesCount
  private d: number; // dimension

  private A: Array<Float32Array>; // 
  private B: Array<Float32Array>; // 
  private AB: Array<Array<Float32Array>>; //

  constructor(samplesCount: number = 1, dimension: number = 1) {
    this.N = samplesCount;
    this.d = dimension;

    this.A = this.generateMatrix(samplesCount, dimension);
    this.B = this.generateMatrix(samplesCount, dimension);
    this.AB = this.createMixMatrax(this.A, this.B);
  }

  private generateMatrix(rowCount: number, colCount: number): Array<Float32Array> {
    const M = new Array<Float32Array>(colCount);

    for (let j = 0; j < colCount; ++j) {
      M[j] = new Float32Array(rowCount);

      for (let i = 0; i < rowCount; ++i)
        M[j][i] = Math.random();
    }

    return M;
  } // generateMatrix

  private createMixMatrax(A: Array<Float32Array>, B: Array<Float32Array>): Array<Array<Float32Array>> {
    const d = this.d;
    const M = new Array<Array<Float32Array>>(d);

    for (let i = 0; i < d; ++i) {
      M[i] = new Array<Float32Array>(0);

      for (let j = 0; j < d; ++j)
        if (i !== j)
          M[i].push(A[j]);
        else
          M[i].push(B[j]);     
    }

    return M;
  } // createMixMatrax

  private getRow(M: Array<Float32Array>, index: number): Float32Array {
    const row = new Float32Array(this.d);

    for (let j = 0; j < this.d; ++j)
      row[j] = M[j][index];

    return row;
  }

  getRowOfMatrixA(index: number): Float32Array {
    return this.getRow(this.A, index);
  }

  getRowOfMatrixB(index: number): Float32Array {
    return this.getRow(this.B, index);
  }

  getRowOfMatrixAB(matrixIndex: number, rowIndex: number): Float32Array {
    return this.getRow(this.AB[matrixIndex], rowIndex);
  }

  print(): void {
    console.log(`N = ${this.N}, d = ${this.d}`);

    console.log('A:');
    console.log(this.A);

    console.log('B:');
    console.log(this.B);

    for (let i = 0; i < this.d; ++i) {
      console.log(`AB[${i}]:`);
      console.log(this.AB[i]);
    }
  }

  printRows(): void {
    console.log('A:');
    for (let i = 0; i < this.N; ++i) 
      console.log(this.getRowOfMatrixA(i));

    console.log('B:');
    for (let i = 0; i < this.N; ++i) 
      console.log(this.getRowOfMatrixB(i));

    for (let j = 0; j < this.d; ++j) {
      console.log(`AB[${j}]:`);
      for (let i = 0; i < this.N; ++i)
        console.log(this.getRowOfMatrixAB(j, i));
    }
  }
} // InputsForVarianceBasedSensitivityAnalysis

const inputs = new InputsForVarianceBasedSensitivityAnalysis(3, 2);

inputs.print();

inputs.printRows();

//console.log(inputs);

