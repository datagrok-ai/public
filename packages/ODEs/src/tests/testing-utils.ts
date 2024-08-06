import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {ODEs} from '../solver-tools/solver-defs';

enum NAMES {
  T = 't',
  X = 'x',
  Y = 'y',
};

type Solution2D = {
  x: number,
  y: number,
};

/** Stiff test problem definition */
const stiffTest = {
  name: 'Stiff 2',
  arg: {name: NAMES.T, start: 0, finish: 4, step: 0.01},
  initial: [52.29, 83.82],
  func: (t: number, _y: Float64Array, _output: Float64Array) => {
    _output[0] = -5 * _y[0] + 3 * _y[1];
    _output[1] = 100 * _y[0] - 301 * _y[1];
  },
  tolerance: 0.00005,
  solutionColNames: [NAMES.X, NAMES.Y],
};

/** Exact solution of the stiff test problem */
const exactSolution = (t: number) => {
  return {
    x: 52.96 * Math.exp(-3.9899 * t) - 0.67 * Math.exp(-302.0101 * t),
    y: 17.83 * Math.exp(-3.9899 * t) + 65.99 * Math.exp(-302.0101 * t),
  };
};

export function evaluateMethod(method: (odes: ODEs) => DG.DataFrame): number {
  let error = 0;

  // Get numerical solution
  const approxSolution = method(stiffTest);
  const rows = approxSolution.rowCount;
  const arg = approxSolution.col(NAMES.T)!.getRawData();
  const xApprox = approxSolution.col(NAMES.X)!.getRawData();
  const yApprox = approxSolution.col(NAMES.Y)!.getRawData();

  let exact: Solution2D;

  // Compute error
  for (let i = 0; i < rows; ++i) {
    exact = exactSolution(arg[i]);

    error = Math.max(
      error,
      Math.abs(xApprox[i] - exact.x),
      Math.abs(yApprox[i] - exact.y),
    );
  }
  return error;
}
