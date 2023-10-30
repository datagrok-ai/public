import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/** Right-hand side of IVP */
type Func = (t: number, y: Float64Array, p: Float64Array) => Float64Array;

/** Initial Value Problem (IVP) type */
type ODEs = {
  name: string,
  arg: {
    name: string,
    start: number,
    finish: number,
    step: number
  },
  initial: Float64Array,
  params: Float64Array,
  func: Func,
  solutionColNames: string[],
};

const ex1 = (t: number, _y: Float64Array, _p: Float64Array) => {
    // constants
    const const1 = 1.0;
    const const2 = 3.0;
  
    // parameters
    const param1 = _p[0];
    const param2 = _p[1];
  
    // function values
    const x = _y[0];
    const y = _y[1];
  
    // expressions
    const coef1 = const1 + param1;
    const coef2 = const2 + param2 + 0.0;
  
    return new Float64Array([
      coef1 * y, // d(x)/dt
      coef2 * x // dy/d(t)
    ]);
} // ex1

  
const ex2 = (t: number, _y: Float64Array, _p: Float64Array) => {  
    // parameters
    const alpha = _p[0];
    const beta = _p[1];
  
    // function values
    const x = _y[0];
    const y = _y[1];  
  
    return new Float64Array([
      alpha * x, // d(x)/dt
      beta * y // dy/d(t)
    ]);
} // ex2

export const example2: ODEs = {
  name: 'Example 2',
  arg: {name: 't', start: 0, finish: 5, step: 0.01},
  initial: new Float64Array([1, 1]),
  params: new Float64Array([2, 3]),
  func: ex2,
  solutionColNames: ['x(t)', 'y(t)']
};

/** Returns derivative with respect to t. */
function tDerivative(t: number, y: Float64Array, p: Float64Array, f: Func, eps: number): Float64Array {
  const size = y.length;
  const res = new Float64Array(size);
  const f0 = f(t, y, p);
  const f1 = f(t + eps, y, p);

  for (let i = 0; i < size; ++i)
    res[i] = (f1[i] - f0[i]) / eps;

  return res;
}

/** Returns Jacobian. */
function Jacobian(t: number, y: Float64Array, p: Float64Array, f: Func, eps: number): Float64Array {
  const size = y.length;
  const res = new Float64Array(size * size);
  const f0 = f(t, y, p);

  for (let j = 0; j < size; ++j) {
    y[j] += eps;

    const f1 = f(t, y, p);

    for (let i = 0; i < size; ++i)
      res[j + i * size] = (f1[i] - f0[i]) / eps;

    y[j] -= eps;
  }

  return res;
}

/** Solve initial value problem. */
export function solveODEs(problem: ODEs): DG.DataFrame {
  const start = problem.arg.start;
  const finish = problem.arg.finish;
  const step = problem.arg.step;
  const initial = problem.initial;

  /** number of solution dataframe rows */  
  const rowCount = Math.trunc((finish - start) / step) + 1;

  /** dimension of the problem */
  const dim = initial.length;

  /**  */
  const tArr = new Float32Array(rowCount);

  /**  */
  const yArrs = Array<Float32Array>(dim);

  for (let i = 0; i < dim; ++i)
    yArrs[i] = new Float32Array(rowCount);

  const solutionDf = DG.DataFrame.fromColumns([DG.Column.fromFloat32Array(problem.arg.name, tArr)]);
  yArrs.forEach((arr, idx) => solutionDf.columns.add(DG.Column.fromFloat32Array(problem.solutionColNames[idx], arr)));
  solutionDf.name = problem.name;

  return solutionDf;
}

/*const y = new Float64Array([1, 12]);
const p = new Float64Array([15, 7]);

//console.log(T(2, y, p, ex2, 0.00000001));

console.log(Jac(0, y, p, ex2, 0.00000001));*/
