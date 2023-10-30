// This is supposed to be generated automatically

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

export function solveExapmle() {
  
}