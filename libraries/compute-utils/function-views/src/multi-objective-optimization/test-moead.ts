/*import {Moead} from './moead';

// Objectives
const func = (x: Float32Array): Float32Array => {
  return new Float32Array([
    x[0]**2 + x[1]**2, // f1
    (x[0] - 2)**2 + (x[1] - 1)**2, // f2
  ]);
};

// Input options
const inputOpts = {
  dim: 2,
  mins: new Float32Array([0, 0]),
  maxs: new Float32Array([2, 2]),
};

// Output dimension
const dim = 2;

// Optimizer
const solution = new Moead(func, inputOpts, dim).perform();

// Print results
for (const sol of solution)
  console.log(sol.point.toString(), ',', sol.objective.toString());*/
