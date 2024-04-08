/** Monte Carlo optimizer */

import {Extremum, OptimizationTask} from './definitions';

export async function getMinimum(task: OptimizationTask): Promise<Extremum> {
  // 1. Check dimensions
  if (task.minVals.length !== task.maxVals.length)
    throw new Error('Non-equal size of arrays with min & max values');

  if (task.minVals.length === 0)
    throw new Error('Empty arrays with min & max values');

  // 2. Init optimization
  const dim = task.minVals.length;
  const argMin = new Float32Array(dim);
  const arg = new Float32Array(dim);
  const step = new Float32Array(dim);

  task.minVals.forEach((val, idx) => {
    argMin[idx] = val;
    step[idx] = task.maxVals[idx] - val;
  });

  let minCost = await task.costFunc(argMin);
  let cost: number;

  // 3. Search for min
  for (let i = 0; i < task.samplesCount; ++i) {
    for (let j = 0; j < dim; ++j)
      arg[j] = task.minVals[j] + step[j] * Math.random();

    cost = await task.costFunc(arg);

    if (cost < minCost) {
      minCost = cost;
      arg.forEach((val, idx) => argMin[idx] = val);
    }
  }

  return {
    arg: argMin,
    cost: minCost,
  };
}
