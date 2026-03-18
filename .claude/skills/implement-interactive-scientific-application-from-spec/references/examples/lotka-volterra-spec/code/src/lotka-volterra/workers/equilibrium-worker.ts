// Lotka-Volterra — Web Worker for equilibrium point exploration

import {EquilibriumTask, EquilibriumResult, LotkaVolterraParams} from '../core';

const ctx: Worker = self as unknown as Worker;

ctx.onmessage = (event: MessageEvent<EquilibriumTask>) => {
  const {paramName, paramMin, paramMax, steps, baseParams} = event.data;

  // DO NOT USE this artificial delay in production code! It's here just to simulate a long-running computation for demonstration purposes.
  const start = Date.now();
  while (Date.now() - start < 5000) { /* busy wait */ }

  try {
    const paramValues: number[] = [];
    const xStar: number[] = [];
    const yStar: number[] = [];

    for (let i = 0; i < steps; i++) {
      const val = paramMin + i * (paramMax - paramMin) / (steps - 1);
      paramValues.push(val);

      const params: LotkaVolterraParams = {...baseParams, [paramName]: val};
      xStar.push(params.gamma / params.delta);
      yStar.push(params.alpha / params.beta);
    }

    const result: EquilibriumResult = {paramValues, xStar, yStar};
    ctx.postMessage(result);
  } catch (err) {
    const result: EquilibriumResult = {
      paramValues: [],
      xStar: [],
      yStar: [],
      error: err instanceof Error ? err.message : 'Unknown error',
    };
    ctx.postMessage(result);
  }
};
