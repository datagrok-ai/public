// Lotka-Volterra — Web Worker for optimization grid search

import {mrt} from 'diff-grok';

import {createLotkaVolterraODE} from '../model';
import {WorkerTask, WorkerResult} from '../core';

const ctx: Worker = self as unknown as Worker;

ctx.onmessage = (event: MessageEvent<WorkerTask[]>) => {
  const tasks = event.data;
  const results: WorkerResult[] = [];

  for (const task of tasks) {
    try {
      const ode = createLotkaVolterraODE({
        alpha: task.alpha,
        beta: task.beta,
        delta: task.delta,
        gamma: task.gamma,
        x0: task.x0,
        y0: task.y0,
        T: task.T,
      });

      const solution = mrt(ode);
      const xValues = solution[1];

      let maxPrey = 0;
      for (let i = 0; i < xValues.length; i++) {
        if (xValues[i] > maxPrey) maxPrey = xValues[i];
      }

      results.push({alpha: task.alpha, beta: task.beta, delta: task.delta, gamma: task.gamma, maxPrey});
    } catch (err) {
      results.push({
        alpha: task.alpha,
        beta: task.beta,
        delta: task.delta,
        gamma: task.gamma,
        maxPrey: -1,
        error: err instanceof Error ? err.message : 'Unknown error',
      });
    }
  }

  ctx.postMessage(results);
};
