// Golf Ball Flight — Web Worker for sensitivity analysis

import {mrt} from 'diff-grok';

import {createGolfBallODE, GolfBallParams} from '../model';
import {SensitivityTask, SensitivityResult} from '../core';

const ctx: Worker = self as unknown as Worker;

ctx.onmessage = (event: MessageEvent<SensitivityTask>) => {
  const {paramName, paramMin, paramMax, steps, baseParams} = event.data;

  try {
    const paramValues: number[] = [];
    const distances: number[] = [];

    for (let i = 0; i < steps; i++) {
      const val = paramMin + i * (paramMax - paramMin) / (steps - 1);
      paramValues.push(val);

      const params: GolfBallParams = {...baseParams, [paramName]: val};
      const ode = createGolfBallODE(params);
      const solution = mrt(ode);
      const xArr = solution[1];
      const yArr = solution[2];

      // Find flight distance (x at landing)
      let dist = xArr[xArr.length - 1];
      for (let j = 1; j < yArr.length; j++) {
        if (yArr[j] <= 0) {
          const y0 = yArr[j - 1];
          const y1 = yArr[j];
          const frac = y0 / (y0 - y1);
          dist = xArr[j - 1] + frac * (xArr[j] - xArr[j - 1]);
          break;
        }
      }

      distances.push(dist);
    }

    const result: SensitivityResult = {paramName, paramValues, distances};
    ctx.postMessage(result);
  } catch (err) {
    const result: SensitivityResult = {
      paramName,
      paramValues: [],
      distances: [],
      error: err instanceof Error ? err.message : 'Unknown error',
    };
    ctx.postMessage(result);
  }
};
