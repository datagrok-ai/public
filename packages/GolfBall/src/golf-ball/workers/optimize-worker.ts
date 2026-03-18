// Golf Ball Flight — Web Worker for height optimization grid search

import {mrt} from 'diff-grok';

import {createGolfBallODE, GolfBallParams} from '../model';
import {OptimizeTask, OptimizeResult} from '../core';

const ctx: Worker = self as unknown as Worker;

ctx.onmessage = (event: MessageEvent<OptimizeTask>) => {
  const {paramName, values, baseParams} = event.data;
  const results: OptimizeResult[] = [];

  for (const val of values) {
    try {
      const params: GolfBallParams = {...baseParams, [paramName]: val};
      const ode = createGolfBallODE(params);
      const solution = mrt(ode);
      const yArr = solution[2];

      let maxHeight = 0;
      for (let i = 0; i < yArr.length; i++) {
        if (yArr[i] > maxHeight) maxHeight = yArr[i];
      }

      results.push({value: val, maxHeight});
    } catch (err) {
      results.push({
        value: val,
        maxHeight: -1,
        error: err instanceof Error ? err.message : 'Unknown error',
      });
    }
  }

  ctx.postMessage(results);
};
