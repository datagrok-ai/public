// Levins Metapopulation Model — Web Worker for optimization task
// Uses mrt from diff-grok via the shared model definition

import {mrt} from 'diff-grok';

import {createLevinsODE} from './model';
import {WorkerTask, WorkerResult} from './core';

const ctx: Worker = self as unknown as Worker;

ctx.onmessage = (event: MessageEvent<WorkerTask>) => {
  const task = event.data;

  try {
    const ode = createLevinsODE({
      p0: task.p0,
      m: task.m_i,
      e0: task.e0,
      rescueEffect: task.rescueEffect,
      t_start: task.t_start,
      t_end: task.t_end,
      t_step: task.t_step,
      tolerance: task.tolerance,
    });

    const solution = mrt(ode);
    const pValues = solution[1];
    const p_end = pValues[pValues.length - 1];

    const result: WorkerResult = {m_i: task.m_i, p_end};
    ctx.postMessage(result);
  } catch (err) {
    const result: WorkerResult = {
      m_i: task.m_i,
      p_end: -1,
      error: err instanceof Error ? err.message : 'Unknown error',
    };
    ctx.postMessage(result);
  }
};
