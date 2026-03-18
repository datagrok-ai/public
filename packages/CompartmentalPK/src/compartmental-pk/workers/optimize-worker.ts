// Multi-Compartment PK — Web Worker for brute-force grid search optimization

import {mrt} from 'diff-grok';

import {
  PkParams, OptimizeWorkerInput, OptimizeWorkerOutput, GridPointResult,
  WeightType, computeMicroConstants, createPkODE, buildSegmentSchedule,
} from '../model';

const ctx: Worker = self as unknown as Worker;

/** Linear interpolation */
function interpolate(t: Float64Array, y: Float64Array, target: number): number {
  if (target <= t[0]) return y[0];
  if (target >= t[t.length - 1]) return y[t.length - 1];

  let lo = 0;
  let hi = t.length - 1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (t[mid] <= target) lo = mid;
    else hi = mid;
  }

  const frac = (target - t[lo]) / (t[hi] - t[lo]);
  return y[lo] + frac * (y[hi] - y[lo]);
}

/** Compute WSSR */
function computeWssr(cObs: number[], cPred: number[], weightType: WeightType): number {
  let wssr = 0;
  for (let i = 0; i < cObs.length; i++) {
    let w = 1;
    if (weightType === '1/C_obs\u00B2' && cObs[i] > 0)
      w = 1 / (cObs[i] * cObs[i]);
    else if (weightType === '1/C_pred\u00B2' && cPred[i] > 0)
      w = 1 / (cPred[i] * cPred[i]);
    const residual = cObs[i] - cPred[i];
    wssr += w * residual * residual;
  }
  return wssr;
}

/** Solve PK and return concentration at observation times */
function solvePkAtTimes(params: PkParams, tObs: number[]): number[] | null {
  try {
    const micro = computeMicroConstants(params);
    const is3comp = params.modelType === '3-Compartment';
    const isOral = params.inputMode === 'Oral';
    const nVars = (is3comp ? 3 : 2) + (isOral ? 1 : 0);

    const segments = buildSegmentSchedule(params);
    const allT: number[] = [];
    const allConc: number[] = [];

    let state = new Array(nVars).fill(0);
    let activeInfusions = 0;

    for (const seg of segments) {
      if (seg.doseEvent) {
        switch (params.inputMode) {
        case 'IV Bolus':
          state[0] += params.dose;
          break;
        case 'IV Infusion':
          activeInfusions++;
          break;
        case 'Oral':
          state[is3comp ? 3 : 2] += params.f * params.dose;
          break;
        }
      }

      if (seg.infusionEnd)
        activeInfusions = Math.max(0, activeInfusions - 1);

      const rIn = params.inputMode === 'IV Infusion' ?
        activeInfusions * params.dose / params.tInf : 0;

      const ode = createPkODE(params, micro, state, seg.start, seg.end, rIn);
      const solution = mrt(ode);

      const tSeg = solution[0];
      const acSeg = solution[1];

      const startIdx = allT.length > 0 ? 1 : 0;
      for (let i = startIdx; i < tSeg.length; i++) {
        allT.push(tSeg[i]);
        allConc.push(acSeg[i] / params.vc);
      }

      const lastIdx = tSeg.length - 1;
      state[0] = solution[1][lastIdx];
      state[1] = solution[2][lastIdx];
      if (is3comp) state[2] = solution[3][lastIdx];
      if (isOral) state[is3comp ? 3 : 2] = solution[is3comp ? 4 : 3][lastIdx];
    }

    // Interpolate at observation times
    const t64 = new Float64Array(allT);
    const c64 = new Float64Array(allConc);
    const cPred: number[] = [];
    for (const to of tObs)
      cPred.push(interpolate(t64, c64, to));

    return cPred;
  } catch {
    return null;
  }
}

ctx.onmessage = (event: MessageEvent<OptimizeWorkerInput>) => {
  try {
    const {gridPoints, paramNames, fixedParams, observedData, weightType} = event.data;
    const results: GridPointResult[] = [];

    for (const point of gridPoints) {
      // Build params by overriding selected parameters
      const params: PkParams = {...fixedParams};
      for (let k = 0; k < paramNames.length; k++)
        (params as any)[paramNames[k]] = point[k];

      const cPred = solvePkAtTimes(params, observedData.tObs);
      if (cPred === null) {
        results.push({paramValues: point, wssr: Infinity});
        continue;
      }

      const wssr = computeWssr(observedData.cObs, cPred, weightType);
      results.push({paramValues: point, wssr});
    }

    const output: OptimizeWorkerOutput = {success: true, results};
    ctx.postMessage(output);
  } catch (err) {
    const output: OptimizeWorkerOutput = {
      success: false,
      error: err instanceof Error ? err.message : 'Unknown error',
    };
    ctx.postMessage(output);
  }
};
