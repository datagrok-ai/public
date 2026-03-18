// Multi-Compartment Pharmacokinetic Model — Computational Core

import {mrt} from 'diff-grok';

import {
  PkParams, PkSolution, ObservedData, DataGenParams,
  ModelType, InputMode, SamplingSchedule, WeightType,
  computeMicroConstants, buildSegmentSchedule, createPkODE,
  SAMPLING_SCHEDULES,
} from './model';

// --- Re-exports ---

export type {
  PkParams, PkSolution, ObservedData, DataGenParams,
  ModelType, InputMode, SamplingSchedule, WeightType,
  SelectedParam, OptimizeParams, OptimizeResult, LandscapeData,
  GridPointResult, OptimizeWorkerInput, OptimizeWorkerOutput,
} from './model';
export {computeMicroConstants, createPkODE, SAMPLING_SCHEDULES} from './model';

// --- Input IDs ---

export type InputId =
  | 'ctrl_model_type' | 'ctrl_input_mode' | 'ctrl_y_scale'
  | 'ctrl_cl' | 'ctrl_vc' | 'ctrl_vp' | 'ctrl_q' | 'ctrl_vd' | 'ctrl_q2'
  | 'ctrl_ka' | 'ctrl_f'
  | 'ctrl_dose' | 'ctrl_repeated' | 'ctrl_tau' | 'ctrl_n_doses' | 'ctrl_t_inf'
  | 'ctrl_t_end' | 'ctrl_mec' | 'ctrl_mtc';

export type ValidationErrors = Map<InputId, string>;

// --- Defaults ---

export const DEFAULTS: PkParams = {
  modelType: '2-Compartment',
  inputMode: 'IV Bolus',
  cl: 5,
  vc: 20,
  vp: 40,
  q: 10,
  vd: 60,
  q2: 2,
  ka: 1.0,
  f: 0.8,
  dose: 100,
  repeatedDosing: false,
  tau: 12,
  nDoses: 5,
  tInf: 1,
  tEnd: 72,
  mec: 1.0,
  mtc: 20.0,
};

// --- Slider ranges ---

export const RANGES: Record<string, {min: number; max: number}> = {
  cl: {min: 0.1, max: 50},
  vc: {min: 1, max: 100},
  vp: {min: 1, max: 200},
  q: {min: 0.1, max: 50},
  vd: {min: 1, max: 200},
  q2: {min: 0.01, max: 5},
  ka: {min: 0.1, max: 5},
  f: {min: 0.01, max: 1.0},
  dose: {min: 1, max: 1000},
  tau: {min: 1, max: 48},
  nDoses: {min: 2, max: 20},
  tInf: {min: 0.5, max: 24},
  tEnd: {min: 1, max: 240},
  mec: {min: 0.01, max: 50},
  mtc: {min: 0.1, max: 100},
};

// --- Validation ---

/** Validate all PK inputs; returns map of input ID -> error message */
export function validate(inputs: PkParams): ValidationErrors {
  const errors: ValidationErrors = new Map();

  // val_01 through val_04: always checked
  if (inputs.cl <= 0)
    errors.set('ctrl_cl', 'Clearance must be positive');
  if (inputs.vc <= 0)
    errors.set('ctrl_vc', 'Central volume must be positive');
  if (inputs.vp <= 0)
    errors.set('ctrl_vp', 'Peripheral volume must be positive');
  if (inputs.q <= 0)
    errors.set('ctrl_q', 'Intercompartmental clearance must be positive');

  // val_05, val_06: 3-compartment only
  if (inputs.modelType === '3-Compartment') {
    if (inputs.vd <= 0)
      errors.set('ctrl_vd', 'Deep compartment volume must be positive');
    if (inputs.q2 <= 0)
      errors.set('ctrl_q2', 'Deep intercompartmental clearance must be positive');
  }

  // val_07, val_08: oral only
  if (inputs.inputMode === 'Oral') {
    if (inputs.ka <= 0)
      errors.set('ctrl_ka', 'Absorption rate must be positive');
    if (inputs.f <= 0 || inputs.f > 1)
      errors.set('ctrl_f', 'Bioavailability must be between 0 (exclusive) and 1 (inclusive)');
  }

  // val_09: dose
  if (inputs.dose <= 0)
    errors.set('ctrl_dose', 'Dose must be positive');

  // val_10, val_11: repeated dosing
  if (inputs.repeatedDosing) {
    if (inputs.tau <= 0)
      errors.set('ctrl_tau', 'Dosing interval must be positive');
    if (inputs.nDoses < 1)
      errors.set('ctrl_n_doses', 'Number of doses must be at least 1');
  }

  // val_12: infusion
  if (inputs.inputMode === 'IV Infusion') {
    if (inputs.tInf <= 0)
      errors.set('ctrl_t_inf', 'Infusion duration must be positive');
  }

  // val_13: infusion + repeated (only if val_10 and val_12 passed)
  if (inputs.inputMode === 'IV Infusion' && inputs.repeatedDosing) {
    if (!errors.has('ctrl_tau') && !errors.has('ctrl_t_inf')) {
      if (inputs.tInf > inputs.tau)
        errors.set('ctrl_t_inf', 'Infusion duration cannot exceed dosing interval');
    }
  }

  // val_14: tEnd
  if (inputs.tEnd <= 0)
    errors.set('ctrl_t_end', 'Simulation time must be positive');

  // val_15, val_16: therapeutic window
  if (inputs.mec < 0)
    errors.set('ctrl_mec', 'MEC cannot be negative');
  if (inputs.mtc <= inputs.mec)
    errors.set('ctrl_mtc', 'MTC must be greater than MEC');

  return errors;
}

// --- PK Solver ---

/** Solve the PK ODE system with segmented integration for multiple dosing */
export function solvePk(params: PkParams): PkSolution {
  const micro = computeMicroConstants(params);
  const is3comp = params.modelType === '3-Compartment';
  const isOral = params.inputMode === 'Oral';
  const nVars = (is3comp ? 3 : 2) + (isOral ? 1 : 0);

  const segments = buildSegmentSchedule(params);

  // Collect all solution arrays
  const allT: number[] = [];
  const allAc: number[] = [];
  const allAp: number[] = [];
  const allAd: number[] = is3comp ? [] : [];
  const allAgut: number[] = isOral ? [] : [];

  // Current state
  let state = new Array(nVars).fill(0);

  // Track active infusions for overlapping infusion doses
  let activeInfusions = 0;

  for (const seg of segments) {
    // Apply dose event
    if (seg.doseEvent) {
      switch (params.inputMode) {
      case 'IV Bolus':
        state[0] += params.dose; // A_c += Dose
        break;
      case 'IV Infusion':
        activeInfusions++;
        break;
      case 'Oral':
        state[is3comp ? 3 : 2] += params.f * params.dose; // A_gut += F * Dose
        break;
      }
    }

    // Handle infusion end
    if (seg.infusionEnd)
      activeInfusions = Math.max(0, activeInfusions - 1);

    // Compute R_in for this segment
    const rIn = params.inputMode === 'IV Infusion' ?
      activeInfusions * params.dose / params.tInf : 0;

    // Solve ODE for this segment
    const ode = createPkODE(params, micro, state, seg.start, seg.end, rIn);
    const solution = mrt(ode);

    // Extract arrays (solution[0] = t, solution[1] = A_c, etc.)
    const tSeg = solution[0];
    const acSeg = solution[1];
    const apSeg = solution[2];
    const adSeg = is3comp ? solution[3] : null;
    const agutSeg = isOral ? solution[is3comp ? 4 : 3] : null;

    // Skip first point of subsequent segments to avoid duplication
    const startIdx = allT.length > 0 ? 1 : 0;

    for (let i = startIdx; i < tSeg.length; i++) {
      allT.push(tSeg[i]);
      allAc.push(acSeg[i]);
      allAp.push(apSeg[i]);
      if (is3comp && adSeg) allAd.push(adSeg[i]);
      if (isOral && agutSeg) allAgut.push(agutSeg[i]);
    }

    // Update state for next segment
    const lastIdx = tSeg.length - 1;
    state[0] = acSeg[lastIdx];
    state[1] = apSeg[lastIdx];
    if (is3comp && adSeg) state[2] = adSeg[lastIdx];
    if (isOral && agutSeg) state[is3comp ? 3 : 2] = agutSeg[lastIdx];
  }

  // Convert to typed arrays
  const t = new Float64Array(allT);
  const ac = new Float64Array(allAc);
  const ap = new Float64Array(allAp);
  const ad = is3comp ? new Float64Array(allAd) : null;
  const agut = isOral ? new Float64Array(allAgut) : null;
  const conc = new Float64Array(allAc.length);
  for (let i = 0; i < allAc.length; i++)
    conc[i] = allAc[i] / params.vc;

  // Compute derived metrics
  const {cMax, tMax, auc} = computePkMetrics(t, conc);

  // Steady-state metrics (only for repeated dosing)
  let cssMax: number | null = null;
  let cssMin: number | null = null;
  let t90ss: number | null = null;

  if (params.repeatedDosing && params.nDoses > 1) {
    const lastDoseTime = (params.nDoses - 1) * params.tau;
    const ssMetrics = computeSteadyStateMetrics(t, conc, lastDoseTime, params.tau);
    cssMax = ssMetrics.cssMax;
    cssMin = ssMetrics.cssMin;
    t90ss = computeT90ss(t, conc, params.tau, cssMax);
  }

  return {t, ac, ap, ad, agut, conc, cMax, tMax, auc, cssMax, cssMin, t90ss};
}

/** Compute C_max, t_max, and AUC from concentration-time arrays */
function computePkMetrics(t: Float64Array, conc: Float64Array): {cMax: number; tMax: number; auc: number} {
  let cMax = 0;
  let tMax = 0;
  let auc = 0;

  for (let i = 0; i < conc.length; i++) {
    if (conc[i] > cMax) {
      cMax = conc[i];
      tMax = t[i];
    }
    // Trapezoidal integration
    if (i > 0)
      auc += 0.5 * (conc[i - 1] + conc[i]) * (t[i] - t[i - 1]);
  }

  return {cMax, tMax, auc};
}

/** Compute steady-state metrics from the last dosing interval */
function computeSteadyStateMetrics(
  t: Float64Array, conc: Float64Array, lastDoseTime: number, tau: number,
): {cssMax: number; cssMin: number} {
  let cssMax = -Infinity;
  let cssMin = Infinity;
  const intervalEnd = lastDoseTime + tau;

  for (let i = 0; i < t.length; i++) {
    if (t[i] >= lastDoseTime && t[i] <= intervalEnd) {
      if (conc[i] > cssMax) cssMax = conc[i];
      if (conc[i] < cssMin) cssMin = conc[i];
    }
  }

  if (cssMax === -Infinity) cssMax = 0;
  if (cssMin === Infinity) cssMin = 0;

  return {cssMax, cssMin};
}

/** Compute time to reach 90% of steady-state maximum */
function computeT90ss(
  t: Float64Array, conc: Float64Array, tau: number, cssMax: number | null,
): number | null {
  if (cssMax === null || cssMax === 0) return null;

  const threshold = 0.9 * cssMax;
  for (let i = 0; i < t.length; i++) {
    // Look at the peak of each interval
    if (i > 0 && conc[i] < conc[i - 1] && conc[i - 1] >= threshold)
      return t[i - 1];
  }

  // Check last point
  if (conc.length > 0 && conc[conc.length - 1] >= threshold)
    return t[t.length - 1];

  return null;
}

// --- Data generation ---

/** Generate synthetic noisy observations from PK solution */
export function generateObservedData(params: DataGenParams): ObservedData {
  const {primarySolution, cv, schedule, tEnd, repeatedDosing, tau, nDoses} = params;
  const baseTimes = SAMPLING_SCHEDULES[schedule];

  // Build sampling times
  let samplingTimes: number[];
  if (repeatedDosing && nDoses > 1) {
    // Extend to cover all dosing intervals, relative to last dose
    samplingTimes = [];
    const lastDoseTime = (nDoses - 1) * tau;
    for (const bt of baseTimes) {
      const t = lastDoseTime + bt;
      if (t <= tEnd) samplingTimes.push(t);
    }
    // Also include some early time points
    for (const bt of baseTimes) {
      if (bt <= tEnd && !samplingTimes.includes(bt))
        samplingTimes.push(bt);
    }
    samplingTimes.sort((a, b) => a - b);
  } else {
    samplingTimes = baseTimes.filter((t) => t <= tEnd);
  }

  // Interpolate primary solution at sampling times
  const cTrue: number[] = [];
  const solT = primarySolution.t;
  const solC = primarySolution.conc;

  for (const ts of samplingTimes)
    cTrue.push(interpolate(solT, solC, ts));

  // Add log-normal noise
  const sigma = Math.sqrt(Math.log(1 + (cv / 100) ** 2));
  const cObs: number[] = [];

  for (const ct of cTrue) {
    if (ct <= 0) {
      cObs.push(0);
      continue;
    }
    const z = boxMullerRandom();
    const noisy = ct * Math.exp(z * sigma - sigma * sigma / 2);
    cObs.push(Math.max(0, noisy));
  }

  return {tObs: samplingTimes, cObs, cTrue};
}

/** Linear interpolation between nearest solution points */
function interpolate(t: Float64Array, y: Float64Array, target: number): number {
  if (target <= t[0]) return y[0];
  if (target >= t[t.length - 1]) return y[t.length - 1];

  // Binary search for bracket
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

/** Box-Muller transform for standard normal random variable */
function boxMullerRandom(): number {
  let u1 = 0;
  let u2 = 0;
  while (u1 === 0) u1 = Math.random();
  while (u2 === 0) u2 = Math.random();
  return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
}

// --- WSSR computation ---

/** Compute weighted sum of squared residuals */
export function computeWssr(
  cObs: number[], cPred: number[], weightType: WeightType,
): number {
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

/** Compute R-squared */
export function computeR2(
  cObs: number[], wssr: number, weightType: WeightType,
): number {
  const mean = cObs.reduce((s, v) => s + v, 0) / cObs.length;
  let ssTot = 0;
  for (let i = 0; i < cObs.length; i++) {
    let w = 1;
    if (weightType === '1/C_obs\u00B2' && cObs[i] > 0)
      w = 1 / (cObs[i] * cObs[i]);
    ssTot += w * (cObs[i] - mean) * (cObs[i] - mean);
  }
  return ssTot > 0 ? 1 - wssr / ssTot : 0;
}

/** Compute AIC */
export function computeAic(n: number, p: number, wssr: number): number {
  return n * Math.log(wssr / n) + 2 * p;
}

// --- Utility ---

/** Generate linearly spaced array */
export function linspace(min: number, max: number, n: number): number[] {
  const result: number[] = [];
  for (let i = 0; i < n; i++)
    result.push(min + i * (max - min) / (n - 1));
  return result;
}

/** Solve PK for a given parameter set and return predicted concentrations at observation times */
export function solvePkAtTimes(
  params: PkParams, tObs: number[],
): number[] {
  const solution = solvePk(params);
  const cPred: number[] = [];
  for (const t of tObs)
    cPred.push(interpolate(solution.t, solution.conc, t));
  return cPred;
}
