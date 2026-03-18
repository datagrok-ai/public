// Multi-Compartment Pharmacokinetic Model — ODE specification

import {ODEs} from 'diff-grok';

// --- Types ---

export type ModelType = '2-Compartment' | '3-Compartment';
export type InputMode = 'IV Bolus' | 'IV Infusion' | 'Oral';
export type SamplingSchedule = 'Standard PK' | 'Rich' | 'Sparse';
export type WeightType = 'Uniform' | '1/C_obs\u00B2' | '1/C_pred\u00B2';

/** All parameters needed for a PK ODE solve */
export interface PkParams {
  modelType: ModelType;
  inputMode: InputMode;
  cl: number;
  vc: number;
  vp: number;
  q: number;
  vd: number;
  q2: number;
  ka: number;
  f: number;
  dose: number;
  repeatedDosing: boolean;
  tau: number;
  nDoses: number;
  tInf: number;
  tEnd: number;
  mec: number;
  mtc: number;
}

/** Solution from the primary PK ODE solve */
export interface PkSolution {
  t: Float64Array;
  ac: Float64Array;
  ap: Float64Array;
  ad: Float64Array | null;
  agut: Float64Array | null;
  conc: Float64Array;
  cMax: number;
  tMax: number;
  auc: number;
  cssMax: number | null;
  cssMin: number | null;
  t90ss: number | null;
}

/** Observed data from synthetic data generation */
export interface ObservedData {
  tObs: number[];
  cObs: number[];
  cTrue: number[];
}

/** Parameters for data generation */
export interface DataGenParams {
  primarySolution: PkSolution;
  cv: number;
  schedule: SamplingSchedule;
  tEnd: number;
  repeatedDosing: boolean;
  tau: number;
  nDoses: number;
}

/** Selected parameter for optimization */
export interface SelectedParam {
  name: string;
  min: number;
  max: number;
}

/** Parameters for optimization */
export interface OptimizeParams {
  observedData: ObservedData;
  selectedParams: SelectedParam[];
  fixedParams: PkParams;
  gridResolution: number;
  weightType: WeightType;
}

/** Result of a single grid point evaluation */
export interface GridPointResult {
  paramValues: number[];
  wssr: number;
}

/** Full optimization result */
export interface OptimizeResult {
  bestParams: Map<string, number>;
  bestWssr: number;
  bestR2: number;
  bestAic: number;
  bestConc: {t: number[]; c: number[]};
  landscape: LandscapeData;
  residuals: {t: number[]; res: number[]};
}

/** WSSR landscape data for visualization */
export interface LandscapeData {
  paramNames: string[];
  paramValues: number[][];
  wssrValues: number[];
  gridResolution: number;
}

/** Message sent to optimization worker */
export interface OptimizeWorkerInput {
  gridPoints: number[][];
  paramNames: string[];
  fixedParams: PkParams;
  tEnd: number;
  observedData: {tObs: number[]; cObs: number[]};
  weightType: WeightType;
}

/** Result from optimization worker */
export interface OptimizeWorkerOutput {
  success: boolean;
  results?: GridPointResult[];
  error?: string;
}

// --- Micro-constants ---

/** Compute micro-rate constants from physiological PK parameters */
export function computeMicroConstants(params: PkParams): {
  k10: number; k12: number; k21: number;
  k13: number; k31: number;
} {
  const k10 = params.cl / params.vc;
  const k12 = params.q / params.vc;
  const k21 = params.q / params.vp;
  const k13 = params.modelType === '3-Compartment' ? params.q2 / params.vc : 0;
  const k31 = params.modelType === '3-Compartment' ? params.q2 / params.vd : 0;
  return {k10, k12, k21, k13, k31};
}

// --- Segment schedule ---

/** Build integration segment schedule for multiple dosing */
export function buildSegmentSchedule(params: PkParams): {start: number; end: number; doseEvent: boolean; infusionEnd: boolean}[] {
  const segments: {start: number; end: number; doseEvent: boolean; infusionEnd: boolean}[] = [];
  const nDoses = params.repeatedDosing ? params.nDoses : 1;
  const tau = params.tau;
  const tEnd = params.tEnd;

  // Collect all event times
  const events: {time: number; doseEvent: boolean; infusionEnd: boolean}[] = [];

  for (let i = 0; i < nDoses; i++) {
    const doseTime = i * tau;
    if (doseTime < tEnd)
      events.push({time: doseTime, doseEvent: true, infusionEnd: false});

    if (params.inputMode === 'IV Infusion') {
      const infEnd = doseTime + params.tInf;
      if (infEnd < tEnd && infEnd > doseTime)
        events.push({time: infEnd, doseEvent: false, infusionEnd: true});
    }
  }

  // Sort by time
  events.sort((a, b) => a.time - b.time);

  // Remove duplicates
  const uniqueEvents: typeof events = [];
  for (const e of events) {
    if (uniqueEvents.length === 0 || Math.abs(e.time - uniqueEvents[uniqueEvents.length - 1].time) > 1e-10)
      uniqueEvents.push(e);
    else {
      // Merge flags
      const last = uniqueEvents[uniqueEvents.length - 1];
      last.doseEvent = last.doseEvent || e.doseEvent;
      last.infusionEnd = last.infusionEnd || e.infusionEnd;
    }
  }

  // Build segments
  for (let i = 0; i < uniqueEvents.length; i++) {
    const start = uniqueEvents[i].time;
    const end = i < uniqueEvents.length - 1 ? uniqueEvents[i + 1].time : tEnd;
    if (end > start) {
      segments.push({
        start,
        end,
        doseEvent: uniqueEvents[i].doseEvent,
        infusionEnd: uniqueEvents[i].infusionEnd,
      });
    }
  }

  return segments;
}

// --- ODE creation ---

/** Create the ODE specification for a single integration segment */
export function createPkODE(
  params: PkParams,
  micro: {k10: number; k12: number; k21: number; k13: number; k31: number},
  initial: number[],
  tStart: number,
  tEnd: number,
  rIn: number,
): ODEs {
  const is3comp = params.modelType === '3-Compartment';
  const isOral = params.inputMode === 'Oral';
  const nVars = (is3comp ? 3 : 2) + (isOral ? 1 : 0);

  return {
    name: 'PK-ODE',
    arg: {name: 't', start: tStart, finish: tEnd, step: 0.05},
    initial: initial,
    func: (_t: number, y: Float64Array, out: Float64Array) => {
      const ac = y[0];
      const ap = y[1];

      let inputRate = rIn;
      if (isOral) {
        const agut = y[is3comp ? 3 : 2];
        inputRate = params.ka * agut;
        out[is3comp ? 3 : 2] = -params.ka * agut; // dA_gut/dt
      }

      // dA_c/dt
      out[0] = -micro.k10 * ac - micro.k12 * ac + micro.k21 * ap + inputRate;
      if (is3comp) {
        const ad = y[2];
        out[0] += -micro.k13 * ac + micro.k31 * ad;
        out[2] = micro.k13 * ac - micro.k31 * ad; // dA_d/dt
      }

      // dA_p/dt
      out[1] = micro.k12 * ac - micro.k21 * ap;
    },
    tolerance: 1e-6,
    solutionColNames: buildSolutionColNames(is3comp, isOral),
  };
}

/** Build solution column names based on model configuration */
function buildSolutionColNames(is3comp: boolean, isOral: boolean): string[] {
  const names = ['A_c', 'A_p'];
  if (is3comp) names.push('A_d');
  if (isOral) names.push('A_gut');
  return names;
}

// --- Sampling schedules ---

export const SAMPLING_SCHEDULES: Record<SamplingSchedule, number[]> = {
  'Standard PK': [0.5, 1, 2, 4, 8, 12, 24],
  'Rich': [0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 10, 12, 16, 20, 24],
  'Sparse': [1, 4, 12, 24],
};
