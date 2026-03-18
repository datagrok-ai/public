// Golf Ball Flight Simulator — Computational Core

import {mrt} from 'diff-grok';

import {GolfBallParams, createGolfBallODE} from './model';

// --- Re-exports ---

export type {GolfBallParams} from './model';
export {createGolfBallODE} from './model';

// --- Types ---

export interface GolfBallSolution {
  t: Float64Array;
  x: Float64Array;
  y: Float64Array;
  v: Float64Array;
  maxHeight: number;
  maxHeightTime: number;
  flightDistance: number;
  flightTime: number;
  landingSpeed: number;
  landingAngle: number;
}

export type InputId = 'ctrl_v0' | 'ctrl_theta' | 'ctrl_m' | 'ctrl_d' |
  'ctrl_Cd' | 'ctrl_rho' | 'ctrl_g';

export type ValidationErrors = Map<InputId, string>;

// --- Defaults ---

export const DEFAULTS: GolfBallParams = {
  v0: 70,
  theta: 45,
  m: 0.0459,
  d: 0.0427,
  Cd: 0.25,
  rho: 1.225,
  g: 9.81,
};

// --- Slider ranges ---

export const RANGES: Record<string, {min: number; max: number}> = {
  v0: {min: 10, max: 120},
  theta: {min: 1, max: 89},
  m: {min: 0.01, max: 0.1},
  d: {min: 0.02, max: 0.08},
  Cd: {min: 0.05, max: 0.6},
  rho: {min: 0.5, max: 1.5},
  g: {min: 1.0, max: 25.0},
};

// --- Validation ---

export function validate(inputs: GolfBallParams): ValidationErrors {
  const errors: ValidationErrors = new Map();
  const {v0, theta, m, d, Cd, rho, g} = inputs;

  if (v0 <= 0)
    errors.set('ctrl_v0', 'Initial speed must be positive');

  if (theta <= 0 || theta >= 90)
    errors.set('ctrl_theta', 'Launch angle must be between 0\u00B0 and 90\u00B0 (exclusive)');

  if (m <= 0)
    errors.set('ctrl_m', 'Ball mass must be positive');

  if (d <= 0)
    errors.set('ctrl_d', 'Ball diameter must be positive');

  if (Cd <= 0)
    errors.set('ctrl_Cd', 'Drag coefficient must be positive');

  if (rho <= 0)
    errors.set('ctrl_rho', 'Air density must be positive');

  if (g <= 0)
    errors.set('ctrl_g', 'Gravitational acceleration must be positive');

  return errors;
}

// --- Solver ---

export function solve(inputs: GolfBallParams): GolfBallSolution {
  const task = createGolfBallODE(inputs);
  const solution = mrt(task);

  const tRaw = solution[0];
  const xRaw = solution[1];
  const yRaw = solution[2];
  const vxRaw = solution[3];
  const vyRaw = solution[4];

  // Find landing index: first index where y < 0 after t > 0
  let landIdx = tRaw.length;
  for (let i = 1; i < tRaw.length; i++) {
    if (yRaw[i] <= 0) {
      landIdx = i;
      break;
    }
  }

  // Interpolate exact landing point if we found a ground crossing
  let landT: number;
  let landX: number;
  let landVx: number;
  let landVy: number;

  if (landIdx < tRaw.length && landIdx > 0) {
    const y0 = yRaw[landIdx - 1];
    const y1 = yRaw[landIdx];
    const frac = y0 / (y0 - y1);

    landT = tRaw[landIdx - 1] + frac * (tRaw[landIdx] - tRaw[landIdx - 1]);
    landX = xRaw[landIdx - 1] + frac * (xRaw[landIdx] - xRaw[landIdx - 1]);
    landVx = vxRaw[landIdx - 1] + frac * (vxRaw[landIdx] - vxRaw[landIdx - 1]);
    landVy = vyRaw[landIdx - 1] + frac * (vyRaw[landIdx] - vyRaw[landIdx - 1]);
  } else {
    // Ball never came back down (shouldn't happen with reasonable params)
    landIdx = tRaw.length;
    landT = tRaw[tRaw.length - 1];
    landX = xRaw[xRaw.length - 1];
    landVx = vxRaw[vxRaw.length - 1];
    landVy = vyRaw[vyRaw.length - 1];
  }

  // Build truncated arrays up to landing point
  const len = landIdx; // use points before ground crossing
  const t = new Float64Array(len + 1);
  const x = new Float64Array(len + 1);
  const y = new Float64Array(len + 1);
  const v = new Float64Array(len + 1);

  let maxHeight = 0;
  let maxHeightTime = 0;

  for (let i = 0; i < len; i++) {
    t[i] = tRaw[i];
    x[i] = xRaw[i];
    y[i] = yRaw[i];
    v[i] = Math.sqrt(vxRaw[i] * vxRaw[i] + vyRaw[i] * vyRaw[i]);

    if (yRaw[i] > maxHeight) {
      maxHeight = yRaw[i];
      maxHeightTime = tRaw[i];
    }
  }

  // Add interpolated landing point
  t[len] = landT;
  x[len] = landX;
  y[len] = 0;
  v[len] = Math.sqrt(landVx * landVx + landVy * landVy);

  const landingSpeed = v[len];
  const landingAngle = Math.abs(Math.atan2(landVy, landVx) * 180 / Math.PI);

  return {
    t, x, y, v,
    maxHeight,
    maxHeightTime,
    flightDistance: landX,
    flightTime: landT,
    landingSpeed,
    landingAngle,
  };
}

// --- Worker message types for optimization ---

export type OptimizeParam = 'theta' | 'v0' | 'Cd';

export const OPTIMIZE_PARAMS: OptimizeParam[] = ['theta', 'v0', 'Cd'];

export interface OptimizeTask {
  paramName: OptimizeParam;
  values: number[];
  baseParams: GolfBallParams;
}

export interface OptimizeResult {
  value: number;
  maxHeight: number;
  error?: string;
}

// --- Worker message types for sensitivity analysis ---

export type SweepParamName = 'v0' | 'theta' | 'Cd' | 'rho' | 'm' | 'd';

export const SWEEP_PARAMS: SweepParamName[] = ['v0', 'theta', 'Cd', 'rho', 'm', 'd'];

export interface SensitivityTask {
  paramName: SweepParamName;
  paramMin: number;
  paramMax: number;
  steps: number;
  baseParams: GolfBallParams;
}

export interface SensitivityResult {
  paramName: string;
  paramValues: number[];
  distances: number[];
  error?: string;
}

// --- Utility ---

export function linspace(min: number, max: number, steps: number): number[] {
  const arr: number[] = [];
  for (let i = 0; i < steps; i++)
    arr.push(min + i * (max - min) / (steps - 1));
  return arr;
}

/** Solve and return only maxHeight (used in optimization workers) */
export function solveMaxHeight(params: GolfBallParams): number {
  const task = createGolfBallODE(params);
  const solution = mrt(task);
  const yArr = solution[2];
  let maxH = 0;
  for (let i = 0; i < yArr.length; i++) {
    if (yArr[i] > maxH) maxH = yArr[i];
  }
  return maxH;
}

/** Solve and return flight distance (used in sensitivity workers) */
export function solveFlightDistance(params: GolfBallParams): number {
  const task = createGolfBallODE(params);
  const solution = mrt(task);
  const xArr = solution[1];
  const yArr = solution[2];

  for (let i = 1; i < yArr.length; i++) {
    if (yArr[i] <= 0) {
      // Interpolate landing x
      const y0 = yArr[i - 1];
      const y1 = yArr[i];
      const frac = y0 / (y0 - y1);
      return xArr[i - 1] + frac * (xArr[i] - xArr[i - 1]);
    }
  }
  return xArr[xArr.length - 1];
}
