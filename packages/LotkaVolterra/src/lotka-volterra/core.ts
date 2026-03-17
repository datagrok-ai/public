// Lotka-Volterra Predator-Prey Model — Computational Core

import {mrt} from 'diff-grok';

import {LotkaVolterraParams, createLotkaVolterraODE, getEquilibrium} from './model';

// --- Re-exports ---

export type {LotkaVolterraParams} from './model';
export {createLotkaVolterraODE, getEquilibrium} from './model';

// --- Types ---

export interface LotkaVolterraSolution {
  t: Float64Array;
  x: Float64Array;
  y: Float64Array;
  xStar: number;
  yStar: number;
  maxPrey: number;
  maxPredators: number;
  stepCount: number;
}

export type InputId = 'ctrl_alpha' | 'ctrl_beta' | 'ctrl_delta' | 'ctrl_gamma' |
  'ctrl_x0' | 'ctrl_y0' | 'ctrl_T';

export type ValidationErrors = Map<InputId, string>;

// --- Defaults ---

export const DEFAULTS: LotkaVolterraParams = {
  alpha: 1.0,
  beta: 0.1,
  delta: 0.075,
  gamma: 1.5,
  x0: 10,
  y0: 5,
  T: 100,
};

// --- Slider ranges ---

export const RANGES: Record<string, {min: number; max: number}> = {
  alpha: {min: 0.1, max: 3.0},
  beta: {min: 0.01, max: 0.5},
  delta: {min: 0.01, max: 0.5},
  gamma: {min: 0.1, max: 3.0},
  x0: {min: 1, max: 200},
  y0: {min: 1, max: 100},
  T: {min: 10, max: 500},
};

// --- Validation ---

export function validate(inputs: LotkaVolterraParams): ValidationErrors {
  const errors: ValidationErrors = new Map();
  const {alpha, beta, delta, gamma, x0, y0, T} = inputs;

  if (alpha <= 0)
    errors.set('ctrl_alpha', 'Prey birth rate must be positive');

  if (beta <= 0)
    errors.set('ctrl_beta', 'Predation rate must be positive');

  if (delta <= 0)
    errors.set('ctrl_delta', 'Predator efficiency must be positive');

  if (gamma <= 0)
    errors.set('ctrl_gamma', 'Predator death rate must be positive');

  if (x0 <= 0)
    errors.set('ctrl_x0', 'Initial prey population must be positive');

  if (y0 <= 0)
    errors.set('ctrl_y0', 'Initial predator population must be positive');

  if (T <= 0)
    errors.set('ctrl_T', 'Simulation time must be positive');

  return errors;
}

// --- Solver ---

export function solve(inputs: LotkaVolterraParams): LotkaVolterraSolution {
  const task = createLotkaVolterraODE(inputs);
  const solution = mrt(task);

  const t = solution[0];
  const x = solution[1];
  const y = solution[2];

  const eq = getEquilibrium(inputs.alpha, inputs.beta, inputs.delta, inputs.gamma);

  let maxPrey = 0;
  let maxPredators = 0;
  for (let i = 0; i < x.length; i++) {
    if (x[i] > maxPrey) maxPrey = x[i];
    if (y[i] > maxPredators) maxPredators = y[i];
  }

  return {
    t, x, y,
    xStar: eq.xStar,
    yStar: eq.yStar,
    maxPrey,
    maxPredators,
    stepCount: t.length,
  };
}

// --- Worker message types ---

export interface WorkerTask {
  alpha: number;
  beta: number;
  delta: number;
  gamma: number;
  x0: number;
  y0: number;
  T: number;
}

export interface WorkerResult {
  alpha: number;
  beta: number;
  delta: number;
  gamma: number;
  maxPrey: number;
  error?: string;
}

// --- Equilibrium exploration types ---

export type ParamName = 'alpha' | 'beta' | 'delta' | 'gamma' | 'x0' | 'y0' | 'T';

export const PARAM_NAMES: ParamName[] = ['alpha', 'beta', 'delta', 'gamma', 'x0', 'y0', 'T'];

export interface EquilibriumTask {
  paramName: ParamName;
  paramMin: number;
  paramMax: number;
  steps: number;
  baseParams: LotkaVolterraParams;
}

export interface EquilibriumResult {
  paramValues: number[];
  xStar: number[];
  yStar: number[];
  error?: string;
}
