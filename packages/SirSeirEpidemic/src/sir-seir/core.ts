// SIR/SEIR Epidemic Model — Computational Core

import {mrt} from 'diff-grok';

import {
  EpidemicParams, EpidemicSolution,
  createSIR_ODE, createSEIR_ODE, POPULATION,
} from './model';

// --- Re-exports ---

export type {EpidemicParams, EpidemicSolution} from './model';
export {createSIR_ODE, createSEIR_ODE, POPULATION, INITIAL_INFECTED} from './model';

// --- Types ---

export type InputId = 'ctrl_r0' | 'ctrl_gamma' | 'ctrl_sigma' | 'ctrl_vaccination';

export type ValidationErrors = Map<InputId, string>;

// --- Defaults ---

export const DEFAULTS: EpidemicParams = {
  modelType: 'SIR',
  r0: 3.0,
  gamma: 0.1,
  sigma: 0.2,
  vaccination: 0,
};

// --- Slider ranges (UI concern, but used by validators for range messages) ---

export const RANGES: Record<string, {min: number; max: number}> = {
  r0: {min: 0.5, max: 8.0},
  gamma: {min: 0.01, max: 0.5},
  sigma: {min: 0.05, max: 1.0},
  vaccination: {min: 0, max: 90},
};

// --- Validation ---

/** Validates epidemic parameters, returns map of input ID → error message */
export function validate(inputs: EpidemicParams): ValidationErrors {
  const errors: ValidationErrors = new Map();
  const {modelType, r0, gamma, sigma, vaccination} = inputs;

  // val_01: R₀ range
  if (r0 < 0.5 || r0 > 8.0)
    errors.set('ctrl_r0', 'R\u2080 must be between 0.5 and 8.0');

  // val_02: γ range
  if (gamma <= 0 || gamma > 0.5)
    errors.set('ctrl_gamma', 'Recovery rate must be in (0, 0.5]');

  // val_03: σ range (SEIR only)
  if (modelType === 'SEIR' && (sigma <= 0 || sigma > 1.0))
    errors.set('ctrl_sigma', 'Incubation rate must be in (0, 1.0]');

  // val_04: vaccination range
  if (vaccination < 0 || vaccination > 90)
    errors.set('ctrl_vaccination', 'Vaccination coverage must be between 0% and 90%');

  // val_05: vaccination feasibility (only if val_04 passed)
  if (!errors.has('ctrl_vaccination')) {
    const v = vaccination / 100;
    if (v >= 1 - 1 / POPULATION)
      errors.set('ctrl_vaccination', 'Vaccination coverage too high \u2014 no susceptible individuals remain');
  }

  return errors;
}

// --- Solver ---

/** Solves the SIR or SEIR ODE system and returns the full solution */
export function solve(inputs: EpidemicParams): EpidemicSolution {
  const isSIR = inputs.modelType === 'SIR';
  const task = isSIR ? createSIR_ODE(inputs) : createSEIR_ODE(inputs);
  const solution = mrt(task);

  const t = solution[0];
  let S: Float64Array;
  let E: Float64Array | null;
  let I: Float64Array;
  let R: Float64Array;

  if (isSIR) {
    S = solution[1];
    E = null;
    I = solution[2];
    R = solution[3];
  } else {
    S = solution[1];
    E = solution[2];
    I = solution[3];
    R = solution[4];
  }

  // Compute R_eff(t) = R₀ · S(t) / N
  const rEff = new Float64Array(t.length);
  for (let i = 0; i < t.length; i++)
    rEff[i] = inputs.r0 * S[i] / POPULATION;

  // Find peak infection
  let peakIdx = 0;
  for (let i = 1; i < I.length; i++) {
    if (I[i] > I[peakIdx])
      peakIdx = i;
  }

  const peakDay = Math.round(t[peakIdx]);
  const peakCount = Math.round(I[peakIdx]);

  // Summary statistics
  const finalRecoveredPct = R[R.length - 1] / POPULATION * 100;
  const herdImmunityThreshold = (1 - 1 / inputs.r0) * 100;

  return {
    t, S, E, I, R, rEff,
    peakDay,
    peakCount,
    finalRecoveredPct,
    herdImmunityThreshold,
  };
}
