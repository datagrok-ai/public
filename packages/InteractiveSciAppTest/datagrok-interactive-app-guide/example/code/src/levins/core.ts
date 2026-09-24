// Levins Metapopulation Model — Computational Core

import {mrt} from 'diff-grok';

import {LevinsParams, createLevinsODE, getEquilibrium} from './model';

// --- Re-exports ---

export type {LevinsParams} from './model';
export {createLevinsODE, getEquilibrium} from './model';

// --- Types ---

export interface LevinsSolution {
  t: Float64Array;
  p: Float64Array;
  p_star: number;
}

export type InputId = 'ctrl_p0' | 'ctrl_m' | 'ctrl_e0' | 'ctrl_rescue' |
  'ctrl_t_start' | 'ctrl_t_end' | 'ctrl_t_step' | 'ctrl_tolerance';

export type ValidationErrors = Map<InputId, string>;

// --- Defaults ---

export const DEFAULTS: LevinsParams = {
  p0: 0.5,
  m: 0.5,
  e0: 0.2,
  rescueEffect: false,
  t_start: 0,
  t_end: 50,
  t_step: 0.1,
  tolerance: 1e-7,
};

// --- Validation ---

export function validate(inputs: LevinsParams): ValidationErrors {
  const errors: ValidationErrors = new Map();
  const {p0, m, e0, rescueEffect, t_start, t_end, t_step, tolerance} = inputs;

  // val_01, val_02
  if (p0 <= 0)
    errors.set('ctrl_p0', 'Initial patch fraction must be greater than 0');
  else if (p0 > 1)
    errors.set('ctrl_p0', 'Initial patch fraction cannot exceed 1');

  // val_03
  if (m <= 0)
    errors.set('ctrl_m', 'Colonization rate must be positive');

  // val_04
  if (e0 <= 0)
    errors.set('ctrl_e0', 'Extinction rate must be positive');

  // val_05 — only if val_03 and val_04 passed
  if (!errors.has('ctrl_m') && !errors.has('ctrl_e0') && !rescueEffect && m <= e0)
    errors.set('ctrl_m', 'Colonization rate must exceed extinction rate (m > e₀). With current values the metapopulation tends to extinction');

  // val_06
  if (t_end <= t_start) {
    errors.set('ctrl_t_end', 'End of interval must be greater than start');
    errors.set('ctrl_t_start', 'End of interval must be greater than start');
  }

  // val_07
  if (t_step <= 0)
    errors.set('ctrl_t_step', 'Step must be positive');

  // val_08 — only if val_06 and val_07 passed
  if (!errors.has('ctrl_t_end') && !errors.has('ctrl_t_step') && t_step >= t_end - t_start)
    errors.set('ctrl_t_step', 'Step must be less than interval length');

  // val_09
  if (tolerance <= 0)
    errors.set('ctrl_tolerance', 'Tolerance must be positive');

  return errors;
}

// --- Solver ---

export function solve(inputs: LevinsParams): LevinsSolution {
  const task = createLevinsODE(inputs);
  const solution = mrt(task);

  return {
    t: solution[0],
    p: solution[1],
    p_star: getEquilibrium(inputs.m, inputs.e0, inputs.rescueEffect),
  };
}

// --- Optimization validation ---

export interface OptimizeInputs {
  m_min: number;
  m_max: number;
}

export type OptInputId = 'dlg_m_min' | 'dlg_m_max';
export type OptValidationErrors = Map<OptInputId, string>;

export function validateOptimize(
  opt: OptimizeInputs, e0: number, rescueEffect: boolean,
): {errors: OptValidationErrors; warning: string | null} {
  const errors: OptValidationErrors = new Map();
  let warning: string | null = null;

  if (opt.m_min <= 0)
    errors.set('dlg_m_min', 'Colonization rate must be positive');

  if (opt.m_max <= 0)
    errors.set('dlg_m_max', 'Colonization rate must be positive');

  if (!errors.has('dlg_m_min') && !errors.has('dlg_m_max') && opt.m_min >= opt.m_max)
    errors.set('dlg_m_min', 'Minimum value must be less than maximum');

  if (errors.size === 0 && !rescueEffect && opt.m_max <= e0)
    warning = 'With current e₀ the entire m range leads to extinction (m ≤ e₀). Increase the maximum or decrease e₀';

  return {errors, warning};
}

// --- Worker message types ---

export interface WorkerTask {
  m_i: number;
  p0: number;
  e0: number;
  rescueEffect: boolean;
  t_start: number;
  t_end: number;
  t_step: number;
  tolerance: number;
}

export interface WorkerResult {
  m_i: number;
  p_end: number;
  error?: string;
}
