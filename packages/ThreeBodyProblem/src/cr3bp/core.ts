// CR3BP — Computational Core

import {mrt} from 'diff-grok';

import {
  CR3BPParams, LagrangePoint, ZVCGrid,
  createCR3BPODE, jacobiConstant, computeLagrangePoints, computeZVCGrid,
} from './model';

// --- Re-exports ---

export type {CR3BPParams, LagrangePoint, ZVCGrid} from './model';
export {createCR3BPODE, jacobiConstant, effectivePotential, computeLagrangePoints, evalRHS} from './model';

// --- Types ---

export interface CR3BPSolution {
  t: Float64Array;
  x: Float64Array;
  y: Float64Array;
  vx: Float64Array;
  vy: Float64Array;
  vmag: Float64Array;
  cj: Float64Array;
  lagrangePoints: LagrangePoint[];
  zvcGrid: ZVCGrid;
}

export type InputId = 'ctrl_mu' | 'ctrl_x0' | 'ctrl_y0' | 'ctrl_vx0' | 'ctrl_vy0' | 'ctrl_T';

export type ValidationErrors = Map<InputId, string>;

// --- Defaults (Free-return trajectory for Earth-Moon) ---

export const DEFAULTS: CR3BPParams = {
  mu: 0.01215,
  x0: 0.994,
  y0: 0.0,
  vx0: 0.0,
  vy0: -2.0016,
  T: 20.0,
};

// --- Slider ranges ---

export const RANGES: Record<string, {min: number; max: number}> = {
  mu: {min: 0.001, max: 0.5},
  x0: {min: -1.5, max: 1.5},
  y0: {min: -1.5, max: 1.5},
  vx0: {min: -2.0, max: 2.0},
  vy0: {min: -2.0, max: 2.0},
  T: {min: 1.0, max: 100.0},
};

// --- Validation ---

export function validate(inputs: CR3BPParams): ValidationErrors {
  const errors: ValidationErrors = new Map();
  const {mu, T} = inputs;

  if (mu <= 0)
    errors.set('ctrl_mu', 'Mass parameter \u03BC must be positive');
  else if (mu > 0.5)
    errors.set('ctrl_mu', 'Mass parameter \u03BC must not exceed 0.5 (by convention m\u2082 \u2264 m\u2081)');

  if (T <= 0)
    errors.set('ctrl_T', 'Integration time must be positive');

  return errors;
}

// --- Solver ---

export function solve(inputs: CR3BPParams): CR3BPSolution {
  const {mu} = inputs;
  const task = createCR3BPODE(inputs);
  const solution = mrt(task);

  const t = solution[0];
  const x = solution[1];
  const y = solution[2];
  const vx = solution[3];
  const vy = solution[4];

  // Velocity magnitude
  const vmag = new Float64Array(t.length);
  for (let i = 0; i < t.length; i++)
    vmag[i] = Math.sqrt(vx[i] ** 2 + vy[i] ** 2);

  // Jacobi constant
  const cj = new Float64Array(t.length);
  for (let i = 0; i < t.length; i++)
    cj[i] = jacobiConstant(x[i], y[i], vx[i], vy[i], mu);

  // Lagrange points
  const lagrangePoints = computeLagrangePoints(mu);

  // ZVC grid based on initial Jacobi constant
  const cj0 = cj[0];
  const zvcGrid = computeZVCGrid(mu, cj0);

  return {t, x, y, vx, vy, vmag, cj, lagrangePoints, zvcGrid};
}
