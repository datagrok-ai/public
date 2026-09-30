// sci-comp single-objective optimizers adapted to IOptimizer, with per-method settings.
//
// Worker-safe on purpose: the Diff Studio fitting worker runs these optimizers,
// so only type imports from `optimizer-misc.ts` (it has a value-level DG import).

import {singleObjective} from '@datagrok-libraries/sci-comp';
import type {IOptimizer, OptimizerBounds, Setting} from './optimizer-misc';

/** Translates the flat settings of a method into the typed settings of its sci-comp optimizer */
type SettingsBuilder<S> = (settings: Map<string, number>, bounds: OptimizerBounds | undefined, x0: Float64Array) => S;

export function adaptOptimizer<S extends singleObjective.CommonSettings>(
  makeOpt: () => singleObjective.Optimizer<S>,
  buildSettings: SettingsBuilder<S>,
): IOptimizer {
  return async (objectiveFunc, x0, settings, threshold, bounds, isCanceled) => {
    // Out-of-bounds points (undefined cost) are penalized as in Nelder-Mead
    const costOutside = 2 * ((await objectiveFunc(x0)) ?? Infinity);
    const maxEvals = settings.get('maxFunctionEvaluations') ?? Infinity;
    let evals = 0;
    const fn = async (x: Float64Array) => {
      ++evals;
      return (await objectiveFunc(x)) ?? costOutside;
    };

    const s = buildSettings(settings, bounds, x0);
    s.onIteration = (state) => (evals >= maxEvals) || ((threshold != null) && (state.bestValue <= threshold)) ||
      (isCanceled?.() === true);

    const res = await makeOpt().minimizeAsync(fn, x0, s);

    return {
      point: res.point,
      cost: res.value,
      iterCosts: Array.from(res.costHistory),
      // The history starts with the initial cost, so it may be one longer than `res.iterations`
      iterCount: res.costHistory.length,
    };
  };
}

const MAX_ITERATIONS: Setting = {
  default: 50,
  min: 1,
  max: 10000,
  caption: 'max iterations',
  tooltipText: 'Maximum number of iterations. Higher = better fit, but longer computation',
  inputType: 'Int',
};

const MAX_EVALUATIONS: Setting = {
  default: 500,
  min: 1,
  max: 1000000,
  caption: 'max evaluations',
  tooltipText: 'Maximum number of model runs per initial point. A numerical gradient takes two runs per fitted input',
  inputType: 'Int',
};

const GRAD_TOLERANCE: Setting = {
  default: 1e-5,
  min: 1e-20,
  max: 1e-1,
  caption: 'gradient tolerance',
  tooltipText: 'Stop once the gradient is smaller. Lower value = more accurate result, but longer computation',
  inputType: 'Float',
};

const HISTORY_SIZE: Setting = {
  default: 10,
  min: 1,
  max: 100,
  caption: 'history size',
  tooltipText: 'Number of previous steps used to approximate the loss function curvature',
  inputType: 'Int',
};

// Larger than the sci-comp default (1e-7): model outputs are often float32, and their rounding
// noise would dominate a smaller step
const FINITE_DIFF_STEP: Setting = {
  default: 1e-5,
  min: 1e-12,
  max: 1e-1,
  caption: 'gradient step',
  tooltipText: 'Input increment used to compute the gradient numerically',
  inputType: 'Float',
};

const PENALTY_COEFFICIENT: Setting = {
  default: 1000,
  min: 1e-3,
  max: 1e12,
  caption: 'penalty',
  tooltipText: 'Penalty coefficient for leaving the input ranges. Higher = stricter ranges, but harder search',
  inputType: 'Float',
};

export const lbfgsbSettingsOpts = new Map<string, Setting>([
  ['maxIterations', MAX_ITERATIONS],
  ['maxFunctionEvaluations', MAX_EVALUATIONS],
  ['gradTolerance', GRAD_TOLERANCE],
  ['historySize', HISTORY_SIZE],
  ['finiteDiffStep', FINITE_DIFF_STEP],
]);

export const psoSettingsOpts = new Map<string, Setting>([
  ['maxIterations', MAX_ITERATIONS],
  ['swarmSize', {
    default: 10,
    min: 2,
    max: 1000,
    caption: 'swarm size',
    tooltipText: 'Number of particles. Higher = wider search, but more model runs per iteration',
    inputType: 'Int',
  }],
  ['inertia', {
    default: 0.7298,
    min: 0,
    max: 1,
    caption: 'inertia',
    tooltipText: 'How strongly particles keep their velocity. Higher = more exploration',
    inputType: 'Float',
  }],
  ['cognitive', {
    default: 1.4962,
    min: 0,
    max: 4,
    caption: 'cognitive',
    tooltipText: 'Pull of each particle towards its own best point',
    inputType: 'Float',
  }],
  ['social', {
    default: 1.4962,
    min: 0,
    max: 4,
    caption: 'social',
    tooltipText: 'Pull of each particle towards the best point of the swarm',
    inputType: 'Float',
  }],
  ['mu', PENALTY_COEFFICIENT],
]);

export const lbfgsSettingsOpts = new Map<string, Setting>([
  ['maxIterations', MAX_ITERATIONS],
  ['maxFunctionEvaluations', MAX_EVALUATIONS],
  ['gradTolerance', GRAD_TOLERANCE],
  ['historySize', HISTORY_SIZE],
  ['finiteDiffStep', FINITE_DIFF_STEP],
  ['mu', PENALTY_COEFFICIENT],
]);

export const adamSettingsOpts = new Map<string, Setting>([
  ['maxIterations', {...MAX_ITERATIONS, default: 200}],
  ['maxFunctionEvaluations', {...MAX_EVALUATIONS, default: 1000}],
  ['learningRate', {
    default: 0.1,
    min: 1e-6,
    max: 10,
    caption: 'learning rate',
    tooltipText: 'Step size of the updates. Higher = faster, but less stable search',
    inputType: 'Float',
  }],
  ['finiteDiffStep', FINITE_DIFF_STEP],
  ['mu', PENALTY_COEFFICIENT],
]);

/** Smooth quadratic penalty for leaving the box, for the methods without native bounds */
function boxPenalty(settings: Map<string, number>, bounds?: OptimizerBounds): singleObjective.CommonSettings {
  if (bounds === undefined)
    return {};

  return {
    constraints: singleObjective.boxConstraints(bounds.lower, bounds.upper),
    penaltyOptions: {mu: settings.get('mu')},
  };
}

/** PSO scatters its swarm with its own PRNG: seeding it from the initial point makes PSO fits
 *  as reproducible as the initial points are */
function seedFromPoint(x0: Float64Array): number {
  return new Int32Array(Float64Array.from(x0).buffer).reduce((hash, word) => (Math.imul(hash, 31) + word) | 0, 17);
}

export function buildLbfgsbSettings(
  settings: Map<string, number>, bounds?: OptimizerBounds,
): singleObjective.LBFGSBSettings {
  return {
    maxIterations: settings.get('maxIterations'),
    maxFunctionEvaluations: settings.get('maxFunctionEvaluations'),
    gradTolerance: settings.get('gradTolerance'),
    historySize: settings.get('historySize'),
    finiteDiffStep: settings.get('finiteDiffStep'),
    bounds,
  };
}

export function buildPsoSettings(
  settings: Map<string, number>, bounds: OptimizerBounds | undefined, x0: Float64Array,
): singleObjective.PSOSettings {
  const isBoxFinite = (bounds !== undefined) && bounds.lower.every(Number.isFinite) &&
    bounds.upper.every(Number.isFinite);

  return {
    maxIterations: settings.get('maxIterations'),
    swarmSize: settings.get('swarmSize'),
    inertia: settings.get('inertia'),
    cognitive: settings.get('cognitive'),
    social: settings.get('social'),
    ...boxPenalty(settings, bounds),
    // Otherwise, PSO scatters the swarm around x0
    searchRange: isBoxFinite ? bounds : undefined,
    seed: seedFromPoint(x0),
  };
}

export function buildLbfgsSettings(
  settings: Map<string, number>, bounds?: OptimizerBounds,
): singleObjective.LBFGSSettings {
  return {
    maxIterations: settings.get('maxIterations'),
    gradTolerance: settings.get('gradTolerance'),
    historySize: settings.get('historySize'),
    finiteDiffStep: settings.get('finiteDiffStep'),
    ...boxPenalty(settings, bounds),
  };
}

export function buildAdamSettings(
  settings: Map<string, number>, bounds?: OptimizerBounds,
): singleObjective.AdamSettings {
  return {
    maxIterations: settings.get('maxIterations'),
    learningRate: settings.get('learningRate'),
    finiteDiffStep: settings.get('finiteDiffStep'),
    ...boxPenalty(settings, bounds),
  };
}
