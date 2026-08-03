// Shared helpers for fitting tests.

import {nelderMeadSettingsOpts} from './imports';
import type {EarlyStoppingSettings, ReproSettings, ValueBoundsData} from './imports';

export function defaultNmSettings(overrides: Record<string, number> = {}): Map<string, number> {
  const m = new Map<string, number>();
  nelderMeadSettingsOpts.forEach((opts, key) => m.set(key, overrides[key] ?? opts.default));
  return m;
}

export function reproSettings(seed = 42): ReproSettings {
  return {reproducible: true, seed};
}

export function noEarlyStopping(): EarlyStoppingSettings {
  return {
    useEarlyStopping: false,
    costFuncThreshold: 0,
    stopAfter: 1,
    useAboveThresholdPoints: false,
  };
}

export function earlyStop(threshold: number, stopAfter: number,
  useAboveThresholdPoints = false): EarlyStoppingSettings {
  return {
    useEarlyStopping: true,
    costFuncThreshold: threshold,
    stopAfter,
    useAboveThresholdPoints,
  };
}

export function constBound(value: number): ValueBoundsData {
  return {type: 'const', value};
}

export function rangeBound(bottom: number, top: number, name: string): ValueBoundsData {
  return {
    type: 'changing',
    bottom: {name, type: 'value', value: bottom},
    top: {name, type: 'value', value: top},
  };
}

export function formulaBound(bottom: string, top: string, name: string): ValueBoundsData {
  return {
    type: 'changing',
    bottom: {name, type: 'formula', formula: bottom},
    top: {name, type: 'formula', formula: top},
  };
}

/** Returns true iff iterCosts[0..iterCount] is non-increasing. */
export function isMonotoneNonIncreasing(iterCosts: number[], iterCount: number, eps = 1e-12): boolean {
  for (let i = 1; i < iterCount; ++i) {
    if (iterCosts[i] - iterCosts[i - 1] > eps)
      return false;
  }
  return true;
}
