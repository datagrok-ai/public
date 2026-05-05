/**
 * Welch's t-test pairwise: each treated group vs the control.
 *
 * Returns raw (uncorrected) p-values rounded to 6 decimals — matching the
 * Python convention. Apply `bonferroniCorrect` for multiplicity adjustment.
 */

import {NumericInput} from '../types';
import {welchTTest} from './welch-t';

export interface TreatedGroup {
  doseLevel: number;
  values: NumericInput;
}

export interface WelchPairwiseResult {
  doseLevel: number;
  pValueWelch: number | null;
}

/**
 * Run Welch's t-test for each treated group vs the control.
 *
 * Returns an empty list when the control has fewer than 2 finite values or
 * when no treated groups are provided.
 */
export function welchPairwise(
  control: NumericInput,
  treated: readonly TreatedGroup[],
): WelchPairwiseResult[] {
  const ctrl = stripNaNCount(control);
  if (ctrl < 2 || treated.length === 0) return [];

  const results: WelchPairwiseResult[] = [];
  for (const t of treated) {
    const r = welchTTest(t.values, control);
    const p = r.pValue;
    results.push({
      doseLevel: t.doseLevel,
      pValueWelch: p === null || Number.isNaN(p) ? null : round6(p),
    });
  }
  return results;
}

function stripNaNCount(x: NumericInput): number {
  let n = 0;
  for (let i = 0; i < x.length; i++)
    if (!Number.isNaN(x[i] as number)) n++;
  return n;
}

function round6(x: number): number {
  return Math.round(x * 1e6) / 1e6;
}
