/**
 * Estimation of c0 — the concentration at the dose time for IV bolus
 * dosing when the data set does not contain an observation at t = dose.
 *
 * Faithful port of PKNCA 0.12 `pk.calc.c0` (see
 * <https://billdenney.github.io/pknca/reference/pk.calc.c0.html>):
 * the default method tries each strategy in order and returns the first
 * non-null result.
 */

import {findCmax} from './cmax';

/** Available c0 estimation strategies. */
export type C0Method = 'c0' | 'logslope' | 'c1' | 'cmin' | 'set0';

/** Default PKNCA-compatible method chain. */
export const C0_DEFAULT_METHODS: ReadonlyArray<C0Method> =
  ['c0', 'logslope', 'c1', 'cmin', 'set0'];

/** Configuration for c0 estimation: dose time and method-chain override. */
export interface C0Options {
  /** Time of dose. Defaults to `0`. */
  readonly timeDose?: number;
  /** Method chain to try in order. Defaults to {@link C0_DEFAULT_METHODS}. */
  readonly methods?: ReadonlyArray<C0Method>;
}

/**
 * Estimate the dose-time concentration for an IV bolus profile.
 *
 * The function never throws; on inputs where every method fails it returns
 * `null`. NaN concentrations are treated as missing.
 *
 * @param time - Time vector, sorted ascending.
 * @param conc - Concentration vector, same length as `time`.
 * @param blqMask - 1 = BLQ, 0 = measurable. Same length as `time`.
 * @param options - Optional dose time and method chain.
 * @returns The estimated c0 (>= 0), or `null` if every strategy fails.
 */
export function estimateC0(
  time: Float64Array, conc: Float64Array, blqMask: Uint8Array,
  options: C0Options = {},
): number | null {
  const timeDose = options.timeDose ?? 0;
  const methods = options.methods ?? C0_DEFAULT_METHODS;
  for (const m of methods) {
    const v = applyMethod(m, time, conc, blqMask, timeDose);
    if (v !== null) return v;
  }
  return null;
}

function applyMethod(
  method: C0Method,
  time: Float64Array, conc: Float64Array, blqMask: Uint8Array,
  timeDose: number,
): number | null {
  switch (method) {
  case 'c0': return c0AtDose(time, conc, blqMask, timeDose);
  case 'logslope': return c0LogSlope(time, conc, blqMask, timeDose);
  case 'c1': return c0FirstPostDose(time, conc, blqMask, timeDose);
  case 'cmin': return c0CMin(conc, blqMask);
  case 'set0': return 0;
  }
}

/** "c0" — return an existing positive observation at the dose time, if any. */
function c0AtDose(
  time: Float64Array, conc: Float64Array, blqMask: Uint8Array,
  timeDose: number,
): number | null {
  for (let i = 0; i < time.length; i++) {
    if (time[i] !== timeDose) continue;
    if (blqMask[i] !== 0) continue;
    const c = conc[i];
    if (Number.isFinite(c) && c > 0) return c;
  }
  return null;
}

/**
 * "logslope" — log-linear back-extrapolation of the first two post-dose
 * measurable concentrations. Used only when c2 < c1 (decay) and c2 != 0.
 *
 * Formula:
 *   slope = (ln c2 − ln c1) / (t2 − t1)
 *   c0    = exp(ln c1 − slope · (t1 − timeDose))
 */
function c0LogSlope(
  time: Float64Array, conc: Float64Array, blqMask: Uint8Array,
  timeDose: number,
): number | null {
  let i1 = -1;
  let i2 = -1;
  for (let i = 0; i < time.length; i++) {
    if (time[i] <= timeDose) continue;
    if (blqMask[i] !== 0) continue;
    const c = conc[i];
    if (!Number.isFinite(c)) continue;
    if (i1 === -1) i1 = i;
    else if (i2 === -1) {i2 = i; break;}
  }
  if (i1 === -1 || i2 === -1) return null;
  const c1 = conc[i1];
  const c2 = conc[i2];
  if (!(c2 < c1) || c2 === 0) return null;
  const t1 = time[i1];
  const t2 = time[i2];
  const slope = (Math.log(c2) - Math.log(c1)) / (t2 - t1);
  return Math.exp(Math.log(c1) - slope * (t1 - timeDose));
}

/** "c1" — use the first post-dose measurable concentration as c0. */
function c0FirstPostDose(
  time: Float64Array, conc: Float64Array, blqMask: Uint8Array,
  timeDose: number,
): number | null {
  for (let i = 0; i < time.length; i++) {
    if (time[i] <= timeDose) continue;
    if (blqMask[i] !== 0) continue;
    const c = conc[i];
    if (Number.isFinite(c)) return c;
  }
  return null;
}

/** "cmin" — minimum non-zero, non-BLQ concentration. */
function c0CMin(conc: Float64Array, blqMask: Uint8Array): number | null {
  let cmin = Infinity;
  let found = false;
  for (let i = 0; i < conc.length; i++) {
    if (blqMask[i] !== 0) continue;
    const c = conc[i];
    if (!Number.isFinite(c) || c <= 0) continue;
    if (c < cmin) {
      cmin = c;
      found = true;
    }
  }
  return found ? cmin : null;
}

/**
 * Insert an estimated c0 at `timeDose` into a profile, returning fresh
 * arrays (the input is not mutated). Useful for IV bolus pre-processing
 * before AUC and lambda_z computation.
 *
 * Re-finds Cmax on the augmented profile; for IV bolus this typically
 * shifts Tmax to `timeDose`.
 *
 * @returns The augmented profile and the new Cmax index, or `null` when
 *          c0 estimation fails entirely.
 */
export function insertC0(
  time: Float64Array, conc: Float64Array, blqMask: Uint8Array,
  options: C0Options = {},
): {time: Float64Array; conc: Float64Array; blqMask: Uint8Array;
    c0: number; cmaxIdx: number} | null {
  const c0 = estimateC0(time, conc, blqMask, options);
  if (c0 === null) return null;
  const timeDose = options.timeDose ?? 0;
  const n = time.length;
  const time2 = new Float64Array(n + 1);
  const conc2 = new Float64Array(n + 1);
  const blqMask2 = new Uint8Array(n + 1);
  time2[0] = timeDose;
  conc2[0] = c0;
  // Original tail
  time2.set(time, 1);
  conc2.set(conc, 1);
  blqMask2.set(blqMask, 1);

  const cmaxR = findCmax(time2, conc2, blqMask2);
  if (cmaxR === null) return null;
  return {time: time2, conc: conc2, blqMask: blqMask2, c0, cmaxIdx: cmaxR.cmaxIdx};
}
