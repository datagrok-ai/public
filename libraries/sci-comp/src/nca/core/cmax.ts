import type {CmaxResult} from './types';

/**
 * Find the observed peak concentration and its time.
 *
 * Returns the **first occurrence** of the maximum when several points share
 * the peak value (PKNCA convention — see
 * <https://billdenney.github.io/pknca/reference/pk.calc.cmax.html>).
 * BLQ points (`blqMask[i] !== 0`) are skipped. `NaN` concentrations (e.g.
 * the result of the `missing` BLQ rule) are also skipped because `NaN > x`
 * is always `false`.
 *
 * @param time - Time vector, sorted ascending.
 * @param conc - Concentration vector, same length as `time`.
 * @param blqMask - 1 = BLQ, 0 = measurable. Same length as `time`.
 * @returns The Cmax/Tmax result, or `null` when no measurable point exists.
 */
export function findCmax(
  time: Float64Array, conc: Float64Array, blqMask: Uint8Array,
): CmaxResult | null {
  let cmax = -Infinity;
  let cmaxIdx = -1;
  for (let i = 0; i < conc.length; i++) {
    if (blqMask[i] !== 0) continue;
    if (conc[i] > cmax) {
      cmax = conc[i];
      cmaxIdx = i;
    }
  }
  if (cmaxIdx === -1) return null;
  return {cmax, tmax: time[cmaxIdx], cmaxIdx};
}
