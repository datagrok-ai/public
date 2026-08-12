/**
 * Derived NCA parameters — simple closed-form formulas that depend on
 * already-computed quantities (lambda_z, AUClast, AUCinf, dose).
 *
 * For extravascular routes the apparent forms `CL/F` and `Vz/F` are
 * obtained by these same formulas — the bioavailability `F` is folded
 * into the quotient implicitly. Caller decides on the physical
 * interpretation by inspecting `ProfileInputs.route`.
 */

/**
 * Terminal half-life: `ln(2) / λz`.
 *
 * Returns `+Infinity` for `λz = 0` and `NaN` for `λz < 0`. Callers should
 * guard upstream — a negative `λz` is a sign that the terminal slope was
 * not estimable.
 */
export function halfLifeFromLambdaZ(lambdaZ: number): number {
  if (lambdaZ < 0) return NaN;
  return Math.LN2 / lambdaZ;
}

/**
 * Clearance: `dose / AUCinf`.
 *
 * For extravascular dosing this is the apparent clearance `CL/F`; for IV
 * dosing it is the absolute clearance `CL`. Returns `+Infinity` for
 * `AUCinf = 0` with positive dose, and `NaN` for the `0/0` case.
 */
export function clearance(dose: number, aucInf: number): number {
  return dose / aucInf;
}

/**
 * Terminal volume of distribution: `dose / (λz · AUCinf)`.
 *
 * For extravascular dosing this is the apparent volume `Vz/F`. Mirrors
 * the divide-by-zero semantics of {@link clearance}.
 */
export function volumeTerminal(
  dose: number, lambdaZ: number, aucInf: number,
): number {
  return dose / (lambdaZ * aucInf);
}

/**
 * Percentage of AUCinf that came from the lambda_z extrapolation tail:
 * `(AUCinf − AUClast) / AUCinf · 100`.
 *
 * Returns `NaN` for the `0/0` case (`AUClast == AUCinf == 0`) and
 * `±Infinity` when `AUCinf == 0` with non-zero `AUClast` (degenerate —
 * caller should not reach this state). A high value (typically > 20%) is
 * a quality flag — see `NcaRules.extrapWarnPct`.
 */
export function pctExtrapolated(aucLast: number, aucInf: number): number {
  return (aucInf - aucLast) / aucInf * 100;
}

/**
 * Mean residence time: `AUMCinf / AUCinf − T_inf/2`.
 *
 * The `T_inf/2` term corrects for the duration of a zero-order IV infusion
 * (Perrier & Mayersohn 1982). For IV-bolus and extravascular routes pass
 * `infusionDuration = 0` (the default) — the correction vanishes and MRT is
 * the plain ratio. For extravascular dosing this is the absorption-inclusive
 * MRT (= MRT_iv + MAT), a legitimate reported quantity.
 *
 * Mirrors the divide-by-zero semantics of the underlying ratio: `+Infinity`
 * for `AUCinf = 0` with positive AUMCinf, `NaN` for the `0/0` case.
 */
export function meanResidenceTime(
  aumcInf: number, aucInf: number, infusionDuration = 0,
): number {
  return aumcInf / aucInf - infusionDuration / 2;
}

/**
 * Steady-state volume of distribution: `Vss = CL · MRT = dose · MRT / AUCinf`.
 *
 * Equivalent to `dose·AUMCinf/AUCinf²` for IV-bolus (`T_inf = 0`) and to
 * `CL·(AUMCinf/AUCinf − T_inf/2)` for IV-infusion (Perrier & Mayersohn 1982
 * Eq. 11). Delegates to {@link meanResidenceTime} so the `T_inf/2` correction
 * is defined in exactly one place.
 *
 * Valid only under linear, time-invariant disposition, and only for IV dosing
 * — the route gate is the caller's job.
 */
export function volumeSteadyState(
  dose: number, aumcInf: number, aucInf: number, infusionDuration = 0,
): number {
  return dose * meanResidenceTime(aumcInf, aucInf, infusionDuration) / aucInf;
}

/**
 * Percentage of `AUMCinf` contributed by the λz extrapolation tail:
 * `(AUMCinf − AUMClast) / AUMCinf · 100`.
 *
 * Because the first-moment tail is weighted by time, this is always ≥ the
 * corresponding %AUCextrap — a high value (no regulatory cutoff exists; see
 * `NcaRules.extrapWarnPctAumc`) flags a fragile MRT/Vss. Same degenerate-case
 * semantics as {@link pctExtrapolated}.
 */
export function pctExtrapolatedAumc(aumcLast: number, aumcInf: number): number {
  return (aumcInf - aumcLast) / aumcInf * 100;
}

/**
 * Absorption lag time: the time of the last `C == 0` sample immediately before
 * the first `C > 0` sample, computed on the BLQ-processed series (BLQ points
 * have been set to 0). Pins Tlag to the same quantifiable boundary AUC uses
 * (see `blq.ts`), so the two never desynchronise on zero-coded exports.
 *
 * Returns `time[0]` when the first sample is already positive (no observed
 * lag), and `NaN` when no sample is positive. Tlag is an absorption concept —
 * the IV route gate is the caller's job.
 */
export function tlag(time: Float64Array, conc: Float64Array): number {
  let firstPos = -1;
  for (let i = 0; i < conc.length; i++) {
    if (conc[i] > 0) {
      firstPos = i;
      break;
    }
  }
  if (firstPos < 0) return NaN;
  if (firstPos === 0) return time[0];
  return time[firstPos - 1];
}
