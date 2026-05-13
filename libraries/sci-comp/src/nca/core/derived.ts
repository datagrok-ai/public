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
 * a quality flag — see PRD FR-222 and `NcaRules.extrapWarnPct`.
 */
export function pctExtrapolated(aucLast: number, aucInf: number): number {
  return (aucInf - aucLast) / aucInf * 100;
}
