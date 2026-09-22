import type {
  ProfileInputs, NcaRules, ComputeResult, ParameterValues,
  ParameterWarning, AucMethod, LambdaZResult, C0Provenance,
} from './types';
import {ROUTE_IV_BOLUS, ROUTE_IV_INFUSION} from './types';
import {applyBlqStrategy} from './blq';
import {augmentProfile} from './augment';
import {lambdaZBestFit, lambdaZManual} from './lambda-z';
import {
  aucLinearNaive, aucLogLinearNaive, aucLinearUpLogDownNaive,
  aucLinearCompensated, aucLogLinearCompensated, aucLinearUpLogDownCompensated,
  aucExtrapolateToInfinity,
} from './auc';
import {
  aumcLinearNaive, aumcLogLinearNaive, aumcLinearUpLogDownNaive,
  aumcLinearCompensated, aumcLogLinearCompensated, aumcLinearUpLogDownCompensated,
  aumcExtrapolateToInfinity,
} from './aumc';
import {
  halfLifeFromLambdaZ, clearance, volumeTerminal, pctExtrapolated,
  meanResidenceTime, volumeSteadyState, pctExtrapolatedAumc, tlag as tlagOf,
} from './derived';

type AucFn = (
  time: Float64Array, conc: Float64Array,
  startIdx: number, endIdx: number,
) => number;

function pickAucFn(method: AucMethod, compensated: boolean): AucFn {
  if (compensated) {
    if (method === 'linear') return aucLinearCompensated;
    if (method === 'log-linear') return aucLogLinearCompensated;
    return aucLinearUpLogDownCompensated;
  }
  if (method === 'linear') return aucLinearNaive;
  if (method === 'log-linear') return aucLogLinearNaive;
  return aucLinearUpLogDownNaive;
}

/** Moment-curve kernel matching the same method/summation choice as AUC. */
function pickAumcFn(method: AucMethod, compensated: boolean): AucFn {
  if (compensated) {
    if (method === 'linear') return aumcLinearCompensated;
    if (method === 'log-linear') return aumcLogLinearCompensated;
    return aumcLinearUpLogDownCompensated;
  }
  if (method === 'linear') return aumcLinearNaive;
  if (method === 'log-linear') return aumcLogLinearNaive;
  return aumcLinearUpLogDownNaive;
}

const NAN_VALUES: ParameterValues = Object.freeze({
  cmax: NaN, tmax: NaN, aucLast: NaN, aucInf: NaN, pctExtrap: NaN,
  lambdaZ: NaN, halfLife: NaN, cl: NaN, vz: NaN,
  aumcLast: NaN, aumcInf: NaN, mrt: NaN, vss: NaN, tlag: NaN,
  pctExtrapAumc: NaN,
});

/**
 * Single-profile NCA pipeline. Stateless. Takes the per-profile inputs
 * (typed arrays + dosing context) and the resolved rules, returns a
 * `ComputeResult` containing parameter values, full provenance, and a
 * status flag.
 *
 * Pipeline:
 * 1–3. `augmentProfile` (one exported kernel, `augment.ts`): BLQ
 *    pre-processing → observed Cmax/Tmax on the un-augmented post-BLQ profile
 *    → route-aware dose-time augmentation. IV bolus: a positive measured t=0
 *    sample is the c0; a missing OR BLQ / non-positive t=0 row is replaced by
 *    `(0, c0)` from the PKNCA c0 chain (`insertC0`). Extravascular and IV
 *    infusion: `(0, 0)` is prepended when no kept t=0 row exists. The kernel
 *    returns TWO masks — the effective BLQ mask (what the observed quantities
 *    read) and the drop set (what integration and λz skip) — plus the
 *    augmented ↔ input index map.
 * 4. Compute AUClast / AUMClast over the augmented profile (skipping the drop
 *    set), using the method and summation strategy in `rules`.
 * 5. Fit lambda_z (auto best-fit or manual) on the augmented profile.
 * 6. Derive AUCinf (= AUClast + cLast/λz), t½, CL, Vz, %AUCextrap; for IV
 *    bolus, the back-extrapolated share of AUCinf (`provenance.c0`).
 * 7. Generate quality warnings.
 *
 * IV-bolus AUC convention: the profile integrates FROM the back-extrapolated
 * c0 (Phoenix WinNonlin `AUC` with `C0`), so the same subject does not lose
 * the dose-time → first-sample area because a pre-dose sample happened to be
 * drawn. Stock PKNCA does not do this — its `auclast` on a raw IV-bolus
 * profile is NA without a t=0 datum and integrates from the observed `(0, 0)`
 * when one is present (measured 2026-09-22) — so the reference fixtures feed
 * PKNCA the augmented profile. A deliberate, documented divergence.
 *
 * Status:
 * - `'failed'`  — no measurable point (every conc was BLQ).
 * - `'partial'` — Cmax/AUClast obtained but lambda_z could not be fit
 *                 (no AUCinf, t½, CL, Vz).
 * - `'ok'`      — every parameter computed.
 */
export function computeNca(inputs: ProfileInputs, rules: NcaRules): ComputeResult {
  // Steps 1–3: the augmentation kernel.
  const aug = augmentProfile(inputs, rules.blq);
  if (aug === null) {
    // No measurable point. The BLQ trace is still reported (the kernel bails
    // after Step 2, so re-run Step 1 for the provenance the auditor expects).
    return {
      values: NAN_VALUES,
      provenance: {
        lambdaZ: null,
        blqApplied: applyBlqStrategy(inputs.conc, inputs.blqMask, inputs.lloq, 0, rules.blq),
        aucMethod: rules.aucMethod,
        compensated: rules.compensatedSummation,
        warnings: [],
        c0: null,
      },
      status: 'failed',
    };
  }
  const {time: augTime, conc: augConc, blqMask: augBlq, dropMask, observedCmax} = aug;
  const cmaxIdxForFit = aug.cmaxIdx;

  // Step 4: AUClast over the augmented profile (skipping the drop set).
  const dense = collectMeasurable(augTime, augConc, dropMask);
  const aucFn = pickAucFn(rules.aucMethod, rules.compensatedSummation);
  let aucLast = NaN;
  let aumcLast = NaN;
  let cLast = NaN;
  let tLast = NaN;
  if (dense.time.length >= 2) {
    const last = dense.time.length - 1;
    const aumcFn = pickAumcFn(rules.aucMethod, rules.compensatedSummation);
    aucLast = aucFn(dense.time, dense.conc, 0, last);
    aumcLast = aumcFn(dense.time, dense.conc, 0, last);
    cLast = dense.conc[last];
    tLast = dense.time[last];
  } else if (dense.time.length === 1) {
    aucLast = 0;
    aumcLast = 0;
    cLast = dense.conc[0];
    tLast = dense.time[0];
  }

  // Tlag is an OBSERVED quantity (independent of lambda_z) and an absorption
  // concept — reported for extravascular routes only; NaN for IV.
  const isIv =
    inputs.route === ROUTE_IV_BOLUS || inputs.route === ROUTE_IV_INFUSION;
  let tlag = NaN;
  if (!isIv && augTime.length >= 1) {
    // Computed on the AUGMENTED series with BLQ/excluded points treated as 0,
    // NOT the dense (BLQ-removed) profile: a BLQ sample before the first
    // measurable point IS the lag boundary, and dropping it would
    // underestimate Tlag. Reads the EFFECTIVE BLQ mask, not the drop set: a
    // substituted BLQ sample is not a quantifiable observation, so it cannot
    // end the lag (documented divergence from PKNCA's tlag on the cleaned
    // profile under exclude / missing / set-half-lloq).
    const tlagConc = new Float64Array(augTime.length);
    for (let i = 0; i < augTime.length; i++) {
      tlagConc[i] =
        (augBlq[i] !== 0 || !Number.isFinite(augConc[i])) ? 0 : augConc[i];
    }
    tlag = tlagOf(augTime, tlagConc);
  }

  // Step 5: lambda_z — on the drop set; the fit itself drops `conc <= 0`.
  const lambdaZRes: LambdaZResult | null =
    (rules.lambdaZ.mode === 'auto-best-fit') ?
      lambdaZBestFit(augTime, augConc, dropMask, cmaxIdxForFit, rules.lambdaZ) :
      (rules.lambdaZ.mode === 'manual-points' && rules.lambdaZ.manualPoints) ?
        lambdaZManual(augTime, augConc, rules.lambdaZ.manualPoints) :
        null;

  // Step 6: AUCinf and derived parameters.
  let aucInf = NaN;
  let aumcInf = NaN;
  let halfLife = NaN;
  let cl = NaN;
  let vz = NaN;
  let vss = NaN;
  let mrt = NaN;
  let pctExtrap = NaN;
  let pctExtrapAumc = NaN;
  let status: 'ok' | 'partial' | 'failed' = 'partial';

  if (
    lambdaZRes !== null && lambdaZRes.lambdaZ > 0 &&
    Number.isFinite(aucLast) && Number.isFinite(cLast) && cLast > 0
  ) {
    // Zero-order infusion correction: only IV-infusion with a known duration
    // carries a non-zero T_inf. Bolus / extravascular → T_inf = 0.
    const tInf = (inputs.route === ROUTE_IV_INFUSION &&
      inputs.infusionDuration !== null && inputs.infusionDuration > 0) ?
      inputs.infusionDuration : 0;

    aucInf = aucLast + aucExtrapolateToInfinity(cLast, lambdaZRes.lambdaZ);
    aumcInf =
      aumcLast + aumcExtrapolateToInfinity(tLast, cLast, lambdaZRes.lambdaZ);
    halfLife = halfLifeFromLambdaZ(lambdaZRes.lambdaZ);
    cl = clearance(inputs.dose, aucInf);
    vz = volumeTerminal(inputs.dose, lambdaZRes.lambdaZ, aucInf);
    mrt = meanResidenceTime(aumcInf, aucInf, tInf);
    // IV-only: for extravascular data Vss would be Vss/F confounded by
    // absorption — a category error, not a high number.
    vss = isIv ? volumeSteadyState(inputs.dose, aumcInf, aucInf, tInf) : NaN;
    pctExtrap = pctExtrapolated(aucLast, aucInf);
    pctExtrapAumc = pctExtrapolatedAumc(aumcLast, aumcInf);
    status = 'ok';
  }

  // c0 provenance (IV bolus only). The back-extrapolated share of AUCinf is
  // the dose-time → first-observation segment — `[0, 1]` of the dense profile,
  // whose slot 0 IS the inserted c0 — integrated with the SAME method and
  // summation the profile used, as a percentage of AUCinf. 0 for an observed
  // dose-time value (nothing was extrapolated); NaN on 'partial' (no AUCinf).
  let c0Prov: C0Provenance | null = null;
  if (aug.c0 !== null) {
    const pct =
      aug.c0.method === 'observed' ? 0 :
        (status === 'ok' && dense.time.length >= 2) ?
          aucFn(dense.time, dense.conc, 0, 1) / aucInf * 100 :
          NaN;
    c0Prov = {...aug.c0, pctAucBackExtrap: pct};
  }

  // Step 7: Quality warnings.
  const warnings: ParameterWarning[] = [];
  // The log-slope was not estimable and c0 rests on a plateau assumption (the
  // first observation, the minimum, or zero). Routine inserts / replacements
  // are recorded in `provenance.c0`, not warned — a warning is for the
  // genuinely fragile anchor.
  if (c0Prov !== null &&
      (c0Prov.method === 'c1' || c0Prov.method === 'cmin' || c0Prov.method === 'set0')) {
    warnings.push({
      code: 'C0_FALLBACK',
      severity: 'warning',
      message:
        `c0 could not be back-extrapolated (log-slope not estimable); ` +
        `fell back to '${c0Prov.method}' (${c0Prov.value}) — the dose-time ` +
        `concentration and AUC from t = 0 rest on a plateau assumption`,
    });
  }
  if (status === 'ok' && pctExtrap > rules.extrapWarnPct) {
    warnings.push({
      code: 'AUC_EXTRAP_HIGH',
      severity: pctExtrap > rules.extrapErrorPct ? 'error' : 'warning',
      message:
        `% AUC extrapolated (${pctExtrap.toFixed(1)}%) exceeds threshold ` +
        `(${rules.extrapWarnPct}%)`,
    });
  }
  if (status === 'ok' && pctExtrapAumc > rules.extrapWarnPctAumc) {
    warnings.push({
      code: 'AUMC_EXTRAP_HIGH',
      severity: 'warning',
      message:
        `% AUMC extrapolated (${pctExtrapAumc.toFixed(1)}%) exceeds threshold ` +
        `(${rules.extrapWarnPctAumc}%) — MRT/Vss are fragile`,
    });
  }
  // A point that is BLQ-flagged but NOT in the drop set is a SUBSTITUTED value
  // (`set-half-lloq` — `set-zero` substitutes are non-positive and the λz filter
  // drops them itself). When one lands inside the accepted window, the terminal
  // slope is fitted partly through a number nobody measured, and the fit
  // statistics cannot say so: a substitute near the trend line scores WELL
  // precisely because it is near the line. Measured on the `06_blq_rules` P4
  // profile — a trailing LLOQ/2 substitute joins a 3-point fit at adj-R² 0.9995
  // and becomes the terminal anchor. Reported, never inferred.
  if (lambdaZRes !== null) {
    const substituted: number[] = [];
    for (let k = 0; k < lambdaZRes.pointsUsed.length; k++) {
      const i = lambdaZRes.pointsUsed[k];
      if (augBlq[i] !== 0 && dropMask[i] === 0) substituted.push(augTime[i]);
    }
    if (substituted.length > 0) {
      warnings.push({
        code: 'LAMBDAZ_SUBSTITUTED_BLQ',
        severity: 'warning',
        message:
          `lambda_z was fitted through ${substituted.length} substituted BLQ ` +
          `value(s) (t = ${substituted.join(', ')}) — the terminal slope and ` +
          `everything derived from it rest partly on unmeasured concentrations`,
      });
    }
  }
  if (lambdaZRes !== null &&
      lambdaZRes.pointsUsed.length <= rules.lambdaZ.minPoints) {
    warnings.push({
      code: 'LAMBDAZ_FEW_POINTS',
      severity: 'info',
      message:
        `lambda_z fit used the minimum allowed number of points ` +
        `(${lambdaZRes.pointsUsed.length})`,
    });
  }
  // Opt-in terminal-phase span diagnostic — see LambdaZResult.spanRatio. Warns
  // but never rejects: a short-span fit is one a scientist should look at, not
  // one the engine should discard for them. No NaN guard needed, since
  // `NaN < threshold` is false.
  if (lambdaZRes !== null && rules.lambdaZ.minSpanRatio !== undefined &&
      lambdaZRes.spanRatio < rules.lambdaZ.minSpanRatio) {
    warnings.push({
      code: 'LAMBDAZ_LOW_SPAN',
      severity: 'warning',
      message:
        `lambda_z window spans ${lambdaZRes.spanRatio.toFixed(2)} half-lives, below ` +
        `the ${rules.lambdaZ.minSpanRatio} threshold — the terminal slope rests on ` +
        `too short a window (adjusted R² cannot detect this at small n)`,
    });
  }
  const blqFraction = countBlq(inputs.blqMask) / inputs.blqMask.length;
  if (blqFraction > 0.5) {
    warnings.push({
      code: 'BLQ_HIGH_FRACTION',
      severity: 'warning',
      message: `BLQ fraction = ${(blqFraction * 100).toFixed(1)}%`,
    });
  }

  return {
    values: {
      cmax: observedCmax.cmax,
      tmax: observedCmax.tmax,
      aucLast,
      aucInf,
      pctExtrap,
      lambdaZ: lambdaZRes !== null ? lambdaZRes.lambdaZ : NaN,
      halfLife,
      cl,
      vz,
      aumcLast,
      aumcInf,
      mrt,
      vss,
      tlag,
      pctExtrapAumc,
    },
    provenance: {
      lambdaZ: lambdaZRes,
      blqApplied: aug.blqApplied,
      aucMethod: rules.aucMethod,
      compensated: rules.compensatedSummation,
      warnings,
      c0: c0Prov,
    },
    status,
  };
}

/** Build dense Float64Arrays of (time, conc) skipping the drop set and NaN entries. */
function collectMeasurable(
  time: Float64Array, conc: Float64Array, dropMask: Uint8Array,
): {time: Float64Array; conc: Float64Array} {
  const tBuf: number[] = [];
  const cBuf: number[] = [];
  for (let i = 0; i < time.length; i++) {
    if (dropMask[i] !== 0) continue;
    if (!Number.isFinite(conc[i])) continue;
    tBuf.push(time[i]);
    cBuf.push(conc[i]);
  }
  // Drop the TRAILING run of non-positive concentrations, mirroring PKNCA's
  // `conc.blq` default. A trailing zero is an unflagged below-LLOQ washout
  // sample — common when a dataset encodes BLQ as conc=0 and carries no
  // LLOQ/BLQ-flag column, leaving `blqMask` all-zeros. It must not anchor the
  // terminal: as `cLast` it is the λz extrapolation base, so `cLast = 0`
  // blocks AUCinf/t½/CL/Vz even for a well-formed λz, and it adds a spurious
  // tail-to-zero area to AUClast. EMBEDDED zeros are deliberately KEPT —
  // `lambdaZBestFit` already excludes `conc ≤ 0` from the regression.
  while (cBuf.length > 0 && cBuf[cBuf.length - 1] <= 0) {
    cBuf.pop();
    tBuf.pop();
  }
  return {time: Float64Array.from(tBuf), conc: Float64Array.from(cBuf)};
}

function countBlq(blqMask: Uint8Array): number {
  let n = 0;
  for (let i = 0; i < blqMask.length; i++) if (blqMask[i] !== 0) n++;
  return n;
}
