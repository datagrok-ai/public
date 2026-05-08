import type {
  ProfileInputs, NcaRules, ComputeResult, ParameterValues,
  ParameterWarning, AucMethod, BlqProcessingResult, LambdaZResult,
} from './types';
import {ROUTE_IV_BOLUS} from './types';
import {applyBlqStrategy} from './blq';
import {findCmax} from './cmax';
import {insertC0} from './c0';
import {lambdaZBestFit, lambdaZManual} from './lambda-z';
import {
  aucLinearNaive, aucLogLinearNaive, aucLinearUpLogDownNaive,
  aucLinearCompensated, aucLogLinearCompensated, aucLinearUpLogDownCompensated,
  aucExtrapolateToInfinity,
} from './auc';
import {
  halfLifeFromLambdaZ, clearance, volumeTerminal, pctExtrapolated,
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

const NAN_VALUES: ParameterValues = Object.freeze({
  cmax: NaN, tmax: NaN, aucLast: NaN, aucInf: NaN, pctExtrap: NaN,
  lambdaZ: NaN, halfLife: NaN, cl: NaN, vz: NaN,
});

/**
 * Single-profile NCA pipeline. Stateless. Takes the per-profile inputs
 * (typed arrays + dosing context) and the resolved rules, returns a
 * `ComputeResult` containing parameter values, full provenance, and a
 * status flag.
 *
 * Pipeline:
 * 1. Apply BLQ pre-processing → modified concentrations + excluded indices.
 * 2. For IV bolus profiles without a t=0 observation, estimate c0 by
 *    log-linear back-extrapolation and insert (0, c0). This mirrors the
 *    PKNCA convention used by every reference NCA tool.
 * 3. Find observed Cmax/Tmax on the original (non-augmented) post-BLQ
 *    profile — for IV bolus the inserted c0 is excluded from the reported
 *    peak.
 * 4. Compute AUClast over the augmented profile (skipping NaN/excluded
 *    points), using the method and summation strategy in `rules`.
 * 5. Fit lambda_z (auto best-fit or manual) on the augmented profile.
 * 6. Derive AUCinf (= AUClast + cLast/λz), t½, CL, Vz, %AUCextrap.
 * 7. Generate quality warnings.
 *
 * Status:
 * - `'failed'`  — no measurable point (every conc was BLQ).
 * - `'partial'` — Cmax/AUClast obtained but lambda_z could not be fit
 *                 (no AUCinf, t½, CL, Vz).
 * - `'ok'`      — every parameter computed.
 */
export function computeNca(inputs: ProfileInputs, rules: NcaRules): ComputeResult {
  // ─────────────────────────────────────────────────────────────────────
  // Step 1: BLQ pre-processing on the raw concentrations.
  // ─────────────────────────────────────────────────────────────────────
  const blqRes: BlqProcessingResult = applyBlqStrategy(
    inputs.conc, inputs.blqMask, inputs.lloq, 0, rules.blq,
  );
  const procConc = blqRes.conc;
  // Effective BLQ mask = original BLQ ∪ excluded by rule.
  const effBlq = new Uint8Array(inputs.blqMask);
  for (let k = 0; k < blqRes.excluded.length; k++)
    effBlq[blqRes.excluded[k]] = 1;

  // ─────────────────────────────────────────────────────────────────────
  // Step 2: Observed Cmax/Tmax on raw post-BLQ profile.
  // ─────────────────────────────────────────────────────────────────────
  const observedCmax = findCmax(inputs.time, procConc, effBlq);
  if (observedCmax === null) {
    return {
      values: NAN_VALUES,
      provenance: {
        lambdaZ: null,
        blqApplied: blqRes,
        aucMethod: rules.aucMethod,
        compensated: rules.compensatedSummation,
        warnings: [],
      },
      status: 'failed',
    };
  }

  // ─────────────────────────────────────────────────────────────────────
  // Step 3: Augment with a t=0 observation when missing, mirroring PKNCA:
  //   - IV bolus      → insert (0, c0) where c0 is back-extrapolated.
  //   - extravascular → insert (0, 0)  (pre-dose conc = 0 by convention).
  // An existing t=0 row counts even when conc(0) = 0 — that's a valid
  // pre-dose observation for extravascular profiles, not a missing value.
  // ─────────────────────────────────────────────────────────────────────
  const hasT0 = (
    inputs.time.length > 0 && inputs.time[0] === 0 && effBlq[0] === 0 &&
    Number.isFinite(procConc[0])
  );
  let augTime = inputs.time;
  let augConc = procConc;
  let augBlq = effBlq;
  let cmaxIdxForFit = observedCmax.cmaxIdx;
  if (!hasT0) {
    if (inputs.route === ROUTE_IV_BOLUS) {
      const aug = insertC0(inputs.time, procConc, effBlq);
      if (aug !== null) {
        augTime = aug.time;
        augConc = aug.conc;
        augBlq = aug.blqMask;
        cmaxIdxForFit = aug.cmaxIdx; // = 0 (inserted c0 is the new peak)
      }
    } else {
      // Extravascular: pre-dose concentration is 0 by convention.
      augTime = prependScalar(inputs.time, 0);
      augConc = prependScalar(procConc, 0);
      augBlq = prependByte(effBlq, 0);
      cmaxIdxForFit = observedCmax.cmaxIdx + 1;
    }
  }

  // ─────────────────────────────────────────────────────────────────────
  // Step 4: AUClast over the augmented profile (skipping NaN/excluded).
  // ─────────────────────────────────────────────────────────────────────
  const dense = collectMeasurable(augTime, augConc, augBlq);
  let aucLast = NaN;
  let cLast = NaN;
  if (dense.time.length >= 2) {
    const aucFn = pickAucFn(rules.aucMethod, rules.compensatedSummation);
    aucLast = aucFn(dense.time, dense.conc, 0, dense.time.length - 1);
    cLast = dense.conc[dense.conc.length - 1];
  } else if (dense.time.length === 1) {
    aucLast = 0;
    cLast = dense.conc[0];
  }

  // ─────────────────────────────────────────────────────────────────────
  // Step 5: lambda_z.
  // ─────────────────────────────────────────────────────────────────────
  const lambdaZRes: LambdaZResult | null =
    (rules.lambdaZ.mode === 'auto-best-fit') ?
      lambdaZBestFit(augTime, augConc, augBlq, cmaxIdxForFit, rules.lambdaZ) :
      (rules.lambdaZ.mode === 'manual-points' && rules.lambdaZ.manualPoints) ?
        lambdaZManual(augTime, augConc, rules.lambdaZ.manualPoints) :
        null;

  // ─────────────────────────────────────────────────────────────────────
  // Step 6: AUCinf and derived parameters.
  // ─────────────────────────────────────────────────────────────────────
  let aucInf = NaN;
  let halfLife = NaN;
  let cl = NaN;
  let vz = NaN;
  let pctExtrap = NaN;
  let status: 'ok' | 'partial' | 'failed' = 'partial';

  if (
    lambdaZRes !== null && lambdaZRes.lambdaZ > 0 &&
    Number.isFinite(aucLast) && Number.isFinite(cLast) && cLast > 0
  ) {
    aucInf = aucLast + aucExtrapolateToInfinity(cLast, lambdaZRes.lambdaZ);
    halfLife = halfLifeFromLambdaZ(lambdaZRes.lambdaZ);
    cl = clearance(inputs.dose, aucInf);
    vz = volumeTerminal(inputs.dose, lambdaZRes.lambdaZ, aucInf);
    pctExtrap = pctExtrapolated(aucLast, aucInf);
    status = 'ok';
  }

  // ─────────────────────────────────────────────────────────────────────
  // Step 7: Quality warnings.
  // ─────────────────────────────────────────────────────────────────────
  const warnings: ParameterWarning[] = [];
  if (status === 'ok' && pctExtrap > rules.extrapWarnPct) {
    warnings.push({
      code: 'AUC_EXTRAP_HIGH',
      severity: pctExtrap > rules.extrapErrorPct ? 'error' : 'warning',
      message:
        `% AUC extrapolated (${pctExtrap.toFixed(1)}%) exceeds threshold ` +
        `(${rules.extrapWarnPct}%)`,
    });
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
    },
    provenance: {
      lambdaZ: lambdaZRes,
      blqApplied: blqRes,
      aucMethod: rules.aucMethod,
      compensated: rules.compensatedSummation,
      warnings,
    },
    status,
  };
}

/** Build dense Float64Arrays of (time, conc) skipping BLQ and NaN entries. */
function collectMeasurable(
  time: Float64Array, conc: Float64Array, blqMask: Uint8Array,
): {time: Float64Array; conc: Float64Array} {
  const tBuf: number[] = [];
  const cBuf: number[] = [];
  for (let i = 0; i < time.length; i++) {
    if (blqMask[i] !== 0) continue;
    if (!Number.isFinite(conc[i])) continue;
    tBuf.push(time[i]);
    cBuf.push(conc[i]);
  }
  return {time: Float64Array.from(tBuf), conc: Float64Array.from(cBuf)};
}

function countBlq(blqMask: Uint8Array): number {
  let n = 0;
  for (let i = 0; i < blqMask.length; i++) if (blqMask[i] !== 0) n++;
  return n;
}

function prependScalar(src: Float64Array, value: number): Float64Array {
  const out = new Float64Array(src.length + 1);
  out[0] = value;
  out.set(src, 1);
  return out;
}

function prependByte(src: Uint8Array, value: number): Uint8Array {
  const out = new Uint8Array(src.length + 1);
  out[0] = value;
  out.set(src, 1);
  return out;
}
