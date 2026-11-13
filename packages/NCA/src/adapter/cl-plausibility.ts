/**
 * Implied-clearance plausibility cross-check — the P2 dose-normalized arm of the
 * conc-mis-scale guard (the implied-CL band was ratified by a DMPK scientist).
 *
 * **Why a NEW post-compute module (mirror of `bioavailability.ts` /
 * `toxicokinetics.ts`).** This is a band check on the engine's *output* CL — it
 * needs the computed `result.values.cl` plus the per-profile body weight and
 * dose-normalization basis. `unit-conversion.ts` is *pre-compute* unit scaling
 * (P1 lives there, on the raw concentration); a post-compute CL-band check is
 * semantically misplaced upstream and the probe explicitly split P1/P2 on the
 * pre/post-compute boundary. One focused primitive, controller-wired (one
 * caller) — testable in isolation without standing up the full pipeline.
 *
 * **Correct layering (rule 5 — no vendored math).** It REUSES the engine's
 * `result.values.cl` (CL = Dose/AUCinf, computed by sci-comp); it does NOT
 * re-derive Dose or AUC. The only arithmetic here is the body-weight
 * normalization `cl / bodyWeightKg` and a band comparison.
 *
 * **Scientific contract (RATIFIED SF-implied-CL-band; per-kg extension RATIFIED
 * 2026-10-23, nca-cl-bw-normalised-units §0b):**
 *  - **Input contract: the ENGINE-basis CL** (architect B1) — `result.values.cl`
 *    with the engine's dose basis, never the reported mL/min/kg column. It runs at
 *    compute time on every row whose CL can be read per kg:
 *    `'absolute'` + body weight → `cl / bodyWeightKg` (L/h/kg); `'per_kg'` → the
 *    engine value as is (it already is L/h/kg).
 *  - **No double-normalization (unit-conversion.ts R2-F1 invariant).** By
 *    construction `doseNormalization === 'per_kg'` ⟺ the engine ran on a per-kg
 *    dose, so its `cl` is never divided by a weight again.
 *  - A ng/mL→mg/L conc mis-scale inflates AUC ~1000× → drives implied CL
 *    ~÷1000 → the **floor** is the primary trigger. The **ceiling** is a
 *    non-physical-high sanity bound (≤ rat total organ blood flow ~15 L/h/kg,
 *    Davies & Morris 1993, PMID 8378254).
 *  - **Known limitation (R2 N-7):** the floor was calibrated on small-molecule
 *    fixtures; a genuinely low-clearance compound (a monoclonal antibody at
 *    ~2e-4–5e-4 L/h/kg) sits below it. The message names that cause; a
 *    modality-aware band is REGISTRY `cl-plausibility-modality-band`.
 *  - Returns `null` (no fire) when `cl` is NaN / non-finite (failed profile), or
 *    the dose is absolute and there is no body weight. Inspect only.
 */

import {nca} from '@datagrok-libraries/sci-comp';

/**
 * Lower bound of the plausible implied per-kg clearance band (L/h/kg).
 *
 * **Value = 5e-3 L/h/kg, calibrated empirically against the REAL ENGINE output
 * over `files/03_rat_iv.csv` (rule 16 — build-established, not the
 * research approximation).** Research approximated the floor as "~1e-3 ... fires
 * 24/24" using a `dose/AUClast` hand-calc. The real engine (`dose/AUCinf`,
 * DEFAULT_RULES) gives a mis-scaled CL/kg band of **[2.97e-4, 1.011e-3]** — the
 * top end (subject R005, the highest-clearance rat, legit CL/kg ≈ 1.01) lands
 * just ABOVE 1e-3, so a 1e-3 floor fires only 23/24 and **misses R005 — a
 * 1 mg/kg subject that P1's ceiling also misses (a coverage hole exactly where
 * P2 is supposed to rescue P1).** The lowest legitimate band-pass profile is
 * theoph at **1.87e-2** L/h/kg. So 5e-3 sits in the real gap (1.011e-3, 1.87e-2):
 * **4.9× above the highest mis-scale (catches 24/24) and 3.7× below the lowest
 * legit profile (zero in-repo false positives).** **Documented residual:** a
 * genuinely ultra-low-clearance drug (CL/kg < ~5e-3 ≈ 0.15% of rat hepatic
 * blood flow) can collide with the floor — a non-blocking advisory, dismissible.
 * The value was moved 1e-3→5e-3 by calibration against the shipped fixtures. */
export const CL_PLAUSIBILITY_FLOOR_LHKG = 5e-3;

/**
 * Upper bound of the plausible implied per-kg clearance band (L/h/kg).
 *
 * **Value = 15 L/h/kg = rat total organ blood flow** (Davies & Morris 1993,
 * PMID 8378254) — clearance cannot physiologically exceed the blood flow
 * delivering drug to the eliminating organs. A non-physical-high sanity bound
 * (e.g. a reverse mg/L→ng/mL mis-scale, or a dose mislabel); all in-repo legit
 * profiles sit well below it (max ~1 L/h/kg). */
export const CL_PLAUSIBILITY_CEILING_LHKG = 15;

/**
 * Cross-check the engine's clearance against the physiological per-kg band and
 * return a non-blocking advisory when it falls outside, else `null`.
 *
 * @param cl                 the engine's `result.values.cl` — L/h on the
 *                           absolute-dose path, L/h/kg on the per-kg path; NaN for
 *                           a failed profile. NEVER the reported mL/min/kg value.
 * @param bodyWeightKg       canonical body weight (kg) the engine used; `null` when absent.
 * @param doseNormalization  the engine dose basis (`extractProfileInputs` / the composite's).
 */
export function impliedClPerKgWarning(
  cl: number,
  bodyWeightKg: number | null,
  doseNormalization: 'absolute' | 'per_kg',
): nca.ParameterWarning | null {
  if (!Number.isFinite(cl)) return null; // failed profile — never fabricate a warning
  let clPerKg: number;
  if (doseNormalization === 'per_kg')
    clPerKg = cl; // the engine ran per kg: already L/h/kg
  else {
    if (bodyWeightKg === null || !Number.isFinite(bodyWeightKg) || bodyWeightKg <= 0) return null;
    clPerKg = cl / bodyWeightKg;
  }
  if (clPerKg >= CL_PLAUSIBILITY_FLOOR_LHKG && clPerKg <= CL_PLAUSIBILITY_CEILING_LHKG)
    return null;

  const tooLow = clPerKg < CL_PLAUSIBILITY_FLOOR_LHKG;
  const cause = tooLow ?
    `implausibly low — a ng/mL→mg/L concentration mis-scale inflates AUC ~1000× ` +
    `and drives implied CL ~÷1000, silently corrupting CL / Vz / AUC — or a genuinely ` +
    `low-clearance compound, e.g. a biologic` :
    `implausibly high — above rat total organ blood flow (~${CL_PLAUSIBILITY_CEILING_LHKG} ` +
    `L/h/kg), which clearance cannot physiologically exceed`;
  return {
    code: 'CL_IMPLAUSIBLE_HEURISTIC_WARNING',
    severity: 'warning',
    message:
      `Implied clearance ${clPerKg.toPrecision(4)} L/h/kg is ${cause}. ` +
      `Plausible band [${CL_PLAUSIBILITY_FLOOR_LHKG}, ${CL_PLAUSIBILITY_CEILING_LHKG}] ` +
      `L/h/kg. Verify the concentration and dose unit pickers.`,
  };
}
