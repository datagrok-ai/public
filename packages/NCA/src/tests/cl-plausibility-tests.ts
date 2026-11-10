/**
 * Adapter — implied-CL plausibility (P2 conc-mis-scale guard) unit category.
 *
 * @see src/adapter/cl-plausibility.ts (implementation under test)
 *
 * Coverage (AC-9 — pure band logic with edge cases):
 *  - absolute dose + BW absent → null; CL NaN → null.
 *  - per_kg engine CL read as is, in or out of the band (S1 extension, nca-cl-bw-normalised-units).
 *  - In-band → null; below floor → warning (the conc-mis-scale primary trigger);
 *    above ceiling → warning (non-physical sanity bound).
 *  - Boundary: exactly at floor / ceiling does NOT fire (inclusive band).
 *  - Message shape: names the implied CL/kg + the band + verify-pickers prompt.
 */

import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {
  impliedClPerKgWarning,
  CL_PLAUSIBILITY_FLOOR_LHKG,
  CL_PLAUSIBILITY_CEILING_LHKG,
} from '../adapter/cl-plausibility';

category('Adapter: implied-CL plausibility (P2)', () => {
  // --- Gate: paths that must NEVER fire (AC-9) -------------------------

  test('AC-9 — body weight absent → null (cannot normalize to per-kg)', async () => {
    // Even a wildly-low absolute CL cannot fire without BW to per-kg-normalize.
    expect(impliedClPerKgWarning(1e-9, null, 'absolute'), null);
  });

  // S1 extension (nca-cl-bw-normalised-units §0b, ratified 2026-10-23; phase-0 corpus check: 0
  // newly flagged in-repo profiles). A per_kg engine CL is already L/h/kg — read as is, never ÷ BW.
  test('S1 — per_kg engine CL outside the band fires (no BW needed)', async () => {
    const w = impliedClPerKgWarning(3e-4, null, 'per_kg');
    expect(w !== null, true, 'a per-kg 3e-4 L/h/kg is below the 5e-3 floor');
    expect(w?.code, 'CL_IMPLAUSIBLE_HEURISTIC_WARNING', 'same code as the absolute path');
    expect(w?.message.includes('biologic') ?? false, true, 'the low-clearance cause is named (R2 N-7)');
  });

  test('S1 — per_kg engine CL inside the band does not fire, and is NOT divided by a BW', async () => {
    expect(impliedClPerKgWarning(0.30929, null, 'per_kg'), null, 'R001 per kg 0.309 L/h/kg in band');
    // 0.0619 L/h/kg ÷ 0.25 would be 0.25 (still in band) — use a value the division would push out:
    // 4 L/h/kg is in band; ÷ 0.25 = 16 would exceed the 15 ceiling. In band ⇒ not divided.
    expect(impliedClPerKgWarning(4, 0.25, 'per_kg'), null, 'per-kg CL read as is (no double normalisation)');
  });

  test('AC-9 — CL NaN (failed profile) → null (never fabricate a warning)', async () => {
    expect(impliedClPerKgWarning(NaN, 0.25, 'absolute'), null);
    expect(impliedClPerKgWarning(Infinity, 0.25, 'absolute'), null);
  });

  test('AC-9 — non-finite / non-positive body weight → null', async () => {
    expect(impliedClPerKgWarning(0.0003, 0, 'absolute'), null);
    expect(impliedClPerKgWarning(0.0003, NaN, 'absolute'), null);
    expect(impliedClPerKgWarning(0.0003, -0.25, 'absolute'), null);
  });

  // --- In-band: no fire (the legit-profile case) ----------------------

  test('in-band — legit rat CL 0.31 L/h/kg → null (no fire)', async () => {
    // R001 per-kg CL = 0.309 L/h/kg: cl 0.0755 L/h ÷ BW 0.2441 kg.
    expect(impliedClPerKgWarning(0.309 * 0.2441, 0.2441, 'absolute'), null);
  });

  test('in-band — at the floor exactly → null (inclusive band)', async () => {
    // cl_per_kg == floor → no fire.
    const cl = CL_PLAUSIBILITY_FLOOR_LHKG * 0.25;
    expect(impliedClPerKgWarning(cl, 0.25, 'absolute'), null);
  });

  test('in-band — at the ceiling exactly → null (inclusive band)', async () => {
    const cl = CL_PLAUSIBILITY_CEILING_LHKG * 0.25;
    expect(impliedClPerKgWarning(cl, 0.25, 'absolute'), null);
  });

  // --- Out of band: fire ----------------------------------------------

  test('AC-7 — below floor (conc mis-scale signature) fires CL_IMPLAUSIBLE_HEURISTIC_WARNING', async () => {
    // Mis-scaled R001: cl_per_kg ≈ 3.1e-4 < 1e-3 floor.
    const clPerKg = 3.1e-4;
    const w = impliedClPerKgWarning(clPerKg * 0.2441, 0.2441, 'absolute');
    expect(w !== null, true, 'below-floor implied CL must fire');
    if (w === null) return;
    expect(w.code, 'CL_IMPLAUSIBLE_HEURISTIC_WARNING');
    expect(w.severity, 'warning', 'non-blocking (BW precedent)');
    expect(w.message.includes('1000'), true, 'low-CL message names the ~1000× AUC inflation');
    expect(w.message.includes('CL'), true, 'message names the CL/Vz/AUC consequence');
  });

  test('above ceiling (non-physical) fires the warning naming the blood-flow bound', async () => {
    // cl_per_kg = 20 L/h/kg > 15 (rat blood flow).
    const w = impliedClPerKgWarning(20 * 0.25, 0.25, 'absolute');
    expect(w !== null, true, 'above-ceiling implied CL must fire');
    if (w === null) return;
    expect(w.code, 'CL_IMPLAUSIBLE_HEURISTIC_WARNING');
    expect(w.message.includes(String(CL_PLAUSIBILITY_CEILING_LHKG)), true,
      'high-CL message names the blood-flow ceiling');
  });

  test('calibrated band constants are the build-established values (rule 16 pins)', async () => {
    // Floor 5e-3: real-engine calibration (mis-scale max R005=1.011e-3, theoph
    // legit min=1.87e-2) — fires 24/24, 0 legit FP. NOT the 1e-3 research approx.
    expect(CL_PLAUSIBILITY_FLOOR_LHKG, 5e-3, 'floor = 5e-3 L/h/kg (build-calibrated, fires 24/24)');
    expect(CL_PLAUSIBILITY_CEILING_LHKG, 15, 'ceiling = rat total organ blood flow (PMID 8378254)');
  });
});
