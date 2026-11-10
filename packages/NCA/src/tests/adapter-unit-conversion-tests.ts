/**
 * Adapter — unit conversion category.
 *
 * @see src/adapter/unit-conversion.ts (implementation under test)
 *
 * Coverage:
 *  - Time / Concentration / Dose / BW / LLOQ rescaling — table correctness.
 *  - LLOQ same-factor invariant (AD-7).
 *  - Infusion duration same-factor + 2× heuristic (AD-10).
 *  - BW species-absent heuristic (AD-6 / R1-finding-2): mouse, rat, dog
 *    boundary, cat, lb known-limitation.
 *  - Identity at canonical (×1 fast path).
 *  - Float64 round-trip epsilon.
 *  - CANONICAL_UNIT_TAGS completeness (13 tagged + 10 untagged).
 */

import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {
  scaleProfileInputs,
  CANONICAL_UNIT_TAGS,
  TIME_TO_HOUR,
  CONC_TO_MG_L,
  DOSE_TO_MG,
  BW_TO_KG,
  convertConcentrationFactor,
  BW_HEURISTIC_THRESHOLD_KG,
  positiveMedian,
  CONC_PLAUSIBILITY_CEILING_MG_L,
  CONC_GUARD_MIN_POSITIVE_COUNT,
  normaliseClVolumes,
  compositeDoseWeight,
  NO_BODY_WEIGHT,
  L_PER_H_TO_ML_PER_MIN,
  type UnitScalars,
  type BodyWeightContext,
} from '../adapter/unit-conversion';

const EPS = 1e-12;

function scalars(timeUnits: 'h' | 'min' | 's' | 'd',
  concUnits: 'mg/L' | 'µg/mL' | 'µg/L' | 'ng/mL' | 'µmol/L' | 'nmol/L',
  doseUnits: 'mg' | 'µg' | 'ng' | 'mg/kg' | 'µg/kg',
): UnitScalars {
  return {timeUnits, concentrationUnits: concUnits, doseUnits};
}

function close(a: number, b: number, eps = EPS): boolean {
  return Math.abs(a - b) <= eps;
}

category('Adapter: unit conversion', () => {
  // ---------------------------------------------------------------
  // Conversion table correctness (one row per unit per axis)
  // ---------------------------------------------------------------

  test('TIME_TO_HOUR — 4 units pin to expected factors', async () => {
    expect(TIME_TO_HOUR['h'], 1);
    expect(close(TIME_TO_HOUR['min'], 1 / 60), true, `min got ${TIME_TO_HOUR['min']}`);
    expect(close(TIME_TO_HOUR['s'], 1 / 3600), true, `s got ${TIME_TO_HOUR['s']}`);
    expect(TIME_TO_HOUR['d'], 24);
  });

  test('CONC_TO_MG_L — 4 mass units pin to expected factors; molar → NaN', async () => {
    expect(CONC_TO_MG_L['mg/L'], 1);
    expect(CONC_TO_MG_L['µg/mL'], 1);
    expect(close(CONC_TO_MG_L['µg/L'], 1 / 1000), true);
    expect(close(CONC_TO_MG_L['ng/mL'], 1 / 1000), true);
    expect(Number.isNaN(CONC_TO_MG_L['µmol/L']), true, 'molar must not silently rescale');
    expect(Number.isNaN(CONC_TO_MG_L['nmol/L']), true, 'molar must not silently rescale');
  });

  test('DOSE_TO_MG — base factors for absolute / per-kg variants', async () => {
    expect(DOSE_TO_MG['mg'], 1);
    expect(close(DOSE_TO_MG['µg'], 1 / 1000), true);
    expect(close(DOSE_TO_MG['ng'], 1 / 1e6), true);
    expect(DOSE_TO_MG['mg/kg'], 1, 'per-kg base factor (multiplied by BW kg at call site)');
    expect(close(DOSE_TO_MG['µg/kg'], 1 / 1000), true);
  });

  test('BW_TO_KG — kg/g', async () => {
    expect(BW_TO_KG['kg'], 1);
    expect(close(BW_TO_KG['g'], 1 / 1000), true);
  });

  test('convertConcentrationFactor — wrapper agrees with table', async () => {
    expect(convertConcentrationFactor('mg/L'), 1);
    expect(close(convertConcentrationFactor('ng/mL'), 1 / 1000), true);
    expect(Number.isNaN(convertConcentrationFactor('µmol/L')), true);
  });

  // ---------------------------------------------------------------
  // scaleProfileInputs — rescaling correctness
  // ---------------------------------------------------------------

  test('time rescaling — min→h: 60 min → 1 h', async () => {
    const out = scaleProfileInputs(
      scalars('min', 'mg/L', 'mg'),
      Float64Array.from([0, 30, 60, 120]),
      Float64Array.from([10, 5, 2.5, 1.25]),
      100, 0, null, null, 'kg');
    expect(close(out.time[0], 0), true);
    expect(close(out.time[1], 0.5), true);
    expect(close(out.time[2], 1), true);
    expect(close(out.time[3], 2), true);
    expect(out.canonicalScalars.timeUnits, 'h');
  });

  test('concentration rescaling — ng/mL→mg/L: 1000 ng/mL → 1 mg/L', async () => {
    const out = scaleProfileInputs(
      scalars('h', 'ng/mL', 'mg'),
      Float64Array.from([0, 1, 2]),
      Float64Array.from([1000, 500, 250]),
      100, 0, null, null, 'kg');
    expect(close(out.conc[0], 1), true, `conc[0] got ${out.conc[0]}`);
    expect(close(out.conc[1], 0.5), true);
    expect(close(out.conc[2], 0.25), true);
    expect(out.canonicalScalars.concentrationUnits, 'mg/L');
  });

  test('dose rescaling — mg/kg + BW=70 kg → 70 mg per 1 mg/kg', async () => {
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'mg/kg'),
      Float64Array.from([0, 1, 2]),
      Float64Array.from([10, 5, 2.5]),
      1, 0, 70, null, 'kg');
    expect(close(out.dose, 70), true, `mg/kg × 70 kg should yield 70 mg; got ${out.dose}`);
    expect(out.canonicalScalars.doseUnits, 'mg');
  });

  test('dose rescaling — µg/kg + BW=70 kg → 0.07 mg per 1 µg/kg', async () => {
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'µg/kg'),
      Float64Array.from([0, 1, 2]),
      Float64Array.from([10, 5, 2.5]),
      1, 0, 70, null, 'kg');
    expect(close(out.dose, 0.07), true, `µg/kg × 70 kg should yield 0.07 mg; got ${out.dose}`);
  });

  // --------- nca-cl-vz-convention AD-2: per-kg dose path (no body weight) ---------

  test('AC-2: mg/kg + null BW → per-kg dose (no ×BW), doseNormalization per_kg, one PER_KG_NO_BW warning', async () => {
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'mg/kg'),
      Float64Array.from([0, 1, 2]),
      Float64Array.from([10, 5, 2.5]),
      1, 0, null, null, 'kg');
    expect(close(out.dose, 1), true, `mg/kg with no BW should stay per-kg (=1), not ×BW; got ${out.dose}`);
    expect(Number.isNaN(out.dose), false, 'no more NaN for the no-BW per-kg case');
    expect(out.doseNormalization, 'per_kg');
    const perKg = out.warnings.filter((w) => w.code === 'PER_KG_NO_BW');
    expect(perKg.length, 1, `expected exactly one PER_KG_NO_BW; got ${JSON.stringify(out.warnings)}`);
    expect(perKg[0].message.includes('per'), true, 'warning names the per-kg assumption');
  });

  test('AC-2 (µg/kg sibling): 1000 µg/kg + null BW → dose 1 mg, per_kg, one PER_KG_NO_BW', async () => {
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'µg/kg'),
      Float64Array.from([0, 1, 2]),
      Float64Array.from([10, 5, 2.5]),
      1000, 0, null, null, 'kg');
    expect(close(out.dose, 1), true, `1000 µg/kg × (1/1000), no ×BW → 1; got ${out.dose}`);
    expect(out.doseNormalization, 'per_kg');
    const perKg = out.warnings.filter((w) => w.code === 'PER_KG_NO_BW');
    expect(perKg.length, 1, `expected one PER_KG_NO_BW for µg/kg; got ${JSON.stringify(out.warnings)}`);
  });

  test('AC-3: mg/kg + BW=0.2441 kg → absolute dose 0.2441, doseNormalization absolute, NO PER_KG_NO_BW', async () => {
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'mg/kg'),
      Float64Array.from([0, 1, 2]),
      Float64Array.from([10, 5, 2.5]),
      1, 0, 0.2441, null, 'kg');
    expect(close(out.dose, 0.2441), true, `mg/kg × 0.2441 kg → 0.2441 mg; got ${out.dose}`);
    expect(out.doseNormalization, 'absolute');
    const perKg = out.warnings.filter((w) => w.code === 'PER_KG_NO_BW');
    expect(perKg.length, 0, `BW-present path must not warn per-kg; got ${JSON.stringify(out.warnings)}`);
  });

  test('AC-3 sibling: absolute dose unit (mg) → doseNormalization absolute, no PER_KG_NO_BW', async () => {
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'mg'),
      Float64Array.from([0, 1, 2]),
      Float64Array.from([10, 5, 2.5]),
      100, 0, null, null, 'kg');
    expect(out.doseNormalization, 'absolute');
    expect(out.warnings.filter((w) => w.code === 'PER_KG_NO_BW').length, 0);
  });

  test('BW rescaling — g→kg: 250 g → 0.25 kg', async () => {
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'mg'),
      Float64Array.from([0, 1]),
      Float64Array.from([10, 5]),
      100, 0, 250, null, 'g');
    expect(close(out.bodyWeight ?? NaN, 0.25), true, `250 g → 0.25 kg; got ${out.bodyWeight}`);
  });

  test('LLOQ same-factor invariant — conc unit ng/mL, lloq=2 → 0.002 mg/L', async () => {
    const out = scaleProfileInputs(
      scalars('h', 'ng/mL', 'mg'),
      Float64Array.from([0, 1, 2]),
      Float64Array.from([1000, 500, 250]),
      100, 2, null, null, 'kg');
    expect(typeof out.lloq, 'number');
    expect(close(out.lloq as number, 0.002), true, `LLOQ 2 ng/mL → 0.002 mg/L; got ${out.lloq}`);
  });

  test('LLOQ array — per-row rescaled by same factor as concentration', async () => {
    const out = scaleProfileInputs(
      scalars('h', 'µg/L', 'mg'),
      Float64Array.from([0, 1]),
      Float64Array.from([1000, 500]),
      100, Float64Array.from([100, 200]), null, null, 'kg');
    expect(out.lloq instanceof Float64Array, true);
    const arr = out.lloq as Float64Array;
    expect(close(arr[0], 0.1), true);
    expect(close(arr[1], 0.2), true);
  });

  test('infusion duration — same factor as time: min + 60 → 1 h', async () => {
    const out = scaleProfileInputs(
      scalars('min', 'mg/L', 'mg'),
      Float64Array.from([0, 30, 60, 120]),
      Float64Array.from([10, 5, 2.5, 1.25]),
      100, 0, null, 60, 'kg');
    expect(close(out.infusionDuration ?? NaN, 1), true,
      `60 min infusion → 1 h; got ${out.infusionDuration}`);
  });

  // ---------------------------------------------------------------
  // Heuristic warnings
  // ---------------------------------------------------------------

  test('BW heuristic — mouse 30 with bwUnit="kg" fires warning (threshold 15)', async () => {
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'mg/kg'),
      Float64Array.from([0, 1, 2]),
      Float64Array.from([10, 5, 2.5]),
      1, 0, 30, null, 'kg');
    expect(out.warnings.some((w) => w.code === 'BW_UNIT_HEURISTIC_WARNING'), true,
      `expected BW_UNIT_HEURISTIC_WARNING for mouse-as-kg case; got ${JSON.stringify(out.warnings)}`);
  });

  test('BW heuristic — rat 250 with bwUnit="kg" fires warning', async () => {
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'mg/kg'),
      Float64Array.from([0, 1, 2]),
      Float64Array.from([10, 5, 2.5]),
      1, 0, 250, null, 'kg');
    expect(out.warnings.some((w) => w.code === 'BW_UNIT_HEURISTIC_WARNING'), true,
      `expected BW_UNIT_HEURISTIC_WARNING for rat-as-kg case; got ${JSON.stringify(out.warnings)}`);
  });

  test('BW heuristic — beagle 15 kg AT threshold boundary does NOT fire (strict >)', async () => {
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'mg/kg'),
      Float64Array.from([0, 1, 2]),
      Float64Array.from([10, 5, 2.5]),
      1, 0, BW_HEURISTIC_THRESHOLD_KG, null, 'kg');
    expect(out.warnings.some((w) => w.code === 'BW_UNIT_HEURISTIC_WARNING'), false,
      `BW = 15 exactly must NOT fire (boundary); got ${JSON.stringify(out.warnings)}`);
  });

  test('BW heuristic — 15.001 kg just above threshold fires', async () => {
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'mg/kg'),
      Float64Array.from([0, 1, 2]),
      Float64Array.from([10, 5, 2.5]),
      1, 0, BW_HEURISTIC_THRESHOLD_KG + 0.001, null, 'kg');
    expect(out.warnings.some((w) => w.code === 'BW_UNIT_HEURISTIC_WARNING'), true,
      `BW = 15.001 must fire (strict-greater)`);
  });

  test('BW heuristic — cat 4 with bwUnit="kg" does NOT fire (below threshold)', async () => {
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'mg/kg'),
      Float64Array.from([0, 1, 2]),
      Float64Array.from([10, 5, 2.5]),
      1, 0, 4, null, 'kg');
    expect(out.warnings.some((w) => w.code === 'BW_UNIT_HEURISTIC_WARNING'), false,
      `cat BW = 4 must NOT fire`);
  });

  test('BW heuristic — bwUnit="g" with BW=30 does NOT fire (correct unit choice)', async () => {
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'mg/kg'),
      Float64Array.from([0, 1, 2]),
      Float64Array.from([10, 5, 2.5]),
      1, 0, 30, null, 'g');
    expect(out.warnings.some((w) => w.code === 'BW_UNIT_HEURISTIC_WARNING'), false,
      `bwUnit='g' must not trigger the kg-heuristic`);
  });

  test('BW heuristic — lb known limitation: 15-lb beagle with bwUnit="kg" at boundary does NOT fire', async () => {
    // A 15-lb beagle reported as BW=15 with bwUnit='kg' (wrong picker) is at
    // the threshold boundary — heuristic does NOT fire. Known v1 limitation
    // (AD-6 R2-3); resolution requires extending the picker to kg|g|lb,
    // forward-deferred.
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'mg/kg'),
      Float64Array.from([0, 1, 2]),
      Float64Array.from([10, 5, 2.5]),
      1, 0, 15, null, 'kg');
    expect(out.warnings.some((w) => w.code === 'BW_UNIT_HEURISTIC_WARNING'), false,
      `15-lb beagle at threshold boundary documented as v1 limitation`);
  });

  test('infusion-duration heuristic — duration > max(time)×2 fires', async () => {
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'mg'),
      Float64Array.from([0, 1, 2, 3, 4]),
      Float64Array.from([10, 5, 2.5, 1.25, 0.625]),
      100, 0, null, 10, 'kg'); // max(time)=4; 10 > 4×2=8 ⇒ fire
    expect(out.warnings.some((w) => w.code === 'INFUSION_DURATION_HEURISTIC_WARNING'), true,
      `expected INFUSION_DURATION_HEURISTIC_WARNING for 10h infusion vs 4h window`);
  });

  test('infusion-duration heuristic — duration = max(time)×2 exactly does NOT fire (strict >)', async () => {
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'mg'),
      Float64Array.from([0, 1, 2, 3, 4]),
      Float64Array.from([10, 5, 2.5, 1.25, 0.625]),
      100, 0, null, 8, 'kg'); // 8 == 4×2 ⇒ no fire
    expect(out.warnings.some((w) => w.code === 'INFUSION_DURATION_HEURISTIC_WARNING'), false,
      `duration at margin boundary must NOT fire (strict-greater)`);
  });

  test('infusion-duration heuristic — continuous infusion within 2× margin does NOT false-positive', async () => {
    // R2-finding-4 false-positive guard: a 24-hour continuous infusion across
    // a 20-hour observation window stays within the 2× margin (24 < 40).
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'mg'),
      Float64Array.from([0, 5, 10, 15, 20]),
      Float64Array.from([10, 5, 2.5, 1.25, 0.625]),
      100, 0, null, 24, 'kg');
    expect(out.warnings.some((w) => w.code === 'INFUSION_DURATION_HEURISTIC_WARNING'), false,
      `continuous infusion within 2× margin must not false-positive`);
  });

  // ---------------------------------------------------------------
  // Identity / round-trip
  // ---------------------------------------------------------------

  test('identity at canonical — h/mg/L/mg/kg input == numerically equal output', async () => {
    const time = Float64Array.from([0, 1, 2, 4, 8]);
    const conc = Float64Array.from([10, 5, 2.5, 0.625, 0.039]);
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'mg'),
      time, conc, 100, 0.01, 70, null, 'kg');
    for (let i = 0; i < time.length; i++) {
      expect(out.time[i], time[i], `time[${i}] must be identity at canonical`);
      expect(out.conc[i], conc[i], `conc[${i}] must be identity at canonical`);
    }
    expect(out.dose, 100);
    expect(out.bodyWeight, 70);
    expect(out.lloq, 0.01);
  });

  test('round-trip — conc ng/mL → mg/L → ng/mL within Float64 epsilon', async () => {
    const original = Float64Array.from([1, 10, 100, 1000, 10000]);
    const toMgL = scaleProfileInputs(
      scalars('h', 'ng/mL', 'mg'),
      Float64Array.from([0, 1, 2, 3, 4]),
      original, 100, 0, null, null, 'kg');
    // Manually reverse: mg/L * 1000 ≡ ng/mL
    for (let i = 0; i < original.length; i++) {
      const back = toMgL.conc[i] * 1000;
      expect(close(back, original[i], 1e-9), true,
        `round-trip conc[${i}]: ${original[i]} → ${toMgL.conc[i]} → ${back}`);
    }
  });

  // ---------------------------------------------------------------
  // CANONICAL_UNIT_TAGS completeness (R2-finding-3 + C0_EXTRAP + FR-200 moment
  // params + the per-kg family and its companions — 23 tagged; nca-cl-bw-normalised-units
  // AD-1 gave CL / VZ / VSS one basis per column, so they live in the map too).
  // ---------------------------------------------------------------

  test('CANONICAL_UNIT_TAGS has the 23 expected tagged columns', async () => {
    const expectedTagged: ReadonlyArray<readonly [string, string]> = [
      ['CMAX', 'mg/L'],
      ['C0_EXTRAP', 'mg/L'],
      ['TMAX', 'h'],
      ['AUC_LAST', 'mg/L·h'],
      ['AUC_INF', 'mg/L·h'],
      ['PCT_EXTRAP', '%'],
      ['LAMBDA_Z', '1/h'],
      ['HALF_LIFE', 'h'],
      // Per kg on every row (D-b) + the absolute companions + the weight.
      ['CL', 'mL/min/kg'],
      ['VZ', 'L/kg'],
      ['VSS', 'L/kg'],
      ['CL_ABS', 'L/h'],
      ['VZ_ABS', 'L'],
      ['VSS_ABS', 'L'],
      ['BW_KG', 'kg'],
      // FR-200 moment / lag params.
      ['AUMC_LAST', 'mg/L·h²'],
      ['AUMC_INF', 'mg/L·h²'],
      ['MRT', 'h'],
      ['TLAG', 'h'],
      ['PCT_AUMC_EXTRAP', '%'],
      ['LAMBDAZ_TFIRST', 'h'],
      ['LAMBDAZ_TLAST', 'h'],
      ['LAMBDAZ_INTERCEPT', 'mg/L'],
    ];
    for (const [name, tag] of expectedTagged)
      expect(CANONICAL_UNIT_TAGS[name], tag, `${name} tag should be ${tag}`);

    expect(Object.keys(CANONICAL_UNIT_TAGS).length, 23,
      `expected exactly 23 tagged columns; got ${Object.keys(CANONICAL_UNIT_TAGS).length}`);
  });

  // ---------------------------------------------------------------
  // P1 conc-mis-scale guard (nca-conc-dose-misscale-guard AD-1):
  // positiveMedian helper + ceiling check + unit-gate + message shape.
  // ---------------------------------------------------------------

  test('positiveMedian — odd count: middle of the ascending positive subset', async () => {
    expect(positiveMedian(Float64Array.from([5, 1, 3])), 3);
    expect(positiveMedian(Float64Array.from([3, 1, 5, 9, 7])), 5);
  });

  test('positiveMedian — even count: mean of the two middle positives', async () => {
    expect(positiveMedian(Float64Array.from([1, 2, 3, 4])), 2.5);
    expect(close(positiveMedian(Float64Array.from([10, 20, 30, 40, 50, 60]))!, 35), true);
  });

  test('positiveMedian — excludes NaN, ±Inf, 0 and negatives (peak-biased)', async () => {
    // Positives are [2, 4, 6] → median 4; the zeros / NaN / Inf / negative drop out.
    const m = positiveMedian(Float64Array.from([0, 2, NaN, 4, -5, Infinity, 6, -Infinity, 0]));
    expect(m, 4, `median over positive subset {2,4,6} = 4; got ${m}`);
  });

  test('positiveMedian — fewer than the count floor positive values → null', async () => {
    // Two positives (floor is 3) → null; a "median" of 1-2 points is one value (R1-F1).
    expect(positiveMedian(Float64Array.from([100, 200])), null);
    expect(positiveMedian(Float64Array.from([0, 0, 500, NaN, 0])), null, 'single positive → null');
    expect(CONC_GUARD_MIN_POSITIVE_COUNT, 3, 'structural count floor is 3');
  });

  test('positiveMedian — all-zero / all-BLQ / empty → null (no positive values)', async () => {
    expect(positiveMedian(Float64Array.from([0, 0, 0, 0, 0])), null);
    expect(positiveMedian(new Float64Array(0)), null);
    expect(positiveMedian(Float64Array.from([NaN, -1, 0, -2])), null);
  });

  test('positiveMedian — exactly at the count floor (3 positives) computes', async () => {
    expect(positiveMedian(Float64Array.from([600, 700, 800])), 700,
      'three positives is the minimum count where a median is not one value');
  });

  test('P1 (AC-2) — conc median above the ceiling fires CONC_MAGNITUDE_HEURISTIC_WARNING', async () => {
    // Median 700 mg/L (> 500 ceiling), canonical-default mg/L → fire.
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'mg'),
      Float64Array.from([0, 1, 2, 4, 8]),
      Float64Array.from([900, 800, 700, 600, 500]),
      100, 0, null, null, 'kg');
    expect(out.warnings.some((w) => w.code === 'CONC_MAGNITUDE_HEURISTIC_WARNING'), true,
      `expected CONC_MAGNITUDE_HEURISTIC_WARNING for median 700 > ${CONC_PLAUSIBILITY_CEILING_MG_L}; ` +
      `got ${JSON.stringify(out.warnings)}`);
  });

  test('P1 (AC-1) — conc median below the ceiling does NOT fire (legit high-dose band)', async () => {
    // Median 298 mg/L (the in-repo legit max, 06_tk D300) < 500 ceiling → no fire.
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'mg'),
      Float64Array.from([0, 1, 2, 4, 8]),
      Float64Array.from([400, 350, 298, 250, 100]),
      100, 0, null, null, 'kg');
    expect(out.warnings.some((w) => w.code === 'CONC_MAGNITUDE_HEURISTIC_WARNING'), false,
      `median 298 (legit high-dose) must NOT fire; got ${JSON.stringify(out.warnings)}`);
  });

  test('P1 (AC-1) — at the ceiling exactly does NOT fire (strict greater)', async () => {
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'mg'),
      Float64Array.from([0, 1, 2]),
      Float64Array.from([CONC_PLAUSIBILITY_CEILING_MG_L, CONC_PLAUSIBILITY_CEILING_MG_L,
        CONC_PLAUSIBILITY_CEILING_MG_L]),
      100, 0, null, null, 'kg');
    expect(out.warnings.some((w) => w.code === 'CONC_MAGNITUDE_HEURISTIC_WARNING'), false,
      `median == ceiling must NOT fire (strict >)`);
  });

  test('P1 (AC-4) — non-canonical unit (ng/mL) NEVER fires even at huge raw magnitude', async () => {
    // The user CORRECTLY labelled ng/mL: rescaling already ÷1000, the high
    // magnitude is gone (median → 0.7 mg/L). The load-bearing FP-avoidance gate.
    const out = scaleProfileInputs(
      scalars('h', 'ng/mL', 'mg'),
      Float64Array.from([0, 1, 2, 4, 8]),
      Float64Array.from([900, 800, 700, 600, 500]),
      100, 0, null, null, 'kg');
    expect(out.warnings.some((w) => w.code === 'CONC_MAGNITUDE_HEURISTIC_WARNING'), false,
      `correctly-labelled ng/mL must NOT fire P1 (concentrationUnits !== 'mg/L' gate); ` +
      `got ${JSON.stringify(out.warnings)}`);
  });

  test('P1 (AC-5) — all-BLQ / all-zero profile does NOT fire (positiveMedian null)', async () => {
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'mg'),
      Float64Array.from([0, 1, 2, 4]),
      Float64Array.from([0, 0, 0, 0]),
      100, 0, null, null, 'kg');
    expect(out.warnings.some((w) => w.code === 'CONC_MAGNITUDE_HEURISTIC_WARNING'), false,
      `all-zero conc → no positive median → no fire`);
  });

  test('P1 (AC-5) — fewer than 3 positive points does NOT fire (count floor)', async () => {
    // Two huge positives (2000) but below the floor → guard does not fire.
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'mg'),
      Float64Array.from([0, 1, 2, 4]),
      Float64Array.from([2000, 2000, 0, 0]),
      100, 0, null, null, 'kg');
    expect(out.warnings.some((w) => w.code === 'CONC_MAGNITUDE_HEURISTIC_WARNING'), false,
      `2 positive points (< floor 3) must NOT fire even above the ceiling`);
  });

  test('P1 (AC-5) — >50% BLQ but ≥3 positives above ceiling DOES fire (median over positives)', async () => {
    // 6 points, 3 BLQ zeros (>50%) but the 3 measurable positives are 700/800/900
    // → positiveMedian excludes the zeros → median 800 > 500 → fire. The fire-side
    // of the BLQ edge case: a peak-biased median is NOT diluted by sub-LLOQ zeros.
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'mg'),
      Float64Array.from([0, 1, 2, 4, 8, 12]),
      Float64Array.from([0, 0, 0, 700, 800, 900]),
      100, 0, null, null, 'kg');
    expect(out.warnings.some((w) => w.code === 'CONC_MAGNITUDE_HEURISTIC_WARNING'), true,
      `>50% BLQ but 3 positives (median 800 > ceiling) must fire — positive-only median`);
  });

  test('P1 (AC-6) — message names the median, the 1000× error, and the CL/Vz/AUC consequence', async () => {
    const out = scaleProfileInputs(
      scalars('h', 'mg/L', 'mg'),
      Float64Array.from([0, 1, 2, 4, 8]),
      Float64Array.from([900, 800, 700, 600, 500]),
      100, 0, null, null, 'kg');
    const w = out.warnings.find((x) => x.code === 'CONC_MAGNITUDE_HEURISTIC_WARNING');
    expect(w !== undefined, true, 'warning must be present');
    if (w === undefined) return;
    expect(w.severity, 'warning', 'non-blocking severity (BW precedent)');
    expect(w.message.includes('700'), true, `message names the observed median; got "${w.message}"`);
    expect(w.message.includes('1000'), true, 'message names the likely 1000× error');
    expect(w.message.includes('CL'), true, 'message names the CL/Vz/AUC consequence');
  });

  test('CANONICAL_UNIT_TAGS does NOT carry tags for the untagged columns', async () => {
    const untagged = ['SUBJECT', 'ROUTE', 'STATUS', 'LAMBDAZ_R2', 'LAMBDAZ_ADJ_R2',
      'LAMBDAZ_NPOINTS', 'AUC_METHOD', 'COMPENSATED', 'WARNINGS_JSON', 'PROVENANCE_JSON',
      'DOSE_NORM', 'CL_BASIS'];
    for (const name of untagged)
      expect(CANONICAL_UNIT_TAGS[name] === undefined, true, `${name} must NOT have a tag`);
  });
});

// ---------------------------------------------------------------------------
// nca-cl-bw-normalised-units F1 — the ONE per-kg formation rule (AD-1 / AD-2) and the shared
// composite dose/weight helper (AD-5b). Expected values are hand constants or computed from the
// inputs in the test — never through the function under test (AD-10 REGRESSION_PIN / IDENTITY).
// ---------------------------------------------------------------------------

const bwCtx = (kg: number | null, range: [number, number] | null = null, perKgDose = false): BodyWeightContext =>
  ({bodyWeightKg: kg, bodyWeightRangeKg: range === null ? null : {min: range[0], max: range[1]}, perKgDose});

category('Adapter: per-kg normalisation (normaliseClVolumes)', () => {
  // R001 (03_rat_iv, PKNCA 0.12.1): absolute CL 0.0754969 L/h on 0.2441 kg; Vz 0.919500 L.
  const r001 = {cl: 0.0754969, vz: 0.9195, vss: 0.6};

  test('AC-1.2: the scale factor is 1000/60 and volumes are unscaled', async () => {
    expect(close(L_PER_H_TO_ML_PER_MIN, 1000 / 60), true, `factor ${L_PER_H_TO_ML_PER_MIN}`);
    const n = normaliseClVolumes({cl: 1, vz: 2, vss: 3}, 'per_kg', NO_BODY_WEIGHT);
    expect(close(n.cl, 1000 / 60), true, 'per-kg L/h/kg → mL/min/kg');
    expect(n.vz, 2, 'Vz per kg unscaled'); expect(n.vss, 3, 'Vss per kg unscaled');
  });

  test('AD-1 row 1: per_kg engine, no BW → engine value per kg; companions + BW NaN; no warning', async () => {
    const n = normaliseClVolumes({cl: 0.30929, vz: 3.7669, vss: 2}, 'per_kg', NO_BODY_WEIGHT);
    expect(close(n.cl, 0.30929 * 1000 / 60, 1e-12), true, `CL ${n.cl}`);
    expect(n.vz, 3.7669); expect(n.vss, 2);
    expect(Number.isNaN(n.clAbs) && Number.isNaN(n.vzAbs) && Number.isNaN(n.vssAbs), true, 'companions NaN');
    expect(Number.isNaN(n.bwKg), true, 'BW NaN');
    expect(n.warnings.length, 0, 'no per-kg note on the per-kg path (PER_KG_NO_BW is the scaler\'s)');
  });

  test('AD-1 row 2: absolute (mg/kg × BW) → abs / BW; R001 5.1548 mL/min/kg (REGRESSION_PIN)', async () => {
    const n = normaliseClVolumes(r001, 'absolute', bwCtx(0.2441, [0.2441, 0.2441], true));
    expect(close(n.cl, 5.1548, 5e-4), true, `R001 CL per kg = 5.1548 mL/min/kg; got ${n.cl}`);
    expect(close(n.vz, 3.7669, 5e-4), true, `R001 Vz per kg = 3.7669 L/kg; got ${n.vz}`);
    expect(n.clAbs, r001.cl, 'CL_ABS = engine'); expect(n.vzAbs, r001.vz); expect(n.vssAbs, r001.vss);
    expect(n.bwKg, 0.2441, 'BW_KG = the weight used'); expect(n.warnings.length, 0, 'no note');
  });

  test('AD-1 row 3: absolute-mg dose, constant BW → abs / BW', async () => {
    const n = normaliseClVolumes(r001, 'absolute', bwCtx(0.25, [0.25, 0.25], false));
    expect(close(n.cl, r001.cl / 0.25 * 1000 / 60, 1e-12), true, `CL ${n.cl}`);
    expect(n.warnings.length, 0, 'constant weight → no note');
  });

  test('AD-2: absolute-mg dose, BW varies within the profile → per-kg refused + PER_KG_NOT_FORMED', async () => {
    const n = normaliseClVolumes(r001, 'absolute', bwCtx(0.25, [0.25, 0.27], false));
    expect(Number.isNaN(n.cl) && Number.isNaN(n.vz) && Number.isNaN(n.vss), true, 'per kg NaN');
    expect(n.clAbs, r001.cl, 'absolute kept'); expect(n.bwKg, 0.25, 'first-row BW');
    expect(n.warnings.length, 1, 'one note');
    expect(n.warnings[0].code, 'PER_KG_NOT_FORMED', 'refusal code');
    expect(n.warnings[0].message.includes('0.25–0.27 kg'), true, `names the range: ${n.warnings[0].message}`);
  });

  test('AD-2 (R2 N-4): mg/kg dose, BW varies → per kg exact + BW_VARIES_WITHIN_PROFILE info', async () => {
    const n = normaliseClVolumes(r001, 'absolute', bwCtx(0.25, [0.25, 0.27], true));
    expect(close(n.cl, r001.cl / 0.25 * 1000 / 60, 1e-12), true, 'per kg formed (weight cancels)');
    expect(n.warnings.length, 1, 'one note');
    expect(n.warnings[0].code, 'BW_VARIES_WITHIN_PROFILE', 'companion note');
    expect(n.warnings[0].severity, 'info', 'information, not a warning');
  });

  test('AD-1 row 5: absolute dose, no BW → per kg NaN + PER_KG_NOT_FORMED; absolute kept', async () => {
    const n = normaliseClVolumes(r001, 'absolute', NO_BODY_WEIGHT);
    expect(Number.isNaN(n.cl), true, 'per kg NaN'); expect(n.clAbs, r001.cl, 'absolute kept');
    expect(Number.isNaN(n.bwKg), true, 'no BW');
    expect(n.warnings.map((w) => w.code).join(), 'PER_KG_NOT_FORMED', 'one refusal');
  });

  test('AC-1.1: BW ≤ 0 / non-finite is treated as absent; NaN engine values stay NaN', async () => {
    for (const bad of [0, -0.25, NaN, Infinity]) {
      const n = normaliseClVolumes(r001, 'absolute', bwCtx(bad));
      expect(Number.isNaN(n.cl), true, `BW ${bad} → per kg NaN`);
      expect(n.warnings[0]?.code, 'PER_KG_NOT_FORMED', `BW ${bad} → refusal`);
    }
    const nan = normaliseClVolumes({cl: NaN, vz: NaN, vss: NaN}, 'absolute', bwCtx(0.25, [0.25, 0.25]));
    expect(Number.isNaN(nan.cl) && Number.isNaN(nan.clAbs), true, 'failed profile → both NaN');
  });

  test('AD-5: per_kg engine with a group-mean BW (sparse /kg composite) → companions = per kg × BW', async () => {
    const n = normaliseClVolumes({cl: 0.3, vz: 2, vss: 1.5}, 'per_kg', bwCtx(0.025, null, true));
    expect(close(n.clAbs, 0.3 * 0.025, 1e-15), true, `CL_ABS ${n.clAbs}`);
    expect(close(n.vzAbs, 2 * 0.025, 1e-15), true, 'VZ_ABS');
    expect(n.bwKg, 0.025, 'BW_KG = group mean');
    expect(close(n.cl, 0.3 * 1000 / 60, 1e-12), true, 'per kg = engine × scale');
  });
});

category('Adapter: composite dose / weight (compositeDoseWeight, AD-5b)', () => {
  test('weight-based cell: mean dose / mean BW = dose per kg exactly; no note', async () => {
    const bw = [0.2, 0.25, 0.3];
    const r = compositeDoseWeight(bw.map((b) => ({dose: 10 * b, bw: b})));
    expect(close(r.meanDose, 2.5, 1e-12), true, `mean dose ${r.meanDose}`);
    expect(close(r.perKgDivisor!, 0.25, 1e-12), true, 'divisor = mean BW');
    expect(close(r.meanDose / r.perKgDivisor!, 10, 1e-12), true, 'mean dose / divisor = 10 mg/kg');
    expect(r.warnings.length, 0, 'weight-based, every member has a BW');
  });

  test('rounded weight-based doses (0.01 mg on a 0.223 kg rat) are not flagged (R2 N-3)', async () => {
    // 10 mg/kg × {0.2234, 0.2468, 0.2615} kg = 2.234 / 2.468 / 2.615 mg, recorded to 0.01 mg:
    // 2.23 / 2.47 / 2.62 → dose/BW 9.982 / 10.008 / 10.019, a spread only the rounding explains.
    const r = compositeDoseWeight([{dose: 2.23, bw: 0.2234}, {dose: 2.47, bw: 0.2468}, {dose: 2.62, bw: 0.2615}]);
    expect(r.warnings.some((w) => w.code === 'COMPOSITE_NOT_WEIGHT_BASED'), false, 'within recorded precision');
  });

  test('fixed doses on different weights → COMPOSITE_NOT_WEIGHT_BASED', async () => {
    // 2.5 mg to every animal, weights recorded to 0.1 g: dose/BW 11.7 vs 8.0, far beyond rounding.
    const r = compositeDoseWeight([{dose: 2.5, bw: 0.2134}, {dose: 2.5, bw: 0.3122}]);
    expect(r.warnings.some((w) => w.code === 'COMPOSITE_NOT_WEIGHT_BASED'), true, 'dose not ∝ BW');
  });

  test('partial BW → COMPOSITE_BW_PARTIAL; divisor gives per kg over the SAME members (F-5)', async () => {
    // Members with a BW: doses 2 and 3 on 0.2 and 0.3 kg (10 mg/kg); a third dose 4 has no BW.
    const r = compositeDoseWeight([{dose: 2, bw: 0.2}, {dose: 3, bw: 0.3}, {dose: 4, bw: NaN}]);
    expect(r.nBw, 2, 'two members with a weight'); expect(r.n, 3, 'three with a dose');
    expect(r.warnings.some((w) => w.code === 'COMPOSITE_BW_PARTIAL' && w.message.includes('2 of 3')), true,
      'names n_BW of n');
    // The engine's dose is the mean over ALL (3); per kg must be mean(2,3)/mean(0.2,0.3) = 10.
    expect(close(r.meanDose / r.perKgDivisor!, 10, 1e-12), true, `per-kg dose ${r.meanDose / r.perKgDivisor!}`);
  });

  test('no member with a BW → divisor null (per kg not formed); zero-dose cell → no weight-based note', async () => {
    expect(compositeDoseWeight([{dose: 2, bw: NaN}]).perKgDivisor, null, 'no BW → null');
    const zero = compositeDoseWeight([{dose: 0, bw: 0.2}, {dose: 0, bw: 0.3}]);
    expect(zero.warnings.length, 0, 'a control cell is not judged');
  });
});
