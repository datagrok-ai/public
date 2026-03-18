// Multi-Compartment PK — Validation tests

import {category, test, expect} from '@datagrok-libraries/utils/src/test';

import {DEFAULTS, validate} from '../compartmental-pk/core';

category('API: Validation', () => {
  // --- val_01: cl <= 0 ---
  test('val_01: cl = 0', async () => {
    const errors = validate({...DEFAULTS, cl: 0});
    expect(errors.has('ctrl_cl'), true, 'Should reject cl = 0');
  });

  test('val_01: cl = -1', async () => {
    const errors = validate({...DEFAULTS, cl: -1});
    expect(errors.has('ctrl_cl'), true, 'Should reject cl < 0');
  });

  // --- val_02: vc <= 0 ---
  test('val_02: vc = 0', async () => {
    const errors = validate({...DEFAULTS, vc: 0});
    expect(errors.has('ctrl_vc'), true, 'Should reject vc = 0');
  });

  test('val_02: vc = -5', async () => {
    const errors = validate({...DEFAULTS, vc: -5});
    expect(errors.has('ctrl_vc'), true, 'Should reject vc < 0');
  });

  // --- val_03: vp <= 0 ---
  test('val_03: vp = 0', async () => {
    const errors = validate({...DEFAULTS, vp: 0});
    expect(errors.has('ctrl_vp'), true, 'Should reject vp = 0');
  });

  // --- val_04: q <= 0 ---
  test('val_04: q = 0', async () => {
    const errors = validate({...DEFAULTS, q: 0});
    expect(errors.has('ctrl_q'), true, 'Should reject q = 0');
  });

  // --- val_05: vd <= 0 (3-comp only) ---
  test('val_05: vd = 0 (3-comp)', async () => {
    const errors = validate({...DEFAULTS, modelType: '3-Compartment', vd: 0});
    expect(errors.has('ctrl_vd'), true, 'Should reject vd = 0 for 3-compartment');
  });

  test('val_05: vd = -1 (3-comp)', async () => {
    const errors = validate({...DEFAULTS, modelType: '3-Compartment', vd: -1});
    expect(errors.has('ctrl_vd'), true, 'Should reject vd < 0 for 3-compartment');
  });

  test('val_05: vd = 0 ignored for 2-comp', async () => {
    const errors = validate({...DEFAULTS, modelType: '2-Compartment', vd: 0});
    expect(errors.has('ctrl_vd'), false, 'Should not check vd for 2-compartment');
  });

  // --- val_06: q2 <= 0 (3-comp only) ---
  test('val_06: q2 = 0 (3-comp)', async () => {
    const errors = validate({...DEFAULTS, modelType: '3-Compartment', q2: 0});
    expect(errors.has('ctrl_q2'), true, 'Should reject q2 = 0 for 3-compartment');
  });

  // --- val_07: ka <= 0 (oral only) ---
  test('val_07: ka = 0 (oral)', async () => {
    const errors = validate({...DEFAULTS, inputMode: 'Oral', ka: 0});
    expect(errors.has('ctrl_ka'), true, 'Should reject ka = 0 for oral');
  });

  test('val_07: ka = 0 ignored for IV Bolus', async () => {
    const errors = validate({...DEFAULTS, inputMode: 'IV Bolus', ka: 0});
    expect(errors.has('ctrl_ka'), false, 'Should not check ka for IV Bolus');
  });

  // --- val_08: f <= 0 or f > 1 (oral only) ---
  test('val_08: f = 0 (oral)', async () => {
    const errors = validate({...DEFAULTS, inputMode: 'Oral', f: 0});
    expect(errors.has('ctrl_f'), true, 'Should reject f = 0');
  });

  test('val_08: f = 1.1 (oral)', async () => {
    const errors = validate({...DEFAULTS, inputMode: 'Oral', f: 1.1});
    expect(errors.has('ctrl_f'), true, 'Should reject f > 1');
  });

  test('val_08: f = 1.0 (oral, valid)', async () => {
    const errors = validate({...DEFAULTS, inputMode: 'Oral', f: 1.0});
    expect(errors.has('ctrl_f'), false, 'Should accept f = 1.0');
  });

  // --- val_09: dose <= 0 ---
  test('val_09: dose = 0', async () => {
    const errors = validate({...DEFAULTS, dose: 0});
    expect(errors.has('ctrl_dose'), true, 'Should reject dose = 0');
  });

  test('val_09: dose = -10', async () => {
    const errors = validate({...DEFAULTS, dose: -10});
    expect(errors.has('ctrl_dose'), true, 'Should reject dose < 0');
  });

  // --- val_10: tau <= 0 (repeated dosing) ---
  test('val_10: tau = 0 (repeated)', async () => {
    const errors = validate({...DEFAULTS, repeatedDosing: true, tau: 0});
    expect(errors.has('ctrl_tau'), true, 'Should reject tau = 0 for repeated dosing');
  });

  test('val_10: tau = 0 ignored when not repeated', async () => {
    const errors = validate({...DEFAULTS, repeatedDosing: false, tau: 0});
    expect(errors.has('ctrl_tau'), false, 'Should not check tau when not repeated');
  });

  // --- val_11: nDoses < 1 (repeated dosing) ---
  test('val_11: nDoses = 0 (repeated)', async () => {
    const errors = validate({...DEFAULTS, repeatedDosing: true, nDoses: 0});
    expect(errors.has('ctrl_n_doses'), true, 'Should reject nDoses = 0');
  });

  // --- val_12: tInf <= 0 (infusion) ---
  test('val_12: tInf = 0 (infusion)', async () => {
    const errors = validate({...DEFAULTS, inputMode: 'IV Infusion', tInf: 0});
    expect(errors.has('ctrl_t_inf'), true, 'Should reject tInf = 0 for infusion');
  });

  test('val_12: tInf = 0 ignored for bolus', async () => {
    const errors = validate({...DEFAULTS, inputMode: 'IV Bolus', tInf: 0});
    expect(errors.has('ctrl_t_inf'), false, 'Should not check tInf for bolus');
  });

  // --- val_13: tInf > tau (infusion + repeated) ---
  test('val_13: tInf > tau (infusion + repeated)', async () => {
    const errors = validate({
      ...DEFAULTS, inputMode: 'IV Infusion', repeatedDosing: true,
      tInf: 15, tau: 12,
    });
    expect(errors.has('ctrl_t_inf'), true, 'Should reject tInf > tau');
  });

  test('val_13: tInf = tau (valid)', async () => {
    const errors = validate({
      ...DEFAULTS, inputMode: 'IV Infusion', repeatedDosing: true,
      tInf: 12, tau: 12,
    });
    expect(errors.has('ctrl_t_inf'), false, 'Should accept tInf = tau');
  });

  // --- val_14: tEnd <= 0 ---
  test('val_14: tEnd = 0', async () => {
    const errors = validate({...DEFAULTS, tEnd: 0});
    expect(errors.has('ctrl_t_end'), true, 'Should reject tEnd = 0');
  });

  // --- val_15: mtc <= mec ---
  test('val_15: mtc = mec', async () => {
    const errors = validate({...DEFAULTS, mec: 5, mtc: 5});
    expect(errors.has('ctrl_mtc'), true, 'Should reject mtc = mec');
  });

  test('val_15: mtc < mec', async () => {
    const errors = validate({...DEFAULTS, mec: 10, mtc: 5});
    expect(errors.has('ctrl_mtc'), true, 'Should reject mtc < mec');
  });

  // --- val_16: mec < 0 ---
  test('val_16: mec = -1', async () => {
    const errors = validate({...DEFAULTS, mec: -1});
    expect(errors.has('ctrl_mec'), true, 'Should reject mec < 0');
  });

  // --- valid defaults ---
  test('Defaults pass validation', async () => {
    const errors = validate(DEFAULTS);
    expect(errors.size, 0, 'Default parameters should be valid');
  });

  // --- multiple simultaneous errors ---
  test('Multiple simultaneous errors', async () => {
    const errors = validate({...DEFAULTS, cl: 0, vc: 0, dose: 0});
    expect(errors.size >= 3, true, 'Should report multiple errors');
    expect(errors.has('ctrl_cl'), true);
    expect(errors.has('ctrl_vc'), true);
    expect(errors.has('ctrl_dose'), true);
  });
});
