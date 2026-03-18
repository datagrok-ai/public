// Golf Ball Flight — API / Validation tests

import {category, test, expect} from '@datagrok-libraries/utils/src/test';

import {DEFAULTS, validate} from '../golf-ball/core';

category('API: Validation', () => {
  // --- val_01: v0 <= 0 ---
  test('val_01: v0 = 0', async () => {
    const errors = validate({...DEFAULTS, v0: 0});
    expect(errors.has('ctrl_v0'), true, 'Should reject v0 = 0');
  });

  test('val_01: v0 = -10', async () => {
    const errors = validate({...DEFAULTS, v0: -10});
    expect(errors.has('ctrl_v0'), true, 'Should reject v0 < 0');
  });

  // --- val_02: theta <= 0 or theta >= 90 ---
  test('val_02: theta = 0', async () => {
    const errors = validate({...DEFAULTS, theta: 0});
    expect(errors.has('ctrl_theta'), true, 'Should reject theta = 0');
  });

  test('val_02: theta = 90', async () => {
    const errors = validate({...DEFAULTS, theta: 90});
    expect(errors.has('ctrl_theta'), true, 'Should reject theta = 90');
  });

  test('val_02: theta = -5', async () => {
    const errors = validate({...DEFAULTS, theta: -5});
    expect(errors.has('ctrl_theta'), true, 'Should reject theta < 0');
  });

  // --- val_03: m <= 0 ---
  test('val_03: m = 0', async () => {
    const errors = validate({...DEFAULTS, m: 0});
    expect(errors.has('ctrl_m'), true, 'Should reject m = 0');
  });

  test('val_03: m = -0.01', async () => {
    const errors = validate({...DEFAULTS, m: -0.01});
    expect(errors.has('ctrl_m'), true, 'Should reject m < 0');
  });

  // --- val_04: d <= 0 ---
  test('val_04: d = 0', async () => {
    const errors = validate({...DEFAULTS, d: 0});
    expect(errors.has('ctrl_d'), true, 'Should reject d = 0');
  });

  test('val_04: d = -0.02', async () => {
    const errors = validate({...DEFAULTS, d: -0.02});
    expect(errors.has('ctrl_d'), true, 'Should reject d < 0');
  });

  // --- val_05: Cd <= 0 ---
  test('val_05: Cd = 0', async () => {
    const errors = validate({...DEFAULTS, Cd: 0});
    expect(errors.has('ctrl_Cd'), true, 'Should reject Cd = 0');
  });

  test('val_05: Cd = -0.1', async () => {
    const errors = validate({...DEFAULTS, Cd: -0.1});
    expect(errors.has('ctrl_Cd'), true, 'Should reject Cd < 0');
  });

  // --- val_06: rho <= 0 ---
  test('val_06: rho = 0', async () => {
    const errors = validate({...DEFAULTS, rho: 0});
    expect(errors.has('ctrl_rho'), true, 'Should reject rho = 0');
  });

  test('val_06: rho = -1', async () => {
    const errors = validate({...DEFAULTS, rho: -1});
    expect(errors.has('ctrl_rho'), true, 'Should reject rho < 0');
  });

  // --- val_07: g <= 0 ---
  test('val_07: g = 0', async () => {
    const errors = validate({...DEFAULTS, g: 0});
    expect(errors.has('ctrl_g'), true, 'Should reject g = 0');
  });

  test('val_07: g = -9.81', async () => {
    const errors = validate({...DEFAULTS, g: -9.81});
    expect(errors.has('ctrl_g'), true, 'Should reject g < 0');
  });

  // --- valid defaults ---
  test('Defaults pass validation', async () => {
    const errors = validate(DEFAULTS);
    expect(errors.size, 0, 'Default parameters should be valid');
  });

  // --- valid boundary values ---
  test('v0 = 10 (valid boundary)', async () => {
    const errors = validate({...DEFAULTS, v0: 10});
    expect(errors.has('ctrl_v0'), false, 'Should accept v0 = 10');
  });

  test('theta = 1 (valid boundary)', async () => {
    const errors = validate({...DEFAULTS, theta: 1});
    expect(errors.has('ctrl_theta'), false, 'Should accept theta = 1');
  });

  test('theta = 89 (valid boundary)', async () => {
    const errors = validate({...DEFAULTS, theta: 89});
    expect(errors.has('ctrl_theta'), false, 'Should accept theta = 89');
  });

  // --- multiple errors ---
  test('Multiple simultaneous errors', async () => {
    const errors = validate({...DEFAULTS, v0: 0, theta: 0, m: 0, Cd: 0});
    expect(errors.size >= 4, true, 'Should report multiple errors');
    expect(errors.has('ctrl_v0'), true);
    expect(errors.has('ctrl_theta'), true);
    expect(errors.has('ctrl_m'), true);
    expect(errors.has('ctrl_Cd'), true);
  });
});
