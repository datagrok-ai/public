// CR3BP — Validation tests

import {category, test, expect} from '@datagrok-libraries/utils/src/test';

import {DEFAULTS, validate} from '../cr3bp/core';

category('API: Validation', () => {
  // --- val_01: mu <= 0 ---
  test('val_01: mu = 0', async () => {
    const errors = validate({...DEFAULTS, mu: 0});
    expect(errors.has('ctrl_mu'), true, 'Should reject mu = 0');
  });

  test('val_01: mu = -0.1', async () => {
    const errors = validate({...DEFAULTS, mu: -0.1});
    expect(errors.has('ctrl_mu'), true, 'Should reject mu < 0');
  });

  // --- val_02: mu > 0.5 ---
  test('val_02: mu = 0.6', async () => {
    const errors = validate({...DEFAULTS, mu: 0.6});
    expect(errors.has('ctrl_mu'), true, 'Should reject mu > 0.5');
  });

  test('val_02: mu = 1.0', async () => {
    const errors = validate({...DEFAULTS, mu: 1.0});
    expect(errors.has('ctrl_mu'), true, 'Should reject mu = 1.0');
  });

  // --- val_03: T <= 0 ---
  test('val_03: T = 0', async () => {
    const errors = validate({...DEFAULTS, T: 0});
    expect(errors.has('ctrl_T'), true, 'Should reject T = 0');
  });

  test('val_03: T = -10', async () => {
    const errors = validate({...DEFAULTS, T: -10});
    expect(errors.has('ctrl_T'), true, 'Should reject T < 0');
  });

  // --- valid defaults ---
  test('Defaults pass validation', async () => {
    const errors = validate(DEFAULTS);
    expect(errors.size, 0, 'Default parameters should be valid');
  });

  // --- valid boundary values ---
  test('mu = 0.001 (valid boundary)', async () => {
    const errors = validate({...DEFAULTS, mu: 0.001});
    expect(errors.has('ctrl_mu'), false, 'Should accept mu = 0.001');
  });

  test('mu = 0.5 (valid boundary)', async () => {
    const errors = validate({...DEFAULTS, mu: 0.5});
    expect(errors.has('ctrl_mu'), false, 'Should accept mu = 0.5');
  });

  test('T = 1 (valid boundary)', async () => {
    const errors = validate({...DEFAULTS, T: 1});
    expect(errors.has('ctrl_T'), false, 'Should accept T = 1');
  });

  // --- multiple errors ---
  test('Multiple simultaneous errors: mu = 0, T = 0', async () => {
    const errors = validate({...DEFAULTS, mu: 0, T: 0});
    expect(errors.size, 2, 'Should report 2 errors');
    expect(errors.has('ctrl_mu'), true);
    expect(errors.has('ctrl_T'), true);
  });
});
