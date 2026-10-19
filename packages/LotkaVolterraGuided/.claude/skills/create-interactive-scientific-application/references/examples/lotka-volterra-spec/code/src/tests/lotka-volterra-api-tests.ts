// Lotka-Volterra — API / Validation tests

import {category, test, expect} from '@datagrok-libraries/utils/src/test';

import {DEFAULTS, validate} from '../lotka-volterra/core';

category('API: Validation', () => {
  // --- val_01: alpha <= 0 ---
  test('val_01: alpha = 0', async () => {
    const errors = validate({...DEFAULTS, alpha: 0});
    expect(errors.has('ctrl_alpha'), true, 'Should reject alpha = 0');
  });

  test('val_01: alpha = -1', async () => {
    const errors = validate({...DEFAULTS, alpha: -1});
    expect(errors.has('ctrl_alpha'), true, 'Should reject alpha < 0');
  });

  // --- val_02: beta <= 0 ---
  test('val_02: beta = 0', async () => {
    const errors = validate({...DEFAULTS, beta: 0});
    expect(errors.has('ctrl_beta'), true, 'Should reject beta = 0');
  });

  test('val_02: beta = -0.1', async () => {
    const errors = validate({...DEFAULTS, beta: -0.1});
    expect(errors.has('ctrl_beta'), true, 'Should reject beta < 0');
  });

  // --- val_03: delta <= 0 ---
  test('val_03: delta = 0', async () => {
    const errors = validate({...DEFAULTS, delta: 0});
    expect(errors.has('ctrl_delta'), true, 'Should reject delta = 0');
  });

  test('val_03: delta = -0.05', async () => {
    const errors = validate({...DEFAULTS, delta: -0.05});
    expect(errors.has('ctrl_delta'), true, 'Should reject delta < 0');
  });

  // --- val_04: gamma <= 0 ---
  test('val_04: gamma = 0', async () => {
    const errors = validate({...DEFAULTS, gamma: 0});
    expect(errors.has('ctrl_gamma'), true, 'Should reject gamma = 0');
  });

  test('val_04: gamma = -1', async () => {
    const errors = validate({...DEFAULTS, gamma: -1});
    expect(errors.has('ctrl_gamma'), true, 'Should reject gamma < 0');
  });

  // --- val_05: x0 <= 0 ---
  test('val_05: x0 = 0', async () => {
    const errors = validate({...DEFAULTS, x0: 0});
    expect(errors.has('ctrl_x0'), true, 'Should reject x0 = 0');
  });

  test('val_05: x0 = -5', async () => {
    const errors = validate({...DEFAULTS, x0: -5});
    expect(errors.has('ctrl_x0'), true, 'Should reject x0 < 0');
  });

  // --- val_06: y0 <= 0 ---
  test('val_06: y0 = 0', async () => {
    const errors = validate({...DEFAULTS, y0: 0});
    expect(errors.has('ctrl_y0'), true, 'Should reject y0 = 0');
  });

  test('val_06: y0 = -5', async () => {
    const errors = validate({...DEFAULTS, y0: -5});
    expect(errors.has('ctrl_y0'), true, 'Should reject y0 < 0');
  });

  // --- val_07: T <= 0 ---
  test('val_07: T = 0', async () => {
    const errors = validate({...DEFAULTS, T: 0});
    expect(errors.has('ctrl_T'), true, 'Should reject T = 0');
  });

  test('val_07: T = -10', async () => {
    const errors = validate({...DEFAULTS, T: -10});
    expect(errors.has('ctrl_T'), true, 'Should reject T < 0');
  });

  // --- valid defaults ---
  test('Defaults pass validation', async () => {
    const errors = validate(DEFAULTS);
    expect(errors.size, 0, 'Default parameters should be valid');
  });

  // --- valid boundary values ---
  test('alpha = 0.1 (valid boundary)', async () => {
    const errors = validate({...DEFAULTS, alpha: 0.1});
    expect(errors.has('ctrl_alpha'), false, 'Should accept alpha = 0.1');
  });

  test('x0 = 1 (valid boundary)', async () => {
    const errors = validate({...DEFAULTS, x0: 1});
    expect(errors.has('ctrl_x0'), false, 'Should accept x0 = 1');
  });

  // --- multiple errors ---
  test('Multiple simultaneous errors', async () => {
    const errors = validate({...DEFAULTS, alpha: 0, beta: 0, x0: 0, T: 0});
    expect(errors.size >= 4, true, 'Should report multiple errors');
    expect(errors.has('ctrl_alpha'), true);
    expect(errors.has('ctrl_beta'), true);
    expect(errors.has('ctrl_x0'), true);
    expect(errors.has('ctrl_T'), true);
  });
});
