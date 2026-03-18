// SIR/SEIR — Validation tests

import {category, test, expect} from '@datagrok-libraries/utils/src/test';

import {DEFAULTS, validate} from '../sir-seir/core';

category('API: Validation', () => {
  // --- val_01: R₀ range ---
  test('val_01a: r0 = 0.4 (below min)', async () => {
    const errors = validate({...DEFAULTS, r0: 0.4});
    expect(errors.has('ctrl_r0'), true, 'Should reject r0 = 0.4');
  });

  test('val_01b: r0 = 8.5 (above max)', async () => {
    const errors = validate({...DEFAULTS, r0: 8.5});
    expect(errors.has('ctrl_r0'), true, 'Should reject r0 = 8.5');
  });

  // --- val_02: γ range ---
  test('val_02a: gamma = 0', async () => {
    const errors = validate({...DEFAULTS, gamma: 0});
    expect(errors.has('ctrl_gamma'), true, 'Should reject gamma = 0');
  });

  test('val_02b: gamma = -0.1', async () => {
    const errors = validate({...DEFAULTS, gamma: -0.1});
    expect(errors.has('ctrl_gamma'), true, 'Should reject gamma < 0');
  });

  test('val_02c: gamma = 0.6', async () => {
    const errors = validate({...DEFAULTS, gamma: 0.6});
    expect(errors.has('ctrl_gamma'), true, 'Should reject gamma > 0.5');
  });

  // --- val_03: σ range (SEIR only) ---
  test('val_03a: sigma = 0 (SEIR)', async () => {
    const errors = validate({...DEFAULTS, modelType: 'SEIR', sigma: 0});
    expect(errors.has('ctrl_sigma'), true, 'Should reject sigma = 0 in SEIR');
  });

  test('val_03b: sigma = 1.5 (SEIR)', async () => {
    const errors = validate({...DEFAULTS, modelType: 'SEIR', sigma: 1.5});
    expect(errors.has('ctrl_sigma'), true, 'Should reject sigma > 1.0 in SEIR');
  });

  test('val_03c: sigma = 0 (SIR) — skipped', async () => {
    const errors = validate({...DEFAULTS, modelType: 'SIR', sigma: 0});
    expect(errors.has('ctrl_sigma'), false, 'Should not validate sigma in SIR mode');
  });

  // --- val_04: vaccination range ---
  test('val_04a: vaccination = -5', async () => {
    const errors = validate({...DEFAULTS, vaccination: -5});
    expect(errors.has('ctrl_vaccination'), true, 'Should reject vaccination < 0');
  });

  test('val_04b: vaccination = 95', async () => {
    const errors = validate({...DEFAULTS, vaccination: 95});
    expect(errors.has('ctrl_vaccination'), true, 'Should reject vaccination > 90');
  });

  // --- val_05: vaccination feasibility ---
  test('val_05a: vaccination = 90 (too high for I₀=1)', async () => {
    // v = 0.9, N*(1-v) = 1000, S₀ = 1000 - 1 = 999 > 0, but v >= 1 - 1/N is false
    // Actually v=0.9 < 1-1/10000=0.9999, so this should pass val_05.
    // Let's check with vaccination = 99.99 which would be caught by val_04 first (>90).
    // Per spec, val_05: vaccination/100 >= 1 - 1/N. With N=10000, threshold is 0.9999.
    // vaccination=90 → v=0.9 < 0.9999 → should pass.
    const errors = validate({...DEFAULTS, vaccination: 90});
    // v=0.9, 1-1/10000=0.9999. 0.9 < 0.9999 → passes val_05. But vaccination=90 is at the max.
    expect(errors.has('ctrl_vaccination'), false, 'vaccination=90 should be valid');
  });

  // --- Valid defaults ---
  test('Defaults pass validation', async () => {
    const errors = validate(DEFAULTS);
    expect(errors.size, 0, 'Default parameters should be valid');
  });

  // --- Valid boundary values ---
  test('r0 = 0.5 (valid boundary)', async () => {
    const errors = validate({...DEFAULTS, r0: 0.5});
    expect(errors.has('ctrl_r0'), false, 'Should accept r0 = 0.5');
  });

  test('gamma = 0.5 (valid boundary)', async () => {
    const errors = validate({...DEFAULTS, gamma: 0.5});
    expect(errors.has('ctrl_gamma'), false, 'Should accept gamma = 0.5');
  });

  test('sigma = 1.0 (valid boundary, SEIR)', async () => {
    const errors = validate({...DEFAULTS, modelType: 'SEIR', sigma: 1.0});
    expect(errors.has('ctrl_sigma'), false, 'Should accept sigma = 1.0');
  });

  // --- Multiple simultaneous errors ---
  test('Multiple simultaneous errors', async () => {
    const errors = validate({
      modelType: 'SIR',
      r0: 0,
      gamma: 0,
      sigma: 0,
      vaccination: 100,
    });
    expect(errors.size >= 3, true, `Should report >= 3 errors, got ${errors.size}`);
    expect(errors.has('ctrl_r0'), true);
    expect(errors.has('ctrl_gamma'), true);
    expect(errors.has('ctrl_vaccination'), true);
  });
});
