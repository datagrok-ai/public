// Levins Metapopulation Model — API tests

import {category, test, expect} from '@datagrok-libraries/utils/src/test';

import {DEFAULTS, validate, validateOptimize} from '../levins/core';

category('API: Validation', () => {
  // --- val_01: p0 <= 0 ---
  test('val_01: p0 = 0', async () => {
    const errors = validate({...DEFAULTS, p0: 0});
    expect(errors.has('ctrl_p0'), true, 'Should reject p0 = 0');
  });

  test('val_01: p0 = -1', async () => {
    const errors = validate({...DEFAULTS, p0: -1});
    expect(errors.has('ctrl_p0'), true, 'Should reject p0 < 0');
  });

  // --- val_02: p0 > 1 ---
  test('val_02: p0 = 1.1', async () => {
    const errors = validate({...DEFAULTS, p0: 1.1});
    expect(errors.has('ctrl_p0'), true, 'Should reject p0 > 1');
  });

  // --- p0 valid boundary ---
  test('p0 = 0.001 (valid boundary)', async () => {
    const errors = validate({...DEFAULTS, p0: 0.001});
    expect(errors.has('ctrl_p0'), false, 'Should accept p0 = 0.001');
  });

  test('p0 = 1 (valid boundary)', async () => {
    const errors = validate({...DEFAULTS, p0: 1});
    expect(errors.has('ctrl_p0'), false, 'Should accept p0 = 1');
  });

  // --- val_03: m <= 0 ---
  test('val_03: m = 0', async () => {
    const errors = validate({...DEFAULTS, m: 0});
    expect(errors.has('ctrl_m'), true, 'Should reject m = 0');
  });

  test('val_03: m = -0.5', async () => {
    const errors = validate({...DEFAULTS, m: -0.5});
    expect(errors.has('ctrl_m'), true, 'Should reject m < 0');
  });

  // --- val_04: e0 <= 0 ---
  test('val_04: e0 = 0', async () => {
    const errors = validate({...DEFAULTS, e0: 0});
    expect(errors.has('ctrl_e0'), true, 'Should reject e0 = 0');
  });

  test('val_04: e0 = -0.1', async () => {
    const errors = validate({...DEFAULTS, e0: -0.1});
    expect(errors.has('ctrl_e0'), true, 'Should reject e0 < 0');
  });

  // --- val_05: m <= e0 (no rescue) ---
  test('val_05: m = e0, no rescue', async () => {
    const errors = validate({...DEFAULTS, m: 0.5, e0: 0.5, rescueEffect: false});
    expect(errors.has('ctrl_m'), true, 'Should reject m = e0 without rescue');
  });

  test('val_05: m < e0, no rescue', async () => {
    const errors = validate({...DEFAULTS, m: 0.3, e0: 0.5, rescueEffect: false});
    expect(errors.has('ctrl_m'), true, 'Should reject m < e0 without rescue');
  });

  test('val_05: m <= e0 with rescue (allowed)', async () => {
    const errors = validate({...DEFAULTS, m: 0.3, e0: 0.5, rescueEffect: true});
    expect(errors.has('ctrl_m'), false, 'Should allow m <= e0 with rescue');
  });

  test('val_05: skipped when val_03 fails', async () => {
    const errors = validate({...DEFAULTS, m: 0, e0: 0.5, rescueEffect: false});
    expect(errors.get('ctrl_m'), 'Colonization rate must be positive');
  });

  // --- val_06: t_end <= t_start ---
  test('val_06: t_end = t_start', async () => {
    const errors = validate({...DEFAULTS, t_start: 10, t_end: 10});
    expect(errors.has('ctrl_t_end'), true, 'Should reject t_end = t_start');
    expect(errors.has('ctrl_t_start'), true, 'Should set error on t_start too');
  });

  test('val_06: t_end < t_start', async () => {
    const errors = validate({...DEFAULTS, t_start: 10, t_end: 5});
    expect(errors.has('ctrl_t_end'), true, 'Should reject t_end < t_start');
  });

  // --- val_07: t_step <= 0 ---
  test('val_07: t_step = 0', async () => {
    const errors = validate({...DEFAULTS, t_step: 0});
    expect(errors.has('ctrl_t_step'), true, 'Should reject t_step = 0');
  });

  test('val_07: t_step = -0.1', async () => {
    const errors = validate({...DEFAULTS, t_step: -0.1});
    expect(errors.has('ctrl_t_step'), true, 'Should reject t_step < 0');
  });

  // --- val_08: t_step >= t_end - t_start ---
  test('val_08: t_step = interval length', async () => {
    const errors = validate({...DEFAULTS, t_start: 0, t_end: 50, t_step: 50});
    expect(errors.has('ctrl_t_step'), true, 'Should reject t_step = interval length');
  });

  test('val_08: t_step > interval length', async () => {
    const errors = validate({...DEFAULTS, t_start: 0, t_end: 50, t_step: 100});
    expect(errors.has('ctrl_t_step'), true, 'Should reject t_step > interval length');
  });

  test('val_08: skipped when val_07 fails', async () => {
    const errors = validate({...DEFAULTS, t_step: -1});
    expect(errors.get('ctrl_t_step'), 'Step must be positive');
  });

  // --- val_09: tolerance <= 0 ---
  test('val_09: tolerance = 0', async () => {
    const errors = validate({...DEFAULTS, tolerance: 0});
    expect(errors.has('ctrl_tolerance'), true, 'Should reject tolerance = 0');
  });

  test('val_09: tolerance = -1e-7', async () => {
    const errors = validate({...DEFAULTS, tolerance: -1e-7});
    expect(errors.has('ctrl_tolerance'), true, 'Should reject tolerance < 0');
  });

  // --- valid defaults ---
  test('Defaults pass validation', async () => {
    const errors = validate(DEFAULTS);
    expect(errors.size, 0, 'Default parameters should be valid');
  });

  // --- multiple errors ---
  test('Multiple simultaneous errors', async () => {
    const errors = validate({...DEFAULTS, p0: 0, m: 0, e0: 0, t_step: 0, tolerance: 0});
    expect(errors.size >= 4, true, 'Should report multiple errors');
    expect(errors.has('ctrl_p0'), true);
    expect(errors.has('ctrl_m'), true);
    expect(errors.has('ctrl_e0'), true);
    expect(errors.has('ctrl_t_step'), true);
    expect(errors.has('ctrl_tolerance'), true);
  });
});

category('API: Optimization Validation', () => {
  // --- opt_val_01: m_min <= 0 ---
  test('opt_val_01: m_min = 0', async () => {
    const {errors} = validateOptimize({m_min: 0, m_max: 1}, 0.2, false);
    expect(errors.has('dlg_m_min'), true, 'Should reject m_min = 0');
  });

  // --- opt_val_02: m_max <= 0 ---
  test('opt_val_02: m_max = -1', async () => {
    const {errors} = validateOptimize({m_min: 0.1, m_max: -1}, 0.2, false);
    expect(errors.has('dlg_m_max'), true, 'Should reject m_max < 0');
  });

  // --- opt_val_03: m_min >= m_max ---
  test('opt_val_03: m_min = m_max', async () => {
    const {errors} = validateOptimize({m_min: 0.5, m_max: 0.5}, 0.2, false);
    expect(errors.has('dlg_m_min'), true, 'Should reject m_min = m_max');
  });

  test('opt_val_03: m_min > m_max', async () => {
    const {errors} = validateOptimize({m_min: 1.0, m_max: 0.5}, 0.2, false);
    expect(errors.has('dlg_m_min'), true, 'Should reject m_min > m_max');
  });

  // --- opt_val_04: warning when m_max <= e0 ---
  test('opt_val_04: m_max <= e0, no rescue — warning', async () => {
    const {errors, warning} = validateOptimize({m_min: 0.05, m_max: 0.1}, 0.2, false);
    expect(errors.size, 0, 'Should not block');
    expect(warning !== null, true, 'Should produce warning');
  });

  test('opt_val_04: m_max <= e0 with rescue — no warning', async () => {
    const {errors, warning} = validateOptimize({m_min: 0.05, m_max: 0.1}, 0.2, true);
    expect(errors.size, 0);
    expect(warning, null, 'No warning with rescue effect');
  });

  // --- valid ---
  test('Valid optimization inputs', async () => {
    const {errors, warning} = validateOptimize({m_min: 0.1, m_max: 1.0}, 0.2, false);
    expect(errors.size, 0, 'Should pass');
    expect(warning, null, 'No warning');
  });
});

