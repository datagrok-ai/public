// Levins Metapopulation Model — Core tests

import {category, test, expect, expectFloat} from '@datagrok-libraries/utils/src/test';

import {
  DEFAULTS, validate, solve, validateOptimize,
  getEquilibrium, LevinsParams,
} from '../levins/core';

category('Levins: Validation', () => {
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

category('Levins: Optimization Validation', () => {
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

category('Levins: Equilibrium', () => {
  test('p* = 1 - e0/m (base model)', async () => {
    expectFloat(getEquilibrium(0.5, 0.2, false), 0.6, 1e-10);
  });

  test('p* = 0 when m <= e0', async () => {
    expectFloat(getEquilibrium(0.2, 0.5, false), 0, 1e-10);
  });

  test('p* = 0 when m = e0', async () => {
    expectFloat(getEquilibrium(0.5, 0.5, false), 0, 1e-10);
  });

  test('p* = NaN with rescue effect', async () => {
    expect(isNaN(getEquilibrium(0.5, 0.2, true)), true, 'Should be NaN with rescue');
  });
});

category('Levins: Solver', () => {
  test('Default parameters produce valid solution', async () => {
    const result = solve(DEFAULTS);
    expect(result.t.length > 0, true, 'Should produce t array');
    expect(result.p.length > 0, true, 'Should produce p array');
    expect(result.t.length, result.p.length, 't and p must have same length');
  });

  test('Solution values in [0, 1]', async () => {
    const result = solve(DEFAULTS);
    for (let i = 0; i < result.p.length; i++) {
      expect(result.p[i] >= 0 && result.p[i] <= 1, true, `p[${i}] = ${result.p[i]} out of [0, 1]`);
    }
  });

  test('p(0) = p0', async () => {
    const result = solve(DEFAULTS);
    expectFloat(result.p[0], DEFAULTS.p0, 1e-6);
  });

  test('t starts at t_start', async () => {
    const result = solve(DEFAULTS);
    expectFloat(result.t[0], DEFAULTS.t_start, 1e-10);
  });

  test('p converges to p* (base model)', async () => {
    const params: LevinsParams = {...DEFAULTS, m: 0.5, e0: 0.2, t_end: 200};
    const result = solve(params);
    const pEnd = result.p[result.p.length - 1];
    expectFloat(pEnd, result.p_star, 0.01);
  });

  test('Rescue effect: solution stays bounded', async () => {
    const params: LevinsParams = {...DEFAULTS, m: 0.3, e0: 0.5, rescueEffect: true, t_end: 100};
    const result = solve(params);
    for (let i = 0; i < result.p.length; i++) {
      expect(result.p[i] >= 0 && result.p[i] <= 1, true, `p[${i}] = ${result.p[i]} out of [0, 1]`);
    }
  });

  test('Higher m leads to higher p(t_end)', async () => {
    const result1 = solve({...DEFAULTS, m: 0.5, e0: 0.2});
    const result2 = solve({...DEFAULTS, m: 1.0, e0: 0.2});
    const pEnd1 = result1.p[result1.p.length - 1];
    const pEnd2 = result2.p[result2.p.length - 1];
    expect(pEnd2 > pEnd1, true, 'Higher m should yield higher p(t_end)');
  });

  test('Custom p0 is used', async () => {
    const result = solve({...DEFAULTS, p0: 0.9});
    expectFloat(result.p[0], 0.9, 1e-6);
  });
});
