// Unit tests for compileFormula / runFormula in formulas-resolver.ts.
// Pin the contract: simple math + conditionals; identifier-extraction
// safety across string literals, comments, property access, method calls,
// template literals, JS globals, unknown names, pure-numeric formulas;
// error paths (runtime throw, NaN result, syntax error). Plus a regression
// for the legacy runFormula(formula, ctx) signature used by
// fitting-view.ts:1308-1309.

import dayjs from 'dayjs';
import {category, test, expect, expectFloat} from '@datagrok-libraries/test/src/test';
import {compileFormula, runFormula} from
  '@datagrok-libraries/compute-utils/function-views/src/fitting/formulas-resolver';

category('ComputeUtils: Fitting / Formula resolver', () => {
  test('compile_and_run_simple_math', async () => {
    const fn = compileFormula('a + b * 2', ['a', 'b']);
    expectFloat(fn({a: 3, b: 4})!, 11, 1e-12);
    expectFloat(fn({a: -1, b: 0.5})!, 0, 1e-12);
  });

  test('compile_and_run_conditional', async () => {
    const fn = compileFormula('a > 0 ? a * scale + offset : 0', ['a', 'scale', 'offset']);
    expectFloat(fn({a: 4, scale: 2.5, offset: 0.7})!, 10.7, 1e-12);
    expectFloat(fn({a: -1, scale: 2.5, offset: 0.7})!, 0, 1e-12);
  });

  test('compile_intersects_with_known_names', async () => {
    // 'positive' / 'negative' are inside string literals — they look like
    // identifiers to the regex but are not in knownNames, so they drop out
    // of the destructure. Formula still works because the literals stay literal.
    const fn = compileFormula(`a > 0 ? 'positive'.length : 'negative'.length`, ['a']);
    expect(fn({a: 1})!, 8);
    expect(fn({a: -1})!, 8);
  });

  test('compile_handles_property_access', async () => {
    const fn = compileFormula('refDf.rowCount * 0.1', ['refDf']);
    expectFloat(fn({refDf: {rowCount: 50}})!, 5, 1e-12);
  });

  test('compile_handles_dayjs_method', async () => {
    const fn = compileFormula('t0.year() / 10000', ['t0']);
    expectFloat(fn({t0: dayjs('2024-06-15')})!, 0.2024, 1e-9);
  });

  test('compile_returns_null_on_runtime_error', async () => {
    const fn = compileFormula('refDf.rowCount * 0.1', ['refDf']);
    expect(fn({refDf: undefined as any}), null);
  });

  test('compile_returns_null_on_nan', async () => {
    const fn = compileFormula('Math.sqrt(a)', ['a']);
    expect(fn({a: -1}), null);
    expectFloat(fn({a: 4})!, 2, 1e-12);
  });

  test('compile_returns_null_on_syntax_error', async () => {
    const fn = compileFormula('a +* 2', ['a']);
    expect(fn({a: 1}), null);
  });

  test('compile_handles_comment_inside_formula', async () => {
    // Identifier inside a /* */ comment looks like a context-var token to
    // the regex. Intersect keeps it; destructure adds an unused binding.
    // Formula behavior must be unaffected by the comment text.
    const fn = compileFormula('a + 1 /* refDf */', ['a', 'refDf']);
    expectFloat(fn({a: 41, refDf: undefined as any})!, 42, 1e-12);
  });

  test('compile_handles_unknown_name_reference', async () => {
    // The formula references a name not in contextNames AND not in JS
    // globals. Should return null (ReferenceError caught) — same as the
    // legacy runFormula behavior.
    const fn = compileFormula('unknown_name + 1', ['a']);
    expect(fn({a: 1}), null);
  });

  test('compile_handles_pure_numeric_formula', async () => {
    // No identifiers — empty destructure line, formula is a literal.
    const fn = compileFormula('42', []);
    expect(fn({}), 42);
  });

  test('compile_handles_math_global', async () => {
    // Math is not in contextNames; formula uses it via global scope.
    const fn = compileFormula('Math.max(a, 0)', ['a']);
    expectFloat(fn({a: 5})!, 5, 1e-12);
    expectFloat(fn({a: -3})!, 0, 1e-12);
  });

  test('compile_handles_template_literal', async () => {
    // Identifier inside ${} of a template literal. Should still bind.
    const fn = compileFormula('`val=${a}`.length', ['a']);
    expect(fn({a: 7}), 5);
  });

  test('runFormula_legacy_path_works', async () => {
    // The fitting-view.ts UI path uses this signature.
    expectFloat(runFormula('a + 1', {a: 41})!, 42, 1e-12);
    expect(runFormula('bogus_name', {}), null);
  });
});
