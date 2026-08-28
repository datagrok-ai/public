/* The param-rules unit contract (plan-w4.md FP-W4-2/4/8): the regex-literal parser, the
   evalRule / evalValidatorExpression verdict mappings over the feature-detected grok_ScriptSync
   global, warn-once, message plaining, the categoryGroups plan expander — and the sanity of the
   test-only scriptSync stub itself. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {dart, scriptSyncStub} from './platform-doubles.mjs';
import {WarnOnce, evalRule, evalValidatorExpression, expandCategoryGroups, plainMessage,
  regexValidator} from '../src/dg/funcs/param-rules.js';

/** Installs the stub for the body and guarantees teardown: global deleted, knob reset. */
function stubbed(name, body) {
  test(name, () => {
    globalThis.grok_ScriptSync = scriptSyncStub;
    try {
      body();
    } finally {
      delete globalThis.grok_ScriptSync;
      dart.scriptSyncError = null;
    }
  });
}

/** No stub: the missing-global paths. */
function bare(name, body) {
  test(name, body);
}

function countWarns(body) {
  const warn = console.warn;
  let count = 0;
  console.warn = () => count++;
  try {
    body();
  } finally {
    console.warn = warn;
  }
  return count;
}

// --- regexValidator ---

bare('regexValidator parses /pattern/ and answers the Dart non-match message', () => {
  const v = regexValidator('/^[a-z]+$/');
  assert.notEqual(v, null);
  assert.equal(v('abc'), null);
  assert.equal(v('Abc'), 'Value doesn\'t match /^[a-z]+$/');
});

bare('the i flag is honored and kept in the message label', () => {
  const v = regexValidator('/^abc$/i');
  assert.equal(v('ABC'), null);
  assert.equal(v('xyz'), 'Value doesn\'t match /^abc$/i');
});

bare('unsupported flags are dropped from the RegExp but stay in the label (vld:125-131)', () => {
  const v = regexValidator('/^a.c$/gs');
  assert.equal(v('abc'), null, 'the g/s flags are stripped, the pattern still applies');
  assert.equal(v('zzz'), 'Value doesn\'t match /^a.c$/gs');
});

bare('a non-string value fails with the same message a non-match gets', () => {
  const v = regexValidator('/^\\d+$/');
  assert.equal(v(42), 'Value doesn\'t match /^\\d+$/');
  assert.equal(v(null), 'Value doesn\'t match /^\\d+$/');
});

bare('anything that is not a regex literal answers null: expressions, a lone slash, no closer', () => {
  assert.equal(regexValidator('minAge > 18'), null);
  assert.equal(regexValidator('/'), null);
  assert.equal(regexValidator('/abc'), null, 'lastIndexOf finds only the opener');
});

// --- evalRule ---

stubbed('evalRule: anything but false is true — the Dart result != false mapping', () => {
  const w = new WarnOnce();
  assert.equal(evalRule('x > 1', {x: 2}, false, 'k1', w), true);
  assert.equal(evalRule('x > 1', {x: 0}, true, 'k2', w), false);
  assert.equal(evalRule('"text"', {}, false, 'k3', w), true, 'a string result is not false');
  assert.equal(evalRule('2 + 2', {}, false, 'k4', w), true, 'a number result is not false');
  assert.equal(evalRule('null', {}, false, 'k5', w), true, 'null is not false (ib:341)');
});

stubbed('evalRule: a script throw keeps the previous state, warning once per key', () => {
  dart.scriptSyncError = true;
  const w = new WarnOnce();
  const warns = countWarns(() => {
    assert.equal(evalRule('x > 1', {x: 2}, true, 'f:visible', w), true, 'previous true kept');
    assert.equal(evalRule('x > 1', {x: 2}, false, 'f:visible', w), false, 'previous false kept');
  });
  assert.equal(warns, 1, 'one warn for two failures under the same key');
});

stubbed('clear() re-arms the once-per-key warn', () => {
  dart.scriptSyncError = true;
  const w = new WarnOnce();
  const warns = countWarns(() => {
    evalRule('x', {}, true, 'k', w);
    w.clear();
    evalRule('x', {}, true, 'k', w);
  });
  assert.equal(warns, 2);
});

bare('evalRule without the global keeps the previous state and warns once', () => {
  const w = new WarnOnce();
  const warns = countWarns(() => {
    assert.equal(evalRule('x > 1', {x: 2}, false, 'g:visible', w), false);
    assert.equal(evalRule('x > 1', {x: 2}, false, 'g:visible', w), false);
  });
  assert.equal(warns, 1);
});

// --- evalValidatorExpression ---

stubbed('evalValidatorExpression: false answers the expression text itself (ib:230-237 parity)', () => {
  assert.equal(evalValidatorExpression('minAge > 18', {minAge: 15}, 'k', new WarnOnce()), 'minAge > 18');
});

stubbed('evalValidatorExpression: a string result is the message, anything else is valid', () => {
  const w = new WarnOnce();
  assert.equal(evalValidatorExpression('x > 2 or "too small"', {x: 1}, 'k', w), 'too small');
  assert.equal(evalValidatorExpression('x > 2', {x: 5}, 'k', w), null);
  assert.equal(evalValidatorExpression('2 + 2', {}, 'k', w), null, 'a number is valid');
});

stubbed('evalValidatorExpression: a throw answers the Dart onError wording', () => {
  dart.scriptSyncError = true;
  assert.equal(evalValidatorExpression('x > 1', {x: 2}, 'k', new WarnOnce()),
    'Error during validation: "x > 1"');
});

bare('evalValidatorExpression without the global is valid and warns once — evalRule symmetry', () => {
  const w = new WarnOnce();
  const warns = countWarns(() => {
    assert.equal(evalValidatorExpression('x > 1', {x: 0}, 'v:validator', w), null);
    assert.equal(evalValidatorExpression('x > 1', {x: 0}, 'v:validator', w), null);
  });
  assert.equal(warns, 1);
});

// --- WarnOnce ---

bare('WarnOnce: one warn per key, independent keys independent', () => {
  const w = new WarnOnce();
  const warns = countWarns(() => {
    w.warn('a', 'm');
    w.warn('a', 'm');
    w.warn('b', 'm');
  });
  assert.equal(warns, 2);
});

bare('WarnOnce: instances neither cross-silence nor cross-reset — one per form', () => {
  const a = new WarnOnce();
  const b = new WarnOnce();
  const warns = countWarns(() => {
    a.warn('k', 'm');
    b.warn('k', 'm');
    a.clear();
    b.warn('k', 'm');
    a.warn('k', 'm');
  });
  assert.equal(warns, 3, 'b warns despite a; a\'s clear does not re-arm b');
});

// --- plainMessage ---

bare('plainMessage strips #{…} command markup and ** emphasis from builtin messages', () => {
  assert.equal(
    plainMessage('Column \'age\' contains missing values. Remove it using **Select | Missing ' +
      'Values...** then **Delete Rows**, or use #{x.MissingValuesImputation.Run}.'),
    'Column \'age\' contains missing values. Remove it using Select | Missing Values... ' +
    'then Delete Rows, or use .');
  assert.equal(plainMessage('plain text'), 'plain text');
});

// --- expandCategoryGroups ---

bare('expandCategoryGroups: a string is one leaf, an array a run of leaves at the depth', () => {
  assert.deepEqual(expandCategoryGroups('Engine', 1), [{isHeader: false, name: 'Engine', level: 1}]);
  assert.deepEqual(expandCategoryGroups(['A', 'B'], 2), [
    {isHeader: false, name: 'A', level: 2},
    {isHeader: false, name: 'B', level: 2}]);
});

bare('expandCategoryGroups: a map emits the header then recurses a level deeper (fpe:379-398)', () => {
  assert.deepEqual(expandCategoryGroups({Power: ['Engine', 'Battery'], Details: {More: 'Contact'}}, 1), [
    {isHeader: true, name: 'Power', level: 1},
    {isHeader: false, name: 'Engine', level: 2},
    {isHeader: false, name: 'Battery', level: 2},
    {isHeader: true, name: 'Details', level: 1},
    {isHeader: true, name: 'More', level: 2},
    {isHeader: false, name: 'Contact', level: 3}]);
});

bare('expandCategoryGroups: non-string list entries and null are dropped', () => {
  assert.deepEqual(expandCategoryGroups(['A', 42, null, 'B'], 1), [
    {isHeader: false, name: 'A', level: 1},
    {isHeader: false, name: 'B', level: 1}]);
  assert.deepEqual(expandCategoryGroups(null, 1), []);
  assert.deepEqual(expandCategoryGroups(7, 1), []);
});

// --- the stub itself ---

bare('scriptSyncStub: ==/!=/and/or/not translate to their JS forms', () => {
  assert.equal(scriptSyncStub('type == "ICE"', {type: 'ICE'}), true);
  assert.equal(scriptSyncStub('type == "ICE"', {type: 'EV'}), false);
  assert.equal(scriptSyncStub('type != "ICE"', {type: 'EV'}), true, 'the != rewrite survives the == pass');
  assert.equal(scriptSyncStub('a > 1 and b > 1', {a: 2, b: 2}), true);
  assert.equal(scriptSyncStub('a > 1 or b > 1', {a: 0, b: 2}), true);
  assert.equal(scriptSyncStub('not (a > 1)', {a: 0}), true);
  assert.equal(typeof scriptSyncStub('2 + 2', {}), 'number');
});

bare('scriptSyncStub: a bad reference throws, and the error knob forces failures', () => {
  assert.throws(() => scriptSyncStub('nosuchvar > 1', {a: 1}));
  dart.scriptSyncError = true;
  try {
    assert.throws(() => scriptSyncStub('a > 1', {a: 2}), /script failed/);
    dart.scriptSyncError = 'boom';
    assert.throws(() => scriptSyncStub('a > 1', {a: 2}), /boom/);
  } finally {
    dart.scriptSyncError = null;
  }
});
