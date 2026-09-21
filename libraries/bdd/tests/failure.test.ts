import assert from 'node:assert/strict';
import {test} from 'node:test';
import {failure, isWaitFailure, journeyFailure, reasonOf, StepFailure} from '../src/runtime/failure.js';

const IN_PAGE = 'locator.evaluate: Error: Box plot has no "statsff" area right now; it has: view, x axis, stats\n' +
  '    at hitArea (eval at evaluate (:311:30), <anonymous>:152:27)\n' +
  '    at Object.menuPoint (eval at evaluate (:311:30), <anonymous>:273:13)\n' +
  '    at UtilityScript.evaluate (<anonymous>:313:16)';

test('an in-page error keeps its sentence and loses the API prefix and the eval stack', () => {
  assert.equal(reasonOf(new Error(IN_PAGE)), 'Box plot has no "statsff" area right now; it has: view, x axis, stats');
});

test('a Playwright timeout reads as seconds; the selector it waited for goes, the rest of the call log stays', () => {
  const e = new Error('locator.hover: Timeout 10000ms exceeded.\nCall log:\n  - waiting for locator(\'.d4-menu-popup\').filter({ visible: true })\n\n');
  e.name = 'TimeoutError';
  assert.equal(reasonOf(e), 'timed out after 10 s');
  const dim = String.fromCharCode(27);
  const covered = new Error(`locator.click: Timeout 10000ms exceeded.\nCall log:\n${dim}[2m  - waiting for locator('.x')${dim}[22m\n${dim}[2m  - <div>…</div> intercepts pointer events${dim}[22m`);
  assert.equal(reasonOf(covered), 'timed out after 10 s\nCall log:\n  - <div>…</div> intercepts pointer events');
  assert.ok(isWaitFailure(e));
  assert.ok(isWaitFailure(new Error('expect(locator).toHaveCount(expected) failed\n\nExpected: 0')));
  assert.ok(!isWaitFailure(new Error(IN_PAGE)));
});

test('a step failure is the feature line, the step, the reason and what the page showed — no frames', () => {
  const f = failure('features/viewers/box-plot.feature:43', 'When user right-clicks on the "statsff" area of box plot viewer',
    new Error(IN_PAGE), 'visible menu items in context menu: Table | Show P Value');
  assert.ok(f instanceof StepFailure);
  assert.equal(f.message, 'features/viewers/box-plot.feature:43\n' +
    '  When user right-clicks on the "statsff" area of box plot viewer\n\n' +
    'Box plot has no "statsff" area right now; it has: view, x axis, stats\n' +
    'visible menu items in context menu: Table | Show P Value');
  assert.equal(f.stack, f.message);
  assert.equal(failure('x:1', 'again', f), f, 'a nested step passes the failure through');
  const located = failure('f.feature:43', 'When', new Error('x'), '', 'C:\\pkg\\bdd\\features\\f.feature:43:1');
  assert.equal(located.stack, `${located.message}\n    at C:\\pkg\\bdd\\features\\f.feature:43:1`, 'the feature line is the one frame');
  assert.match(journeyFailure([{name: 'S', error: located}], 2).stack!, /\n    at C:\\pkg\\bdd\\features\\f\.feature:43:1$/);
});

test('a programming error in a binding keeps its own frames and drops Playwright\'s', () => {
  const e = new TypeError('Cannot read properties of undefined (reading \'x\')');
  e.stack = `TypeError: ${e.message}\n    at hover (C:\\pkg\\bdd\\bindings\\box-plot.ts:12:5)\n` +
    '    at TestTypeImpl._step (C:\\lib\\node_modules\\playwright\\lib\\common\\testType.js:250:9)';
  const f = failure('f.feature:7', 'When user hovers', e);
  assert.match(f.stack!, /f\.feature:7\n  When user hovers\n\nCannot read properties of undefined \(reading 'x'\)\n    at hover \(C:\\pkg\\bdd\\bindings\\box-plot\.ts:12:5\)$/);
});

test('a journey lists the failed scenarios with their step reports', () => {
  const e = journeyFailure([
    {name: 'Menu regions', error: failure('f.feature:43', 'When user right-clicks', new Error(IN_PAGE))},
    {name: 'Coloring', error: new Error('boom')},
  ], 13);
  assert.equal(e.message, '2 of 13 scenarios failed\n\n' +
    'Menu regions\n  f.feature:43\n    When user right-clicks\n\n  Box plot has no "statsff" area right now; it has: view, x axis, stats\n\n' +
    'Coloring\n  boom');
  assert.equal(e.stack, e.message);
});
