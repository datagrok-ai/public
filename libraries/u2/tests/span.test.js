/* core/span.ts: the one span parser — text check, resolution against now, the non-enumerable mark. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {isSpanText, markSpan, resolveSpan, spanOf} from '../src/core/span.js';
import {FilterError} from '../src/core/filter/index.js';

const NOW = new Date(2026, 0, 31, 12);

test('isSpanText: signed unit spans and now; anything else is not a span', () => {
  assert.deepEqual(['-1w', '2d', '0h', 'now', '12m', '-3y'].map(isSpanText), [true, true, true, true, true, true]);
  assert.deepEqual(['', '1x', '- 1w', '1.5d', 'soon', '2026-01-01'].map(isSpanText),
    [false, false, false, false, false, false]);
});

test('resolveSpan: a fresh Date offset from now; now is a copy; a bad span throws FilterError', () => {
  assert.equal(resolveSpan('-1w', NOW).getTime(), NOW.getTime() - 7 * 86400e3);
  assert.equal(resolveSpan('2h', NOW).getTime(), NOW.getTime() + 2 * 3600e3);
  const same = resolveSpan('now', NOW);
  assert.equal(same.getTime(), NOW.getTime());
  assert.notEqual(same, NOW);
  assert.throws(() => resolveSpan('soon', NOW), FilterError);
});

test('markSpan / spanOf: the tag rides the Date but stays out of JSON and enumeration', () => {
  const date = markSpan(resolveSpan('-1w', NOW), '-1w');
  assert.equal(spanOf(date), '-1w');
  assert.equal(spanOf(new Date()), undefined);
  assert.equal(JSON.stringify({date}), JSON.stringify({date: date.toISOString()}));
  assert.deepEqual(Object.keys(date), []);
});
