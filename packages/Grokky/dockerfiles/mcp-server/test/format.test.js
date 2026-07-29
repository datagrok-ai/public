// Result shaping. A tool result stays in the transcript and is re-read on every later hop of the
// turn, so these caps are a per-turn multiplier, not a one-off saving.

import test from 'node:test';
import assert from 'node:assert';
import {formatResult, formatError} from '../dist/format.js';

const text = (r) => r.content[0].text;

test('results are compact JSON, not pretty-printed', () => {
  const out = text(formatResult({a: 1, b: {c: 2}}));
  assert.equal(out, '{"a":1,"b":{"c":2}}');
  assert.ok(!out.includes('\n'), 'pretty-printing spends tokens on whitespace');
});

test('array results are paged with a note pointing at offset/limit', () => {
  const data = Array.from({length: 120}, (_, i) => ({i}));
  const out = JSON.parse(text(formatResult(data, {paged: true, args: {}})));
  assert.equal(out.total, 120);
  assert.equal(out.items.length, 50, 'default page is 50');
  assert.equal(out.showing, '1-50');
  assert.match(out.note, /offset/);
});

test('paging honours offset and limit', () => {
  const data = Array.from({length: 120}, (_, i) => ({i}));
  const out = JSON.parse(text(formatResult(data, {paged: true, args: {offset: 100, limit: 5}})));
  assert.equal(out.items.length, 5);
  assert.equal(out.items[0].i, 100);
  assert.equal(out.showing, '101-105');
});

test('a short array is returned bare — no paging envelope to unwrap', () => {
  const out = JSON.parse(text(formatResult([{a: 1}], {paged: true, args: {}})));
  assert.deepEqual(out, [{a: 1}]);
});

test('oversized results are truncated with an actionable hint', () => {
  const out = text(formatResult('x'.repeat(50_000)));
  assert.ok(out.length < 21_000, 'must be capped');
  assert.match(out, /truncated at 20000 chars/);
});

test('file payloads get the larger cap and are never paged', () => {
  const body = Array.from({length: 200}, (_, i) => `line ${i}`);
  const out = text(formatResult(body, {raw: true, paged: false}));
  assert.deepEqual(JSON.parse(out).length, 200, 'file content must not be silently paged away');
});

test('errors come back as a readable result, not a throw', () => {
  const out = JSON.parse(text(formatError('unknown op', {available: ['list', 'get']})));
  assert.equal(out.error, 'unknown op');
  assert.deepEqual(out.available, ['list', 'get']);
});
