/* StatCard: value formatting, the signed delta, and the loading / error / ready states an
   AsyncValue source drives. Same contract as the other UI suites — every test leaves the
   live-scope count where it found it. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {signal, Scope} from '../src/index.js';
import {AsyncValue} from '../src/core/async-value.js';
import {StatCard} from '../src/components/display/stat-card.js';

function ui(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

function deferred() {
  const box = {};
  box.promise = new Promise((resolve, reject) => {
    box.resolve = resolve;
    box.reject = reject;
  });
  return box;
}

const value = (card) => card.root.querySelector('.u2-stat-value').textContent;
const delta = (card) => card.root.querySelector('.u2-stat-delta');

ui('the value area formats numbers, passes strings through, and honours a custom format', () => {
  const number = new StatCard({label: 'Rows'});
  assert.equal(number.root.dataset.u2, 'stat-card');
  assert.equal(number.root.classList.contains('u2-card'), true, 'it shares the card surface');
  assert.equal(number.root.classList.contains('u2-stat-card'), true);
  assert.equal(value(number), '—', 'no value reads as an em dash');
  assert.equal(number.root.querySelector('.u2-stat-label').textContent, 'Rows');
  assert.equal(delta(number), null, 'no delta option, no delta line');
  number.dispose();

  const local = new StatCard({label: 'Rows', value: 1234567});
  assert.equal(value(local), (1234567).toLocaleString());
  local.dispose();

  const literal = new StatCard({label: 'Rows', value: '1.2M'});
  assert.equal(value(literal), '1.2M');
  literal.dispose();

  const formatted = new StatCard({label: 'Rows', value: 1234567,
    format: (x) => `${(x / 1e6).toFixed(1)}M`});
  assert.equal(value(formatted), '1.2M');
  formatted.dispose();

  const withIcon = new StatCard({label: 'Revenue', value: '1.2M', icon: 'chart-line'});
  assert.equal(withIcon.root.querySelector('.u2-stat-icon').classList.contains('fa-chart-line'), true);
  withIcon.dispose();
});

ui('a signal value updates live and stops when the card is disposed', () => {
  const rows = signal(10);
  const card = new StatCard({label: 'Rows', value: rows});
  assert.equal(value(card), (10).toLocaleString());

  rows.value = 20000;
  assert.equal(value(card), (20000).toLocaleString());

  rows.value = undefined;
  assert.equal(value(card), '—');

  rows.value = 5;
  card.dispose();
  rows.value = 6;
  assert.equal(value(card), (5).toLocaleString(), 'the effect died with the card');
});

ui('the delta colors by direction, flips under deltaInverted, and stays neutral at zero', () => {
  const up = new StatCard({label: 'Revenue', delta: 0.12});
  assert.equal(delta(up).classList.contains('u2-stat-delta-up'), true);
  assert.equal(delta(up).textContent, '+12%');
  up.dispose();

  const down = new StatCard({label: 'Revenue', delta: -0.034});
  assert.equal(delta(down).classList.contains('u2-stat-delta-down'), true);
  assert.equal(delta(down).textContent, '-3%');
  down.dispose();

  const inverted = new StatCard({label: 'Error rate', delta: 0.12, deltaInverted: true});
  assert.equal(delta(inverted).classList.contains('u2-stat-delta-down'), true,
    'growth is the bad direction here');
  assert.equal(delta(inverted).textContent, '+12%', 'the sign still reads as it is');
  inverted.dispose();

  const improved = new StatCard({label: 'Error rate', delta: -0.05, deltaInverted: true});
  assert.equal(delta(improved).classList.contains('u2-stat-delta-up'), true);
  improved.dispose();

  const custom = new StatCard({label: 'Latency', delta: 0.2, deltaFormat: (d) => `${d * 1000} ms`});
  assert.equal(delta(custom).textContent, '200 ms');
  custom.dispose();
});

ui('a signal delta follows the sign, and hides itself when it has nothing to say', () => {
  const change = signal(0.2);
  const card = new StatCard({label: 'Revenue', delta: change});
  assert.equal(delta(card).classList.contains('u2-stat-delta-up'), true);

  change.value = -0.2;
  assert.equal(delta(card).classList.contains('u2-stat-delta-down'), true);
  assert.equal(delta(card).classList.contains('u2-stat-delta-up'), false);

  change.value = 0;
  assert.equal(delta(card).hidden, false, 'zero is a fact worth showing');
  assert.equal(delta(card).classList.contains('u2-stat-delta-up'), false);
  assert.equal(delta(card).classList.contains('u2-stat-delta-down'), false);
  assert.equal(delta(card).textContent, '0%');

  change.value = undefined;
  assert.equal(delta(card).hidden, true);
  assert.equal(delta(card).textContent, '');
  card.dispose();
});

ui('a source drives the value area: skeleton while loading, message on error, value when ready', async () => {
  const run = deferred();
  const source = new AsyncValue(() => run.promise, {auto: false});
  const card = new StatCard({label: 'Queries today', source, value: 'ignored'});
  assert.notEqual(card.root.querySelector('.u2-skeleton'), null, 'idle shows the skeleton');

  const landed = source.refresh();
  assert.notEqual(card.root.querySelector('.u2-skeleton'), null, 'and so does loading');
  run.resolve(1200);
  await landed;
  await flush();
  assert.equal(value(card), (1200).toLocaleString(), 'the source wins over the value option');
  assert.equal(card.root.querySelector('.u2-skeleton'), null);
  card.dispose();
  source.dispose();
});

ui('a failing source shows the message, and a later run recovers', async () => {
  let next = deferred();
  const source = new AsyncValue(() => next.promise, {auto: false});
  const card = new StatCard({label: 'Queries today', source});

  let landed = source.refresh();
  next.reject(new Error('connection refused'));
  await landed;
  await flush();
  assert.equal(card.root.querySelector('.u2-stat-error').textContent, 'connection refused');

  next = deferred();
  landed = source.refresh();
  next.resolve('42');
  await landed;
  await flush();
  assert.equal(value(card), '42');
  assert.equal(card.root.querySelector('.u2-stat-error'), null);
  card.dispose();
  source.dispose();
});
