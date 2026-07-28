const test = require('node:test');
const assert = require('node:assert');
const {raceIdle, IDLE, childProcesses} = require('../dist/watchdog');

const after = (ms, value) => new Promise((r) => setTimeout(() => r(value), ms));

test('passes through a value that settles in time', async () => {
  assert.equal(await raceIdle(after(10, 'event'), 200), 'event');
});

test('reports IDLE when nothing settles in time', async () => {
  assert.equal(await raceIdle(after(200, 'event'), 20), IDLE);
});

test('IDLE is distinguishable from a legitimate falsy value', async () => {
  // A turn can legitimately yield undefined/null; conflating that with a timeout would kill
  // healthy turns, so the sentinel must not be falsy-equal to anything real.
  assert.notEqual(await raceIdle(after(5, undefined), 200), IDLE);
  assert.notEqual(await raceIdle(after(5, null), 200), IDLE);
  assert.notEqual(await raceIdle(after(5, 0), 200), IDLE);
  assert.notEqual(await raceIdle(after(5, false), 200), IDLE);
});

test('rejection propagates rather than being swallowed as a timeout', async () => {
  await assert.rejects(() => raceIdle(Promise.reject(new Error('boom')), 200), /boom/);
});

test('does not hold the event loop open after resolving early', async () => {
  // If the idle timer were left pending, a finished process would hang for the timeout duration.
  const t0 = Date.now();
  await raceIdle(after(5, 'x'), 60000);
  assert.ok(Date.now() - t0 < 1000, 'resolved promptly');
});

test('childProcesses returns an array and never throws off-Linux', () => {
  const kids = childProcesses();
  assert.ok(Array.isArray(kids));
  for (const k of kids) {
    assert.equal(typeof k.pid, 'number');
    assert.equal(typeof k.command, 'string');
  }
});

test('childProcesses reports no children for an implausible parent pid', () => {
  assert.deepEqual(childProcesses(2147483646), []);
});
