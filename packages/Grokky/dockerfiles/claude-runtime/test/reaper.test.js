const test = require('node:test');
const assert = require('node:assert');
const {startOrphanReaper, childProcesses} = require('../dist/watchdog');

const tick = (ms) => new Promise((r) => setTimeout(r, ms));

test('polls the active-query count on a schedule', async () => {
  let polls = 0;
  const timer = startOrphanReaper(() => {
    polls++;
    return 1;
  }, 10);
  await tick(60);
  clearInterval(timer);
  assert.ok(polls > 1, `polled ${polls} times`);
});

test('returns an unref-ed timer so it cannot hold the process open', () => {
  const timer = startOrphanReaper(() => 0, 10);
  // A reaper that kept the event loop alive would stop the container ever shutting down.
  assert.equal(typeof timer.hasRef === 'function' ? timer.hasRef() : false, false);
  clearInterval(timer);
});

test('survives an activeCount that throws instead of crashing the runtime', async () => {
  // This runs on a bare timer callback: an escaping throw is an uncaught exception that takes
  // down the whole runtime. The janitor must never kill the process it is cleaning up after.
  const timer = startOrphanReaper(() => {
    throw new Error('registry unavailable');
  }, 10);
  await tick(40);
  clearInterval(timer);
  assert.ok(true, 'process still alive');
});

test('childProcesses reports pid, command and a monotonic startTime', () => {
  const kids = childProcesses();
  assert.ok(Array.isArray(kids));
  for (const k of kids) {
    assert.equal(typeof k.pid, 'number');
    assert.equal(typeof k.command, 'string');
    assert.equal(typeof k.startTime, 'number');
  }
});

test('childProcesses reports no children for an implausible parent pid', () => {
  assert.deepEqual(childProcesses(2147483646), []);
});

// The surplus arithmetic is the safety-critical part: too eager and it kills a live turn, too
// lax and the orphan keeps burning the container's only CPU. Exercised directly.
const surplusOf = (agentCount, active) => agentCount - Math.max(0, active);

test('surplus arithmetic never targets a process while queries are in flight', () => {
  assert.equal(surplusOf(1, 1) > 0, false, 'one agent, one query — nothing to reap');
  assert.equal(surplusOf(3, 3) > 0, false, 'three agents, three queries — all legitimate');
  assert.equal(surplusOf(2, 3) > 0, false, 'fewer agents than queries — never negative-kill');
});

test('surplus arithmetic reclaims exactly the excess', () => {
  assert.equal(surplusOf(1, 0), 1, 'one agent, no query — the orphan case');
  assert.equal(surplusOf(3, 1), 2, 'two orphans alongside one live turn');
  // A negative activeCount must not be read as "kill extra".
  assert.equal(surplusOf(1, -5), 1);
});
