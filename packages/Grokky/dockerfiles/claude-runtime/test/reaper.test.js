const test = require('node:test');
const assert = require('node:assert');
const {startOrphanReaper, isAgentProcess, descendantPids} = require('../dist/watchdog');

const tick = (ms) => new Promise((r) => setTimeout(r, ms));
// Every process counts as owned, so a sweep on Linux never kills anything of this machine's.
const everyPid = {has: () => true};

test('polls the pool\'s owned pids on a schedule', async () => {
  let polls = 0;
  const timer = startOrphanReaper(() => {
    polls++;
    return everyPid;
  }, 10);
  await tick(60);
  clearInterval(timer);
  assert.ok(polls > 1, `polled ${polls} times`);
});

test('returns an unref-ed timer so it cannot hold the process open', () => {
  const timer = startOrphanReaper(() => everyPid, 10);
  // A reaper that kept the event loop alive would stop the container ever shutting down.
  assert.equal(typeof timer.hasRef === 'function' ? timer.hasRef() : false, false);
  clearInterval(timer);
});

test('survives an ownedPids that throws instead of crashing the runtime', async () => {
  // This runs on a bare timer callback: an escaping throw is an uncaught exception that takes
  // down the whole runtime. The janitor must never kill the process it is cleaning up after.
  const timer = startOrphanReaper(() => {
    throw new Error('registry unavailable');
  }, 10);
  await tick(40);
  clearInterval(timer);
});

test('isAgentProcess matches the agent module on the command line, whatever the interpreter', () => {
  assert.ok(isAgentProcess('/opt/hermes/.venv/bin/python -m grokky_agent '));
  assert.ok(isAgentProcess('/opt/hermes/.venv/bin/python3.14 -m grokky_agent'));
  assert.ok(isAgentProcess('python3 -u -m grokky_agent '));
});

test('isAgentProcess never matches a script a terminal tool call runs', () => {
  assert.ok(!isAgentProcess('python3 analyze.py '));
  assert.ok(!isAgentProcess('python3 -c print("-m grokky_agent is fine") '));
  assert.ok(!isAgentProcess('python3 -m grokky_agent_tools '));
  assert.ok(!isAgentProcess('python3 -m pytest grokky_agent '));
  assert.ok(!isAgentProcess('bash -c ls /opt/grokky_agent-notes '));
  assert.ok(!isAgentProcess('node /app/dist/server.js '));
});

test('descendantPids collects the whole tree under a pid, at any depth, and nothing else', () => {
  const table = [
    {pid: 1, ppid: 0}, {pid: 10, ppid: 1}, {pid: 11, ppid: 10}, {pid: 12, ppid: 11}, {pid: 13, ppid: 11},
    {pid: 14, ppid: 10}, {pid: 20, ppid: 1}, {pid: 21, ppid: 20},
  ];
  assert.deepEqual(descendantPids(10, table).sort((a, b) => a - b), [11, 12, 13, 14]);
  assert.deepEqual(descendantPids(20, table), [21]);
  assert.deepEqual(descendantPids(12, table), []);
});
