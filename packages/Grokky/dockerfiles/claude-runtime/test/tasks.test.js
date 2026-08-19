const test = require('node:test');
const assert = require('node:assert');
const {claimTask, waitForClaim, endTask, releaseTask} = require('../dist/tasks');

let n = 0;
const freshId = () => `task-${++n}`;

test('worker claim resolves when the turn ends', async () => {
  const id = freshId();
  const claim = claimTask(id, 's1');
  endTask(id);
  assert.deepEqual(await claim, {status: 'ended'});
});

test('a claim after the turn already ended resolves immediately', async () => {
  const id = freshId();
  waitForClaim(id, 's1', new AbortController().signal);
  endTask(id);
  assert.deepEqual(await claimTask(id, 's1'), {status: 'ended'});
});

test('turn holds until the worker claims', async () => {
  const id = freshId();
  const held = waitForClaim(id, 's1', new AbortController().signal);
  claimTask(id, 's1');
  assert.equal(await held, true);
  endTask(id);
});

test('abort while holding resolves false and the task still ends cleanly', async () => {
  const id = freshId();
  const ac = new AbortController();
  const held = waitForClaim(id, 's1', ac.signal);
  ac.abort();
  assert.equal(await held, false);
  endTask(id);
  assert.deepEqual(await claimTask(id, 's1'), {status: 'ended'});
});

test('release frees both waiters and reports the session to abort', async () => {
  const id = freshId();
  const held = waitForClaim(id, 's1', new AbortController().signal);
  const claim = claimTask(id, 's1');
  assert.equal(releaseTask(id), 's1');
  assert.equal(await held, true); // was claimed before the release
  assert.deepEqual(await claim, {status: 'released'});
});

test('end after release keeps the released status; release after end is a no-op', async () => {
  const id = freshId();
  claimTask(id, 's1');
  releaseTask(id);
  endTask(id);
  assert.deepEqual(await claimTask(id, 's1'), {status: 'released'});
  const id2 = freshId();
  endTask(id2 /* unknown task */);
  assert.equal(releaseTask(id2), null);
});
