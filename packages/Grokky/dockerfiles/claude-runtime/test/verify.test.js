// Verifier gate — the state machine that decides whether a turn may Stop.
// Tests run against the compiled dist/, so `npm test` builds first.
// Note: the hook callbacks are async; every invocation must be awaited.

const test = require('node:test');
const assert = require('node:assert');
const {Verifier, isActionTool, isReadonlyTool, bareToolName} = require('../dist/verify.js');

const post = (toolName, response) => ({hook_event_name: 'PostToolUse', tool_name: toolName, tool_response: response});
const stop = () => ({hook_event_name: 'Stop'});
const mcpResult = (o) => ({content: [{type: 'text', text: JSON.stringify(o)}]});
const blocked = (r) => r.decision === 'block';

const EXEC = 'mcp__datagrok-browser__datagrok_exec';
const VERIFY = 'mcp__datagrok-browser__datagrok_verify';

test('bareToolName strips the mcp server prefix', () => {
  assert.equal(bareToolName('mcp__datagrok__call_function'), 'call_function');
  assert.equal(bareToolName(EXEC), 'datagrok_exec');
  assert.equal(bareToolName('Read'), 'Read');
});

test('isReadonlyTool recognises the read-only name prefixes', () => {
  for (const n of ['list_functions', 'get_project', 'search_project', 'read_space_file', 'download_file', 'whoami'])
    assert.ok(isReadonlyTool(n), `${n} should be read-only`);
  assert.ok(isReadonlyTool('datagrok_show_entities'), 'show_entities is a render, not a mutation');
  for (const n of ['create_project', 'upload_file', 'delete_space', 'call_function'])
    assert.ok(!isReadonlyTool(n), `${n} should not be read-only`);
});

test('isActionTool: exec and mutating mcp tools count, reads and plain tools do not', () => {
  assert.ok(isActionTool(EXEC));
  assert.ok(isActionTool('mcp__datagrok__call_function'));
  assert.ok(isActionTool('mcp__datagrok__create_script'));
  assert.ok(!isActionTool('mcp__datagrok__list_functions'));
  assert.ok(!isActionTool('mcp__datagrok__whoami'));
  assert.ok(!isActionTool('Read'), 'non-mcp built-ins are never actions');
  assert.ok(!isActionTool('Bash'));
});

test('an unverified action blocks the Stop', async () => {
  const v = new Verifier();
  await v.postToolUse(post(EXEC, mcpResult({success: true})));
  assert.ok(blocked(await v.stop(stop())));
});

test('a turn with no actions never blocks', async () => {
  const v = new Verifier();
  await v.postToolUse(post('Read', {}));
  await v.postToolUse(post('mcp__datagrok__list_functions', mcpResult([])));
  assert.ok(!blocked(await v.stop(stop())));
  assert.ok(!v.hadActions);
});

test('a self-verifying exec satisfies the gate without a separate verify round-trip', async () => {
  const v = new Verifier();
  await v.postToolUse(post(EXEC, mcpResult({success: true, verified: {passed: true, observed: 1}})));
  assert.ok(!blocked(await v.stop(stop())), 'in-round-trip verification should count as a pass');
  assert.equal(JSON.parse(v.statsLine()).passes, 1);
});

test('a self-verifying exec that FAILED still blocks', async () => {
  const v = new Verifier();
  await v.postToolUse(post(EXEC, mcpResult({success: true, verified: {passed: false}})));
  assert.ok(blocked(await v.stop(stop())));
});

test('a passing datagrok_verify clears everything pending', async () => {
  const v = new Verifier();
  await v.postToolUse(post(EXEC, mcpResult({success: true})));
  await v.postToolUse(post(EXEC, mcpResult({success: true})));
  await v.postToolUse(post(VERIFY, mcpResult({passed: true})));
  assert.ok(!blocked(await v.stop(stop())), 'one passing verify covers all pending actions in the turn');
});

test('a failing verify blocks with the stronger "did NOT take effect" reason', async () => {
  const v = new Verifier();
  await v.postToolUse(post(EXEC, mcpResult({success: true})));
  await v.postToolUse(post(VERIFY, mcpResult({passed: false})));
  const r = await v.stop(stop());
  assert.ok(blocked(r));
  assert.match(r.reason, /did NOT take effect/);
});

test('block reasons are framed as internal feedback the model must not quote', async () => {
  const v = new Verifier();
  await v.postToolUse(post(EXEC, mcpResult({success: true})));
  const r = await v.stop(stop());
  assert.match(r.reason, /never mention, quote, or allude to it/);
});

test('the gate gives up after 3 blocks and marks the turn unverified', async () => {
  const v = new Verifier();
  await v.postToolUse(post(EXEC, mcpResult({success: true})));
  for (let i = 0; i < 3; i++)
    assert.ok(blocked(await v.stop(stop())), `block ${i + 1} should still gate`);
  assert.ok(!blocked(await v.stop(stop())), 'fourth Stop is allowed through');
  assert.ok(v.exhausted, 'and the turn is flagged unverified for the UI');
});

test('onBlock fires once per block so the runtime can start a hidden revision', async () => {
  let calls = 0;
  const v = new Verifier(() => calls++);
  await v.postToolUse(post(EXEC, mcpResult({success: true})));
  await v.stop(stop());
  await v.stop(stop());
  assert.equal(calls, 2);
});

test('a malformed tool response does not crash the gate', async () => {
  const v = new Verifier();
  await v.postToolUse(post(EXEC, {content: [{text: 'not json'}]}));
  await v.postToolUse(post(EXEC, undefined));
  assert.ok(blocked(await v.stop(stop())), 'unparseable means unverified, not verified');
});

test('non-Stop / non-PostToolUse events are ignored', async () => {
  const v = new Verifier();
  assert.ok(!blocked(await v.stop({hook_event_name: 'PreToolUse'})));
  await v.postToolUse({hook_event_name: 'Stop', tool_name: EXEC});
  assert.ok(!v.hadActions, 'a Stop routed to postToolUse must not register an action');
});
