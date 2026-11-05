// Verifier gate — the state machine that decides whether a turn may Stop.
// Tests run against the compiled dist/, so `npm test` builds first.

const test = require('node:test');
const assert = require('node:assert');
const {Verifier, isActionTool, isReadonlyTool, bareToolName} = require('../dist/verify.js');

const blocked = (feedback) => feedback !== undefined;

const EXEC = 'datagrok_exec';
const VERIFY = 'datagrok_verify';

test('bareToolName strips the mcp server prefix', () => {
  assert.equal(bareToolName('mcp__datagrok__call_function'), 'call_function');
  assert.equal(bareToolName('mcp__datagrok-browser__datagrok_exec'), 'datagrok_exec');
  assert.equal(bareToolName('mcp__datagrok__datagrok_spaces'), 'datagrok_spaces');
  assert.equal(bareToolName(EXEC), 'datagrok_exec');
  assert.equal(bareToolName('read_file'), 'read_file');
});

test('isReadonlyTool recognises the read-only name prefixes', () => {
  for (const n of ['list_functions', 'get_project', 'search_project', 'read_space_file', 'download_file', 'whoami'])
    assert.ok(isReadonlyTool(n), `${n} should be read-only`);
  assert.ok(isReadonlyTool('datagrok_show_entities'), 'show_entities is a render, not a mutation');
  for (const n of ['create_project', 'upload_file', 'delete_space', 'call_function'])
    assert.ok(!isReadonlyTool(n), `${n} should not be read-only`);
});

// Domain tools are one tool name covering both reads and mutations, so the gate has to read the
// op. Getting this wrong is silent: a mutating call classified read-only skips verification.
test('isReadonlyTool reads the op for domain tools, not the tool name', () => {
  for (const op of ['list', 'get', 'list_children', 'read_file', 'search', 'whoami'])
    assert.ok(isReadonlyTool('datagrok_spaces', {op}), `op ${op} should be read-only`);
  for (const op of ['create', 'delete', 'write_file', 'add_entity', 'remove_entity', 'share'])
    assert.ok(!isReadonlyTool('datagrok_spaces', {op}), `op ${op} should count as an action`);
  assert.ok(isReadonlyTool('datagrok_spaces', {}), 'no op = asking for the catalog');
  assert.ok(isReadonlyTool('datagrok_projects', {op: 'list', args: {filter: 'x'}}));
  assert.ok(!isReadonlyTool('datagrok_projects', {op: 'create', args: {name: 'x'}}));
});

test('isActionTool: exec and mutating mcp tools count, mcp reads do not', () => {
  assert.ok(isActionTool(EXEC));
  assert.ok(isActionTool('mcp__datagrok__datagrok_functions', {op: 'call', args: {name: 'X:y'}}));
  assert.ok(isActionTool('mcp__datagrok__datagrok_functions', {op: 'create_script', args: {script: 'x'}}));
  assert.ok(!isActionTool('mcp__datagrok__datagrok_functions', {op: 'list'}));
  assert.ok(!isActionTool('mcp__datagrok__datagrok_platform', {op: 'whoami'}));
});

test('the agent\'s own tools and read-only browser and view tools are not actions', () => {
  for (const n of ['read_file', 'write_file', 'patch', 'search_files', 'terminal', 'process_manage', 'skills_list',
    'skill_view', 'skill_manage', 'web_search', 'web_extract', 'AskUserQuestion', 'datagrok_show_entities',
    'list_view_functions', 'list_view_widgets', 'get_view_function_result'])
    assert.ok(!isActionTool(n, {}), `${n} should not be an action`);
  for (const n of [EXEC, 'call_view_function', 'mcp__datagrok-view__call_view_function'])
    assert.ok(isActionTool(n, {}), `${n} should be an action`);
});

// Fail-closed: a domain tool whose op we do not recognise must count as an action so the verify
// gate over-asks, rather than letting an unverified mutation through.
test('isActionTool treats an unrecognised domain op as an action', () => {
  assert.ok(isActionTool('mcp__datagrok__datagrok_spaces', {op: 'reticulate'}));
});

test('an unverified action blocks the Stop', () => {
  const v = new Verifier();
  v.afterTool(EXEC, {}, JSON.stringify({success: true}));
  assert.ok(blocked(v.atStop()));
});

test('a turn with no actions never blocks', () => {
  const v = new Verifier();
  v.afterTool('read_file', {}, '');
  v.afterTool('mcp__datagrok__list_functions', {}, JSON.stringify([]));
  assert.ok(!blocked(v.atStop()));
  assert.ok(!v.hadActions);
});

test('a self-verifying exec satisfies the gate without a separate verify round-trip', () => {
  const v = new Verifier();
  v.afterTool(EXEC, {}, JSON.stringify({success: true, verified: {passed: true, observed: 1}}));
  assert.ok(!blocked(v.atStop()), 'in-round-trip verification should count as a pass');
  assert.equal(JSON.parse(v.statsLine()).passes, 1);
});

test('a self-verifying exec that FAILED still blocks', () => {
  const v = new Verifier();
  v.afterTool(EXEC, {}, JSON.stringify({success: true, verified: {passed: false}}));
  assert.ok(blocked(v.atStop()));
});

test('a passing datagrok_verify clears everything pending', () => {
  const v = new Verifier();
  v.afterTool(EXEC, {}, JSON.stringify({success: true}));
  v.afterTool(EXEC, {}, JSON.stringify({success: true}));
  v.afterTool(VERIFY, {}, JSON.stringify({passed: true}));
  assert.ok(!blocked(v.atStop()), 'one passing verify covers all pending actions in the turn');
});

// Regression: an earlier unverified exec used to keep pendingActions > 0 through a later
// self-verified exec, so the turn paid a hidden revision + verify round-trip for nothing.
test('a self-verified exec pass clears earlier pending actions like a datagrok_verify pass', () => {
  const v = new Verifier();
  v.afterTool(EXEC, {}, JSON.stringify({success: true, verified: {passed: false}}));
  v.afterTool(EXEC, {}, JSON.stringify({success: true, verified: {passed: true, observed: 4}}));
  assert.ok(!blocked(v.atStop()), 'the passing self-verify covers the turn');
});

test('a failing verify blocks with a reason that forbids redoing the action', () => {
  const v = new Verifier();
  v.afterTool(EXEC, {}, JSON.stringify({success: true}));
  v.afterTool(VERIFY, {}, JSON.stringify({passed: false}));
  const r = v.atStop();
  assert.ok(blocked(r));
  assert.match(r, /Verification has failed/);
  assert.match(r, /NEVER redo a state-changing action/);
  assert.match(r, /assertion itself was faulty/);
});

test('the first (no-failure) block also forbids redoing the actions', () => {
  const v = new Verifier();
  v.afterTool(EXEC, {}, JSON.stringify({success: true}));
  const r = v.atStop();
  assert.ok(blocked(r));
  assert.match(r, /Do NOT redo the actions/);
});

test('block reasons are framed as internal feedback the model must not quote', () => {
  const v = new Verifier();
  v.afterTool(EXEC, {}, JSON.stringify({success: true}));
  const r = v.atStop();
  assert.match(r, /never mention, quote, or allude to it/);
});

test('the gate gives up after 3 blocks and marks the turn unverified', () => {
  const v = new Verifier();
  v.afterTool(EXEC, {}, JSON.stringify({success: true}));
  for (let i = 0; i < 3; i++)
    assert.ok(blocked(v.atStop()), `block ${i + 1} should still gate`);
  assert.ok(!blocked(v.atStop()), 'fourth Stop is allowed through');
  assert.ok(v.exhausted, 'and the turn is flagged unverified for the UI');
});

test('onBlock fires once per block so the runtime can start a hidden revision', () => {
  let calls = 0;
  const v = new Verifier(() => calls++);
  v.afterTool(EXEC, {}, JSON.stringify({success: true}));
  v.atStop();
  v.atStop();
  assert.equal(calls, 2);
});

test('a malformed tool response does not crash the gate', () => {
  const v = new Verifier();
  v.afterTool(EXEC, {}, 'not json');
  v.afterTool(EXEC, {}, '');
  assert.ok(blocked(v.atStop()), 'unparseable means unverified, not verified');
});

test('a view function call blocks the Stop until a verify passes; listing and reading do not', () => {
  const v = new Verifier();
  v.afterTool('list_view_functions', {query: 'range'}, JSON.stringify({success: true, functions: []}));
  v.afterTool('get_view_function_result', {name: 'getRange', parameters: {}}, JSON.stringify({success: true}));
  assert.ok(!blocked(v.atStop()));
  v.afterTool('call_view_function', {name: 'setRange', parameters: {min: 1}}, JSON.stringify({success: true}));
  assert.ok(blocked(v.atStop()));
  v.afterTool(VERIFY, {assertion: 'return 1', description: 'range set'}, JSON.stringify({passed: true}));
  assert.ok(!blocked(v.atStop()));
});

test('a Datagrok domain tool blocks the Stop by its op', () => {
  const r = new Verifier();
  r.afterTool('mcp__datagrok__datagrok_spaces', {op: 'list'}, JSON.stringify({ok: true}));
  r.afterTool('mcp__datagrok__datagrok_functions', {op: 'get', args: {name: 'Chem:smilesToMw'}}, '{}');
  r.afterTool('mcp__datagrok__datagrok_functions', {}, '{"ops": []}');
  assert.ok(!blocked(r.atStop()));
  const w = new Verifier();
  w.afterTool('mcp__datagrok__datagrok_spaces', {op: 'delete', args: {id: 's'}}, JSON.stringify({ok: true}));
  assert.ok(blocked(w.atStop()));
});

// The agent's browser tools report a call that failed or was aborted as {passed: false, error}.
test('a verify call the browser never answered counts as a failed verify', () => {
  const v = new Verifier();
  v.afterTool(EXEC, {code: 'x'}, JSON.stringify({success: false, error: 'aborted'}));
  v.afterTool(VERIFY, {assertion: 'a', description: 'd'}, JSON.stringify({passed: false, error: 'aborted'}));
  assert.match(v.atStop(), /Verification has failed 1 time/);
});
