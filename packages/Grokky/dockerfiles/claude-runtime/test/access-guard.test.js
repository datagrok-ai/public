// makeAccessGuard — PreToolUse hook: blocks absolute /workspace paths, `..` path
// segments, and everything under /users except the session user's own directory (audit F-16).

const test = require('node:test');
const assert = require('node:assert');
const {makeAccessGuard} = require('../dist/query-options.js');

const ME = '/users/a1b2-c3d4';
const g = makeAccessGuard(ME);
const run = (guard, toolInput) => guard({hook_event_name: 'PreToolUse', tool_name: 'x', tool_input: toolInput});
const blocked = (r) => r.decision === 'block';

test('allows own dir and cwd-relative workspace paths', async () => {
  assert.ok(!blocked(await run(g, {file_path: `${ME}/agents/skills/x.md`})));
  assert.ok(!blocked(await run(g, {pattern: '**/*.ts', path: `${ME}/workspace/packages`})));
  assert.ok(!blocked(await run(g, {file_path: 'workspace/packages/Chem/x.ts'})));
  assert.ok(!blocked(await run(g, {command: 'echo hi'})));
  assert.ok(!blocked(await run(g, {command: 'git log HEAD..main && echo {1..5} ...done'})));
  assert.ok(!blocked(await run(makeAccessGuard(undefined), {file_path: '/workspace/packages/Chem/x.ts'})));
});

test('blocks other users, bare /users, and absolute /workspace', async () => {
  assert.ok(blocked(await run(g, {file_path: '/users/other-user/agents/secret.md'})));
  assert.ok(blocked(await run(g, {file_path: `${ME}x/f`})));
  assert.ok(blocked(await run(g, {command: 'ls /users'})));
  assert.ok(blocked(await run(g, {command: 'grep -r foo /workspace'})));
  assert.ok(blocked(await run(g, {file_path: '/data/foo/workspace/x.ts'})));
  assert.ok(blocked(await run(g, {command: 'ls /workspace;true'})));
  assert.ok(blocked(await run(g, {command: `cp ${ME}/workspace/a /workspace/b`})));
  assert.ok(blocked(await run(g, {command: 'cat ../other-user/agents/f'})));
  assert.ok(blocked(await run(g, {pattern: 'secret', path: '..'})));
  assert.ok(blocked(await run(makeAccessGuard(undefined), {file_path: '/users/anyone/x'})));
});

test('non-PreToolUse events pass through', async () => {
  const r = await g({hook_event_name: 'PostToolUse', tool_name: 'x', tool_input: {command: 'ls /users'}});
  assert.ok(!blocked(r));
});
