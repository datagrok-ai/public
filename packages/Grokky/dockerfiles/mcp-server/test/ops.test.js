// The op registry behind the domain tools. Tests run against the compiled dist/.
//
// The important invariant here is the naming convention: the claude-runtime verify gate decides
// whether a call mutated anything by matching the OP NAME against the same read-verb prefixes
// (claude-runtime/src/verify.ts READONLY_NAME_RE). The two regexes live in different containers
// and cannot import each other, so an op named outside the convention silently changes whether
// the runtime demands verification. These tests are the guard against that drift.

import test from 'node:test';
import assert from 'node:assert';
import {DOMAINS, catalog, opMenu, missingParams, isReadonlyOp} from '../dist/ops.js';

const allOps = () => DOMAINS.flatMap((d) => Object.entries(d.ops).map(([name, op]) => ({d, name, op})));

const READ_VERBS = ['whoami', 'list', 'get', 'search', 'read', 'download', 'find'];
const WRITE_VERBS = ['create', 'delete', 'add', 'remove', 'write', 'upload', 'share', 'attach', 'call', 'update'];

test('every op is runnable and documented', () => {
  for (const {d, name, op} of allOps()) {
    assert.equal(typeof op.run, 'function', `${d.tool}.${name} has no run()`);
    assert.ok(op.desc?.length > 10, `${d.tool}.${name} needs a real description`);
  }
});

test('op names follow the read/write verb convention the runtime verify gate depends on', () => {
  for (const {d, name} of allOps()) {
    const verb = [...READ_VERBS, ...WRITE_VERBS].find((v) => name.startsWith(v));
    assert.ok(verb, `${d.tool}.${name} starts with no known verb — the runtime cannot classify it`);
    assert.equal(isReadonlyOp(name), READ_VERBS.some((v) => name.startsWith(v)),
      `${d.tool}.${name} is classified inconsistently with its verb`);
  }
});

test('tool names are unique and namespaced', () => {
  const names = DOMAINS.map((d) => d.tool);
  assert.equal(new Set(names).size, names.length, 'duplicate tool name');
  for (const n of names)
    assert.match(n, /^datagrok_[a-z]+$/);
});

test('catalog exposes every op with its parameters', () => {
  const spaces = DOMAINS.find((d) => d.tool === 'datagrok_spaces');
  const c = catalog(spaces);
  assert.equal(Object.keys(c.operations).length, Object.keys(spaces.ops).length);
  assert.match(c.operations.create_subspace.params.spaceId, /required/);
  assert.ok(!/required/.test(c.operations.create_subspace.params.link), 'link is optional');
  assert.equal(c.operations.list.description, spaces.ops.list.desc);
});

test('opMenu lists every op name — it is what the model sees in the tool description', () => {
  for (const d of DOMAINS) {
    const menu = opMenu(d);
    for (const name of Object.keys(d.ops))
      assert.ok(menu.includes(name), `${d.tool} menu is missing ${name}`);
  }
});

test('missingParams reports only absent required params', () => {
  const spaces = DOMAINS.find((d) => d.tool === 'datagrok_spaces');
  assert.deepEqual(missingParams(spaces.ops.create_subspace, {}), ['spaceId', 'name']);
  assert.deepEqual(missingParams(spaces.ops.create_subspace, {spaceId: 'a'}), ['name']);
  assert.deepEqual(missingParams(spaces.ops.create_subspace, {spaceId: 'a', name: 'b'}), []);
  assert.deepEqual(missingParams(spaces.ops.list, {}), [], 'filter is optional');
});

test('the consolidated surface still covers every operation the split tools had', () => {
  // The 32 tools this replaced, by the domain+op they map onto. A regression here means an
  // operation was dropped in the consolidation rather than renamed.
  const expected = {
    datagrok_functions: ['list', 'get', 'call', 'list_scripts', 'list_queries', 'create_script', 'create_query'],
    datagrok_files: ['list', 'download', 'upload'],
    datagrok_projects: ['list', 'get', 'create', 'delete', 'search', 'list_recent', 'attach_entity',
      'share', 'list_shares'],
    datagrok_spaces: ['list', 'get', 'create', 'delete', 'create_subspace', 'list_children',
      'add_entity', 'remove_entity', 'read_file', 'write_file', 'delete_file'],
    datagrok_platform: ['whoami', 'list_users', 'list_groups', 'list_connections'],
  };
  for (const [tool, ops] of Object.entries(expected)) {
    const d = DOMAINS.find((x) => x.tool === tool);
    assert.ok(d, `missing domain tool ${tool}`);
    for (const op of ops)
      assert.ok(d.ops[op], `${tool} lost operation ${op}`);
  }
});
