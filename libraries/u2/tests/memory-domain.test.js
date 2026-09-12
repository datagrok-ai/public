/* MemoryDomainBackend (WO-7): a schema.json as properties, info and access; queries evaluated by
   the filter feature's own `toMask`, sorted, paged, projected, with the per-row access columns on
   demand; transactions all-or-nothing with `$ref` resolution, the schema's rules and version
   conflicts; the rows as a frame with the writer attached. The cases every backend must agree on
   are listed in docs/domain-backend-contract.md. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {MemoryDomainBackend} from '../src/sources/memory-domain.js';
import {Filters} from '../src/core/filter/index.js';

import {SCHEMA, backend} from './domain-fixtures.mjs';

test('schema: columns become properties, system columns first; info from the declaration', async () => {
  const issue = await backend().table('grit.issue');
  const byName = Object.fromEntries(issue.properties.map((p) => [p.name, p]));
  assert.deepEqual(issue.properties.map((p) => p.name).slice(0, 5),
    ['id', 'version', 'created_on', 'updated_on', 'author_id']);
  assert.equal(byName.author_id.semType, 'User');
  assert.equal(byName.number.propertyType, 'int');
  assert.equal(byName.weight.propertyType, 'double');
  assert.equal(byName.done.propertyType, 'bool');
  assert.equal(byName.due.propertyType, 'datetime');
  assert.equal(byName.tags.propertyType, 'string_list');
  assert.equal(byName.project_id.propertyType, 'string');
  assert.equal(byName.project_id.semType, 'grit.project', 'a ref carries its target as the semType');
  assert.equal(byName.reporter.semType, 'User');
  assert.equal(byName.title.nullable, false);
  assert.equal(byName.description.nullable, true);
  assert.equal(byName.description.editor, 'textarea');
  assert.deepEqual(byName.priority.choices, ['low', 'high']);
  assert.equal(byName.priority.friendlyName, 'Priority');
  assert.equal(byName.id.set, undefined, 'system columns have no setter');
  assert.equal(typeof byName.title.set, 'function');
  assert.deepEqual(issue.info, {nameColumn: 'title', businessKey: ['project_id', 'number'],
    singularName: 'Issue', pluralName: 'Issues'});
  const project = await backend().table('grit.project');
  assert.equal(project.info.nameColumn, 'name', 'the convention fallback');

  const access = await issue.access();
  assert.equal(access.can.delete, true);
  assert.equal(access.fields.id, 'readonly');
  assert.equal(access.fields.author_id, 'readonly');
  assert.equal(access.fields.title, 'editable');
  await assert.rejects(backend().table('grit.nope'), (e) => e.code === 'not-found');
});

test('query: copies of the store, filtered by the grammar or a tree, sorted, paged, projected', async () => {
  const be = backend();
  const issue = await be.table('grit.issue');
  const all = await issue.query();
  assert.deepEqual(all.map((r) => r.title), ['Aspirin', 'Ibuprofen', 'Naproxen']);
  assert.equal(all[0].version, 1);
  all[0].title = 'mutated';
  assert.equal((await issue.query())[0].title, 'Aspirin', 'a query answers copies');

  const titles = async (spec) => (await issue.query(spec)).map((r) => r.title);
  assert.deepEqual(await titles({filter: 'title starts "A"'}), ['Aspirin']);
  assert.deepEqual(await titles({filter: 'done = true'}), ['Aspirin']);
  assert.deepEqual(await titles({filter: 'weight > 1'}), ['Aspirin', 'Naproxen']);
  assert.deepEqual(await titles({filter: 'number = 1 and project_id = "p2"'}), ['Naproxen']);
  assert.deepEqual(await titles({filter: 'priority = null'}), ['Aspirin', 'Ibuprofen']);
  assert.deepEqual(await titles({filter: 'due < "2026-06-01"'}), ['Aspirin']);
  assert.deepEqual(await titles({filter: [{property: 'project_id', operator: '=', value: 'p1'}]}), ['Aspirin', 'Ibuprofen']);
  assert.deepEqual(await titles({filter: Filters.toDomainTree(Filters.group('and', [Filters.cond('done', '=', false)]))}),
    ['Ibuprofen', 'Naproxen']);
  assert.deepEqual(await titles({sort: '!weight'}), ['Naproxen', 'Aspirin', 'Ibuprofen']);
  assert.deepEqual(await titles({sort: 'project_id,!number'}), ['Ibuprofen', 'Aspirin', 'Naproxen']);
  assert.deepEqual(await titles({sort: 'priority'}), ['Naproxen', 'Aspirin', 'Ibuprofen'], 'nulls last ascending');
  assert.deepEqual(await titles({sort: '!priority'}), ['Aspirin', 'Ibuprofen', 'Naproxen'], 'first descending, as Postgres');
  assert.deepEqual(await titles({limit: 2}), ['Aspirin', 'Ibuprofen']);
  assert.deepEqual(await titles({limit: 2, offset: 2}), ['Naproxen']);
  assert.deepEqual(await issue.query({columns: ['id', 'title'], limit: 1}), [{id: 'i1', title: 'Aspirin'}]);
  assert.equal(await issue.count(), 3);
  assert.equal(await issue.count('done = false'), 2);
  await assert.rejects(issue.query({filter: 'title starts'}), (e) => e.code === 'validation');
  await assert.rejects(issue.query({filter: 'nope = 1'}), /Unknown column "nope"/);
});

test('query withAccess: edit and delete are the table\'s answer, share is null, unless the row carries its own', async () => {
  const be = new MemoryDomainBackend(SCHEMA, {
    rows: {issue: [{id: 'i1', project_id: 'p1', title: 'a'}, {id: 'i2', project_id: 'p1', title: 'b', '~can_edit': false},
      {id: 'i3', project_id: 'p1', title: 'c', '~can_share': false}]},
    access: {can: {view: true, insert: false, edit: true, delete: false, share: true},
      fields: {title: 'editable', number: 'readonly'}},
  });
  const issue = await be.table('grit.issue');
  const rows = await issue.query({withAccess: true});
  assert.deepEqual(rows.map((r) => [r['~can_edit'], r['~can_delete'], r['~can_share']]),
    [[true, false, null], [false, false, null], [true, false, false]], 'the server\'s JSON shape off row mode');
  assert.equal('~can_edit' in (await issue.query())[0], false, 'not asked, not answered');
  assert.deepEqual(await issue.access(), {can: {view: true, insert: false, edit: true, delete: false, share: true},
    fields: {title: 'editable', number: 'readonly'}});
});

test('transaction: insert, update with expectedVersion, delete; results per op', async () => {
  const be = backend();
  const issue = await be.table('grit.issue');
  const [inserted, updated, deleted] = await issue.transaction([
    {op: 'insert', table: 'issue', values: {project_id: 'p1', title: 'New'}},
    {op: 'update', table: 'issue', id: 'i2', values: {title: 'Ibuprofen 400'}, expectedVersion: 1},
    {op: 'delete', table: 'issue', id: 'i3'},
  ]);
  assert.match(inserted.id, /^[0-9a-f-]{36}$/);
  assert.equal(inserted.version, 1);
  assert.deepEqual(updated, {id: 'i2', version: 2});
  assert.deepEqual(deleted, {id: 'i3'});
  const rows = await issue.query();
  assert.deepEqual(rows.map((r) => r.title), ['Aspirin', 'Ibuprofen 400', 'New']);
  assert.equal(rows[1].version, 2);
  assert.equal(rows[2].created_on !== undefined && rows[2].updated_on !== undefined, true);
});

test('transaction: a version conflict, a missing row or a required column refuses the WHOLE batch', async () => {
  const be = backend();
  const issue = await be.table('grit.issue');
  await assert.rejects(issue.transaction([
    {op: 'insert', table: 'issue', values: {project_id: 'p1', title: 'New'}},
    {op: 'update', table: 'issue', id: 'i2', values: {title: 'x'}, expectedVersion: 7},
  ]), (e) => e.code === 'version-conflict' && /expected 7/.test(e.message));
  assert.equal(await issue.count(), 3, 'the insert before the conflict did not land');
  await assert.rejects(issue.transaction([{op: 'update', table: 'issue', id: 'nope', values: {}}]),
    (e) => e.code === 'not-found');
  await assert.rejects(issue.transaction([{op: 'insert', table: 'issue', values: {project_id: 'p1'}}]),
    (e) => e.code === 'validation' && /column "title": Value can't be empty/.test(e.message));
  await assert.rejects(issue.transaction([{op: 'update', table: 'issue', id: 'i1', values: {title: null}}]),
    (e) => e.code === 'validation');
  assert.equal((await issue.query({filter: 'id = "i1"'}))[0].version, 1);
});

test('transaction: choices, min and max refuse the batch; author_id and the timestamps are stamped on insert', async () => {
  const issue = await backend({author: 'u7'}).table('grit.issue');
  await assert.rejects(issue.transaction([{op: 'insert', table: 'issue', values: {project_id: 'p1', title: 'x', priority: 'medium'}}]),
    (e) => e.code === 'validation' && /column "priority": Must be one of: low, high/.test(e.message));
  await assert.rejects(issue.transaction([{op: 'update', table: 'issue', id: 'i1', values: {number: 0}}]),
    (e) => e.code === 'validation' && /column "number": Must be at least 1/.test(e.message));
  const [inserted] = await issue.transaction([{op: 'insert', table: 'issue', values: {project_id: 'p1', title: 'x', number: 1}}]);
  const [row] = await issue.query({filter: `id = "${inserted.id}"`});
  assert.equal(row.author_id, 'u7');
  assert.equal(typeof row.created_on, 'string');
  assert.equal(row.updated_on, row.created_on);
  assert.equal((await backend().table('grit.issue').then((t) => t.query()))[0].author_id, 'me', 'the default author');
});

test('frame: the rows as a MemoryFrame with the writer attached; a page appends into the same frame', async () => {
  const issue = await backend().table('grit.issue');
  const frame = await issue.frame({limit: 2, withAccess: true});
  assert.equal(frame.df.rowCount, 2);
  assert.deepEqual(frame.df.columns.names().slice(-4), ['~state', '~can_edit', '~can_delete', '~can_share']);
  assert.equal(frame.df.get('title', 1), 'Ibuprofen');
  assert.equal(frame.df.columns.byName('priority').isNone(0), true, 'a none cell is known as one');
  assert.equal(frame.edit.isDirty.value, false);
  assert.equal(await frame.append({limit: 2, offset: 2}), 1);
  assert.equal(frame.df.rowCount, 3);
  frame.edit.setValue('i3', 'title', 'Naproxen 250');
  assert.equal(frame.df.get('title', 2), 'Naproxen 250');
  assert.equal(frame.edit.changeCount.value, 1);
  assert.equal(await frame.edit.save(), true);
  assert.equal(issue.rows[2].title, 'Naproxen 250');
  frame.dispose();
  const bare = await issue.frame({withAccess: false});
  assert.equal(bare.df.columns.byName('~can_edit'), null, 'no access columns unless asked');
  bare.dispose();
});

test('transaction: $ref names an earlier insert, $$ escapes, a forward ref is refused', async () => {
  const be = backend();
  const project = await be.table('grit.project');
  const issue = await be.table('grit.issue');
  const [p] = await project.transaction([{op: 'insert', table: 'project', ref: 'p', values: {key: 'X', name: '$$100'}}]);
  assert.equal((await project.query({filter: 'key = "X"'}))[0].name, '$100');
  const [, child] = await issue.transaction([
    {op: 'insert', table: 'issue', ref: 'parent', values: {project_id: p.id, title: 'parent'}},
    {op: 'insert', table: 'issue', values: {project_id: p.id, title: 'child', tags: ['$parent']}},
  ]);
  const rows = await issue.query({filter: `id = "${child.id}"`});
  assert.equal(rows[0].tags[0].length, 36, 'resolved element-wise inside a list');
  await assert.rejects(issue.transaction([
    {op: 'insert', table: 'issue', values: {project_id: '$later', title: 'a'}},
    {op: 'insert', table: 'issue', ref: 'later', values: {project_id: 'p1', title: 'b'}},
  ]), (e) => e.code === 'bad-ref' && /forward reference/.test(e.message));
  await assert.rejects(issue.transaction([{op: 'insert', table: 'issue', values: {project_id: '$nope', title: 'a'}}]),
    (e) => e.code === 'bad-ref' && /unknown reference/i.test(e.message));
});
