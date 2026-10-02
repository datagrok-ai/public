/* MemoryDomainBackend (WO-7): a schema.json as properties, info and access; queries evaluated by
   the filter feature's own `toMask`, sorted, paged, projected, with the per-row access columns on
   demand; transactions all-or-nothing with `$ref` resolution, the schema's rules and version
   conflicts; the rows as a frame with the writer attached. The cases every backend must agree on
   are listed in docs/domain-backend-contract.md. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {MemoryDomainBackend, MemoryTable} from '../src/sources/memory-domain.js';
import {Filters} from '../src/core/filter/index.js';

import {SCHEMA, backend, LOCATIONS, TREE, locations} from './domain-fixtures.mjs';

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
    singularName: 'Issue', pluralName: 'Issues', searchableColumns: ['title', 'description'],
    constraints: [{name: 'weight_positive', expr: 'weight >= 0', message: 'Weight is positive'}], refFilters: {},
    permissions: ['escalate'], childTables: [], hierarchy: false, parentColumn: null});
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
  assert.equal(await issue.count({filter: 'done = false'}), 2);
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
  assert.deepEqual(await issue.access(), {can: {escalate: true, view: true, insert: false, edit: true, delete: false, share: true},
    fields: {title: 'editable', number: 'readonly'}}, 'a declared permission is granted unless the option says otherwise');
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
  // a write answers with the system columns it stamped: the platform re-reads the row for them
  assert.equal(typeof inserted.created_on === 'string' && typeof inserted.author_id === 'string', true);
  assert.equal(updated.id, 'i2');
  assert.equal(updated.version, 2);
  assert.equal(typeof updated.updated_on, 'string');
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

test('transaction: $ref names an insert anywhere in the batch, $$ escapes, an unknown ref is refused', async () => {
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

  const [first, later] = await issue.transaction([
    {op: 'insert', table: 'issue', values: {project_id: '$later', title: 'a'}},
    {op: 'insert', table: 'project', ref: 'later', values: {key: 'L', name: 'Later'}},
  ]);
  assert.equal((await issue.query({filter: `id = "${first.id}"`}))[0].project_id, later.id,
    'a forward ref runs after the insert it names; results keep the request order');
  assert.equal((await project.query({filter: 'key = "L"'}))[0].id, later.id);
  await assert.rejects(issue.transaction([{op: 'insert', table: 'issue', values: {project_id: '$nope', title: 'a'}}]),
    (e) => e.code === 'bad-ref' && /unknown reference/i.test(e.message));
  await assert.rejects(issue.transaction([
    {op: 'insert', table: 'issue', ref: 'a', values: {project_id: 'p1', title: '$b'}},
    {op: 'insert', table: 'issue', ref: 'b', values: {project_id: 'p1', title: '$a'}},
  ]), (e) => e.code === 'bad-ref' && /Operation 0: circular reference among refs "a", "b"/.test(e.message));
  await assert.rejects(issue.transaction([{op: 'insert', table: 'nope', values: {}}]), (e) => e.code === 'not-found');
});

test('transaction: any table in one batch, all or nothing across tables; deletes run child-first', async () => {
  const be = backend();
  const project = await be.table('grit.project');
  const issue = await be.table('grit.issue');
  await assert.rejects(issue.transaction([
    {op: 'insert', table: 'grit.project', ref: 'p', values: {key: 'N', name: 'New'}},
    {op: 'insert', table: 'issue', values: {project_id: '$p', title: 'ok'}},
    {op: 'insert', table: 'issue', values: {project_id: '$p'}},
  ]), (e) => e.code === 'validation' && /Operation 2/.test(e.message));
  assert.equal(await project.count(), 2, 'the other table rolled back too');
  assert.equal(await issue.count(), 3);

  await assert.rejects(project.transaction([{op: 'delete', table: 'project', id: 'p1'}]),
    (e) => e.code === 'validation' && /referenced by issue.project_id/.test(e.message), 'a referenced row stays');
  const results = await project.transaction([
    {op: 'delete', table: 'project', id: 'p2'},
    {op: 'delete', table: 'grit.issue', id: 'i3'},
  ]);
  assert.deepEqual(results, [{id: 'p2'}, {id: 'i3'}], 'parent first as requested, child first as run; request order');
  assert.equal(await project.count(), 1);
  assert.equal(issue.history.at(-1).tx_id, project.history.at(-1).tx_id, 'one transaction');
});

test('search: case-insensitive over the searchable columns, the name column by default; count agrees', async () => {
  const be = backend();
  const issue = await be.table('grit.issue');
  assert.deepEqual(issue.info.searchableColumns, ['title', 'description']);
  const titles = async (spec) => (await issue.query(spec)).map((r) => r.title);
  assert.deepEqual(await titles({search: 'PRO'}), ['Ibuprofen', 'Naproxen']);
  assert.deepEqual(await titles({search: 'pro', filter: 'done = false', sort: '!title'}), ['Naproxen', 'Ibuprofen']);
  assert.equal(await issue.count({search: 'pro'}), 2);
  assert.equal(await issue.count({filter: 'done = true', search: 'pro'}), 0);
  const project = await be.table('grit.project');
  assert.deepEqual(project.info.searchableColumns, ['name']);
  assert.deepEqual((await project.query({search: 'grok'})).map((r) => r.key), ['DG']);
  const bare = new MemoryDomainBackend({name: 's', tables: {t: {columns: {x: {type: 'int'}}}}});
  await assert.rejects(bare.tableSync('s.t').query({search: '1'}),
    (e) => e.code === 'validation' && /no searchable column/.test(e.message));
});

test('info: constraints, refFilters, permissions and childTables from the schema; access grants a declared permission', async () => {
  const be = backend();
  const issue = await be.table('grit.issue');
  const project = await be.table('grit.project');
  assert.deepEqual(issue.info.constraints, [{name: 'weight_positive', expr: 'weight >= 0', message: 'Weight is positive'}]);
  assert.deepEqual(issue.info.permissions, ['escalate']);
  assert.deepEqual(project.info.childTables, [{schema: 'grit', table: 'issue', fkColumn: 'project_id', label: 'project_id'}]);
  assert.deepEqual(issue.info.childTables, []);
  assert.equal((await issue.access()).can.escalate, true);
  const denied = backend({access: {can: {view: true, escalate: false}, fields: {}}});
  assert.equal((await denied.table('grit.issue').then((t) => t.access())).can.escalate, false);
  const filtered = new MemoryDomainBackend({name: 's', tables: {
    country: {columns: {name: {type: 'string', isName: true}}},
    city: {columns: {country_id: {type: 'ref', ref: 'country'}, name: {type: 'string'}}},
    person: {columns: {country_id: {type: 'ref', ref: 'country'},
      city_id: {type: 'ref', ref: 'city', filter: 'country_id = $country_id', friendlyName: 'City'}}},
  }});
  assert.deepEqual(filtered.tableSync('s.person').info.refFilters, {city_id: 'country_id = $country_id'});
  assert.deepEqual(filtered.tableSync('s.city').info.childTables, [{schema: 's', table: 'person', fkColumn: 'city_id', label: 'City'}]);
  assert.deepEqual(filtered.tableSync('s.country').info.childTables.map((c) => c.table), ['city', 'person']);
});

test('audit: one line per op under one tx_id, before and after the row', async () => {
  const be = backend();
  const issue = await be.table('grit.issue');
  assert.deepEqual(await issue.audit('i1'), [], 'seeding is not history');
  const [inserted] = await issue.transaction([
    {op: 'insert', table: 'issue', values: {project_id: 'p1', title: 'New'}},
    {op: 'update', table: 'issue', id: 'i1', values: {title: 'Aspirin 100'}},
    {op: 'delete', table: 'issue', id: 'i2'},
  ]);
  const [update] = await issue.audit('i1');
  assert.deepEqual([update.op, update.actor_id, update.before.title, update.after.title, update.after.version],
    ['update', null, 'Aspirin', 'Aspirin 100', 2]);
  const [insert] = await issue.audit(inserted.id);
  assert.deepEqual([insert.op, insert.before, insert.after.title], ['insert', null, 'New']);
  const [del] = await issue.audit('i2');
  assert.deepEqual([del.op, del.before.title, del.after], ['delete', 'Ibuprofen', null]);
  assert.equal(new Set([update, insert, del].map((h) => h.tx_id)).size, 1);
  assert.equal(typeof update.ts, 'string');
});

test('saveAll: every writer\'s batch as one transaction, the results sliced back to each', async () => {
  const be = backend();
  const project = await be.table('grit.project');
  const issue = await be.table('grit.issue');
  const projects = await project.frame({});
  const issues = await issue.frame({});
  const p = projects.edit.newRow({key: 'S', name: 'Session'});
  issues.edit.newRow({project_id: p, title: 'Child of a draft'});
  issues.edit.setValue('i1', 'title', 'Aspirin 100');
  const ops = issues.edit.buildOps();
  assert.deepEqual(ops.map((x) => x.op.op), ['update', 'insert']);
  assert.deepEqual(ops[1].op, {op: 'insert', table: 'grit.issue', ref: ops[1].row.id,
    values: {project_id: `$${p}`, title: 'Child of a draft'}}, 'the draft reference is the server\'s $ref');
  assert.equal(await be.saveAll([projects.edit, issues.edit]), true);
  assert.equal(projects.edit.isDirty.value, false);
  assert.equal(issues.edit.isDirty.value, false);
  const saved = projects.df.rows.at(-1);
  assert.match(saved.id, /^[0-9a-f-]{36}$/);
  assert.equal(issues.df.rows.at(-1).project_id, saved.id, 'the child\'s FK is the parent\'s real id');
  assert.equal(issue.rows.at(-1).project_id, saved.id);
  assert.equal(issue.history.at(-1).tx_id, project.history.at(-1).tx_id, 'one transaction');
  projects.dispose();
  issues.dispose();
});

test('saveAll: every writer is closed for the whole transaction — an edit made meanwhile is refused', async () => {
  const be = backend();
  const issue = await be.table('grit.issue');
  const issues = await issue.frame({});
  issues.edit.setValue('i1', 'title', 'First');
  const saving = be.saveAll([issues.edit]);
  assert.equal(issues.edit.isSaving.value, true, 'closed before the first await');
  issues.edit.setValue('i1', 'title', 'Second');
  issues.edit.markDeleted('i2');
  assert.throws(() => issues.edit.newRow({project_id: 'p1', title: 'Third'}),
    /cannot add a row while the batch is being saved/);
  assert.equal(await saving, true);
  assert.equal(issues.edit.isSaving.value, false);
  assert.equal(issue.rows[0].title, 'First', 'what was sent is what landed');
  assert.equal(issues.df.get('title', 0), 'First', 'and the frame says so — no edit was lost silently');
  assert.equal(issues.df.rowCount, 3, 'the delete was refused too');
  assert.equal(issues.edit.isDirty.value, false);
  issues.edit.setValue('i1', 'title', 'Second');
  assert.equal(issues.edit.isDirty.value, true, 'and the writer takes edits again once it is open');
  issues.dispose();
});

test('transaction: a restore undoes a landed delete — one version, one `undelete` audit line', async () => {
  const be = backend();
  const issue = await be.table('grit.issue');
  await issue.transaction([{op: 'delete', table: 'issue', id: 'i2'}]);
  const [result] = await issue.transaction([{op: 'restore', table: 'issue', id: 'i2'}]);
  assert.deepEqual([result.id, result.version], ['i2', 3]);
  assert.deepEqual((await issue.query()).map((r) => r.title), ['Aspirin', 'Ibuprofen', 'Naproxen']);
  const [undelete] = (await issue.audit('i2')).slice(-1);
  assert.deepEqual([undelete.op, undelete.before.is_deleted, undelete.after.is_deleted],
    ['undelete', true, false]);
  await assert.rejects(issue.transaction([{op: 'restore', table: 'issue', id: 'i2'}]),
    (e) => e.code === 'not-found' && /no deleted row "i2"/.test(e.message));
});

test('transaction: a parent and its child restored in one batch run parent-first, whatever the order', async () => {
  const be = backend();
  const project = await be.table('grit.project');
  const issue = await be.table('grit.issue');
  await issue.transaction([{op: 'delete', table: 'issue', id: 'i3'}]);
  await project.transaction([{op: 'delete', table: 'project', id: 'p2'}]);
  await be.transaction([{op: 'restore', table: 'issue', id: 'i3'},
    {op: 'restore', table: 'project', id: 'p2'}]);
  assert.equal(await issue.count(), 3, 'the child listed first still landed');
  assert.equal(await project.count(), 2);
});

test('frame: the mask columns a spec projects — ~state always, ~can_* and ~is_deleted on demand', async () => {
  const issue = await backend().table('grit.issue');
  const plain = await issue.frame({});
  assert.deepEqual(plain.df.columns.names().slice(-1), ['~state']);
  plain.dispose();
  const masked = await issue.frame({deleted: 'only', withAccess: true});
  assert.deepEqual(masked.df.columns.names().slice(-5),
    ['~state', '~can_edit', '~can_delete', '~can_share', '~is_deleted']);
  masked.dispose();
});

test('captions: a ref column projects the target\'s display name; the refusals are the server\'s', async () => {
  const be = backend({rows: {project: [{id: 'p1', key: 'GRIT', name: 'Grit'}],
    issue: [{id: 'i1', project_id: 'p1', number: 1, title: 'Aspirin'},
      {id: 'i2', project_id: null, number: 2, title: 'Orphan'}]}});
  const issue = await be.table('grit.issue');
  const rows = await issue.query({captions: ['project_id']});
  assert.deepEqual(rows.map((r) => r['~caption_project_id']), ['Grit', null],
    'the target\'s name, and null where there is no target');
  assert.equal(rows[0]['~caption_project_id'], 'Grit');
  const frame = await issue.frame({captions: ['project_id']});
  assert.equal(frame.df.columns.names().at(-1), '~caption_project_id');
  frame.dispose();

  const refused = async (captions, message) => await assert.rejects(issue.query({captions}),
    (e) => e.code === 'filter' && e.message === message);
  await refused(['title'], 'Unknown or inaccessible caption column "title"');
  await refused(['nope'], 'Unknown or inaccessible caption column "nope"');
  await refused(['project_id.name'], 'Nested caption "project_id.name" is not supported');
  await refused(['project_id', 'project_id'], 'Duplicate caption "project_id"');
});

test('seq: one bump per transaction per table, none when the transaction throws', async () => {
  const be = backend();
  const issue = be.tableSync('grit.issue');
  const project = be.tableSync('grit.project');
  assert.deepEqual([issue.seq, project.seq], [0, 0]);
  assert.deepEqual(await issue.probe(), {count: -1, last: '0'}, 'an unscoped read is the token alone');
  await issue.transaction([{op: 'update', table: 'issue', id: 'i1', values: {title: 'A'}},
    {op: 'update', table: 'issue', id: 'i2', values: {title: 'B'}}]);
  assert.deepEqual([issue.seq, project.seq], [1, 0], 'one bump for two ops on one table');
  assert.deepEqual(await issue.probe(), {count: -1, last: '1'});
  assert.equal((await issue.probe({deleted: 'include'})).count, 3, 'a scoped read still counts');
  await assert.rejects(issue.transaction([{op: 'update', table: 'issue', id: 'nope', values: {title: 'C'}}]));
  assert.equal(issue.seq, 1, 'a transaction that threw bumped nothing');

  await issue.transaction([{op: 'delete', table: 'issue', id: 'i3'}]);
  await issue.transaction([{op: 'delete', table: 'project', id: 'p2'}]);
  assert.deepEqual([issue.seq, project.seq], [2, 1],
    'the delete scanned the referencing child — a read is not a write');
  await issue.transaction([{op: 'restore', table: 'project', id: 'p2'}]);
  await issue.transaction([{op: 'restore', table: 'issue', id: 'i3'}]);
  assert.deepEqual([issue.seq, project.seq], [3, 2], 'the restore read the parent — a read is not a write');
});

test('batch: an invalid row is reported once — as an error, never also as a duplicate', async () => {
  const issue = await backend().table('grit.issue');
  // the second row has no title and repeats the first one's business key: it is one row of the
  // report, not two
  const payload = [{project_id: 'p1', number: 9, title: 'Fine'}, {project_id: 'p1', number: 9}];
  const dry = await issue.validate(payload, {allOrNothing: false});
  assert.deepEqual(dry.rows.map((r) => r.predicted), ['insert', 'error']);
  const report = await issue.batch(payload, {allOrNothing: false});
  assert.equal(report.errorCount, 1);
  assert.equal(report.skipped, 0, 'a row that lands nowhere does not clash with the one that does');
  assert.deepEqual(report.rows.map((r) => [r.index, r.status]), [[1, 'error'], [0, 'inserted']]);
});

test('updateWhere: the filter selects, the cap limits, hasMore says the filter matched past it', async () => {
  const issue = await backend().table('grit.issue');
  assert.deepEqual(await issue.updateWhere('id in ("i1", "i3")', {priority: 'low'}),
    {updated: 2, hasMore: false});
  const rows = await issue.query();
  assert.deepEqual(rows.map((r) => r.priority), ['low', undefined, 'low']);
  assert.deepEqual(rows.map((r) => r.version), [2, 1, 2], 'only the matched rows were written');
  const audit = await issue.audit('i1');
  assert.deepEqual(audit.map((a) => a.op), ['update'], 'one audit line per row');

  assert.deepEqual(await issue.updateWhere('done = false', {priority: 'high'}, {limit: 1}),
    {updated: 1, hasMore: true}, 'the cap stops at limit and reports the rest');
  assert.deepEqual(await issue.updateWhere('done = false', {priority: 'low'}, {limit: 9999}),
    {updated: 2, hasMore: false}, 'and a limit past the cap is clamped down to it');
  assert.deepEqual(await issue.updateWhere('title starts "zzz"', {priority: 'high'}),
    {updated: 0, hasMore: false});
});

test('updateWhere: the Edit predicate narrows the selection silently', async () => {
  const rows = {project: [], issue: [
    {id: 'i1', project_id: 'p1', number: 1, title: 'Mine', '~can_edit': true},
    {id: 'i2', project_id: 'p1', number: 2, title: 'Theirs', '~can_edit': false},
  ]};
  const issue = new MemoryDomainBackend(SCHEMA, {rows}).tableSync('grit.issue');
  assert.deepEqual(await issue.updateWhere('number >= 1', {priority: 'low'}), {updated: 1, hasMore: false},
    'the row the caller may not edit is left out, not refused');
  assert.deepEqual((await issue.query()).map((r) => r.priority), ['low', undefined]);

  const readOnly = new MemoryDomainBackend(SCHEMA,
    {rows, access: {can: {view: true, edit: false}, fields: {priority: 'editable'}}}).tableSync('grit.issue');
  assert.deepEqual(await readOnly.updateWhere('number >= 1', {priority: 'low'}), {updated: 1, hasMore: false},
    'a row\'s own ~can_edit is the truth the table-level denial cannot override');
});

test('ancestors: root-first, without the row itself; a non-hierarchy table does not answer it at all', async () => {
  const location = locations();
  assert.equal(location.info.hierarchy, true);
  assert.equal(location.info.parentColumn, 'parent_id');
  assert.equal(location.support.ancestors, true);
  assert.deepEqual(await location.ancestors('l1'), [], 'a root has none');
  assert.deepEqual(await location.ancestors('l9'), []);
  assert.deepEqual(await location.ancestors('nope'), [], 'a row the caller cannot see answers no path');

  // declared, never guessed: the member is absent and `support` says why, so a caller refuses by
  // name instead of calling something that would throw
  const issue = await backend().table('grit.issue');
  assert.equal(issue.support.ancestors, false);
  assert.equal(issue.ancestors, undefined);
  assert.throws(() => new MemoryDomainBackend({name: 's', tables: {t: {hierarchy: true,
    columns: {x: {type: 'int'}}}}}), (e) => /invalid-hierarchy/.test(e.message));
});

test('support: what the memory backend declares, and every optional member installed to match', async () => {
  const issue = await backend().table('grit.issue');
  assert.deepEqual(issue.support, {systemColumns: ['id', 'version', 'created_on', 'updated_on', 'author_id'],
    writes: true, deleted: true, restore: true, audit: true, ancestors: false, probe: true, version: true,
    watch: false});
  for (const [member, flag] of [['restore', 'restore'], ['ancestors', 'ancestors'], ['updateWhere', 'writes'],
    ['batch', 'writes'], ['probe', 'probe'], ['audit', 'audit']])
    assert.equal(issue[member] !== undefined, issue.support[flag], `${member} follows support.${flag}`);
  assert.equal(locations().support.ancestors, true, 'a hierarchy table declares the walk it can do');
});

test('ancestors: the chain stops at an ancestor out of sight, at a cycle and at the depth cap', async () => {
  // no-oracle: an ancestor the caller cannot see is simply not there, and the chain ends where
  // it stops resolving
  const orphan = locations({rows: {location: [...TREE.location, {id: 'l5', name: 'Orphan', parent_id: 'hidden'}]}});
  assert.deepEqual(await orphan.ancestors('l5'), []);

  const cyclic = locations();
  await cyclic.transaction([{op: 'update', table: 'location', id: 'l1', values: {parent_id: 'l4'}}]);
  assert.deepEqual((await cyclic.ancestors('l4')).map((a) => a.id), ['l1', 'l2', 'l3'],
    'a cycle terminates at the row it came back to');

  const deep = new MemoryDomainBackend(LOCATIONS, {rows: {location:
    Array.from({length: 80}, (_, i) => ({id: `d${i}`, name: `L${i}`, ...(i === 0 ? {} : {parent_id: `d${i - 1}`})}))}})
    .tableSync('stock.location');
  const chain = await deep.ancestors('d79');
  assert.equal(chain.length, MemoryTable.maxPathDepth);
  assert.equal(chain.at(-1).id, 'd78', 'root-first, truncated at the far end');
});
