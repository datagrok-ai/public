import {emsBackends} from '../src/core/domain-backend.js';
import {test, testOwner} from '@datagrok-libraries/u2/tests/owned-test.js';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {backends, Scope} from '@datagrok-libraries/u2';
import {MemoryDomainBackend} from '../src/core/memory-domain.js';
import {backend} from './domain-fixtures.mjs';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const {DomainPick} = await import('../src/dg/pick.js');

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const saved = {...backends};
    const savedEms = {...emsBackends};
    try {
      await body();
    } finally {
      for (const key of Object.keys(backends))
        delete backends[key];
      Object.assign(backends, saved);
      for (const key of Object.keys(emsBackends))
        delete emsBackends[key];
      Object.assign(emsBackends, savedEms);
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

function type(input, value) {
  input.value = value;
  fire(input, 'input');
}

function options() {
  return [...document.body.querySelectorAll('.u2-typeahead-option')].map((el) => el.textContent);
}

function mount(pick) {
  document.body.append(pick.root);
  return pick.root.querySelector('input');
}

scoped('candidates come from the target\'s name column; a pick writes the id, only ✕ clears it', async () => {
  emsBackends.domain = backend();
  const pick = new DomainPick('grit.project', {label: 'Project', debounceMs: 0});
  const input = mount(pick);
  assert.equal(pick.root.dataset.u2, 'domain-pick');
  assert.equal(pick.value, null);
  input.focus();
  type(input, 'Gri');
  await flush();
  await flush();
  assert.deepEqual(options(), ['Grit'], 'contains, so "Gr" would list Datagrok too');
  fire(input, 'keydown', {key: 'Enter'});
  await flush();
  assert.equal(pick.value, 'p1');
  assert.equal(input.value, 'Grit');
  type(input, 'Gri');
  await flush();
  assert.equal(pick.value, 'p1', 'typed text is not a pick: the last pick holds');
  fire(input, 'blur');
  await flush();
  assert.equal(input.value, 'Grit', 'and its name comes back on blur');
  type(input, '');
  await flush();
  assert.equal(pick.value, 'p1', 'clearing the text is not a clear: the box only looks empty while focused');
  fire(input, 'blur');
  await flush();
  assert.equal(input.value, 'Grit', 'and the held name comes back on blur, with nothing to say');
  fire(pick.root.querySelector('.u2-input-clear'), 'click');
  await flush();
  assert.equal(pick.value, null, '✕ clears');
  pick.dispose();
});

scoped('a written id resolves to its name; an unknown id shows itself', async () => {
  emsBackends.domain = backend();
  const pick = new DomainPick('grit.project', {label: 'Project', value: 'p2'});
  const input = mount(pick);
  await flush();
  assert.equal(input.value, 'Datagrok');
  assert.equal(pick.value, 'p2', 'resolving did not clear the value');
  pick.value = 'nope';
  await flush();
  assert.equal(input.value, 'nope');
  pick.value = null;
  await flush();
  assert.equal(input.value, '');
  pick.dispose();
});

scoped('a backend that declares resolveNames answers an id there; one it names nothing reads the rows once per tick',
  async () => {
    const calls = [];
    const queries = [];
    const rows = {CHEMBL1: {id: 'CHEMBL1', pref_name: 'ASPIRIN'}, CHEMBL2: {id: 'CHEMBL2', pref_name: 'IBRUTINIB'}};
    emsBackends.domain = {
      table: async () => ({info: {nameColumn: 'pref_name', businessKey: ['chembl_id'], searchableColumns: []}, query: async ({filter}) => {
        queries.push(filter);
        return filter.value.map((id) => rows[id]).filter((r) => r !== undefined);
      }}),
      saveAll: () => Promise.resolve(true),
      resolveNames: (table, ids) => {
        calls.push({table, ids: [...ids]});
        return Promise.resolve({p1: 'Grit'});
      },
    };
    assert.deepEqual(await DomainPick.resolve('grit.project', 'p1'), {id: 'p1', name: 'Grit'});
    assert.deepEqual(calls, [{table: 'grit.project', ids: ['p1']}]);
    const named = await Promise.all(['CHEMBL1', 'CHEMBL2', 'CHEMBL1', 'gone']
      .map((id) => DomainPick.resolve('chembl.molecule', id)));
    assert.deepEqual(named.map((item) => item.name), ['ASPIRIN', 'IBRUTINIB', 'ASPIRIN', 'gone'],
      'an id nobody names is shown as itself');
    assert.deepEqual(queries, [{property: 'id', operator: 'in', value: ['CHEMBL1', 'CHEMBL2', 'gone']}],
      'the ids asked for together are read in one query');
  });

scoped('an extra filter narrows the candidates; an empty query lists them all', async () => {
  emsBackends.domain = backend();
  const pick = new DomainPick('grit.project', {filter: 'key = "DG"', debounceMs: 0});
  const input = mount(pick);
  input.focus();
  fire(input, 'keydown', {key: 'ArrowDown'});
  await flush();
  await flush();
  assert.deepEqual(options(), ['Datagrok']);
  pick.dispose();

  const all = testOwner().build(() => new DomainPick('grit.project', {debounceMs: 0}));
  const box = mount(all);
  box.focus();
  fire(box, 'keydown', {key: 'ArrowDown'});
  await flush();
  await flush();
  assert.deepEqual(options(), ['Datagrok', 'Grit'], 'sorted by what they show');
  all.dispose();
});

scoped('the name column falls back to the business key, then to the id', async () => {
  emsBackends.domain = new MemoryDomainBackend({name: 'x', tables: {
    thing: {businessKey: ['code'], columns: {code: {type: 'string', required: true}}},
    blob: {columns: {size: {type: 'int'}}},
  }}, {rows: {thing: [{id: 't1', code: 'AB'}, {id: 't2', code: 'CD'}], blob: [{id: 'b1', size: 1}]}});
  const table = await emsBackends.domain.table('x.thing');
  assert.equal(DomainPick.nameColumn(table.info), 'code');
  assert.deepEqual(await DomainPick.search('x.thing', 'c'), [{id: 't2', name: 'CD'}]);
  assert.deepEqual(await DomainPick.resolve('x.thing', 't1'), {id: 't1', name: 'AB'});
  const blob = await emsBackends.domain.table('x.blob');
  assert.equal(DomainPick.nameColumn(blob.info), 'id');
  assert.deepEqual(await DomainPick.search('x.blob', ''), [{id: 'b1', name: 'b1'}]);
});

scoped('DomainPick.filter: the like condition, AND-ed with the extra filter as a tree', () => {
  assert.equal(DomainPick.filter('name', '  ', undefined), undefined);
  assert.deepEqual(DomainPick.filter('name', ' ab ', undefined),
    [{property: 'name', operator: 'like', value: '%ab%'}]);
  const extra = {property: 'key', operator: '=', value: 'DG'};
  assert.deepEqual(DomainPick.filter('name', '', extra), [extra]);
  assert.deepEqual(DomainPick.filter('name', 'a', extra),
    [{property: 'name', operator: 'like', value: '%a%'}, 'and', extra]);
  const parsed = DomainPick.filter('name', 'a', 'key = "DG"');
  assert.equal(parsed.length, 3);
  assert.equal(parsed[1], 'and');
  assert.equal(Array.isArray(parsed[2]), true, 'a string filter becomes a nested tree');
  assert.equal(DomainPick.filter('name', '', '   '), undefined);
});

scoped('without a backend the search reports the failure in the popup', async () => {
  delete emsBackends.domain;
  const pick = new DomainPick('grit.project', {debounceMs: 0});
  const input = mount(pick);
  input.focus();
  type(input, 'a');
  await flush();
  await flush();
  assert.match(document.body.querySelector('.u2-typeahead-error').textContent, /no platform backend/);
  pick.dispose();
});

scoped('ergonomics: the candidates open on focus sorted by name with nothing highlighted, no placeholder, ' +
  '✕ clears', async () => {
    emsBackends.domain = backend();
    const pick = new DomainPick('grit.project', {label: 'Project', debounceMs: 0, value: 'p2'});
    const input = mount(pick);
    // the caption is read with the batch of its tick
    await flush();
    await flush();
    assert.equal(input.placeholder, '', 'no placeholder repeating the caption');
    assert.equal(input.value, 'Datagrok');
    input.focus();
    await flush();
    await flush();
    assert.equal(pick.typeAhead.isOpen.value, true, 'open on focus');
    assert.deepEqual(options(), ['Datagrok', 'Grit'], 'the empty query lists the first candidates by name');
    assert.equal(document.body.querySelector('.u2-typeahead-option-active'), null, 'nothing highlighted over a held value');
    const clear = pick.root.querySelector('.u2-input-clear');
    assert.equal(clear.hidden, false);
    fire(clear, 'click');
    await flush();
    assert.equal(pick.value, null);
    assert.equal(input.value, '');
    assert.equal(clear.hidden, true);
    pick.dispose();
  });

scoped('stray text never survives a blur: cleared, or the held pick\'s name put back', async () => {
  emsBackends.domain = backend();
  const pick = new DomainPick('grit.project', {label: 'Project', debounceMs: 0});
  const input = mount(pick);
  input.focus();
  type(input, 'zzz');
  await flush();
  await flush();
  fire(input, 'blur');
  await flush();
  assert.equal(input.value, '', 'text that looks chosen is worse than none');
  assert.equal(pick.value, null);
  assert.equal(pick.validity.value, null);
  input.focus();
  type(input, 'Grit');
  await flush();
  await flush();
  fire(input, 'keydown', {key: 'Enter'});
  await flush();
  assert.equal(pick.value, 'p1');
  fire(input, 'blur');
  await flush();
  assert.equal(input.value, 'Grit', 'a pick stays');
  pick.dispose();

  const held = testOwner().build(() => new DomainPick('grit.project', {label: 'Project', debounceMs: 0, value: 'p2'}));
  const box = mount(held);
  await flush();
  assert.equal(box.value, 'Datagrok');
  box.focus();
  type(box, 'Dat');
  await flush();
  fire(box, 'blur');
  await flush();
  assert.equal(held.value, 'p2', 'typing over a held value does not drop it');
  assert.equal(box.value, 'Datagrok', 'the text is put back on blur');
  assert.equal(held.validity.value, null);
  held.dispose();
});

scoped('opened on focus, Enter alone commits the highlighted first candidate', async () => {
  emsBackends.domain = backend();
  const pick = new DomainPick('grit.project', {label: 'Project', debounceMs: 0});
  const input = mount(pick);
  input.focus();
  await flush();
  await flush();
  fire(input, 'keydown', {key: 'Enter'});
  await flush();
  assert.equal(pick.value, 'p2', 'Datagrok, the first by name');
  assert.equal(input.value, 'Datagrok');
  pick.dispose();
});

scoped('a blur that puts the held pick back says what was typed and what was kept; the next focus clears it', async () => {
  emsBackends.domain = backend();
  const pick = new DomainPick('grit.project', {label: 'Project', debounceMs: 0, value: 'p2'});
  const input = mount(pick);
  await flush();
  const hint = pick.root.querySelector('[data-u2-part="hint"]');
  assert.equal(hint.hidden, true);
  input.focus();
  type(input, 'zzz');
  await flush();
  await flush();
  fire(input, 'blur');
  await flush();
  assert.equal(pick.value, 'p2');
  assert.equal(input.value, 'Datagrok');
  assert.equal(hint.hidden, false);
  assert.equal(hint.textContent, 'No match for "zzz" — kept Datagrok');
  fire(input, 'focus');
  assert.equal(hint.hidden, true, 'cleared on the next focus');
  fire(input, 'blur');
  await flush();
  assert.equal(hint.hidden, true, 'nothing to say when nothing was typed');
  pick.dispose();
});

scoped('a fast Enter: the pick is made when the candidates land — the name typed, else the first', async () => {
  emsBackends.domain = backend();
  const pick = new DomainPick('grit.project', {label: 'Project', debounceMs: 0});
  const input = mount(pick);
  input.focus();
  type(input, 'gr');
  fire(input, 'keydown', {key: 'Enter'});
  assert.equal(pick.typeAhead.isPickPending, true, 'nothing to pick from yet');
  assert.equal(pick.value, null);
  await flush();
  await flush();
  assert.equal(pick.value, 'p2', 'the first candidate: Datagrok, by name');
  assert.equal(input.value, 'Datagrok');
  assert.equal(pick.typeAhead.isPickPending, false);

  type(input, 'grit');
  fire(input, 'keydown', {key: 'Enter'});
  fire(input, 'blur');
  await flush();
  await flush();
  assert.equal(pick.value, 'p1', 'the exact name wins over the first candidate, and a blur meanwhile does not clear it');
  assert.equal(input.value, 'Grit');
  pick.dispose();
});

scoped('an empty table says "No <rows> yet" for the empty query; a query that finds nothing says "No matches"',
  async () => {
    emsBackends.domain = new MemoryDomainBackend(
      {name: 'grit', tables: {project: {friendlyName: 'Projects', columns: {name: {type: 'string', isName: true}}}}},
      {rows: {project: []}});
    const pick = new DomainPick('grit.project', {label: 'Project', debounceMs: 0});
    const input = mount(pick);
    input.focus();
    await flush();
    await flush();
    assert.equal(document.body.querySelector('.u2-typeahead-empty')?.textContent, 'No projects yet');
    type(input, 'x');
    await flush();
    await flush();
    assert.equal(document.body.querySelector('.u2-typeahead-empty')?.textContent, 'No matches');
    pick.dispose();
  });

scoped('a dependent picker: `$name` binds to the sibling and AND-s with the search; unbound → no candidates, a hint',
  async () => {
    emsBackends.domain = backend();
    let row = {project_id: 'p1'};
    const pick = new DomainPick('grit.issue', {label: 'Parent', debounceMs: 0, filter: 'project_id = $project_id',
      params: () => row, siblings: {properties: [{name: 'project_id', type: 'string', friendlyName: 'Project'}]}});
    const input = mount(pick);
    input.focus();
    await flush();
    await flush();
    assert.deepEqual(options(), ['Aspirin', 'Ibuprofen']);
    type(input, 'ibu');
    await flush();
    await flush();
    assert.deepEqual(options(), ['Ibuprofen']);
    row = {project_id: 'p2'};
    type(input, '');
    await flush();
    await flush();
    assert.deepEqual(options(), ['Naproxen'], 'bound afresh at every search');
    row = {project_id: ''};
    type(input, 'a');
    await flush();
    await flush();
    assert.deepEqual(options(), []);
    assert.equal(document.body.querySelector('.u2-typeahead-empty')?.textContent, 'Pick a project first');
    pick.dispose();
  });

scoped('DomainPick.bind: a bound filter is a tree, an empty sibling leaves the name unbound, a tree passes through', () => {
  assert.deepEqual(DomainPick.bind('project_id = $project_id', {project_id: 'p1'}),
    {filter: [{property: 'project_id', operator: '=', value: 'p1'}], unbound: []});
  assert.deepEqual(DomainPick.bind('project_id = $project_id', {project_id: ''}),
    {filter: undefined, unbound: ['project_id']});
  assert.deepEqual(DomainPick.bind('project_id = $project_id and done = $done', {project_id: 'p1'}).unbound, ['done']);
  const tree = {property: 'key', operator: '=', value: 'DG'};
  assert.deepEqual(DomainPick.bind(tree), {filter: tree, unbound: []});
  assert.deepEqual(DomainPick.bind(undefined), {unbound: []});
});
