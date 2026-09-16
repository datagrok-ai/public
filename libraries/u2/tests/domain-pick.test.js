/* `domains.pick` (WO-11) over the memory backend: candidates by the target's name column, a pick
   writes the id, a written id resolves to its name, an extra filter narrows, the name column falls
   back to the business key and then the id — and `PickInput`'s face keeps value and box in step. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {backends} from '../src/sources/backends.js';
import {MemoryDomainBackend} from '../src/sources/memory-domain.js';
import {backend} from './domain-fixtures.mjs';

register('./dg-stub.mjs', import.meta.url);
const {DomainPick} = await import('../src/dg/domain/pick.js');

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const saved = {...backends};
    try {
      await body();
    } finally {
      for (const key of Object.keys(backends))
        delete backends[key];
      Object.assign(backends, saved);
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

scoped('candidates come from the target\'s name column; a pick writes the id, typing clears it', async () => {
  backends.domain = backend();
  const pick = new DomainPick('grit.project', {label: 'Project', debounceMs: 0});
  const input = mount(pick);
  assert.equal(pick.root.dataset.u2, 'domain-pick');
  assert.equal(pick.value.value, null);
  input.focus();
  type(input, 'Gri');
  await flush();
  await flush();
  assert.deepEqual(options(), ['Grit'], 'contains, so "Gr" would list Datagrok too');
  fire(input, 'keydown', {key: 'Enter'});
  await flush();
  assert.equal(pick.value.value, 'p1');
  assert.equal(input.value, 'Grit');
  type(input, 'Gri');
  await flush();
  assert.equal(pick.value.value, 'p1', 'typed text is not a pick: the last pick holds');
  fire(input, 'blur');
  await flush();
  assert.equal(input.value, 'Grit', 'and its name comes back on blur');
  type(input, '');
  await flush();
  assert.equal(pick.value.value, null, 'clearing the text clears the value');
  pick.dispose();
});

scoped('a written id resolves to its name; an unknown id shows itself', async () => {
  backends.domain = backend();
  const pick = new DomainPick('grit.project', {label: 'Project', value: 'p2'});
  const input = mount(pick);
  await flush();
  assert.equal(input.value, 'Datagrok');
  assert.equal(pick.value.value, 'p2', 'resolving did not clear the value');
  pick.value.value = 'nope';
  await flush();
  assert.equal(input.value, 'nope');
  pick.value.value = null;
  await flush();
  assert.equal(input.value, '');
  pick.dispose();
});

scoped('a backend that declares resolveNames answers an id there: one call, no query', async () => {
  const calls = [];
  backends.domain = {
    table: () => assert.fail('the resolver path reads no row'),
    saveAll: () => Promise.resolve(true),
    resolveNames: (table, ids) => {
      calls.push({table, ids: [...ids]});
      return Promise.resolve({p1: 'Grit'});
    },
  };
  assert.deepEqual(await DomainPick.resolve('grit.project', 'p1'), {id: 'p1', name: 'Grit'});
  assert.deepEqual(calls, [{table: 'grit.project', ids: ['p1']}]);
  assert.deepEqual(await DomainPick.resolve('grit.project', 'gone'), {id: 'gone', name: 'gone'},
    'an id the resolver answers null for is shown as itself');
});

scoped('an extra filter narrows the candidates; an empty query lists them all', async () => {
  backends.domain = backend();
  const pick = new DomainPick('grit.project', {filter: 'key = "DG"', debounceMs: 0});
  const input = mount(pick);
  input.focus();
  fire(input, 'keydown', {key: 'ArrowDown'});
  await flush();
  await flush();
  assert.deepEqual(options(), ['Datagrok']);
  pick.dispose();

  const all = new DomainPick('grit.project', {debounceMs: 0});
  const box = mount(all);
  box.focus();
  fire(box, 'keydown', {key: 'ArrowDown'});
  await flush();
  await flush();
  assert.deepEqual(options(), ['Grit', 'Datagrok'], 'in the table\'s own order, not by name');
  all.dispose();
});

scoped('the name column falls back to the business key, then to the id', async () => {
  backends.domain = new MemoryDomainBackend({name: 'x', tables: {
    thing: {businessKey: ['code'], columns: {code: {type: 'string', required: true}}},
    blob: {columns: {size: {type: 'int'}}},
  }}, {rows: {thing: [{id: 't1', code: 'AB'}, {id: 't2', code: 'CD'}], blob: [{id: 'b1', size: 1}]}});
  const table = await backends.domain.table('x.thing');
  assert.equal(DomainPick.nameColumn(table.info), 'code');
  assert.deepEqual(await DomainPick.search('x.thing', 'c'), [{id: 't2', name: 'CD'}]);
  assert.deepEqual(await DomainPick.resolve('x.thing', 't1'), {id: 't1', name: 'AB'});
  const blob = await backends.domain.table('x.blob');
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
  delete backends.domain;
  const pick = new DomainPick('grit.project', {debounceMs: 0});
  const input = mount(pick);
  input.focus();
  type(input, 'a');
  await flush();
  await flush();
  assert.match(document.body.querySelector('.u2-typeahead-error').textContent, /no platform backend/);
  pick.dispose();
});

scoped('ergonomics: the candidates open on focus in the table\'s order with the first highlighted, the placeholder ' +
  'is the caption, ✕ clears', async () => {
    backends.domain = backend();
    const pick = new DomainPick('grit.project', {label: 'Project', debounceMs: 0, value: 'p1'});
    const input = mount(pick);
    await flush();
    assert.equal(input.placeholder, 'Project…');
    assert.equal(input.value, 'Grit');
    input.focus();
    await flush();
    await flush();
    assert.equal(pick.typeAhead.isOpen.value, true, 'open on focus');
    assert.deepEqual(options(), ['Grit', 'Datagrok'], 'the empty query lists the first candidates as the table orders them');
    assert.equal(document.body.querySelector('.u2-typeahead-option-active')?.textContent, 'Grit', 'the first is highlighted');
    const clear = pick.root.querySelector('.u2-input-clear');
    assert.equal(clear.hidden, false);
    fire(clear, 'click');
    await flush();
    assert.equal(pick.value.value, null);
    assert.equal(input.value, '');
    assert.equal(clear.hidden, true);
    pick.dispose();
  });

scoped('stray text never survives a blur: cleared, or the held pick\'s name put back', async () => {
  backends.domain = backend();
  const pick = new DomainPick('grit.project', {label: 'Project', debounceMs: 0});
  const input = mount(pick);
  input.focus();
  type(input, 'zzz');
  await flush();
  await flush();
  fire(input, 'blur');
  await flush();
  assert.equal(input.value, '', 'text that looks chosen is worse than none');
  assert.equal(pick.value.value, null);
  assert.equal(pick.validity.value, null);
  input.focus();
  type(input, 'Grit');
  await flush();
  await flush();
  fire(input, 'keydown', {key: 'Enter'});
  await flush();
  assert.equal(pick.value.value, 'p1');
  fire(input, 'blur');
  await flush();
  assert.equal(input.value, 'Grit', 'a pick stays');
  pick.dispose();

  const held = new DomainPick('grit.project', {label: 'Project', debounceMs: 0, value: 'p2'});
  const box = mount(held);
  await flush();
  assert.equal(box.value, 'Datagrok');
  box.focus();
  type(box, 'Dat');
  await flush();
  fire(box, 'blur');
  await flush();
  assert.equal(held.value.value, 'p2', 'typing over a held value does not drop it');
  assert.equal(box.value, 'Datagrok', 'the text is put back on blur');
  assert.equal(held.validity.value, null);
  held.dispose();
});

scoped('opened on focus, Enter alone commits the highlighted first candidate', async () => {
  backends.domain = backend();
  const pick = new DomainPick('grit.project', {label: 'Project', debounceMs: 0});
  const input = mount(pick);
  input.focus();
  await flush();
  await flush();
  fire(input, 'keydown', {key: 'Enter'});
  await flush();
  assert.equal(pick.value.value, 'p1', 'Grit, the first in the table\'s order');
  assert.equal(input.value, 'Grit');
  pick.dispose();
});

scoped('a blur that puts the held pick back says what was typed and what was kept; the next focus clears it', async () => {
  backends.domain = backend();
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
  assert.equal(pick.value.value, 'p2');
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
  backends.domain = backend();
  const pick = new DomainPick('grit.project', {label: 'Project', debounceMs: 0});
  const input = mount(pick);
  input.focus();
  type(input, 'gr');
  fire(input, 'keydown', {key: 'Enter'});
  assert.equal(pick.typeAhead.isPickPending, true, 'nothing to pick from yet');
  assert.equal(pick.value.value, null);
  await flush();
  await flush();
  assert.equal(pick.value.value, 'p1', 'the first candidate: Grit, in the table\'s order');
  assert.equal(input.value, 'Grit');
  assert.equal(pick.typeAhead.isPickPending, false);

  type(input, 'grit');
  fire(input, 'keydown', {key: 'Enter'});
  fire(input, 'blur');
  await flush();
  await flush();
  assert.equal(pick.value.value, 'p1', 'the exact name wins over the first candidate, and a blur meanwhile does not clear it');
  assert.equal(input.value, 'Grit');
  pick.dispose();
});

scoped('an empty table says "No <rows> yet" for the empty query; a query that finds nothing says "No matches"',
  async () => {
    backends.domain = new MemoryDomainBackend(
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
    backends.domain = backend();
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
