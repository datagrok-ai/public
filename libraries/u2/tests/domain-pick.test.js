/* `domainPick` (WO-11) over the memory backend: candidates by the target's name column, a pick
   writes the id, a written id resolves to its name, an extra filter narrows, the name column falls
   back to the business key and then the id — and `PickInput`'s face keeps value and box in step. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {backends} from '../src/sources/backends.js';
import {MemoryDomainBackend} from '../src/sources/memory-domain.js';
import {domainPick, DomainPick} from '../src/dg/domain/pick.js';
import {backend} from './domain-fixtures.mjs';

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
  const pick = domainPick('grit.project', {label: 'Project', debounceMs: 0});
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
  const pick = domainPick('grit.project', {label: 'Project', value: 'p2'});
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

scoped('an extra filter narrows the candidates; an empty query lists them all', async () => {
  backends.domain = backend();
  const pick = domainPick('grit.project', {filter: 'key = "DG"', debounceMs: 0});
  const input = mount(pick);
  input.focus();
  fire(input, 'keydown', {key: 'ArrowDown'});
  await flush();
  await flush();
  assert.deepEqual(options(), ['Datagrok']);
  pick.dispose();

  const all = domainPick('grit.project', {debounceMs: 0});
  const box = mount(all);
  box.focus();
  fire(box, 'keydown', {key: 'ArrowDown'});
  await flush();
  await flush();
  assert.deepEqual(options(), ['Datagrok', 'Grit'], 'by name');
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
  const pick = domainPick('grit.project', {debounceMs: 0});
  const input = mount(pick);
  input.focus();
  type(input, 'a');
  await flush();
  await flush();
  assert.match(document.body.querySelector('.u2-typeahead-error').textContent, /no platform backend/);
  pick.dispose();
});

scoped('ergonomics: the candidates open on focus with the first highlighted, the placeholder is the caption, ✕ clears',
  async () => {
    backends.domain = backend();
    const pick = domainPick('grit.project', {label: 'Project', debounceMs: 0, value: 'p1'});
    const input = mount(pick);
    await flush();
    assert.equal(input.placeholder, 'Project…');
    assert.equal(input.value, 'Grit');
    input.focus();
    await flush();
    await flush();
    assert.equal(pick.typeAhead.isOpen.value, true, 'open on focus');
    assert.deepEqual(options(), ['Datagrok', 'Grit'], 'the empty query lists the first candidates');
    assert.equal(document.body.querySelector('.u2-typeahead-option-active')?.textContent, 'Datagrok', 'the first is highlighted');
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
  const pick = domainPick('grit.project', {label: 'Project', debounceMs: 0});
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

  const held = domainPick('grit.project', {label: 'Project', debounceMs: 0, value: 'p2'});
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
  const pick = domainPick('grit.project', {label: 'Project', debounceMs: 0});
  const input = mount(pick);
  input.focus();
  await flush();
  await flush();
  fire(input, 'keydown', {key: 'Enter'});
  await flush();
  assert.equal(pick.value.value, 'p2', 'Datagrok, first by name');
  assert.equal(input.value, 'Datagrok');
  pick.dispose();
});

scoped('a fast Enter: the pick is made when the candidates land — the name typed, else the first', async () => {
  backends.domain = backend();
  const pick = domainPick('grit.project', {label: 'Project', debounceMs: 0});
  const input = mount(pick);
  input.focus();
  type(input, 'gr');
  fire(input, 'keydown', {key: 'Enter'});
  assert.equal(pick.typeAhead.isPickPending, true, 'nothing to pick from yet');
  assert.equal(pick.value.value, null);
  await flush();
  await flush();
  assert.equal(pick.value.value, 'p2', 'the first candidate: Datagrok, by name');
  assert.equal(input.value, 'Datagrok');
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
    const pick = domainPick('grit.project', {label: 'Project', debounceMs: 0});
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
