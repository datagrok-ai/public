/* Popup and picker inputs on the u2 ↔ bdd contract (WO-C10): popups are keyed — a keystroke keeps
   the option elements that still match and builds only those it brings in — options carry
   role/data-u2-item/aria-selected/item-label, the state block reports open/busy/items, revealItem
   opens the popup over everything, and the events pass the control. */

import {test, testOwner} from './owned-test.js';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {state} from 'datagrok-api/u2core';
import {Scope, Combobox, DateInput} from '../src/index.js';
import {MultiSelect} from '../src/components/inputs/multi-select.js';
import {TypeAhead} from '../src/components/inputs/typeahead.js';
import {SuggestInput} from '../src/components/inputs/suggest-input.js';
import {TagsInput} from '../src/components/inputs/tags-input.js';
import {MapInput} from '../src/components/inputs/map-input.js';
import {RangeSlider} from '../src/components/inputs/range-slider.js';
import {MessageInput} from '../src/components/inputs/message-input.js';

function smoke(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

function mount(component) {
  document.body.append(component.root);
  return component;
}

function type(input, text) {
  input.value = text;
  fire(input, 'input');
}

const options = (popup, cls) => [...popup.querySelectorAll(cls)];

smoke('MultiSelect: a filter keystroke over 500 items builds no option, clearing it builds only the hidden ones', async () => {
  const items = Array.from({length: 500}, (_, i) => `item ${i}`);
  const ms = mount(new MultiSelect({label: 'Items', items}));
  fire(ms.root.querySelector('[role="combobox"]'), 'click');
  await flush();
  const popup = document.body.querySelector('.u2-multi-select-popup');
  const before = new Set(options(popup, '.u2-multi-select-option'));
  assert.equal(before.size, 500);

  const filter = popup.querySelector('.u2-multi-select-filter-input');
  type(filter, 'item 4');
  await flush();
  const narrowed = options(popup, '.u2-multi-select-option');
  assert.equal(narrowed.length, 111, 'item 4, item 40-49, item 400-499');
  assert.equal(narrowed.filter((el) => !before.has(el)).length, 0, 'a narrowing keystroke builds nothing');

  type(filter, 'item 4');
  type(filter, '');
  await flush();
  const all = options(popup, '.u2-multi-select-option');
  assert.equal(all.length, 500);
  assert.equal(all.filter((el) => !narrowed.includes(el)).length, 389, 'only the options that came back are built');
  assert.equal(all.indexOf(narrowed[0]), 4, 'kept options move to their place');
  ms.dispose();
});

smoke('MultiSelect: options carry role, key, selection and an item-label; state reports open, items and selection', async () => {
  const ms = mount(new MultiSelect({label: 'Props', items: ['MW', {value: 'hbd', label: 'H-bond donors'}],
    value: ['hbd']}));
  assert.equal(state(ms.root).open, false);
  assert.deepEqual(state(ms.root).selection, ['hbd']);
  assert.equal(state(ms.root).text, 'H-bond donors');
  assert.deepEqual(state(ms.root).items.list.map((i) => [i.key, i.label, i.selected, i.rendered]),
    [['MW', 'MW', false, false], ['hbd', 'H-bond donors', true, false]]);

  const el = await ms.revealItem('H-bond donors');
  assert.equal(state(ms.root).open, true);
  assert.equal(el.dataset.u2Item, 'hbd');
  assert.equal(el.getAttribute('role'), 'option');
  assert.equal(el.getAttribute('aria-selected'), 'true');
  assert.equal(el.querySelector('[data-u2-part="item-label"]').textContent, 'H-bond donors');
  assert.equal(state(ms.root).items.list[1].rendered, true);
  ms.dispose();
});

smoke('Combobox: keyed options, stable ids for aria-activedescendant, state and revealItem over everything', async () => {
  const combobox = mount(new Combobox({items: ['alpha', 'beta', 'gamma', 'delta'], label: 'Greek'}));
  const input = combobox.root.querySelector('input');
  input.focus();
  assert.equal(state(combobox.root).label, 'Greek');
  assert.equal(state(combobox.root).open, false);
  assert.equal(state(combobox.root).busy, false);

  fire(input, 'keydown', {key: 'ArrowDown'});
  await flush();
  assert.equal(state(combobox.root).open, true);
  const popup = document.body.querySelector('.u2-combobox-popup');
  const beta = popup.querySelector('[data-u2-item="beta"]');
  assert.equal(beta.getAttribute('role'), 'option');
  assert.equal(beta.querySelector('[data-u2-part="item-label"]').textContent, 'beta');

  type(input, 'ta');
  await flush();
  assert.deepEqual(options(popup, '.u2-combobox-option').map((el) => el.dataset.u2Item), ['beta', 'delta']);
  assert.equal(popup.querySelector('[data-u2-item="beta"]'), beta, 'the matching option is kept, not rebuilt');
  fire(input, 'keydown', {key: 'ArrowDown'});
  assert.equal(input.getAttribute('aria-activedescendant'), beta.id);
  assert.equal(beta.getAttribute('aria-selected'), 'true');

  const items = state(combobox.root).items;
  assert.equal(items.count, 4, 'items are every option, not the filtered popup');
  assert.deepEqual(items.list.map((i) => i.rendered), [false, true, false, true]);

  type(input, 'zz');
  await flush();
  assert.equal(popup.querySelector('.u2-combobox-empty')?.textContent, 'No matches');
  type(input, 'ta');
  await flush();
  assert.equal(popup.querySelector('.u2-combobox-empty'), null, 'a hidden status row keeps no row class');

  fire(input, 'keydown', {key: 'Escape'});
  const gamma = await combobox.revealItem('Gamma');
  assert.equal(gamma?.dataset.u2Item, 'gamma', 'revealItem opens the popup unfiltered and matches the label');
  assert.equal(options(popup, '.u2-combobox-option').length, 4);
  combobox.dispose();
});

smoke('Combobox: change passes the value and the control; Enter and Escape go to onEnter/onEscape while closed', async () => {
  const seen = [];
  const keys = [];
  const combobox = mount(new Combobox({items: ['alpha', 'beta'],
    onEnter: (e, c) => keys.push(['enter', e.key, c === combobox]),
    onEscape: (e, c) => keys.push(['escape', e.key, c === combobox])}));
  combobox.on('change', (value, source) => seen.push([value, source === combobox]));
  const input = combobox.root.querySelector('input');
  input.focus();
  fire(input, 'keydown', {key: 'ArrowDown'});
  fire(input, 'keydown', {key: 'ArrowDown'});
  fire(input, 'keydown', {key: 'Enter'});
  assert.deepEqual(seen, [['alpha', true]]);
  assert.deepEqual(keys, [], 'Enter in the open popup picks');

  const enter = fire(input, 'keydown', {key: 'Enter', cancelable: true});
  fire(input, 'keydown', {key: 'Escape'});
  assert.deepEqual(keys, [['enter', 'Enter', true], ['escape', 'Escape', true]]);
  assert.equal(enter, false, 'a handled Enter is kept from the form around');
  assert.equal(combobox.value, 'alpha', 'an Escape with onEscape does not clear');
  combobox.dispose();
});

smoke('TypeAhead: an async source marks the root busy while it fetches; keystrokes keep matching rows', async () => {
  const calls = [];
  const typeahead = mount(new TypeAhead({
    source: (query) => new Promise((resolve) => calls.push({query, resolve})),
    itemText: (x) => x, placeholder: 'Pick', debounceMs: 0}));
  const input = typeahead.root.querySelector('input');
  input.focus();
  assert.equal(state(typeahead.root).busy, false);
  assert.equal(state(typeahead.root).label, 'Pick', 'the placeholder names an unlabelled type-ahead');
  type(input, 'a');
  assert.equal(typeahead.root.getAttribute('aria-busy'), 'true');
  assert.equal(state(typeahead.root).busy, true);
  await flush();
  calls[0].resolve(['alpha', 'gamma', 'beta']);
  await flush();
  assert.equal(typeahead.root.hasAttribute('aria-busy'), false);
  assert.equal(state(typeahead.root).busy, false);
  assert.equal(state(typeahead.root).open, true);
  assert.deepEqual(state(typeahead.root).items.list.map((i) => i.key), ['alpha', 'gamma', 'beta']);
  typeahead.dispose();
});

smoke('TypeAhead: sync rows are kept by key; revealItem lists everything; change passes the item and the control', async () => {
  const seen = [];
  const typeahead = mount(new TypeAhead({source: ['alpha', 'beta', 'gamma'], itemText: (x) => x}));
  typeahead.on('change', (item, source) => seen.push([item, source === typeahead]));
  const input = typeahead.root.querySelector('input');
  input.focus();
  fire(input, 'keydown', {key: 'ArrowDown'});
  await flush();
  const popup = document.body.querySelector('.u2-typeahead-popup');
  const gamma = popup.querySelector('[data-u2-item="gamma"]');
  assert.equal(gamma.getAttribute('role'), 'option');
  assert.ok(gamma.querySelector('[data-u2-part="item-label"]'));
  type(input, 'mm');
  await flush();
  assert.deepEqual(options(popup, '.u2-typeahead-option').map((el) => el.dataset.u2Item), ['gamma']);
  assert.equal(popup.querySelector('[data-u2-item="gamma"]'), gamma);

  fire(input, 'keydown', {key: 'Escape'});
  const beta = await typeahead.revealItem('beta');
  assert.equal(beta?.dataset.u2Item, 'beta');
  fire(beta, 'pointerdown');
  assert.deepEqual(seen, [['beta', true]]);
  assert.equal(state(typeahead.root).value, 'beta');
  assert.equal(state(typeahead.root).text, 'beta');
  typeahead.dispose();
});

smoke('TypeAhead: revealItem lists everything past minChars; the next keystroke brings the hint back (P109)', async () => {
  const calls = [];
  const typeahead = mount(new TypeAhead({
    source: (query) => new Promise((resolve) => calls.push({query, resolve})),
    itemText: (x) => x, minChars: 2, debounceMs: 0}));
  const input = typeahead.root.querySelector('input');
  input.focus();
  const revealed = typeahead.revealItem('beta');
  await flush();
  assert.deepEqual(calls.map((c) => c.query), [''], 'the reveal asks for the empty query');
  calls[0].resolve(['alpha', 'beta']);
  const beta = await revealed;
  assert.equal(beta?.dataset.u2Item, 'beta');
  const popup = document.body.querySelector('.u2-typeahead-popup');
  assert.equal(popup.querySelector('.u2-typeahead-hint'), null, 'no hint row during the reveal');

  type(input, 'a');
  await flush();
  assert.ok(popup.querySelector('.u2-typeahead-hint'), 'one character is under minChars again');
  assert.equal(options(popup, '.u2-typeahead-option').length, 0);
  typeahead.dispose();

  const sync = mount(testOwner().build(() => new TypeAhead({source: ['alpha', 'beta'], itemText: (x) => x, minChars: 1})));
  assert.equal((await sync.revealItem('alpha'))?.dataset.u2Item, 'alpha', 'a sync source too');
  sync.dispose();
});

smoke('SuggestInput: busy while the source fetches, open with the popup', async () => {
  const calls = [];
  const input = mount(new SuggestInput({label: 'City',
    source: (query) => new Promise((resolve) => calls.push({query, resolve})), debounceMs: 0}));
  const box = input.root.querySelector('input');
  box.focus();
  type(box, 'Ber');
  assert.equal(input.root.getAttribute('aria-busy'), 'true');
  assert.equal(state(input.root).busy, true);
  assert.equal(state(input.root).open, true);
  await flush();
  calls[0].resolve(['Berlin', 'Bern']);
  await flush();
  assert.equal(state(input.root).busy, false);
  assert.equal(input.root.hasAttribute('aria-busy'), false);
  assert.deepEqual(state(input.root).items.list.map((i) => [i.key, i.rendered]), [['Berlin', true], ['Bern', true]]);
  input.dispose();
});

smoke('TagsInput: chips are list items kept by text; the state lists them', async () => {
  const tags = mount(new TagsInput({label: 'Tags', items: ['red', 'green', 'blue'], value: ['red']}));
  const red = tags.root.querySelector('[data-u2-item="red"]');
  assert.equal(red.getAttribute('role'), 'listitem');
  assert.equal(red.querySelector('[data-u2-part="item-label"]').textContent, 'red');
  tags.add('blue');
  assert.equal(tags.root.querySelector('[data-u2-item="red"]'), red, 'adding a tag keeps the others');
  assert.deepEqual(state(tags.root).items.list.map((i) => i.key), ['red', 'blue']);
  assert.equal(state(tags.root).open, false);
  tags.dispose();
});

smoke('DateInput: open while the calendar shows; the value reads as ISO', async () => {
  const date = mount(new DateInput({label: 'Due', value: new Date(Date.UTC(2026, 9, 15))}));
  assert.equal(state(date.root).open, false);
  assert.equal(state(date.root).value, '2026-10-15T00:00:00.000Z');
  fire(date.root.querySelector('.u2-date-toggle'), 'click');
  await flush();
  assert.equal(state(date.root).open, true);
  date.dispose();
});

smoke('MapInput: rows are items keyed by key; typing keeps the cell, a write from outside reloads', async () => {
  const map = mount(new MapInput({label: 'Headers', value: {a: '1', b: '2'}}));
  assert.deepEqual(state(map.root).items.list.map((i) => i.key), ['a', 'b']);
  const cell = map.root.querySelector('[data-u2-item="a"] .u2-map-value');
  cell.value = '10';
  fire(cell, 'input');
  await flush();
  assert.deepEqual(map.value, {a: '10', b: '2'});
  assert.equal(map.root.querySelector('[data-u2-item="a"] .u2-map-value'), cell, 'an edit never rebuilds the rows');
  map.value = {c: '3'};
  assert.deepEqual(state(map.root).items.list.map((i) => i.key), ['c']);
  map.dispose();
});

smoke('RangeSlider: a label names it; the value is [lo, hi]; change passes both ends and the control', async () => {
  const seen = [];
  const changed = [];
  const slider = mount(new RangeSlider({label: 'Price', lo: 20, hi: 60,
    onChanged: (lo, hi, source) => changed.push([lo, hi, source === slider])}));
  slider.on('change', (value, source) => seen.push([value, source === slider]));
  assert.equal(state(slider.root).label, 'Price');
  assert.deepEqual(state(slider.root).value, [20, 60]);
  slider.setValues(30, 50);
  assert.deepEqual(seen, [[[30, 50], true]]);
  assert.deepEqual(changed, [[30, 50, true]]);
  slider.dispose();
});

smoke('MessageInput: the placeholder names it; the state carries the text', async () => {
  const message = mount(new MessageInput({placeholder: 'Write a message'}));
  message.value = 'hello';
  assert.equal(state(message.root).label, 'Write a message');
  assert.equal(state(message.root).text, 'hello');
  assert.equal(state(message.root).open, false);
  message.dispose();
});
