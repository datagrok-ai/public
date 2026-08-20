/* TypeAhead: sync filtering with a custom row renderer, keyboard picking, selection/text
   coupling, and the async loading/error/retry states driven by hand-resolved promises.
   `userInput` (src/dg/user-input.ts) is not covered here — it imports datagrok-api, which does
   not load headless; it is exercised in the platform gallery instead. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope, Combobox} from '../src/index.js';
import {TypeAhead} from '../src/components/typeahead.js';

const USERS = [
  {name: 'Ada Almeida', login: 'adaalmeida'},
  {name: 'Bruno Bauer', login: 'brunobauer'},
  {name: 'Chen Costa', login: 'chencosta'},
];

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

function texts(root, selector) {
  return root.querySelectorAll(selector).map((el) => el.textContent);
}

function userRow(user) {
  const row = document.createElement('div');
  row.className = 'test-user';
  row.textContent = `${user.name} <${user.login}>`;
  return row;
}

function users(options = {}) {
  return mount(new TypeAhead({
    source: USERS,
    itemText: (u) => u.name,
    render: userRow,
    placeholder: 'Assignee',
    ...options,
  }));
}

function type(input, value) {
  input.value = value;
  fire(input, 'input');
}

smoke('sync: filters by itemText and renders every option through the callback', async () => {
  const typeahead = users();
  const input = typeahead.root.querySelector('input');
  input.focus();
  fire(input, 'keydown', {key: 'ArrowDown'});
  await flush();

  const popup = document.body.querySelector('.u2-typeahead-popup');
  assert.equal(typeahead.isOpen.value, true);
  assert.equal(typeahead.machineState.value, 'open');
  assert.equal(input.getAttribute('aria-expanded'), 'true');
  assert.equal(popup.getAttribute('role'), 'listbox');
  assert.equal(popup.querySelectorAll('.u2-typeahead-option').length, 3);
  assert.deepEqual(texts(popup, '.u2-typeahead-option .test-user'),
    ['Ada Almeida <adaalmeida>', 'Bruno Bauer <brunobauer>', 'Chen Costa <chencosta>']);
  assert.equal(popup.querySelector('.u2-typeahead-option').getAttribute('role'), 'option');

  type(input, 'bau');
  await flush();
  assert.deepEqual(texts(popup, '.u2-typeahead-option .test-user'), ['Bruno Bauer <brunobauer>']);

  type(input, 'zzz');
  await flush();
  assert.equal(popup.querySelectorAll('.u2-typeahead-option').length, 0);
  assert.equal(popup.querySelector('.u2-typeahead-empty').textContent, 'No matches');
  typeahead.dispose();
});

smoke('keyboard: ArrowDown/Enter picks the active option and fills the input', async () => {
  const typeahead = users();
  const input = typeahead.root.querySelector('input');
  input.focus();
  fire(input, 'keydown', {key: 'ArrowDown'});
  await flush();

  const popup = document.body.querySelector('.u2-typeahead-popup');
  fire(input, 'keydown', {key: 'ArrowDown'});
  fire(input, 'keydown', {key: 'ArrowDown'});
  const active = popup.querySelectorAll('.u2-typeahead-option')[1];
  assert.equal(active.getAttribute('aria-selected'), 'true');
  assert.equal(active.classList.contains('u2-typeahead-option-active'), true);
  assert.equal(input.getAttribute('aria-activedescendant'), active.id);

  fire(input, 'keydown', {key: 'ArrowUp'});
  fire(input, 'keydown', {key: 'Enter'});
  assert.equal(typeahead.selected.value, USERS[0]);
  assert.equal(typeahead.text.value, 'Ada Almeida');
  assert.equal(input.value, 'Ada Almeida');
  assert.equal(typeahead.isOpen.value, false);
  assert.equal(popup.isConnected, false);
  typeahead.dispose();
});

smoke('a typed query highlights the first match, so plain Enter accepts it', async () => {
  const typeahead = users();
  const input = typeahead.root.querySelector('input');
  input.focus();
  type(input, 'bau');
  await flush();

  const popup = document.body.querySelector('.u2-typeahead-popup');
  const first = popup.querySelector('.u2-typeahead-option');
  assert.equal(first.classList.contains('u2-typeahead-option-active'), true);
  assert.equal(input.getAttribute('aria-activedescendant'), first.id);

  fire(input, 'keydown', {key: 'Enter'});
  assert.equal(typeahead.selected.value, USERS[1]);
  assert.equal(input.value, 'Bruno Bauer');

  fire(input, 'keydown', {key: 'ArrowDown'});
  await flush();
  type(input, '');
  await flush();
  assert.equal(document.body.querySelector('.u2-typeahead-option-active'), null,
    'an empty box offers everything and pre-selects nothing');
  typeahead.dispose();
});

smoke('pointerdown on a row picks it; typing afterwards drops the stale selection', async () => {
  const typeahead = users();
  const input = typeahead.root.querySelector('input');
  input.focus();
  fire(input, 'keydown', {key: 'ArrowDown'});
  await flush();

  const popup = document.body.querySelector('.u2-typeahead-popup');
  fire(popup.querySelectorAll('.u2-typeahead-option .test-user')[2], 'pointerdown');
  assert.equal(typeahead.selected.value, USERS[2]);
  assert.equal(input.value, 'Chen Costa');

  type(input, 'Chen Cost');
  await flush();
  assert.equal(typeahead.selected.value, null, 'an edited text never keeps the old object');
  assert.equal(typeahead.text.value, 'Chen Cost', 'clearing the selection leaves the typed text');
  typeahead.dispose();
});

smoke('programmatic selection writes the text; Esc closes, second Esc clears', async () => {
  const typeahead = users();
  const input = typeahead.root.querySelector('input');
  typeahead.selected.value = USERS[1];
  assert.equal(input.value, 'Bruno Bauer');
  assert.equal(typeahead.text.value, 'Bruno Bauer');

  input.focus();
  fire(input, 'keydown', {key: 'ArrowDown'});
  await flush();
  assert.equal(typeahead.isOpen.value, true);

  fire(input, 'keydown', {key: 'Escape'});
  assert.equal(typeahead.isOpen.value, false);
  assert.equal(typeahead.machineState.value, 'focused');
  assert.equal(typeahead.selected.value, USERS[1], 'the first Esc only closes');

  fire(input, 'keydown', {key: 'Escape'});
  assert.equal(typeahead.selected.value, null);
  assert.equal(typeahead.text.value, '');
  assert.equal(input.value, '');
  typeahead.dispose();
});

smoke('Tab closes and keeps the pick', async () => {
  const typeahead = users();
  const input = typeahead.root.querySelector('input');
  input.focus();
  fire(input, 'keydown', {key: 'ArrowDown'});
  await flush();
  fire(input, 'keydown', {key: 'ArrowDown'});
  fire(input, 'keydown', {key: 'Enter'});
  fire(input, 'keydown', {key: 'ArrowDown'});
  await flush();

  fire(input, 'keydown', {key: 'Tab'});
  assert.equal(typeahead.isOpen.value, false);
  assert.equal(typeahead.selected.value, USERS[0]);
  assert.equal(input.value, 'Ada Almeida');
  typeahead.dispose();
});

smoke('minChars: the popup hints until the query is long enough', async () => {
  let calls = 0;
  const typeahead = mount(new TypeAhead({
    source: (query) => {
      calls++;
      return Promise.resolve(USERS.filter((u) => u.name.toLowerCase().includes(query)));
    },
    itemText: (u) => u.name,
    minChars: 2,
    debounceMs: 0,
  }));
  const input = typeahead.root.querySelector('input');
  input.focus();
  type(input, 'a');
  await flush();

  const popup = document.body.querySelector('.u2-typeahead-popup');
  assert.equal(calls, 0, 'a short query never reaches the source');
  assert.equal(popup.querySelector('.u2-typeahead-hint').textContent, 'Type 2 or more characters');

  type(input, 'al');
  await flush();
  assert.equal(calls, 1);
  assert.deepEqual(texts(popup, '.u2-typeahead-option'), ['Ada Almeida']);

  type(input, 'a');
  await flush();
  assert.equal(popup.querySelectorAll('.u2-typeahead-option').length, 0, 'stale options are dropped');
  assert.equal(popup.querySelector('.u2-typeahead-hint').isConnected, true);
  typeahead.dispose();
});

smoke('async: debounced loading, error row, Retry, then a pick', async () => {
  const calls = [];
  const typeahead = mount(new TypeAhead({
    source: (query, abort) => {
      const call = {query, abort};
      const promise = new Promise((resolve, reject) => Object.assign(call, {resolve, reject}));
      calls.push(call);
      return promise;
    },
    itemText: (u) => u.name,
    render: userRow,
    debounceMs: 0,
  }));
  const input = typeahead.root.querySelector('input');
  input.focus();
  type(input, 'b');
  await flush();

  const popup = document.body.querySelector('.u2-typeahead-popup');
  assert.equal(calls.length, 1);
  assert.equal(popup.querySelector('.u2-typeahead-loading').textContent, 'Loading…');

  calls[0].reject(new Error('Directory unavailable'));
  await flush();
  assert.equal(popup.querySelector('.u2-typeahead-error').textContent, 'Directory unavailableRetry');

  fire(popup.querySelector('.u2-typeahead-retry'), 'pointerdown');
  await flush();
  assert.equal(calls.length, 2);
  assert.equal(calls[1].query, 'b');

  calls[1].resolve([USERS[1]]);
  await flush();
  assert.deepEqual(texts(popup, '.u2-typeahead-option .test-user'), ['Bruno Bauer <brunobauer>']);

  fire(input, 'keydown', {key: 'ArrowDown'});
  fire(input, 'keydown', {key: 'Enter'});
  assert.equal(typeahead.selected.value, USERS[1]);
  assert.equal(input.value, 'Bruno Bauer');
  typeahead.dispose();
});

smoke('dispose: closes the popup and releases the option scopes', async () => {
  const before = Scope.liveCount;
  const typeahead = users();
  const input = typeahead.root.querySelector('input');
  input.focus();
  fire(input, 'keydown', {key: 'ArrowDown'});
  await flush();

  const popup = document.body.querySelector('.u2-typeahead-popup');
  assert.ok(Scope.liveCount > before);
  typeahead.dispose();
  assert.equal(popup.isConnected, false);
  assert.equal(Scope.liveCount, before);

  fire(input, 'keydown', {key: 'ArrowDown'});
  await flush();
  assert.equal(document.body.querySelector('.u2-typeahead-popup'), null, 'listeners died with the scope');
});

/* The `EditableChoiceInput` exemption (plan WO-14.3): the platform's select-or-type hybrid keeps
   unmatched text AS the value (`editable_choice_input.dart:42-48`). Combobox must do the same,
   or the exemption would be wrong and it would need a `freeText` option. */
smoke('Combobox: typed text that matches nothing IS the value', async () => {
  const combobox = mount(new Combobox({items: ['acid', 'base']}));
  const input = combobox.root.querySelector('input');

  input.value = 'ester';
  fire(input, 'input');
  await flush();
  assert.equal(combobox.value.value, 'ester', 'free text is the value, no item needed');
  assert.equal(document.body.querySelector('.u2-combobox-empty').textContent, 'No matches');

  fire(input, 'keydown', {key: 'Enter'});
  assert.equal(combobox.value.value, 'ester', 'Enter with nothing highlighted leaves it alone');

  input.value = 'aci';
  fire(input, 'input');
  await flush();
  fire(input, 'keydown', {key: 'ArrowDown'});
  fire(input, 'keydown', {key: 'Enter'});
  assert.equal(combobox.value.value, 'acid', 'and picking an item still commits the item');
  combobox.dispose();
});
