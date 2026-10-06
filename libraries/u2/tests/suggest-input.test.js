/* SuggestInput (plan-w2.md FP-W2-6): a free-text input whose value IS the typed text, with a
   SuggestionList popup over an async fetch — the `u2-typeahead` row vocabulary reused. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {SuggestInput} from '../src/components/inputs/suggest-input.js';

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

function box(input) {
  return input.root.querySelector('input');
}

function popup() {
  return document.body.querySelector('.u2-typeahead-popup');
}

function rows() {
  return popup()?.querySelectorAll('.u2-typeahead-option') ?? [];
}

function type(input, text) {
  box(input).value = text;
  fire(box(input), 'input');
}

function key(input, name, init = {}) {
  return fire(box(input), 'keydown', {key: name, ...init});
}

async function until(cond, what, timeoutMs = 1000) {
  const start = Date.now();
  while (!cond()) {
    if (Date.now() - start > timeoutMs)
      assert.fail(`timed out waiting for ${what}`);
    await flush();
  }
}

smoke('every keystroke writes the value signal, and typing opens the popup with the fetched rows', async () => {
  const calls = [];
  const input = mount(new SuggestInput({label: 'Drug', debounceMs: 1,
    source: async (q) => {
      calls.push(q);
      return [`${q}-1`, `${q}-2`];
    }}));
  assert.equal(input.root.dataset.u2, 'suggest-input');
  assert.ok(input.root.classList.contains('u2-suggest-input'));
  assert.equal(box(input).type, 'text');

  type(input, 'as');
  assert.equal(input.value.value, 'as', 'the value is the typed text, on every keystroke');
  await until(() => rows().length === 2, 'the popup rows');
  assert.deepEqual(calls, ['as'], 'queried with the typed text');
  assert.deepEqual(rows().map((r) => r.textContent), ['as-1', 'as-2']);
  input.dispose();
});

smoke('a pointer pick commits the row value and dismisses the popup', async () => {
  const changes = [];
  const input = mount(new SuggestInput({label: 'Drug', debounceMs: 1,
    onChanged: (v) => changes.push(v),
    source: async (q) => [`${q}-1`]}));
  type(input, 'b');
  await until(() => rows().length === 1, 'the popup rows');
  fire(rows()[0], 'pointerdown');
  await flush();
  assert.equal(input.value.value, 'b-1');
  assert.equal(box(input).value, 'b-1', 'the box shows the pick');
  assert.equal(popup(), null, 'the popup is dismissed');
  assert.equal(changes[changes.length - 1], 'b-1');
  input.dispose();
});

smoke('keyboard: arrows move the active row, Enter commits it and is contained, a bare Enter passes', async () => {
  const input = mount(new SuggestInput({label: 'Drug', debounceMs: 1,
    source: async () => ['aa', 'ab']}));
  const seen = [];
  const host = (e) => {
    if (e.key === 'Enter')
      seen.push(e.key);
  };
  document.body.addEventListener('keydown', host);

  type(input, 'a');
  await until(() => rows().length === 2, 'the popup rows');
  key(input, 'ArrowDown');
  key(input, 'ArrowDown');
  assert.equal(rows()[1].classList.contains('u2-typeahead-option-active'), true);
  key(input, 'ArrowUp');
  assert.equal(rows()[0].classList.contains('u2-typeahead-option-active'), true);

  const passed = key(input, 'Enter');
  await flush();
  assert.equal(input.value.value, 'aa', 'Enter commits the active row');
  assert.equal(popup(), null);
  assert.equal(passed, false, 'preventDefault while the popup is open — dialog-OK protection');
  assert.deepEqual(seen, [], 'the committing Enter never reaches the host');

  key(input, 'Enter');
  assert.deepEqual(seen, ['Enter'], 'with the popup closed Enter passes through');
  document.body.removeEventListener('keydown', host);
  input.dispose();
});

smoke('Escape dismisses the popup; closed, it reverts to the focus-time value; blur and Tab dismiss', async () => {
  const input = mount(new SuggestInput({label: 'Drug', debounceMs: 1,
    source: async () => ['x1']}));
  fire(box(input), 'focus');
  type(input, 'x');
  await until(() => rows().length === 1, 'the popup rows');
  key(input, 'Escape');
  assert.equal(popup(), null, 'Esc closes the popup');
  assert.equal(input.value.value, 'x', 'and keeps the text');
  key(input, 'Escape');
  assert.equal(input.value.value, '', 'Esc with the popup closed reverts to the value at focus');
  key(input, 'Escape');
  assert.equal(input.value.value, '', 'and is a no-op once nothing changed');

  type(input, 'x');
  await until(() => rows().length === 1, 'the popup again');
  key(input, 'Tab');
  assert.equal(popup(), null, 'Tab dismisses');

  type(input, 'x');
  await until(() => rows().length === 1, 'and once more');
  fire(box(input), 'blur');
  assert.equal(popup(), null, 'blur dismisses');
  input.dispose();
});

smoke('a committed pick is the revert point: a type-over then Escape restores it', async () => {
  const input = mount(new SuggestInput({label: 'Drug', debounceMs: 1,
    source: async (q) => [`${q}-1`]}));
  fire(box(input), 'focus');
  type(input, 'bo');
  await until(() => rows().length === 1, 'the popup rows');
  fire(rows()[0], 'pointerdown');
  await flush();
  assert.equal(input.value.value, 'bo-1');

  type(input, 'cl');
  assert.equal(input.value.value, 'cl', 'the type-over still writes through on every keystroke');
  await until(() => rows().length === 1, 'the popup reopens');
  key(input, 'Escape');
  assert.equal(popup(), null, 'the first Esc only dismisses the popup');
  assert.equal(input.value.value, 'cl');
  key(input, 'Escape');
  assert.equal(input.value.value, 'bo-1', 'the second Esc reverts to the committed pick');
  assert.equal(box(input).value, 'bo-1', 'the box shows the revert');
  input.dispose();
});

smoke('an external write while focused moves the revert point — Escape keeps the external value', async () => {
  const input = mount(new SuggestInput({label: 'Drug', debounceMs: 1,
    source: async (q) => [`${q}-1`]}));
  fire(box(input), 'focus');
  type(input, 'ty');
  await until(() => rows().length === 1, 'the popup rows');
  input.value.value = 'ext';
  key(input, 'Escape');
  assert.equal(input.value.value, 'ext', 'the first Esc only dismisses the popup');
  key(input, 'Escape');
  assert.equal(input.value.value, 'ext', 'the external write is the new revert point, not the focus-time value');
  assert.equal(box(input).value, 'ext');
  input.dispose();
});

smoke('minChars: below the threshold the popup hints and the source is never queried', async () => {
  const calls = [];
  const input = mount(new SuggestInput({label: 'Drug', minChars: 2, debounceMs: 1,
    source: async (q) => {
      calls.push(q);
      return [q];
    }}));
  type(input, 'a');
  await flush();
  await flush();
  assert.deepEqual(calls, [], 'below minChars nothing is queried');
  assert.equal(popup().querySelector('.u2-typeahead-hint').textContent, 'Type 2 or more characters');

  type(input, 'ab');
  await until(() => rows().length === 1, 'queried at the threshold');
  assert.deepEqual(calls, ['ab']);
  input.dispose();
});

smoke('dispose closes the popup and releases every scope', async () => {
  const before = Scope.liveCount;
  const input = mount(new SuggestInput({label: 'Drug', debounceMs: 1,
    source: async () => ['a']}));
  type(input, 'a');
  await until(() => rows().length === 1, 'the popup rows');
  const open = popup();
  input.dispose();
  assert.equal(open.isConnected, false);
  assert.equal(Scope.liveCount, before);
});
