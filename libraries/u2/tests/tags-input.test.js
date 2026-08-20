/* TagsInput: chips over T, suggestion popup (sync and async), free creation, the keyboard
   contract, and the per-render scopes that keep rendered rows and chips leak-free. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/index.js';
import {TagsInput} from '../src/components/tags-input.js';

const KEYWORDS = ['acid', 'base', 'ester', 'amine'];
const USERS = [
  {id: 1, name: 'Alice'},
  {id: 2, name: 'Bob'},
  {id: 3, name: 'Carol'},
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

function box(input) {
  return input.root.querySelector('.u2-tags-input');
}

function chips(input) {
  return input.root.querySelectorAll('.u2-tags-chip .u2-tag-text').map((el) => el.textContent);
}

function popup() {
  return document.body.querySelector('.u2-tags-popup');
}

function rows() {
  return popup().querySelectorAll('.u2-tags-option');
}

function type(input, text) {
  box(input).value = text;
  fire(box(input), 'input');
}

function key(input, name, init = {}) {
  return fire(box(input), 'keydown', {key: name, ...init});
}

/** Debounce timer plus the fetch continuation behind it. */
async function settle() {
  await flush();
  await flush();
}

smoke('typing opens the popup and filters the suggestions', async () => {
  const input = mount(new TagsInput({label: 'Tags', items: KEYWORDS}));
  assert.equal(input.root.dataset.u2, 'tags-input');
  assert.equal(input.isOpen.value, false);

  type(input, 'a');
  await flush();
  assert.equal(input.isOpen.value, true);
  assert.equal(box(input).getAttribute('aria-expanded'), 'true');
  assert.deepEqual(rows().map((r) => r.textContent), ['acid', 'base', 'amine']);
  assert.equal(popup().querySelector('.u2-tags-option-active'), null,
    'rows are offered but none is highlighted — Enter belongs to the typed text (allowNew)');

  type(input, 'am');
  await flush();
  assert.deepEqual(rows().map((r) => r.textContent), ['amine']);

  type(input, 'zzz');
  await flush();
  assert.equal(rows().length, 0);
  assert.equal(popup().querySelector('.u2-tags-empty').textContent, 'No matches');
  input.dispose();
});

smoke('picking a suggestion adds a chip, clears the box and drops the item from the list',
  async () => {
    const changes = [];
    const input = mount(new TagsInput({label: 'Tags', items: KEYWORDS,
      onChanged: (v) => changes.push(v)}));
    type(input, 'a');
    await flush();

    fire(rows()[0], 'pointerdown');
    await flush();
    assert.deepEqual(input.value.value, ['acid']);
    assert.deepEqual(chips(input), ['acid']);
    assert.equal(changes.length, 1);
    assert.equal(box(input).value, '', 'the box is cleared for the next tag');
    assert.equal(input.isOpen.value, false, 'and the popup closes, as the platform does');

    type(input, 'a');
    await flush();
    assert.deepEqual(rows().map((r) => r.textContent), ['base', 'amine'],
      'an item already picked is not suggested again');

    fire(input.root.querySelectorAll('.u2-tag-remove')[0], 'click');
    await flush();
    assert.deepEqual(input.value.value, []);
    assert.deepEqual(chips(input), []);
    input.dispose();
  });

smoke('keyboard: arrows and Home/End move, Enter accepts, Esc closes, Backspace pops the last chip',
  async () => {
    const input = mount(new TagsInput({label: 'Tags', items: KEYWORDS}));
    type(input, '');
    await flush();
    assert.equal(rows().length, 4);

    key(input, 'ArrowDown');
    key(input, 'ArrowDown');
    assert.equal(rows()[1].classList.contains('u2-tags-option-active'), true);
    assert.equal(box(input).getAttribute('aria-activedescendant'), rows()[1].id);

    key(input, 'ArrowUp');
    assert.equal(rows()[0].classList.contains('u2-tags-option-active'), true);
    key(input, 'End');
    assert.equal(rows()[3].classList.contains('u2-tags-option-active'), true);
    key(input, 'Home');
    assert.equal(rows()[0].classList.contains('u2-tags-option-active'), true);

    key(input, 'Enter');
    await flush();
    assert.deepEqual(input.value.value, ['acid']);
    assert.equal(input.isOpen.value, false);

    key(input, 'ArrowDown');
    await flush();
    assert.equal(input.isOpen.value, true, 'ArrowDown reopens the popup');
    key(input, 'Escape');
    assert.equal(input.isOpen.value, false);
    assert.equal(popup(), null, 'the overlay is gone');

    key(input, 'Backspace');
    assert.deepEqual(input.value.value, [], 'Backspace on an empty box drops the last chip');

    input.value.value = ['acid'];
    type(input, 'ba');
    key(input, 'Backspace');
    assert.deepEqual(input.value.value, ['acid'], 'but never while something is typed');
    input.dispose();
  });

smoke('an Enter that adds a chip is consumed; one that does nothing reaches the host', async () => {
  const input = mount(new TagsInput({label: 'Tags', items: KEYWORDS}));
  const seen = [];
  const host = (e) => {
    if (e.key === 'Enter')
      seen.push(e.key);
  };
  document.body.addEventListener('keydown', host);

  type(input, '');
  await flush();
  key(input, 'ArrowDown');
  key(input, 'Enter');
  await flush();
  assert.deepEqual(input.value.value, ['acid']);
  assert.deepEqual(seen, [], 'the accepting Enter never bubbles to the dialog around it');

  key(input, 'Enter');
  assert.deepEqual(seen, ['Enter'], 'an Enter with nothing to add still OKs the host');
  document.body.removeEventListener('keydown', host);
  input.dispose();
});

smoke('Ctrl+Enter belongs to the hosting dialog, not to the box', async () => {
  const input = mount(new TagsInput({label: 'Tags', items: KEYWORDS, allowNew: true}));
  const seen = [];
  const host = (e) => {
    if (e.key === 'Enter')
      seen.push(e.ctrlKey ? 'ctrl' : 'plain');
  };
  document.body.addEventListener('keydown', host);

  type(input, 'ketone');
  await flush();
  key(input, 'Enter', {ctrlKey: true});
  assert.deepEqual(input.value.value, [], 'Ctrl+Enter never creates a tag');
  assert.deepEqual(seen, ['ctrl'], 'and reaches the dialog OK it is meant for');

  key(input, 'Enter', {metaKey: true});
  assert.deepEqual(input.value.value, [], 'the same on macOS');

  key(input, 'Enter');
  assert.deepEqual(input.value.value, ['ketone'], 'a plain Enter still creates');
  document.body.removeEventListener('keydown', host);
  input.dispose();
});

smoke('allowNew: Enter turns the typed text into an item, deduplicating by text', async () => {
  const plain = mount(new TagsInput({label: 'Tags', items: KEYWORDS}));
  type(plain, 'ketone');
  await flush();
  key(plain, 'Enter');
  assert.deepEqual(plain.value.value, [], 'creating is off by default');
  plain.dispose();

  const input = mount(new TagsInput({label: 'Tags', items: KEYWORDS, allowNew: true}));
  type(input, 'ketone');
  await flush();
  key(input, 'Enter');
  await flush();
  assert.deepEqual(input.value.value, ['ketone']);
  assert.deepEqual(chips(input), ['ketone']);
  assert.equal(box(input).value, '');

  type(input, 'ketone');
  await flush();
  key(input, 'Enter');
  assert.deepEqual(input.value.value, ['ketone'], 'the same text never doubles up');
  input.dispose();

  const typed = mount(new TagsInput({label: 'Users', items: USERS, allowNew: true,
    itemText: (u) => u.name, createItem: (text) => ({id: 0, name: text})}));
  type(typed, 'Dave');
  await flush();
  key(typed, 'Enter');
  assert.deepEqual(typed.value.value, [{id: 0, name: 'Dave'}], 'createItem builds the object');
  typed.dispose();
});

smoke('renderer: chips use markup, rows use listItem, and caption is the fallback text',
  async () => {
    const renderer = {
      caption: (u) => u.name,
      markup: (u) => {
        const el = document.createElement('span');
        el.className = 'demo-chip';
        el.textContent = `@${u.name}`;
        return el;
      },
      listItem: (u) => {
        const el = document.createElement('div');
        el.className = 'demo-row';
        el.textContent = `${u.name} (${u.id})`;
        return el;
      },
    };
    const input = mount(new TagsInput({label: 'Users', items: USERS, renderer, value: [USERS[0]]}));
    assert.deepEqual(input.root.querySelectorAll('.demo-chip').map((el) => el.textContent), ['@Alice']);

    type(input, 'o');
    await flush();
    assert.deepEqual(popup().querySelectorAll('.demo-row').map((el) => el.textContent),
      ['Bob (2)', 'Carol (3)']);

    fire(rows()[0], 'pointerdown');
    await flush();
    assert.deepEqual(input.value.value.map((u) => u.name), ['Alice', 'Bob']);
    assert.deepEqual(input.root.querySelectorAll('.demo-chip').map((el) => el.textContent),
      ['@Alice', '@Bob']);
    input.dispose();
  });

smoke('deleteEnabled hides the ✕ and protects the chip from Backspace', async () => {
  const input = mount(new TagsInput({label: 'Tags', items: KEYWORDS, value: ['acid', 'base'],
    deleteEnabled: (item) => item !== 'acid'}));
  assert.deepEqual(input.root.querySelectorAll('.u2-tag-remove').map((b) => b.getAttribute('aria-label')),
    ['Remove base']);

  key(input, 'Backspace');
  assert.deepEqual(input.value.value, ['acid'], 'the deletable one goes');
  key(input, 'Backspace');
  assert.deepEqual(input.value.value, ['acid'], 'the protected one stays');
  input.dispose();
});

smoke('async source: minChars hint, debounced query, error row and retry', async () => {
  const calls = [];
  let fail = false;
  const input = mount(new TagsInput({label: 'Users', minChars: 2, debounceMs: 1,
    itemText: (u) => u.name,
    source: async (query) => {
      calls.push(query);
      if (fail)
        throw new Error('offline');
      return USERS.filter((u) => u.name.toLowerCase().includes(query.toLowerCase()));
    }}));

  type(input, 'a');
  await flush();
  assert.deepEqual(calls, [], 'below minChars nothing is queried');
  assert.equal(popup().querySelector('.u2-tags-hint').textContent, 'Type 2 or more characters');

  type(input, 'al');
  await settle();
  assert.deepEqual(calls, ['al']);
  assert.deepEqual(rows().map((r) => r.textContent), ['Alice']);

  fail = true;
  type(input, 'bo');
  await settle();
  assert.equal(popup().querySelector('.u2-tags-error').textContent, 'offlineRetry');
  fail = false;
  fire(popup().querySelector('.u2-tags-retry'), 'pointerdown');
  await settle();
  assert.deepEqual(rows().map((r) => r.textContent), ['Bob']);
  input.dispose();
});

smoke('dispose closes the popup and releases the row and chip scopes', async () => {
  const before = Scope.liveCount;
  const input = mount(new TagsInput({label: 'Tags', items: KEYWORDS, value: ['acid']}));
  type(input, 'b');
  await flush();
  const open = popup();
  assert.ok(Scope.liveCount > before);

  input.dispose();
  assert.equal(open.isConnected, false);
  assert.equal(Scope.liveCount, before);

  type(input, 'e');
  await flush();
  assert.equal(popup(), null, 'listeners died with the scope');
});
