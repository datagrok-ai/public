/* focusMark / refocus: the focused element is named by the keys on the way down from the host — one
   per level of keyed nesting — and found again in the rebuilt DOM: itself where focusable, else its
   first enabled focusable descendant. A key no longer there hands focus to a neighbour; an unkeyed
   element or focus outside the host is left alone. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {resetDom} from './dom-shim.js';
import {focusMark, focusPath, refocus} from '../src/core/focus.js';

function el(tag, data = {}, ...children) {
  const node = document.createElement(tag);
  for (const [k, v] of Object.entries(data))
    node.dataset[k] = v;
  if (tag === 'input')
    node.type = 'checkbox';
  node.append(...children);
  return node;
}

/** The same shape a context panel rebuilds: an access grid by name, rows and boxes by key, a
 * bool input by name, a bare keyed button, an unkeyed input. */
function content() {
  return [
    el('div', {u2Name: 'access'}, el('table', {}, el('tbody', {},
      el('tr', {u2Key: 'Sales'}, el('td', {}, el('input', {u2Key: 'view'}), el('input', {u2Key: 'edit'}))),
      el('tr', {u2Key: 'Dev'}, el('td', {}, el('input', {u2Key: 'view'})))))),
    el('div', {u2Name: 'writable'}, el('span', {}, el('input'))),
    el('button', {u2Key: 'plain'}),
    el('input'),
  ];
}

function build() {
  const host = document.createElement('div');
  host.append(...content());
  document.body.append(host);
  return host;
}

const rebuild = (host) => host.replaceChildren(...content());

/** What the designer does around a structural patch: mark, rebuild, refocus. */
function kept(host, body) {
  const mark = focusMark(host);
  body();
  refocus(host, mark);
}
const inputs = (host) => host.querySelectorAll('input');

function ui(name, body) {
  test(name, async () => {
    try {
      await body();
    } finally {
      resetDom();
    }
  });
}

ui('focusPath names the focused element by its keyed ancestors; empty without any, null outside the host', () => {
  const host = build();
  inputs(host)[1].focus();
  assert.deepEqual(focusPath(host), ['access', 'Sales', 'edit']);
  host.querySelector('[data-u2-name="writable"] input').focus();
  assert.deepEqual(focusPath(host), ['writable'], 'an input root by its name');
  inputs(host)[4].focus();
  assert.deepEqual(focusPath(host), []);
  document.body.focus();
  assert.equal(focusPath(host), null);
});

ui('refocus finds the same keys in the rebuilt DOM, an input root resolving to its focusable', () => {
  const host = build();
  inputs(host)[1].focus();
  const before = document.activeElement;
  kept(host, () => rebuild(host));
  assert.equal(document.activeElement === before, false);
  assert.equal(document.activeElement === inputs(host)[1], true);
  host.querySelector('[data-u2-name="writable"] input').focus();
  kept(host, () => rebuild(host));
  assert.equal(document.activeElement === host.querySelector('[data-u2-name="writable"] input'), true);
  host.querySelector('button').focus();
  kept(host, () => rebuild(host));
  assert.equal(document.activeElement === host.querySelector('button'), true);
});

ui('a key is matched among the nearest keyed descendants only, never deeper down another branch', () => {
  const host = build();
  refocus(host, ['access', 'edit']);
  assert.equal(document.activeElement === inputs(host)[0], true,
    'edit is two levels down, under a row: not it, but the level\'s first focusable');
  refocus(host, ['access', 'Sales', 'edit']);
  assert.equal(document.activeElement === inputs(host)[1], true);
});

ui('an unkeyed element, a gone key or a disabled target: nothing is focused in its place', () => {
  const host = build();
  inputs(host)[4].focus();
  kept(host, () => rebuild(host));
  assert.equal(document.body.contains(document.activeElement), false, 'an empty path refocuses nothing');
  const plain = host.querySelector('button');
  plain.focus();
  refocus(host, null);
  assert.equal(document.activeElement === plain, true);
  refocus(host, ['nope']);
  assert.equal(document.activeElement === inputs(host)[0], true, 'a key gone without a neighbour: the host\'s first focusable');
  inputs(host)[2].focus();
  kept(host, () => {
    rebuild(host);
    inputs(host)[2].disabled = true;
  });
  assert.equal(document.activeElement === inputs(host)[2], false, 'a disabled box cannot take focus');
});

ui('a keyed element gone from the rebuilt DOM hands focus to the neighbour after it, else before, else the level', () => {
  const host = build();
  const rows = () => host.querySelector('[data-u2-name="access"] tbody');
  rows().append(el('tr', {u2Key: 'Ops'}, el('td', {}, el('input', {u2Key: 'view'}), el('button', {u2Key: 'remove'}))));
  rows().querySelector('[data-u2-key="Dev"] input').focus();
  kept(host, () => rows().querySelector('[data-u2-key="Dev"]').remove());
  assert.equal(document.activeElement === rows().querySelector('[data-u2-key="Ops"] input'), true,
    'the row after, on the same key');
  rows().querySelector('[data-u2-key="Ops"] button').focus();
  kept(host, () => rows().querySelector('[data-u2-key="Ops"]').remove());
  assert.equal(document.activeElement === rows().querySelector('[data-u2-key="Sales"] input'), true,
    'the last row gone: the one before, its first focusable since it has no remove');
  rows().querySelector('[data-u2-key="Sales"] input').focus();
  kept(host, () => rows().querySelector('[data-u2-key="Sales"]').remove());
  assert.equal(document.activeElement === host.querySelector('[data-u2-name="writable"] input'), true,
    'no neighbour left: the level\'s first focusable');
});
