/* MessageInput core: the markup serialization walk, token → chip restore, the @ accelerator over a
   stubbed caret, the keyboard matrix and the draft. The shim has neither a selection nor a Range,
   so this file installs its own `window.getSelection` (settable position, no ranges to hand back —
   real caret behavior is e2e's job) and a `localStorage` map; dom-shim.js is never edited. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {MessageInput} from '../src/components/inputs/message-input.js';

let selection = null;
window.getSelection = () => selection;

const store = new Map();
window.localStorage = {
  getItem: (key) => store.has(key) ? store.get(key) : null,
  setItem: (key, value) => store.set(key, String(value)),
  removeItem: (key) => store.delete(key),
};

/** Puts the stubbed caret at `offset` inside `node`; `null` means "the host offers no selection". */
function caret(node, offset) {
  selection = node === null ? null :
    {anchorNode: node, anchorOffset: offset, isCollapsed: true, rangeCount: 1};
}

function smoke(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      caret(null);
      store.clear();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

const PEOPLE = [
  {id: 'u-1', name: 'Ada Almeida'},
  {id: 'u-2', name: 'Alan Bauer'},
  {id: 'u-3', name: 'Chen Costa'},
];

const TOKEN = '<span>#\\{x\\.([^."]+)\\."(.*?)"\\}<\\/span>';
const markupOf = (p) => `<span>#{x.${p.id}."${p.name}"}</span>`;

function people(options = {}) {
  return {
    trigger: '@',
    minChars: 0,
    debounceMs: 0,
    fetch: (query) => Promise.resolve(query === '' ? PEOPLE :
      PEOPLE.filter((p) => p.name.toLowerCase().includes(query.toLowerCase()))),
    caption: (p) => p.name,
    toMarkup: markupOf,
    tokenPattern: TOKEN,
    renderToken: (m) => {
      const el = document.createElement('span');
      el.className = 'test-token';
      el.textContent = m[2];
      return el;
    },
    ...options,
  };
}

function compose(options = {}) {
  const control = new MessageInput({mentions: [people()], ...options});
  document.body.append(control.root);
  return control;
}

const editorOf = (control) => control.root.querySelector('.u2-msg-editor');
const popup = () => document.body.querySelector('.u2-msg-popup');
const chips = (control) => editorOf(control).querySelectorAll('.u2-msg-chip');

/** Types at the end of the editor the way the browser does — text into the caret's node, then
 * `input` — creating the text node on the first call. */
function type(control, text) {
  const editor = editorOf(control);
  let node = editor.childNodes[editor.childNodes.length - 1];
  if (!node || node.nodeType !== 3) {
    node = document.createTextNode('');
    editor.append(node);
  }
  node.textContent += text;
  caret(node, node.textContent.length);
  fire(editor, 'input');
}

/** Types into the open query span, as the browser would. */
function typeQuery(control, text) {
  const span = editorOf(control).querySelector('.u2-msg-query');
  span.textContent += text;
  caret(span.firstChild, span.textContent.length);
  fire(editorOf(control), 'input');
}

smoke('serializes text, chips as their token, and br/div boundaries as newlines', async () => {
  const control = compose();
  const editor = editorOf(control);
  control.value = `hello ${markupOf(PEOPLE[0])}`;
  editor.append(document.createElement('br'));
  const block = document.createElement('div');
  block.textContent = 'second line';
  editor.append(block);

  assert.equal(chips(control).length, 1);
  assert.equal(chips(control)[0].dataset.token, markupOf(PEOPLE[0]));
  assert.equal(chips(control)[0].getAttribute('contenteditable'), 'false');
  assert.equal(control.value, `hello ${markupOf(PEOPLE[0])}\nsecond line`,
    'the block boundary emits one newline and the trailing one is trimmed');
  control.dispose();
});

smoke('parses tokens into chips, round-trips byte-exact, leaves unknown text literal', async () => {
  const control = compose();
  const markup = `Hi ${markupOf(PEOPLE[0])} and ${markupOf(PEOPLE[2])}, please review.`;
  control.value = markup;

  assert.equal(chips(control).length, 2);
  assert.deepEqual([...chips(control)].map((c) => c.textContent), ['Ada Almeida', 'Chen Costa']);
  assert.equal(editorOf(control).querySelectorAll('.test-token').length, 2,
    'the restore path renders through renderToken');
  assert.equal(control.value, markup);

  control.value = 'plain <b>text</b> with #{no token} in it';
  assert.equal(chips(control).length, 0);
  assert.equal(control.value, 'plain <b>text</b> with #{no token} in it');
  control.dispose();
});

smoke('a multi-line value restores through <br> and serializes back to \\n', async () => {
  const control = compose();
  control.value = 'first\nsecond';
  assert.equal(editorOf(control).querySelectorAll('br').length, 1);
  assert.equal(control.value, 'first\nsecond');
  control.dispose();
});

smoke('insertMention appends the chip plus a space, and value carries toMarkup', async () => {
  const control = compose();
  type(control, 'ping ');
  control.insertMention(PEOPLE[1]);

  assert.equal(chips(control).length, 1);
  assert.equal(control.value, `ping ${markupOf(PEOPLE[1])} `);
  assert.equal(control.empty.value, false);
  control.dispose();
});

smoke('Backspace next to a chip removes it whole; elsewhere it is left to the browser', async () => {
  const control = compose();
  const editor = editorOf(control);
  control.value = `ping ${markupOf(PEOPLE[0])} x`;
  const chip = chips(control)[0];
  const after = chip.nextSibling;

  caret(after, 1);
  const kept = fire(editor, 'keydown', {key: 'Backspace'});
  assert.equal(chips(control).length, 1, 'mid-text Backspace is the browser\'s own');
  assert.equal(kept, true, 'and is not prevented');

  caret(after, 0);
  const eaten = fire(editor, 'keydown', {key: 'Backspace'});
  assert.equal(chips(control).length, 0);
  assert.equal(eaten, false, 'the chip deletion is prevented');
  assert.equal(control.value, 'ping  x');

  control.value = `${markupOf(PEOPLE[2])}tail`;
  caret(editorOf(control), 0);
  fire(editor, 'keydown', {key: 'Delete'});
  assert.equal(chips(control).length, 0, 'Delete eats the chip in front of the caret');
  control.dispose();
});

smoke('Enter matrix: ctrlEnter sends, plain Enter breaks; enter mode is the mirror', async () => {
  const sent = [];
  const ctrl = compose({onSend: (m) => sent.push(m)});
  type(ctrl, 'draft');

  fire(editorOf(ctrl), 'keydown', {key: 'Enter'});
  assert.deepEqual(sent, [], 'plain Enter never sends in ctrlEnter mode');
  type(ctrl, 'more');
  assert.equal(ctrl.value, 'draft\nmore', 'it broke the line instead');

  fire(editorOf(ctrl), 'keydown', {key: 'Enter', ctrlKey: true});
  assert.deepEqual(sent, ['draft\nmore']);
  await flush();
  assert.equal(ctrl.value, '');
  ctrl.dispose();

  const quick = compose({sendOn: 'enter', onSend: (m) => sent.push(m)});
  type(quick, 'now');
  fire(editorOf(quick), 'keydown', {key: 'Enter', shiftKey: true});
  type(quick, 'later');
  assert.equal(quick.value, 'now\nlater', 'Shift+Enter still breaks the line');
  fire(editorOf(quick), 'keydown', {key: 'Enter'});
  assert.deepEqual(sent, ['draft\nmore', 'now\nlater']);
  await flush();
  assert.equal(quick.value, '');
  quick.dispose();
});

smoke('send clears the box and drops the draft; empty gates it', async () => {
  const sent = [];
  const control = compose({draftKey: 'chat-7', onSend: (m) => sent.push(m)});
  assert.equal(control.empty.value, true);

  control.send();
  assert.deepEqual(sent, [], 'an empty box never sends');

  type(control, 'hello ');
  control.insertMention(PEOPLE[0]);
  assert.equal(control.empty.value, false);
  assert.equal(store.get('u2-msg-draft-chat-7'), `hello ${markupOf(PEOPLE[0])} `);

  control.send();
  assert.deepEqual(sent, [`hello ${markupOf(PEOPLE[0])} `]);
  await flush();
  assert.equal(control.value, '');
  assert.equal(control.empty.value, true);
  assert.equal(store.has('u2-msg-draft-chat-7'), false, 'the draft key is gone');
  control.dispose();
});

smoke('a handler that rejects keeps the message and its draft where they are', async () => {
  const control = compose({draftKey: 'chat-8',
    onSend: () => Promise.reject(new Error('offline'))});
  type(control, 'not going anywhere');
  const errors = [];
  const real = console.error;
  console.error = (...args) => errors.push(args[0]);
  try {
    control.send();
    await flush();
  } finally {
    console.error = real;
  }
  assert.equal(control.value, 'not going anywhere', 'the text survived the failure');
  assert.equal(control.empty.value, false);
  assert.equal(store.get('u2-msg-draft-chat-8'), 'not going anywhere', 'and so did the draft');
  assert.equal(errors.length, 1, 'the failure is reported, not swallowed');
  control.dispose();
});

smoke('a stored draft is restored on construction, chips and all; options.value wins', async () => {
  store.set('u2-msg-draft-chat-9', `left off at ${markupOf(PEOPLE[1])}`);
  const restored = compose({draftKey: 'chat-9'});
  assert.equal(chips(restored).length, 1);
  assert.equal(restored.value, `left off at ${markupOf(PEOPLE[1])}`);
  restored.dispose();

  const explicit = compose({draftKey: 'chat-9', value: 'fresh'});
  assert.equal(explicit.value, 'fresh');
  explicit.dispose();
});

smoke('typing @ opens a caret-anchored popup; ArrowDown + Enter inserts the chip', async () => {
  const control = compose();
  const editor = editorOf(control);
  type(control, 'hey @');
  await flush();

  const span = editor.querySelector('.u2-msg-query');
  assert.ok(span, 'the trigger is wrapped into a query span');
  assert.equal(span.textContent, '@');
  assert.equal(control.isMentionOpen.value, true);
  assert.equal(popup().querySelectorAll('.u2-msg-option').length, 3);

  typeQuery(control, 'al');
  await flush();
  assert.deepEqual([...popup().querySelectorAll('.u2-msg-option')].map((o) => o.textContent),
    ['Ada Almeida', 'Alan Bauer']);

  fire(editor, 'keydown', {key: 'ArrowDown'});
  fire(editor, 'keydown', {key: 'Enter'});
  await flush();

  assert.equal(control.isMentionOpen.value, false);
  assert.equal(editor.querySelector('.u2-msg-query'), null, 'the query span is gone');
  assert.equal(chips(control).length, 1);
  assert.equal(control.value, `hey ${markupOf(PEOPLE[1])} `);
  assert.equal(popup(), null);
  control.dispose();
});

smoke('Esc unwraps the query span back into plain text', async () => {
  const control = compose();
  const editor = editorOf(control);
  type(control, '@ad');
  await flush();
  assert.equal(control.isMentionOpen.value, true);

  fire(editor, 'keydown', {key: 'Escape'});
  await flush();
  assert.equal(control.isMentionOpen.value, false);
  assert.equal(editor.querySelector('.u2-msg-query'), null);
  assert.equal(control.value, '@ad', 'the typed text survives the dismissal');
  assert.equal(chips(control).length, 0);
  control.dispose();
});

smoke('openMentions inserts the trigger and opens; a pick lands after it', async () => {
  const control = compose();
  control.openMentions();
  await flush();

  assert.equal(control.isMentionOpen.value, true);
  assert.equal(editorOf(control).querySelector('.u2-msg-query').textContent, '@');
  fire(popup().querySelectorAll('.u2-msg-option')[2], 'pointerdown');
  await flush();
  assert.equal(control.value, `${markupOf(PEOPLE[2])} `);
  control.dispose();
});

smoke('minChars: the popup hints until the query is worth searching', async () => {
  let calls = 0;
  const control = compose({mentions: [people({
    minChars: 2,
    fetch: (query) => {
      calls++;
      return Promise.resolve(PEOPLE.filter((p) => p.name.toLowerCase().includes(query)));
    },
  })]});
  type(control, '@a');
  await flush();
  assert.equal(calls, 0, 'a short query never reaches the provider');
  assert.equal(popup().querySelector('.u2-msg-hint').textContent, 'Type 2 or more characters');

  typeQuery(control, 'l');
  await flush();
  assert.equal(calls, 1);
  assert.deepEqual([...popup().querySelectorAll('.u2-msg-option')].map((o) => o.textContent),
    ['Ada Almeida', 'Alan Bauer']);
  control.dispose();
});

smoke('dispose with an open popup releases every session scope', async () => {
  const before = Scope.liveCount;
  const control = compose();
  type(control, '@a');
  await flush();
  assert.ok(Scope.liveCount > before);

  const open = popup();
  control.dispose();
  assert.equal(open.isConnected, false);
  assert.equal(Scope.liveCount, before);

  type(control, '@a');
  await flush();
  assert.equal(popup(), null, 'listeners died with the scope');
});
