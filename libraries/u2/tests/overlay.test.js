/* Overlay: the shared popup layer's lifecycle — what closing does, and what it leaves behind on
   the owner's scope. Long-lived inputs (combobox, date, multi-select, typeahead, tags, font) open
   and close the same overlay hundreds of times over their life. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/index.js';
import {Overlay, OVERLAY_CLOSE_EVENT} from '../src/core/overlay.js';

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

function anchored() {
  const anchor = document.createElement('div');
  document.body.append(anchor);
  return anchor;
}

smoke('reopening never piles disposers up on the owner', async () => {
  const scope = new Scope();
  const anchor = anchored();
  const content = document.createElement('div');

  for (let i = 0; i < 10; i++) {
    const close = Overlay.show(anchor, content, scope);
    assert.equal(content.isConnected, true);
    close();
    close();
    assert.equal(content.isConnected, false);
  }
  // the private list IS the leak: one dead closure per open/close round is what accumulates
  assert.equal(scope._disposers?.size ?? 0, 0, 'a closed overlay leaves no cleanup behind');

  Overlay.show(anchor, content, scope);
  assert.equal(scope._disposers?.size ?? 0, 1, 'the open one is still owned');
  scope.dispose();
  assert.equal(content.isConnected, false, 'and disposing the owner closes it');
});

smoke('a cleanup that disowns itself does not skip its neighbours', async () => {
  const scope = new Scope();
  const anchor = anchored();
  const ran = [];

  scope.own(() => ran.push('before'));
  Overlay.show(anchor, document.createElement('div'), scope);
  scope.own(() => ran.push('after'));
  scope.dispose();
  assert.deepEqual(ran, ['before', 'after']);
});

smoke('Escape closes the top popup only, stops there, and returns focus that was inside it to the owner', () => {
  const scope = new Scope();
  const anchorA = anchored();
  const anchorB = document.createElement('button');
  document.body.append(anchorB);
  const a = document.createElement('div');
  const b = document.createElement('div');
  const field = document.createElement('input');
  b.append(field);
  Overlay.show(anchorA, a, scope);
  Overlay.show(anchorB, b, scope);
  field.focus();
  let below = 0;
  anchorB.addEventListener('keydown', () => below++);
  assert.equal(fire(field, 'keydown', {key: 'Escape'}), false, 'taken');
  assert.deepEqual([a.isConnected, b.isConnected], [true, false]);
  assert.equal(document.activeElement, anchorB, 'focus back on the owner');
  fire(anchorB, 'keydown', {key: 'Escape'});
  assert.equal(below, 0, 'the control under the popup never sees the Escape');
  assert.equal(a.isConnected, false);
  scope.dispose();
});

smoke('an outside press closes the popups above the entry it lands in; a press on an anchor keeps its popup', () => {
  const scope = new Scope();
  const anchorA = anchored();
  const a = document.createElement('div');
  const inner = document.createElement('button');
  a.append(inner);
  Overlay.show(anchorA, a, scope);
  const b = document.createElement('div');
  Overlay.show(inner, b, scope);
  fire(anchorA, 'pointerdown');
  assert.deepEqual([a.isConnected, b.isConnected], [true, false], 'a press on A\'s anchor closes only what is above A');
  Overlay.show(inner, b, scope);
  fire(inner, 'pointerdown');
  assert.deepEqual([a.isConnected, b.isConnected], [true, true], 'the anchor of B keeps B');
  fire(document.body, 'pointerdown');
  assert.deepEqual([a.isConnected, b.isConnected], [false, false]);
  scope.dispose();
});

smoke('a removed anchor closes its popup through the one body observer', async () => {
  const scope = new Scope();
  const anchor = anchored();
  const content = document.createElement('div');
  let closed = 0;
  content.addEventListener(OVERLAY_CLOSE_EVENT, () => closed++);
  Overlay.show(anchor, content, scope);
  anchor.remove();
  await flush();
  assert.deepEqual([content.isConnected, closed], [false, 1]);
  scope.dispose();
});

smoke('the owner announces the popup: aria-controls, aria-expanded, aria-haspopup; data-u2-part on the popup', () => {
  const scope = new Scope();
  const anchor = document.createElement('div');
  const editor = document.createElement('input');
  anchor.append(editor);
  document.body.append(anchor);
  const content = document.createElement('div');
  content.setAttribute('role', 'listbox');
  const close = Overlay.show(anchor, content, scope);
  assert.match(content.id, /^u2-\d+-popup$/);
  assert.equal(content.dataset.u2Part, 'popup');
  assert.deepEqual(['aria-controls', 'aria-expanded', 'aria-haspopup'].map((a) => editor.getAttribute(a)),
    [content.id, 'true', 'listbox'], 'on the focus target inside the anchor');
  close();
  assert.deepEqual([editor.getAttribute('aria-controls'), editor.getAttribute('aria-expanded')], [null, 'false']);

  const preset = document.createElement('div');
  preset.id = 'own-popup';
  editor.setAttribute('aria-controls', 'own-popup');
  Overlay.show(anchor, preset, scope)();
  assert.equal(editor.getAttribute('aria-controls'), 'own-popup', 'a link the owner made stays');
  scope.dispose();
});

smoke('positioning follows events, not a frame loop', async () => {
  const scope = new Scope();
  const raf = globalThis.requestAnimationFrame;
  let frames = 0;
  globalThis.requestAnimationFrame = (fn) => {
    frames++;
    return raf(fn);
  };
  try {
    Overlay.show(anchored(), document.createElement('div'), scope);
    await flush();
    await flush();
  } finally {
    globalThis.requestAnimationFrame = raf;
  }
  assert.equal(frames, 0);
  scope.dispose();
});
