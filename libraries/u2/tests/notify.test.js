/* notify: the toast service contract — per-type balloons in the top-right stack, auto-hide with
   a hover pause, one-time and single keys, copy/close affordances, string content as text. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope, notify} from '../src/index.js';

function smoke(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      notify.closeAll();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

function sleep(ms) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

function balloons() {
  return notify.host.querySelectorAll('.u2-notify');
}

smoke('per-type balloons stack in the container; strings render as text', async () => {
  const info = notify.info('Saved');
  const warning = notify.warning('Careful');
  const error = notify.error('<b>Broke</b>');
  assert.equal(balloons().length, 3);
  assert.ok(info.root.classList.contains('u2-notify-info'));
  assert.ok(warning.root.classList.contains('u2-notify-warning'));
  assert.ok(error.root.classList.contains('u2-notify-error'));
  assert.equal(info.root.getAttribute('role'), 'status');
  assert.equal(error.root.getAttribute('role'), 'alert');
  const content = error.root.querySelector('.u2-notify-content');
  assert.equal(content.textContent, '<b>Broke</b>', 'markup in a string stays literal');
  assert.equal(content.querySelector('b'), null, 'no element was created from it');

  notify.closeAll();
  assert.equal(balloons().length, 0);
});

smoke('info auto-hides; warning and error stay; explicit autoHide overrides', async () => {
  const info = notify.info('Bye', {timeout: 20});
  const warning = notify.warning('Stay');
  const sticky = notify.info('Also stay', {autoHide: false});
  assert.ok(info.root.querySelector('.u2-notify-timer'), 'auto-hiding balloon shows the countdown');
  assert.equal(warning.root.querySelector('.u2-notify-timer'), null);
  assert.equal(sticky.root.querySelector('.u2-notify-timer'), null);

  await sleep(50);
  assert.equal(info.root.isConnected, false, 'info closed itself');
  assert.equal(warning.root.isConnected, true);
  assert.equal(sticky.root.isConnected, true);
});

smoke('hovering pauses the countdown; leaving resumes it', async () => {
  const handle = notify.info('Hover me', {timeout: 30});
  fire(handle.root, 'pointerenter');
  await sleep(60);
  assert.equal(handle.root.isConnected, true, 'paused while hovered');
  fire(handle.root, 'pointerleave');
  await sleep(60);
  assert.equal(handle.root.isConnected, false, 'resumed and closed after the pointer left');
});

smoke('a click closes the balloon; the copy button does not', async () => {
  const handle = notify.error('Click to dismiss');
  fire(handle.root.querySelector('.u2-icon-btn'), 'click');
  assert.equal(handle.root.isConnected, true, 'copy keeps the balloon open');
  fire(handle.root, 'click');
  assert.equal(handle.root.isConnected, false);

  const pinned = notify.error('Not by click', {closeOnClick: false});
  fire(pinned.root, 'click');
  assert.equal(pinned.root.isConnected, true);
  const [, close] = pinned.root.querySelectorAll('.u2-icon-btn');
  fire(close, 'click');
  assert.equal(pinned.root.isConnected, false, 'the close button always works');
});

smoke('singleKey allows one balloon at a time and frees the key on close', async () => {
  const first = notify.warning('Busy', {singleKey: 'busy'});
  assert.notEqual(first, null);
  assert.equal(notify.warning('Busy', {singleKey: 'busy'}), null, 'second is suppressed');
  first.close();
  first.close();
  assert.notEqual(notify.warning('Busy', {singleKey: 'busy'}), null, 'free again after close');
});

smoke('oneTimeKey persists through localStorage; without storage every show is a first', async () => {
  const original = Object.getOwnPropertyDescriptor(globalThis, 'localStorage');
  const store = new Map();
  Object.defineProperty(globalThis, 'localStorage', {configurable: true, value: {
    getItem: (k) => store.has(k) ? store.get(k) : null,
    setItem: (k, v) => store.set(k, String(v)),
  }});
  try {
    assert.notEqual(notify.info('Welcome', {oneTimeKey: 'welcome'}), null);
    assert.equal(notify.info('Welcome', {oneTimeKey: 'welcome'}), null, 'remembered');
    assert.ok(store.has('u2-notify-welcome'));
  } finally {
    if (original)
      Object.defineProperty(globalThis, 'localStorage', original);
    else
      delete globalThis.localStorage;
  }
  assert.notEqual(notify.info('Welcome', {oneTimeKey: 'other'}), null,
    'storage-less hosts still show');
});

smoke('element content mounts as given; the container is recreated after a DOM reset', async () => {
  const el = document.createElement('div');
  el.textContent = 'rich';
  const handle = notify.info(el);
  assert.ok(handle.root.querySelector('.u2-notify-content').contains(el));
  handle.close();

  resetDom();
  const after = notify.info('Fresh');
  assert.equal(after.root.isConnected, true, 'a fresh container was minted');
  after.close();
});
