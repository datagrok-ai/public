/* notify on a platform page: routed into the platform's balloon stack through grok.shell (src/dg/shell/
   balloons.ts) — types, text-only strings, timeouts in seconds, u2's own keys; off the platform u2
   keeps its toasts. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from './dom-shim.js';
import {Scope, notify} from '../src/index.js';

register('./dg-stub.mjs', import.meta.url);
const grok = await import('datagrok-api/grok');
await import('../src/dg/shell/balloons.js');

/** The platform's side, as `balloon.dart` does it: a `.d4-balloon` appended to one container. */
function platform() {
  const calls = [];
  const show = (type) => (content, options) => {
    calls.push({type, content, options});
    let container = document.querySelector('.d4-balloon-container');
    if (container === null) {
      container = document.createElement('div');
      container.className = 'd4-balloon-container';
      document.body.append(container);
    }
    const balloon = document.createElement('div');
    balloon.className = `d4-balloon ${type}`;
    balloon.append(typeof content === 'string' ? content : content);
    container.append(balloon);
  };
  const saved = {info: grok.shell.info, warning: grok.shell.warning, error: grok.shell.error};
  Object.assign(grok.shell, {info: show('info'), warning: show('warning'), error: show('error')});
  globalThis.grok_Balloon = () => {};
  return {calls, restore: () => {
    Object.assign(grok.shell, saved);
    delete globalThis.grok_Balloon;
  }};
}

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

smoke('off the platform u2 keeps its own toasts', () => {
  const handle = notify.info('Saved');
  assert.equal(handle.root.classList.contains('u2-notify'), true);
  assert.equal(document.querySelector('.d4-balloon'), null);
});

smoke('on a platform page every type goes into the platform stack, a string as text', async () => {
  const {calls, restore} = platform();
  try {
    const info = notify.info('Saved <b>3</b> rows', {copyText: 'details'});
    notify.warning('Careful', {timeout: 3000});
    notify.error('Refused');
    assert.equal(document.querySelector('.u2-notify'), null, 'no second stack');
    assert.deepEqual(calls.map((c) => c.type), ['info', 'warning', 'error']);
    assert.equal(calls[0].content.textContent, 'Saved <b>3</b> rows', 'never markup');
    assert.equal(calls[0].content.querySelector?.('b') ?? null, null);
    assert.deepEqual(calls[0].options, {autoHide: true, copyText: 'details', timeout: 5});
    assert.deepEqual(calls[1].options, {autoHide: false, copyText: undefined, timeout: 3}, 'a warning stays, seconds');
    assert.equal(calls[2].content, 'Refused', 'an error string the platform shows as text already');
    assert.equal(info.root.classList.contains('d4-balloon'), true, 'the handle is the platform balloon');
    info.close();
    assert.equal(info.root.isConnected, false);
    notify.closeAll();
    assert.equal(document.querySelectorAll('.d4-balloon').length, 0, 'closeAll reaches them too');
  } finally {
    restore();
  }
});

smoke('u2 keeps singleKey and oneTimeKey; a balloon the platform closes releases its key', async () => {
  const {calls, restore} = platform();
  try {
    const first = notify.warning('Busy', {singleKey: 'busy'});
    assert.equal(notify.warning('Busy', {singleKey: 'busy'}), null, 'one at a time');
    first.root.remove();
    await flush();
    assert.notEqual(notify.warning('Busy', {singleKey: 'busy'}), null, 'released once the platform closed it');
    assert.equal(calls.length, 2);
  } finally {
    restore();
  }
});

smoke('a balloon the presenter released before it returned is not held open', () => {
  const saved = notify.presenter;
  let closes = 0;
  notify.presenter = (content, type, options, closed) => {
    closed();
    return {root: document.createElement('div'), close: () => closes++};
  };
  try {
    notify.info('Saved', {singleKey: 'saved'});
    notify.closeAll();
    assert.equal(closes, 0, 'closeAll has nothing of it to close');
    assert.notEqual(notify.info('Saved', {singleKey: 'saved'}), null, 'its key is free');
  } finally {
    notify.presenter = saved;
  }
});
