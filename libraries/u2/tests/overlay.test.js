/* Overlay: the shared popup layer's lifecycle — what closing does, and what it leaves behind on
   the owner's scope. Long-lived inputs (combobox, date, multi-select, typeahead, tags, font) open
   and close the same overlay hundreds of times over their life. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/index.js';
import {Overlay} from '../src/core/overlay.js';

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
  assert.equal(scope._disposers.length, 0, 'a closed overlay leaves no cleanup behind');

  Overlay.show(anchor, content, scope);
  assert.equal(scope._disposers.length, 1, 'the open one is still owned');
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
