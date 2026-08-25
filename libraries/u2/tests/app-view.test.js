/* appView(): the shell chrome over a u2 content component (V-4): the content is released through
   the platform kill channel — `data-kill-on-close` + `Widget.registerCleanup` — and the chrome
   components are disposed with it through `content.own`. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from './dom-shim.js';
import {signal} from '../src/core/signals.js';
import {Scope} from '../src/core/scope.js';
import {Control} from '../src/core/component.js';
import {platform} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const {appView} = await import('../src/dg/shell/app-view.js');
const DG = await import('datagrok-api/dg');

function viewed(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      platform.reset();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

viewed('the view is over the content root, released through the kill channel', () => {
  const content = new Control();
  const view = appView({name: 'App', content});
  assert.ok(view instanceof DG.View);
  assert.equal(view.root, content.root);
  assert.equal(view.name, 'App');
  assert.equal(content.root.getAttribute('data-kill-on-close'), 'true');
  assert.equal(platform.cleanups.length, 1);

  globalThis.grok_Widget_Kill(content.root);
  assert.equal(content.scope.isDisposed, true, 'the registered cleanup disposes the content');
  assert.equal(platform.cleanups.length, 0, 'the cleanup ran once and is gone');
});

viewed('chrome components are disposed with the content; plain elements pass as given', () => {
  const content = new Control();
  const ribbonItem = new Control();
  const el = document.createElement('span');
  const toolbox = new Control();
  const view = appView({name: 'App', content, ribbon: [[ribbonItem, el]], toolbox});
  assert.deepEqual(view.ribbonPanels, [[ribbonItem.root, el]]);
  assert.equal(view.toolbox, toolbox.root);
  content.dispose();
  assert.equal(ribbonItem.scope.isDisposed, true);
  assert.equal(toolbox.scope.isDisposed, true);
});

viewed('a status signal renders as a live text panel whose effect dies with the content', () => {
  const status = signal('loading');
  const content = new Control();
  const view = appView({name: 'App', content, status});
  const [el] = view.statusBarPanels;
  assert.equal(el.textContent, 'loading');
  status.value = 'ready';
  assert.equal(el.textContent, 'ready');
  content.dispose();
  status.value = 'later';
  assert.equal(el.textContent, 'ready', 'the effect died with the content');
});
