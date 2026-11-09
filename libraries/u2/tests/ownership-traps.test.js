/* Ownership traps: what a row, a handler or a field builds must go with what it belongs to. A context menu closes
   with its owner; the async primitives go with whoever built them; a balloon's content is the balloon's. Every test
   leaves the live-scope count where it found it. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {layout} from './layout.mjs';
import {signal, Scope, Control, AsyncSource, AsyncPager, AsyncValue, Tooltip, Toolbar, button, confirmDialog, div,
  each, link, notify, span} from '../src/index.js';
import {Menu} from '../src/components/navigation/menu.js';

function owned(name, body) {
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

owned('a context menu closes when the row it was opened on goes', () => {
  const rows = signal(['a', 'b']);
  const menu = Menu.popup().item('Delete', () => {});
  const host = div([]);
  const list = Control.build(() => {
    each(host, rows, (key) => key, (key) => {
      const row = div([key]);
      menu.bindContext(row);
      return row;
    });
    return host;
  });
  document.body.append(list.root);
  fire(host.children[0], 'contextmenu', {clientX: 10, clientY: 10});
  assert.equal(menu.isOpen.peek(), true);
  rows.value = ['b'];
  assert.equal(menu.isOpen.peek(), false, 'the row gone, its menu goes');
  assert.equal(document.body.querySelector('.u2-menu'), null);

  fire(host.children[0], 'contextmenu', {clientX: 10, clientY: 10});
  assert.equal(menu.isOpen.peek(), true);
  rows.value = ['b', 'c'];
  assert.equal(menu.isOpen.peek(), true, 'another row\'s change leaves it open');
  list.dispose();
  assert.equal(menu.isOpen.peek(), false);
});

owned('bindContext takes the builder as its owner; without one it says so', () => {
  assert.throws(() => Menu.popup().bindContext(div([])), /bindContext\(\) needs an owner/);
});

owned('AsyncSource and AsyncPager go with whoever built them, as AsyncValue does', async () => {
  let resolveSource;
  let resolvePage;
  let source;
  let pager;
  const control = Control.build(() => {
    source = new AsyncSource(() => new Promise((r) => resolveSource = r), {debounceMs: 0});
    pager = new AsyncPager({fetchPage: () => new Promise((r) => resolvePage = r), pageSize: 2});
    return div([]);
  });
  source.query('a');
  pager.reset();
  await flush();
  control.dispose();
  resolveSource(['x']);
  resolvePage(['p1', 'p2']);
  await flush();
  assert.equal(source.state.peek().kind, 'loading', 'a disposed source drops the late answer');
  assert.deepEqual(pager.items.peek(), [], 'a disposed pager drops the late page');
});

owned('a source disposed early leaves nothing on its long-lived owner', () => {
  const control = Control.build(() => div([]));
  const sizes = [];
  for (let i = 0; i < 3; i++) {
    const source = control.build(() => new AsyncSource(async () => [], {debounceMs: 0}));
    source.dispose();
    sizes.push(control.scope._disposers?.size ?? 0);
  }
  assert.equal(new Set(sizes).size, 1, 'no closure piles up per source');
  control.dispose();
});

owned('the owner option: after an await, the async primitives go with the owner they are given', async () => {
  const session = Control.build(() => div([]));
  await flush();
  const value = new AsyncValue(async () => 1, {auto: false, owner: session});
  const source = new AsyncSource(async () => ['x'], {debounceMs: 0, owner: session.scope});
  const pager = new AsyncPager({fetchPage: async () => ['p'], owner: session});
  session.dispose();
  await value.refresh();
  source.query('x');
  pager.reset();
  await flush();
  assert.equal(value.state.peek().kind, 'idle');
  assert.equal(source.state.peek().kind, 'idle');
  assert.deepEqual(pager.items.peek(), []);
});

owned('Menu.bindContext takes a control as its owner too', () => {
  const row = Control.build(() => div([]));
  const menu = Menu.popup().item('Open', () => {});
  menu.bindContext(row.root, {owner: row});
  document.body.append(row.root);
  fire(row.root, 'contextmenu', {clientX: 1, clientY: 1});
  assert.equal(menu.isOpen.peek(), true);
  row.dispose();
  assert.equal(menu.isOpen.peek(), false);
});

test('Tooltip.titleIfCut: the whole text while the line is cut, gone once it fits; a title of its own stays', () => {
  const cell = span('Kinase hit-to-lead follow-up');
  layout(cell, {scrollWidth: 300, clientWidth: 100});
  Tooltip.titleIfCut(cell);
  assert.equal(cell.title, 'Kinase hit-to-lead follow-up');
  layout(cell, {clientWidth: 400});
  Tooltip.titleIfCut(cell);
  assert.equal(cell.getAttribute('title'), null);
  const own = span('Oct 1, 2026');
  own.title = 'Oct 1, 2026, 09:30 AM';
  layout(own, {scrollWidth: 300, clientWidth: 100});
  Tooltip.titleIfCut(own);
  assert.equal(own.title, 'Oct 1, 2026, 09:30 AM');
});

test('a toolbar and a ribbon group fold into the same "…" — More actions', () => {
  const owner = Scope.root('test');
  const bar = owner.build(() => new Toolbar());
  const overflow = bar.root.querySelector('.u2-toolbar-overflow');
  assert.equal(overflow.getAttribute('aria-label'), Toolbar.OVERFLOW_LABEL);
  assert.equal(Toolbar.OVERFLOW_LABEL, 'More actions');
  assert.ok(overflow.querySelector(`.fa-${Toolbar.OVERFLOW_ICON}`));
  owner.dispose();
});

owned('a failed page keeps its message', async () => {
  const owner = Scope.root('test');
  const pager = owner.build(() => new AsyncPager({fetchPage: async () => { throw new Error('Server is down'); }}));
  pager.reset();
  await flush();
  assert.equal(pager.state.peek(), 'error');
  assert.equal(pager.error.peek(), 'Server is down');
  owner.dispose();
});

owned('a balloon builds its content in its own scope, ends it on close and says so', async () => {
  const undone = [];
  let closed = 0;
  const deleted = signal('KHL-3');
  const row = Control.build(() => div([]));
  const handle = row.build(() => notify.info(
    () => div([span(deleted), ' deleted. ', link('Undo', () => undone.push(deleted.peek()))]),
    {onClose: () => closed++, autoHide: false}));
  row.dispose();
  const undo = handle.root.querySelector('.u2-link');
  assert.equal(handle.root.textContent.includes('KHL-3 deleted.'), true);
  deleted.value = 'KHL-4';
  assert.equal(handle.root.textContent.includes('KHL-4 deleted.'), true, 'its live text is its own');
  undo.dispatchEvent(new MouseEvent('click', {bubbles: false}));
  assert.deepEqual(undone, ['KHL-4']);
  handle.close();
  assert.equal(closed, 1);
  undo.dispatchEvent(new MouseEvent('click', {bubbles: false}));
  assert.deepEqual(undone, ['KHL-4'], 'its handlers end with it');
  handle.close();
  assert.equal(closed, 1, 'once');
});

owned('a dialog or a confirm opened from a row stays open and answers through the row\'s rebuild', async () => {
  const rows = signal([{id: 'a', n: 1}]);
  let answer;
  const host = div([]);
  const list = Control.build(() => {
    each(host, rows, (x) => x.id, (x) => button(`Delete ${x.n}`, () => {
      answer = confirmDialog({title: 'Delete', message: `Delete ${x.n}?`, confirmText: 'Delete'});
    }));
    return host;
  });
  document.body.append(list.root);
  fire(host.querySelector('button'), 'click');
  rows.value = [{id: 'a', n: 2}];
  const dialog = document.querySelector('.u2-dialog');
  assert.ok(dialog?.isConnected, 'the row rebuilt, the confirm stays');
  fire([...dialog.querySelectorAll('button')].find((b) => b.textContent.trim() === 'Delete'), 'click');
  assert.equal(await answer, true);
  list.dispose();
});
