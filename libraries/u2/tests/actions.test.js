/* Item actions: the hover block renders the icon-bearing subset, the context menu the full list,
   and VirtualList wires right-click to it. The timestamp factory rides along — same recipe. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {timestamp} from '../src/core/elements.js';
import {rowActions, actionsMenu} from '../src/components/actions.js';
import {VirtualList} from '../src/components/list.js';

function scoped(name, body) {
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

function sampleActions(log) {
  return [
    {name: 'Copy', icon: 'copy', run: () => log.push('copy')},
    {name: 'Open ticket', icon: 'external-link-alt', run: () => log.push('open')},
    {name: 'Copy id', run: () => log.push('id')},
  ];
}

scoped('rowActions: icon buttons for icon-bearing actions only, tooltip carries the name', () => {
  const log = [];
  const block = rowActions(sampleActions(log));
  assert.equal(block.dataset.u2, 'row-actions');
  const buttons = [...block.querySelectorAll('[data-u2="icon-button"]')];
  assert.equal(buttons.length, 2, 'the icon-less action appears only in menus');
  assert.deepEqual(buttons.map((b) => b.getAttribute('aria-label')), ['Copy', 'Open ticket']);
  fire(buttons[0], 'click');
  assert.deepEqual(log, ['copy']);
});

scoped('actionsMenu: every action becomes an item, picking one runs it and closes', () => {
  const log = [];
  const menu = actionsMenu(sampleActions(log));
  menu.show({x: 10, y: 10});
  const items = [...document.querySelectorAll('[role="menuitem"]')];
  assert.deepEqual(items.map((el) => el.querySelector('.u2-menu-label').textContent),
    ['Copy', 'Open ticket', 'Copy id']);
  fire(items[2], 'click');
  assert.deepEqual(log, ['id']);
  assert.equal(menu.isOpen.value, false);
});

scoped('VirtualList contextActions: right-click selects the row and opens the full menu', () => {
  const log = [];
  const list = new VirtualList({
    itemHeight: 20,
    render: (item) => {
      const el = document.createElement('span');
      el.textContent = item;
      return el;
    },
    contextActions: (item) => [{name: `Act on ${item}`, run: () => log.push(item)}],
  });
  document.body.append(list.root);
  list.setItems(['a', 'b', 'c']);

  const row = list.root.querySelector('[data-index="1"]');
  fire(row.firstChild, 'contextmenu');
  assert.equal(list.selectedIndex.value, 1, 'right-click selects');
  const item = document.querySelector('[role="menuitem"]');
  assert.equal(item.textContent.trim(), 'Act on b');
  fire(item, 'click');
  assert.deepEqual(log, ['b']);
  list.dispose();
});

scoped('timestamp: short visual, full title, year only when not current, empty when invalid', () => {
  const now = new Date();
  const thisYear = new Date(now.getFullYear(), 7, 11, 14, 30);
  const el = timestamp(thisYear);
  assert.equal(el.className, 'u2-timestamp');
  assert.equal(el.textContent,
    thisYear.toLocaleDateString(undefined, {month: 'short', day: 'numeric'}));
  assert.equal(el.title, thisYear.toLocaleString(undefined,
    {month: 'short', day: 'numeric', year: 'numeric', hour: '2-digit', minute: '2-digit'}));

  const old = new Date(2019, 2, 5);
  assert.equal(timestamp(old).textContent,
    old.toLocaleDateString(undefined, {month: 'short', day: 'numeric', year: 'numeric'}));

  assert.equal(timestamp({toDate: () => thisYear}).textContent, el.textContent,
    'anything with toDate() — dayjs included — is accepted');
  assert.equal(timestamp('not a date').textContent, '');
  assert.equal(timestamp('not a date').title, '');
});
