/* Menu bar on the node DOM shim: opening and switching by click, hover and keyboard, the openMenu
   signal, disabled items, and disposal. Same contract as smoke.test.js — every test leaves the
   live-scope count where it found it. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/index.js';
import {MenuBar} from '../src/components/navigation/menu-bar.js';

function ui(name, body) {
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

function menuPanel() {
  return document.body.querySelector('.u2-menu');
}

function menuItems() {
  return menuPanel().querySelectorAll('.u2-menu-item');
}

function menuLabels() {
  return menuPanel().querySelectorAll('.u2-menu-label').map((el) => el.textContent);
}

function items(bar) {
  return bar.root.querySelectorAll('.u2-menu-bar-item');
}

function tabStops(bar) {
  return items(bar).filter((el) => el.tabIndex === 0).map((el) => el.textContent);
}

function expanded(bar) {
  return items(bar).filter((el) => el.getAttribute('aria-expanded') === 'true').map((el) => el.textContent);
}

function mount(fired, options) {
  const bar = new MenuBar(options);
  bar
    .menu('File', (m) => m.item('New', () => fired.push('New')).item('Open', () => fired.push('Open')))
    .menu('Edit', (m) => m.item('Undo', () => fired.push('Undo')))
    .menu('View', (m) => m.item('Grid', () => fired.push('Grid')));
  document.body.append(bar.root);
  return bar;
}

ui('MenuBar: click opens a menu anchored to the item, and openMenu names it', async () => {
  const bar = mount([]);
  assert.equal(bar.root.dataset.u2, 'menu-bar');
  assert.equal(bar.root.getAttribute('role'), 'menubar');
  assert.equal(items(bar).length, 3);
  assert.equal(items(bar)[0].getAttribute('role'), 'menuitem');
  assert.equal(items(bar)[0].getAttribute('aria-haspopup'), 'menu');
  assert.equal(items(bar)[0].getAttribute('aria-expanded'), 'false');
  assert.equal(bar.openMenu.value, null);
  assert.deepEqual(tabStops(bar), ['File'], 'a single tab stop');

  fire(items(bar)[0], 'click');
  assert.equal(bar.openMenu.value, 'File');
  assert.deepEqual(expanded(bar), ['File']);
  assert.equal(items(bar)[0].classList.contains('u2-menu-bar-item-open'), true);
  assert.deepEqual(menuLabels(), ['New', 'Open']);

  bar.dispose();
  await flush();
});

ui('MenuBar: clicking the open item closes it, clicking another switches', async () => {
  const bar = mount([]);
  fire(items(bar)[0], 'click');
  const opened = menuPanel();

  fire(items(bar)[0], 'click');
  assert.equal(opened.isConnected, false);
  assert.equal(menuPanel(), null);
  assert.equal(bar.openMenu.value, null);
  assert.deepEqual(expanded(bar), []);

  fire(items(bar)[0], 'click');
  fire(items(bar)[1], 'click');
  assert.equal(bar.openMenu.value, 'Edit');
  assert.deepEqual(expanded(bar), ['Edit'], 'only one item is expanded at a time');
  assert.deepEqual(menuLabels(), ['Undo']);

  bar.dispose();
  await flush();
});

ui('MenuBar: hovering another item switches the open menu, and does nothing while closed', async () => {
  const bar = mount([]);
  fire(items(bar)[1], 'pointerenter');
  assert.equal(bar.openMenu.value, null, 'hover alone never opens the bar');

  fire(items(bar)[0], 'click');
  fire(items(bar)[2], 'pointerenter');
  assert.equal(bar.openMenu.value, 'View');
  assert.deepEqual(menuLabels(), ['Grid']);
  assert.deepEqual(expanded(bar), ['View']);

  fire(items(bar)[2], 'pointerenter');
  assert.equal(bar.openMenu.value, 'View', 'hovering the open item keeps its menu');

  bar.dispose();
  await flush();
});

ui('MenuBar: ←/→ rove the tab stop, Home/End jump, and both switch the open menu', async () => {
  const bar = mount([]);
  fire(bar.root, 'keydown', {key: 'ArrowRight'});
  assert.deepEqual(tabStops(bar), ['Edit']);
  assert.equal(document.activeElement, items(bar)[1]);
  assert.equal(bar.openMenu.value, null, 'roving alone opens nothing');

  fire(bar.root, 'keydown', {key: 'ArrowRight'});
  fire(bar.root, 'keydown', {key: 'ArrowRight'});
  assert.deepEqual(tabStops(bar), ['File'], 'the roving wraps around');
  fire(bar.root, 'keydown', {key: 'ArrowLeft'});
  assert.deepEqual(tabStops(bar), ['View']);
  fire(bar.root, 'keydown', {key: 'Home'});
  assert.deepEqual(tabStops(bar), ['File']);
  fire(bar.root, 'keydown', {key: 'End'});
  assert.deepEqual(tabStops(bar), ['View']);

  fire(bar.root, 'keydown', {key: 'Home'});
  fire(bar.root, 'keydown', {key: 'ArrowDown'});
  assert.equal(bar.openMenu.value, 'File');
  assert.equal(document.activeElement, menuPanel(), 'the focus moves into the menu');

  fire(menuPanel(), 'keydown', {key: 'ArrowRight'});
  assert.equal(bar.openMenu.value, 'Edit', '→ from inside the menu switches the whole menu');
  assert.deepEqual(menuLabels(), ['Undo']);
  fire(menuPanel(), 'keydown', {key: 'ArrowLeft'});
  assert.equal(bar.openMenu.value, 'File');

  bar.dispose();
  await flush();
});

ui('MenuBar: activating an item runs it, closes everything and clears openMenu', async () => {
  const fired = [];
  const bar = mount(fired);
  fire(items(bar)[0], 'click');
  fire(menuItems()[1], 'click');
  assert.deepEqual(fired, ['Open']);
  assert.equal(menuPanel(), null);
  assert.equal(bar.openMenu.value, null);
  assert.deepEqual(expanded(bar), []);
  assert.equal(document.activeElement, items(bar)[0], 'the focus comes back to the bar item');

  fire(items(bar)[0], 'click');
  fire(menuItems()[0], 'click');
  assert.deepEqual(fired, ['Open', 'New'], 'the reopened menu is live too');

  bar.dispose();
  await flush();
});

ui('MenuBar: Esc closes the open menu and returns focus to the item', async () => {
  const bar = mount([]);
  fire(items(bar)[1], 'click');
  assert.equal(bar.openMenu.value, 'Edit');

  fire(document, 'keydown', {key: 'Escape'});
  assert.equal(menuPanel(), null);
  assert.equal(bar.openMenu.value, null);
  assert.equal(document.activeElement, items(bar)[1]);

  fire(bar.root, 'keydown', {key: 'Escape'});
  assert.equal(bar.openMenu.value, null, 'Esc on a closed bar is a no-op');

  bar.dispose();
  await flush();
});

ui('MenuBar: a disabled item never opens, is skipped by roving and by hover switching', async () => {
  const bar = mount([], {items: [{label: 'Tools', build: (m) => m.item('Options', () => {}), enabled: false}]});
  const [tools, file, edit] = items(bar);
  assert.equal(tools.getAttribute('aria-disabled'), 'true');
  assert.equal(file.getAttribute('aria-disabled'), 'false');
  assert.deepEqual(tabStops(bar), ['File'], 'the tab stop skips the disabled item');

  fire(tools, 'click');
  assert.equal(bar.openMenu.value, null);

  fire(file, 'click');
  fire(tools, 'pointerenter');
  assert.equal(bar.openMenu.value, 'File', 'a disabled item is never hovered into');
  fire(bar.root, 'keydown', {key: 'ArrowLeft'});
  assert.equal(bar.openMenu.value, 'View', '← from File wraps past Tools');

  bar.setEnabled('Tools', true);
  assert.equal(tools.getAttribute('aria-disabled'), 'false');
  fire(tools, 'click');
  assert.equal(bar.openMenu.value, 'Tools');
  assert.deepEqual(tabStops(bar), ['Tools']);

  bar.setEnabled('Tools', false);
  assert.equal(bar.openMenu.value, null, 'disabling the open item closes it');
  assert.equal(menuPanel(), null);
  fire(bar.root, 'keydown', {key: 'Home'});
  assert.deepEqual(tabStops(bar), ['File']);

  bar.setEnabled('Missing', true);
  assert.equal(items(bar).length, 4, 'an unknown label is ignored');
  assert.equal(edit.getAttribute('aria-disabled'), 'false');

  bar.dispose();
  await flush();
});

ui('MenuBar: dispose closes the open menu and releases every listener', async () => {
  const fired = [];
  const bar = mount(fired);
  fire(items(bar)[0], 'click');
  assert.equal(bar.openMenu.value, 'File');

  bar.dispose();
  assert.equal(menuPanel(), null, 'disposal closes the open menu');
  assert.equal(bar.openMenu.value, null);
  assert.deepEqual(expanded(bar), []);

  fire(items(bar)[0], 'click');
  fire(items(bar)[1], 'pointerenter');
  fire(bar.root, 'keydown', {key: 'ArrowDown'});
  assert.equal(menuPanel(), null, 'the click, hover and key listeners died with the bar');
  assert.deepEqual(fired, []);
  await flush();
});
