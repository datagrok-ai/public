/* Icon buttons and the button group on the node DOM shim: composition and icon side, the toggle
   signal round trip, radio vs multi vs action semantics, roving focus, auto-tooltips and enabling.
   Same contract as smoke.test.js — every test leaves the live-scope count where it found it. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {signal, Scope, Control} from '../src/index.js';
import {iconButton, buttonWithIcon, dropDownButton} from '../src/components/buttons.js';
import {ButtonGroup} from '../src/components/button-group.js';

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

function cells(group) {
  return group.root.querySelectorAll('.u2-btn-group-item');
}

function cell(group, id) {
  return group.root.querySelector(`[data-id="${id}"]`);
}

function tabStops(group) {
  return [...cells(group)].filter((el) => el.tabIndex === 0).map((el) => el.dataset.id);
}

ui('buttonWithIcon: icon then label, primary skin, and iconRight flips the order', () => {
  const clicks = [];
  const el = buttonWithIcon('Save', 'save', () => clicks.push('save'));
  assert.equal(el.classList.contains('u2-btn'), true);
  assert.equal(el.dataset.u2, 'button');
  assert.equal(el.children.length, 2);
  assert.equal(el.children[0].dataset.u2, 'icon');
  assert.equal(el.children[0].classList.contains('fa-save'), true);
  assert.equal(el.children[1].className, 'u2-btn-label');
  assert.equal(el.textContent, 'Save');
  fire(el, 'click');
  assert.deepEqual(clicks, ['save']);

  const primary = buttonWithIcon('Add', 'plus', () => {}, {primary: true});
  assert.equal(primary.classList.contains('u2-btn-primary'), true);

  const right = buttonWithIcon('Next', 'arrow-right', () => {}, {iconRight: true});
  assert.equal(right.children[0].className, 'u2-btn-label');
  assert.equal(right.children[1].dataset.u2, 'icon');
});

ui('buttonWithIcon and iconButton: a tooltip without an owner falls back to title', () => {
  const el = buttonWithIcon('Delete', 'trash-alt', () => {}, {tooltip: 'Delete the selected rows'});
  assert.equal(el.title, 'Delete the selected rows');

  const square = iconButton('pencil', () => {}, {tooltip: 'Edit'});
  assert.equal(square.title, 'Edit');
  assert.equal(square.getAttribute('aria-label'), 'Edit', 'the tooltip names an icon-only button');
});

ui('iconButton: icon-only square button that fires its click', () => {
  const clicks = [];
  const el = iconButton('sync', () => clicks.push('refresh'), {variant: 'solid'});
  assert.equal(el.dataset.u2, 'icon-button');
  assert.equal(el.classList.contains('u2-icon-btn'), true);
  assert.equal(el.children.length, 1, 'no label — the icon is the whole button');
  assert.equal(el.children[0].classList.contains('fas'), true, 'the variant reaches the icon');
  fire(el, 'click');
  fire(el, 'click');
  assert.deepEqual(clicks, ['refresh', 'refresh']);
});

ui('iconButton: the toggle signal drives aria-pressed, and the button drives the signal', () => {
  const starred = signal(false);
  const clicks = [];
  const owner = Control.build(() => iconButton('star', () => clicks.push('star'), {toggle: starred}));
  const el = owner.root;
  assert.equal(el.getAttribute('aria-pressed'), 'false');

  fire(el, 'click');
  assert.equal(starred.value, true);
  assert.equal(el.getAttribute('aria-pressed'), 'true');
  assert.equal(el.classList.contains('u2-icon-btn-pressed'), true);
  assert.deepEqual(clicks, ['star'], 'onClick fires alongside the toggle');

  fire(el, 'click');
  assert.equal(starred.value, false);
  assert.equal(el.getAttribute('aria-pressed'), 'false');

  starred.value = true;
  assert.equal(el.getAttribute('aria-pressed'), 'true', 'a programmatic write shows up too');

  owner.dispose();
  starred.value = false;
  assert.equal(el.getAttribute('aria-pressed'), 'true', 'the effect died with its owner');
});

ui('iconButton: a toggling button needs an ambient owner', () => {
  assert.throws(() => iconButton('star', () => {}, {toggle: signal(false)}), /needs an owner/);
});

ui('ButtonGroup single: radio semantics — selects, never deselects, one tab stop', () => {
  const group = new ButtonGroup({items: [
    {id: 'left', label: 'Left', icon: 'align-left'},
    {id: 'center', label: 'Center', icon: 'align-center'},
    {id: 'right', label: 'Right', icon: 'align-right'},
  ], toggle: 'single'});
  assert.equal(group.root.dataset.u2, 'button-group');
  assert.equal(group.root.getAttribute('role'), 'radiogroup');
  assert.equal(cells(group).length, 3);
  assert.equal(cell(group, 'left').getAttribute('role'), 'radio');
  assert.equal(cell(group, 'left').getAttribute('aria-checked'), 'false');
  assert.deepEqual(group.selected.value, []);
  assert.deepEqual(tabStops(group), ['left'], 'a single tab stop');

  fire(cell(group, 'center'), 'click');
  assert.deepEqual(group.selected.value, ['center']);
  assert.equal(cell(group, 'center').getAttribute('aria-checked'), 'true');
  assert.equal(cell(group, 'center').classList.contains('u2-btn-group-selected'), true);

  fire(cell(group, 'center'), 'click');
  assert.deepEqual(group.selected.value, ['center'], 'clicking the active cell keeps it selected');

  fire(cell(group, 'right'), 'click');
  assert.deepEqual(group.selected.value, ['right']);
  assert.equal(cell(group, 'center').getAttribute('aria-checked'), 'false');
  group.dispose();
});

ui('ButtonGroup: ←/→ rove the single tab stop and move focus', () => {
  const group = new ButtonGroup({items: [
    {id: 'a', label: 'A'}, {id: 'b', label: 'B'}, {id: 'c', label: 'C'},
  ], toggle: 'single'});
  document.body.append(group.root);

  fire(group.root, 'keydown', {key: 'ArrowRight'});
  assert.deepEqual(tabStops(group), ['b']);
  assert.equal(document.activeElement, cell(group, 'b'));

  fire(group.root, 'keydown', {key: 'ArrowRight'});
  fire(group.root, 'keydown', {key: 'ArrowRight'});
  assert.deepEqual(tabStops(group), ['a'], 'the roving wraps around');

  fire(group.root, 'keydown', {key: 'ArrowLeft'});
  assert.deepEqual(tabStops(group), ['c']);
  assert.equal(cell(group, 'c').getAttribute('aria-checked'), 'false', 'roving alone selects nothing');

  cell(group, 'b').focus();
  assert.deepEqual(tabStops(group), ['b'], 'focus from outside syncs the roving index');
  group.dispose();
});

ui('ButtonGroup multi: independent toggles, aria-pressed, fresh arrays on every write', () => {
  const group = new ButtonGroup({items: [
    {id: 'grid', label: 'Grid'}, {id: 'cards', label: 'Cards'}, {id: 'chart', label: 'Chart'},
  ], toggle: 'multi'});
  assert.equal(group.root.getAttribute('role'), 'group');
  assert.equal(cell(group, 'grid').hasAttribute('role'), false, 'multi stays plain buttons');
  assert.equal(cell(group, 'grid').getAttribute('aria-pressed'), 'false');

  const empty = group.selected.value;
  fire(cell(group, 'grid'), 'click');
  const one = group.selected.value;
  fire(cell(group, 'chart'), 'click');
  const two = group.selected.value;
  assert.deepEqual(one, ['grid']);
  assert.deepEqual(two, ['grid', 'chart']);
  assert.notEqual(one, empty, 'a write never mutates the previous array');
  assert.notEqual(two, one);
  assert.deepEqual(empty, [], 'the first array is still empty');
  assert.equal(cell(group, 'chart').getAttribute('aria-pressed'), 'true');

  fire(cell(group, 'grid'), 'click');
  assert.deepEqual(group.selected.value, ['chart'], 'a second click toggles it back off');
  assert.equal(cell(group, 'grid').getAttribute('aria-pressed'), 'false');
  group.dispose();
});

ui('ButtonGroup none: fires onClick and selects nothing', () => {
  const fired = [];
  const group = new ButtonGroup({items: [
    {id: 'add', label: 'Add row', icon: 'plus', onClick: () => fired.push('add')},
    {id: 'copy', label: 'Copy', icon: 'copy', onClick: () => fired.push('copy')},
  ]});
  fire(cell(group, 'add'), 'click');
  fire(cell(group, 'add'), 'click');
  fire(cell(group, 'copy'), 'click');
  assert.deepEqual(fired, ['add', 'add', 'copy']);
  assert.deepEqual(group.selected.value, []);
  assert.equal(cell(group, 'add').hasAttribute('aria-pressed'), false);
  assert.equal(cell(group, 'add').classList.contains('u2-btn-group-selected'), false);
  group.dispose();
});

ui('ButtonGroup iconOnly: labels become the tooltip and the accessible name', () => {
  const group = new ButtonGroup({items: [
    {id: 'add', label: 'Add row', icon: 'plus'},
    {id: 'delete', label: 'Delete', icon: 'trash-alt', tooltip: 'Delete the selected rows'},
  ], iconOnly: true});
  assert.equal(group.root.classList.contains('u2-btn-group-icons'), true);
  assert.equal(cell(group, 'add').querySelector('.u2-btn-group-label'), null, 'no visible label');
  assert.equal(cell(group, 'add').textContent, '');
  assert.equal(cell(group, 'add').getAttribute('aria-label'), 'Add row');
  assert.equal(cell(group, 'delete').getAttribute('aria-label'), 'Delete the selected rows',
    'an explicit tooltip wins over the label');
  assert.equal(cell(group, 'add').children[0].classList.contains('fa-plus'), true);
  group.dispose();

  const labelled = new ButtonGroup({items: [{id: 'add', label: 'Add row', icon: 'plus'}]});
  assert.equal(cell(labelled, 'add').querySelector('.u2-btn-group-label').textContent, 'Add row');
  assert.equal(cell(labelled, 'add').hasAttribute('aria-label'), false, 'a visible label is the name');
  labelled.dispose();
});

ui('ButtonGroup: densities are class-driven, and setEnabled moves a cell out of the tab order', () => {
  for (const density of ['form', 'toolbar', 'ribbon']) {
    const group = new ButtonGroup({items: [{id: 'a', label: 'A'}], density});
    assert.equal(group.root.classList.contains(`u2-btn-group-${density}`), true);
    group.dispose();
  }

  const group = new ButtonGroup({items: [
    {id: 'run', label: 'Run'},
    {id: 'stop', label: 'Stop', enabled: false},
    {id: 'reset', label: 'Reset'},
  ], density: 'toolbar'});
  assert.equal(cell(group, 'stop').disabled, true);

  document.body.append(group.root);
  fire(group.root, 'keydown', {key: 'ArrowRight'});
  assert.deepEqual(tabStops(group), ['reset'], 'a disabled cell is skipped');

  group.setEnabled('stop', true);
  assert.equal(cell(group, 'stop').disabled, false);
  fire(group.root, 'keydown', {key: 'ArrowLeft'});
  assert.deepEqual(tabStops(group), ['stop'], 'enabling puts it back in the roving order');

  group.setEnabled('run', false);
  assert.equal(cell(group, 'run').disabled, true);
  group.dispose();
});

ui('ButtonGroup: dispose releases the selection effect and the click listeners', () => {
  const fired = [];
  const group = new ButtonGroup({items: [
    {id: 'a', label: 'A', onClick: () => fired.push('a')},
  ], toggle: 'multi'});
  const el = cell(group, 'a');
  group.dispose();

  fire(el, 'click');
  assert.deepEqual(fired, [], 'the click listener died with the group');
  group.selected.value = ['a'];
  assert.equal(el.getAttribute('aria-pressed'), 'false', 'the selection effect died with the group');
});

function mount(el) {
  document.body.append(el);
  return el;
}

function menuPanel() {
  return document.body.querySelector('.u2-menu');
}

function menuLabels() {
  return menuPanel().querySelectorAll('.u2-menu-label').map((el) => el.textContent);
}

function items() {
  return menuPanel().querySelectorAll('.u2-menu-item');
}

ui('dropDownButton: click opens the menu, an item closes it, and it reopens', async () => {
  const fired = [];
  const el = mount(dropDownButton('Export', (m) => m
    .item('CSV', () => fired.push('CSV'))
    .item('JSON', () => fired.push('JSON'))
    .separator()
    .group('More')
    .item('Deep', () => fired.push('Deep'))
    .endGroup()));
  assert.equal(el.dataset.u2, 'dropdown-button');
  assert.equal(el.classList.contains('u2-btn'), true);
  assert.equal(el.children[0].className, 'u2-btn-label');
  assert.equal(el.children[1].classList.contains('u2-dropdown-chevron'), true);
  assert.equal(el.children[1].classList.contains('fa-chevron-down'), true);
  assert.equal(el.getAttribute('aria-haspopup'), 'menu');
  assert.equal(el.getAttribute('aria-expanded'), 'false');

  fire(el, 'click');
  assert.equal(el.getAttribute('aria-expanded'), 'true');
  assert.deepEqual(menuLabels(), ['CSV', 'JSON', 'More']);

  fire(items()[0], 'click');
  assert.deepEqual(fired, ['CSV']);
  assert.equal(el.getAttribute('aria-expanded'), 'false');
  assert.equal(menuPanel(), null);

  fire(el, 'click');
  assert.equal(el.getAttribute('aria-expanded'), 'true');
  fire(items()[1], 'click');
  assert.deepEqual(fired, ['CSV', 'JSON'], 'the reopened menu is live too');
  await flush();
});

ui('dropDownButton: clicking the open button closes the menu instead of reopening it', async () => {
  const el = mount(dropDownButton('Export', (m) => m.item('CSV', () => {})));
  fire(el, 'click');
  const opened = menuPanel();
  fire(el, 'click');
  assert.equal(opened.isConnected, false);
  assert.equal(menuPanel(), null);
  assert.equal(el.getAttribute('aria-expanded'), 'false');

  fire(el, 'click');
  assert.notEqual(menuPanel(), null, 'a third click opens it again');
  fire(el, 'click');
  await flush();
});

ui('dropDownButton split: the main part fires its action, only the arrow opens the menu', async () => {
  const fired = [];
  const el = mount(dropDownButton('Run', (m) => m.item('Run with parameters', () => fired.push('params')),
    {icon: 'play', primary: true, split: () => fired.push('run')}));
  assert.equal(el.dataset.u2, 'dropdown-button');
  assert.equal(el.classList.contains('u2-dropdown-split'), true);
  const [main, arrow] = el.children;
  assert.equal(main.classList.contains('u2-btn-primary'), true);
  assert.equal(arrow.classList.contains('u2-btn-primary'), true, 'primary skins the whole control');
  assert.equal(main.hasAttribute('aria-haspopup'), false);
  assert.equal(arrow.getAttribute('aria-haspopup'), 'menu');
  assert.equal(arrow.getAttribute('aria-label'), 'More');
  assert.equal(main.children[0].classList.contains('fa-play'), true);
  assert.equal(arrow.children[0].classList.contains('u2-dropdown-chevron'), true);

  fire(main, 'click');
  assert.deepEqual(fired, ['run']);
  assert.equal(menuPanel(), null, 'the main part never opens the menu');
  assert.equal(arrow.getAttribute('aria-expanded'), 'false');

  fire(arrow, 'click');
  assert.equal(arrow.getAttribute('aria-expanded'), 'true');
  fire(items()[0], 'click');
  assert.deepEqual(fired, ['run', 'params']);
  assert.equal(arrow.getAttribute('aria-expanded'), 'false');
  await flush();
});

ui('dropDownButton: ArrowDown opens the menu, from the button and from either split segment', async () => {
  const plain = mount(dropDownButton('View', (m) => m.item('Grid', () => {})));
  fire(plain, 'keydown', {key: 'ArrowDown'});
  assert.equal(plain.getAttribute('aria-expanded'), 'true');
  fire(document, 'keydown', {key: 'Escape'});
  assert.equal(plain.getAttribute('aria-expanded'), 'false');

  const el = mount(dropDownButton('Run', (m) => m.item('Run with parameters', () => {}), {split: () => {}}));
  const [main, arrow] = el.children;
  fire(main, 'keydown', {key: 'ArrowDown'});
  assert.equal(arrow.getAttribute('aria-expanded'), 'true', 'ArrowDown anywhere in the control opens');
  fire(document, 'keydown', {key: 'Escape'});
  assert.equal(arrow.getAttribute('aria-expanded'), 'false');
  await flush();
});

ui('dropDownButton: the menu is rebuilt on every open, so checked state is never stale', async () => {
  const wrap = signal(false);
  const el = mount(dropDownButton('Options', (m) => m
    .item('Word wrap', () => wrap.value = !wrap.peek(), {check: wrap.peek()})));
  fire(el, 'click');
  assert.equal(items()[0].getAttribute('aria-checked'), 'false');
  fire(items()[0], 'click');
  assert.equal(wrap.value, true);

  fire(el, 'click');
  assert.equal(items()[0].getAttribute('aria-checked'), 'true', 'the second open reads the new value');
  fire(el, 'click');
  await flush();
});

ui('dropDownButton: disposing the owner closes the menu and releases the listeners', async () => {
  const fired = [];
  const owner = Control.build(() =>
    dropDownButton('Export', (m) => m.item('CSV', () => fired.push('CSV'))));
  const el = mount(owner.root);
  fire(el, 'click');
  assert.equal(el.getAttribute('aria-expanded'), 'true');

  owner.dispose();
  assert.equal(el.getAttribute('aria-expanded'), 'false', 'disposal closes the open menu');
  assert.equal(menuPanel(), null);

  fire(el, 'keydown', {key: 'ArrowDown'});
  fire(el, 'click');
  assert.equal(menuPanel(), null, 'both listeners died with the owner');
  assert.deepEqual(fired, []);
  await flush();
});
