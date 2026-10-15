/* Containers, navigation, actions and display (WO-C11): the state they report (CONTRACT §5.2–5.7),
   their events (one idiom: `onX` options and `control.on('x', h)`), the Wizard's runs (P96), the
   Card's handler and selection (P83), and no page listener while a tab menu or a menu is closed
   (P24). Every test leaves the live-scope count where it found it. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {fire, flush, resetDom} from './dom-shim.js';
import {state} from 'datagrok-api/u2core';
import {signal, Scope, button} from '../src/index.js';
import {debug} from '../src/core/scope.js';
import {Wizard} from '../src/components/containers/wizard.js';
import {TabStrip} from '../src/components/containers/tabs.js';
import {Card} from '../src/components/containers/card.js';
import {Section} from '../src/components/containers/section.js';
import {Accordion} from '../src/components/containers/accordion.js';
import {Splitter} from '../src/components/containers/splitter.js';
import {Menu} from '../src/components/navigation/menu.js';
import {MenuBar} from '../src/components/navigation/menu-bar.js';
import {Toolbar} from '../src/components/navigation/toolbar.js';
import {Breadcrumbs} from '../src/components/navigation/breadcrumbs.js';
import {ButtonGroup} from '../src/components/actions/button-group.js';
import {iconButton, dropDownButton} from '../src/components/actions/buttons.js';
import {ProgressBar} from '../src/components/display/progress.js';
import {StatCard} from '../src/components/display/stat-card.js';
import {AsyncView} from '../src/components/display/async-view.js';
import {AsyncSource} from '../src/core/async-source.js';

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

function tracked(body) {
  debug.enable({orphanSeconds: 0});
  try {
    return body((type) => debug.listeners().filter((l) => l.target === 'document' && l.type === type).length);
  } finally {
    debug.disable();
  }
}

function content(text) {
  const el = document.createElement('div');
  el.textContent = text;
  return el;
}

const footer = (w, text) => w.root.querySelectorAll('.u2-wizard-buttons button').find((b) => b.textContent === text);

// --- tabs ------------------------------------------------------------------------------------------

ui('TabStrip: closing the same tab id twice notifies twice, through the option and the event (P69)', () => {
  const seen = [];
  const strip = new TabStrip({onTabClosed: (id, s) => seen.push(`option ${id} ${s === strip}`)});
  strip.on('tabClosed', (id) => seen.push(`event ${id}`));
  for (let i = 0; i < 2; i++) {
    strip.addTab({id: 'x', label: 'X', closable: true, content: content('x')});
    strip.closeTab('x');
  }
  assert.deepEqual(seen, ['option x true', 'event x', 'option x true', 'event x']);
  assert.equal('onTabClosed' in strip, false, 'no signal stands for the event');
  strip.dispose();
});

ui('TabStrip: tab headers are addressable items; the state names the current tab', async () => {
  const strip = new TabStrip({tabs: [
    {id: 'data', label: 'Data', content: content('d')},
    {id: 'style', label: 'Style', closable: true, content: content('s')},
  ]});
  document.body.append(strip.root);
  const headers = strip.root.querySelectorAll('[role="tab"]');
  assert.deepEqual(headers.map((h) => h.dataset.u2Item), ['data', 'style']);
  assert.equal(headers[1].querySelector('[data-u2-part="item-label"]').textContent, 'Style');
  assert.notEqual(headers[1].querySelector('[data-u2-part="close"]'), null);
  strip.activeTab.value = 'style';
  const s = state(strip.root);
  assert.equal(s.current, 'style');
  assert.deepEqual(s.items.list.map((i) => [i.key, i.label, i.current]), [['data', 'Data', false], ['style', 'Style', true]]);
  assert.deepEqual(s.selection, ['style']);
  assert.equal(await strip.revealItem('Data'), headers[0]);
  strip.dispose();
});

ui('TabStrip: the page is watched only while the overflow menu is open (P24)', () => tracked((listeners) => {
  const strip = new TabStrip({tabs: [{id: 'a', label: 'A', content: content('a')}]});
  document.body.append(strip.root);
  assert.equal(listeners('pointerdown') + listeners('keydown'), 0, 'closed: nothing on the document');
  strip.root.querySelector('.u2-tabs-overflow').click();
  assert.equal(listeners('pointerdown'), 1);
  assert.equal(listeners('keydown'), 1);
  fire(document.body, 'pointerdown');
  assert.equal(strip.root.querySelector('.u2-tabs-overflow-menu-open'), null, 'an outside press closes it');
  assert.equal(listeners('pointerdown') + listeners('keydown'), 0);
  strip.dispose();
}));

// --- wizard ----------------------------------------------------------------------------------------

ui('Wizard: a rejected commit keeps the step and shows the message in the error row (P96)', async () => {
  let fail = true;
  const w = new Wizard({steps: [
    {id: 'one', title: 'One', content: content('1'),
      commit: () => fail ? Promise.reject(new Error('Server said no')) : Promise.resolve(true)},
    {id: 'two', title: 'Two', content: content('2')},
  ]});
  document.body.append(w.root);
  const row = w.root.querySelector('[data-u2-part="error-row"]');
  assert.equal(row.hidden, true);
  w.next();
  assert.equal(w.root.getAttribute('aria-busy'), 'true', 'busy on the root while the commit runs');
  assert.equal(state(w.root).busy, true);
  await flush();
  assert.equal(w.currentStep.value, 'one', 'a rejection is a refusal');
  assert.equal(row.hidden, false);
  assert.equal(row.textContent, 'Server said no');
  assert.equal(w.root.getAttribute('aria-busy'), 'false');
  assert.deepEqual([state(w.root).valid, state(w.root).error], [false, 'Server said no']);
  fail = false;
  w.next();
  await flush();
  assert.equal(w.currentStep.value, 'two');
  assert.equal(row.hidden, true, 'the next run clears the message');
  w.dispose();
});

ui('Wizard: a throwing commit or finish shows its message too, and nothing moves', () => {
  const w = new Wizard({steps: [
    {id: 'one', title: 'One', content: content('1'), commit: () => {
      throw new Error('Bad input');
    }},
    {id: 'two', title: 'Two', content: content('2')},
  ]});
  w.next();
  assert.equal(w.currentStep.value, 'one');
  assert.equal(w.error.value, 'Bad input');
  w.dispose();
});

ui('Wizard: two NEXTs in one task commit once (P96)', () => {
  let commits = 0;
  const w = new Wizard({steps: [
    {id: 'one', title: 'One', content: content('1'), commit: () => {
      commits++;
      return true;
    }},
    {id: 'two', title: 'Two', content: content('2')},
    {id: 'three', title: 'Three', content: content('3')},
  ]});
  w.next();
  w.next();
  assert.equal(commits, 1);
  assert.equal(w.currentStep.value, 'three', 'the second NEXT is the next step\'s');
  w.dispose();
});

ui('Wizard: Esc and ✕ abort the running finish; its signal says so, and the run\'s answer still counts', async () => {
  const scope = new Scope();
  let ctx;
  let settle;
  const w = Scope.runWith(scope, () => new Wizard({
    steps: [{id: 'one', title: 'One', content: content('1')}],
    onFinish: (c) => {
      ctx = c;
      return new Promise((resolve, reject) => {
        settle = {resolve, reject};
        c.signal.addEventListener('abort', () => reject(Object.assign(new Error('aborted'), {name: 'AbortError'})));
      });
    },
  }));
  const dialog = w.openInDialog('Create');
  footer(w, 'FINISH').click();
  assert.equal(ctx.signal.aborted, false);
  fire(dialog.root, 'keydown', {key: 'Escape'});
  assert.equal(ctx.signal.aborted, true, 'Esc aborts');
  assert.equal(dialog.isOpen.value, true, 'and closes nothing while the run is in flight');
  await flush();
  assert.equal(w.busy.value, false, 'the aborted run settled: the footer is back');
  assert.equal(footer(w, 'FINISH').disabled, false);
  assert.equal(w.error.value, null, 'an abort shows nothing');
  assert.equal(w.completed.value, false);
  footer(w, 'FINISH').click();
  dialog.root.querySelector('.u2-dialog-close').click();
  assert.equal(ctx.signal.aborted, true, '✕ aborts too');
  settle.resolve(true);
  await flush();
  scope.dispose();
});

ui('Wizard: finish, cancel and step are events; onActivate content lives until the step\'s next visit', () => {
  const scope = new Scope();
  const events = [];
  const built = [];
  const w = Scope.runWith(scope, () => new Wizard({
    steps: [
      {id: 'one', title: 'One', content: content('1')},
      {id: 'two', title: 'Two', content: content('2'), onActivate: () => built.push(Scope.ambient)},
    ],
  }));
  w.on('step', (id) => events.push(`step ${id}`));
  w.on('finish', () => events.push('finish'));
  w.on('cancel', () => events.push('cancel'));
  const dialog = w.openInDialog('Set up');
  w.next();
  w.back();
  assert.equal(built[0].isDisposed, false, 'BACK keeps what the visit built');
  w.next();
  assert.equal(built[0].isDisposed, true, 'the next visit replaces it');
  assert.equal(built[1].isDisposed, false);
  dialog.close();
  w.openInDialog('Set up');
  w.next();
  assert.deepEqual(events, ['step two', 'step one', 'step two', 'cancel', 'finish']);
  scope.dispose();
  assert.equal(built[1].isDisposed, true, 'and the wizard\'s end ends it');
});

ui('Wizard: the steps are items with the current one, the completed ones and the state of the gate', () => {
  const name = signal('');
  const w = new Wizard({steps: [
    {id: 'source', title: 'Source', content: content('1')},
    {id: 'name', title: 'Name', content: content('2'), canProceed: () => name.value ? null : 'Enter a name'},
    {id: 'review', title: 'Review', content: content('3')},
  ]});
  const markers = w.root.querySelectorAll('.u2-wizard-step');
  assert.deepEqual(markers.map((m) => m.dataset.u2Item), ['source', 'name', 'review']);
  w.next();
  w.next();
  const s = w.getState();
  assert.equal(s.current, 'name');
  assert.deepEqual([s.valid, s.error], [false, 'Enter a name']);
  assert.deepEqual(s.items.list.map((i) => [i.key, i.current, i.completed]),
    [['source', false, true], ['name', true, false], ['review', false, false]]);
  assert.equal(s.busy, false);
  w.dispose();
});

// --- card ------------------------------------------------------------------------------------------

ui('Card: the handler gets the event and the card; select mode never unselects (P83)', () => {
  const got = [];
  const card = new Card({title: 'Aspirin', selectable: true, selectMode: 'select',
    onClick: (e, c) => got.push([e.type, c === card])});
  fire(card.root, 'click');
  fire(card.root, 'keydown', {key: 'Enter'});
  assert.deepEqual(got, [['click', true], ['keydown', true]]);
  assert.equal(card.selected.value, true, 'a second activation keeps it selected');
  assert.equal(card.root.getAttribute('aria-pressed'), 'true');
  assert.deepEqual([state(card.root).label, state(card.root).selected], ['Aspirin', true]);
  card.dispose();
});

ui('Card: a selected card shows it even when a click cannot change it; compact is a class', () => {
  const current = signal(false);
  const card = new Card({title: 'Ibuprofen', selected: current, compact: true});
  assert.equal(card.root.classList.contains('u2-card-compact'), true);
  assert.equal(card.root.getAttribute('aria-current'), null);
  current.value = true;
  assert.equal(card.root.classList.contains('u2-card-selected'), true);
  assert.equal(card.root.getAttribute('aria-current'), 'true');
  assert.equal(card.root.getAttribute('aria-pressed'), null, 'not a toggle');
  assert.equal(card.getState().selected, true);
  assert.equal(card.root.querySelector('[data-u2-part="title"]').textContent, 'Ibuprofen');
  card.dispose();
});

ui('Card: in dev, a click on a card its owner disposed while it stayed on the page warns (P83)', () => tracked(() => {
  const owner = new Scope();
  const card = Scope.runWith(owner, () => new Card({title: 'Orphan', onClick: () => {}}));
  document.body.append(card.root);
  owner.dispose();
  const warnings = [];
  const warn = console.warn;
  console.warn = (message) => warnings.push(String(message));
  try {
    fire(card.root, 'click');
  } finally {
    console.warn = warn;
  }
  assert.deepEqual(warnings, ['u2: a click on card "Orphan", which was disposed with its owner']);
}));

// --- section, accordion, splitter -----------------------------------------------------------------

ui('Section reports expanded and its title; Accordion lists its panes', () => {
  const section = new Section({title: 'Advanced', expanded: false});
  assert.deepEqual([state(section.root).label, state(section.root).expanded], ['Advanced', false]);
  section.header.click();
  assert.equal(section.getState().expanded, true);
  section.dispose();

  const accordion = new Accordion();
  accordion.addPane('General', content('g'), true);
  accordion.addPane('Advanced', content('a'));
  const s = state(accordion.root);
  assert.deepEqual(s.items.list.map((i) => [i.key, i.expanded]), [['General', true], ['Advanced', false]]);
  assert.deepEqual(accordion.root.querySelectorAll('[data-u2-item]').map((h) => h.dataset.u2Item), ['General', 'Advanced']);
  accordion.dispose();

  const splitter = new Splitter([content('l'), content('r')], {direction: 'horizontal', sizes: [0.3, 0.7]});
  assert.equal(splitter.getState().extra.sizes, '30,70');
  assert.equal(splitter.root.querySelectorAll('[data-u2-part="sash"]').length, 1);
  splitter.dispose();
});

// --- menus -----------------------------------------------------------------------------------------

ui('Menu: a context menu at a point is on the overlay stack — no page listener of its own (P24)', () => tracked((listeners) => {
  const picked = [];
  const menu = Menu.popup().item('Copy', () => picked.push('copy')).group('Export').item('CSV', () => picked.push('csv'))
    .endGroup();
  assert.equal(listeners('keydown') + listeners('pointerdown'), 0);
  menu.show({x: 40, y: 50});
  assert.equal(menu.isOpen.value, true);
  assert.equal(listeners('keydown'), 1, 'the overlay stack\'s one');
  const panel = document.body.querySelector('.u2-menu');
  assert.deepEqual(panel.querySelectorAll('[data-u2-item]').map((el) => el.dataset.u2Item), ['Copy', 'Export']);
  const s = state(panel);
  assert.equal(s.kind, 'menu');
  assert.equal(s.open, true);
  assert.deepEqual(s.items.list.map((i) => i.key), ['Copy', 'Export']);
  fire(panel, 'keydown', {key: 'Escape'});
  assert.equal(menu.isOpen.value, false, 'Escape closes it through the stack');
  assert.equal(listeners('keydown') + listeners('pointerdown'), 0);
  assert.equal(document.body.querySelector('.u2-menu-point'), null, 'the point anchor goes with it');

  menu.show({x: 10, y: 10});
  fire(document.body, 'pointerdown');
  assert.equal(menu.isOpen.value, false, 'an outside press closes it');
}));

ui('Menu: an item\'s handler runs owned by whoever added it, after the menu closed', () => {
  const owner = new Scope();
  let ran;
  const menu = Scope.runWith(owner, () => Menu.popup().item('Open', () => ran = Scope.ambient));
  const anchor = button('More', () => {});
  document.body.append(anchor);
  menu.show({anchor});
  document.body.querySelector('.u2-menu-item').click();
  assert.equal(ran, owner);
  assert.equal(menu.isOpen.value, false);
  owner.dispose();
});

ui('MenuBar and Toolbar list their items; a toolbar toggle is selected while pressed', () => {
  const bar = new MenuBar({items: [{label: 'File', build: (m) => m.item('New', () => {})},
    {label: 'Edit', build: () => {}, enabled: false}]});
  assert.deepEqual(bar.getItems().list.map((i) => [i.key, i.disabled]), [['File', false], ['Edit', true]]);
  bar.dispose();

  const on = signal(true);
  const toolbar = new Toolbar().addButton('Run', () => {}).addSeparator().addToggle('Wrap', on);
  const s = state(toolbar.root);
  assert.deepEqual(s.items.list.map((i) => [i.key, i.selected]), [['Run', undefined], ['Wrap', true]]);
  assert.equal(toolbar.root.querySelector('[data-u2-item="Wrap"]').getAttribute('aria-pressed'), 'true');
  toolbar.dispose();
});

// the shim loads no stylesheets: the cascade the » logic depends on is pinned at the source
test('toolbar.css: a hidden item outranks the button display, and the toolbar never sizes to its items', () => {
  const css = readFileSync(new URL('../css/toolbar.css', import.meta.url), 'utf8');
  assert.match(css, /\.u2-toolbar > \.u2-toolbar-hidden \{\s*display: none;/,
    '.u2-toolbar-btn (same weight, declared later) would otherwise keep folded items and the » on screen');
  assert.match(css, /\.u2-toolbar \{[^}]*flex: 1 1 0;/, 'a content-sized toolbar shrinks with every item it folds');
});

ui('Breadcrumbs: items by text, the last one current; the handler gets the event and the breadcrumbs', () => {
  const got = [];
  const crumbs = new Breadcrumbs({items: ['Home', 'Projects', 'Demo'], onClick: (i, item, e, c) =>
    got.push([i, item, e.type, c === crumbs])});
  const s = state(crumbs.root);
  assert.equal(s.current, 'Demo');
  assert.deepEqual(s.items.list.map((i) => i.key), ['Home', 'Projects', 'Demo']);
  crumbs.root.querySelector('[data-u2-item="Projects"]').click();
  assert.deepEqual(got, [[1, 'Projects', 'click', true]]);
  crumbs.dispose();
});

// --- actions ---------------------------------------------------------------------------------------

ui('button() roots carry data-u2; icon and group handlers get the event; a group reports its selection', () => {
  const scope = new Scope();
  const got = [];
  Scope.runWith(scope, () => {
    const plain = button('Save', (e) => got.push(`save ${e.type}`));
    assert.equal(plain.dataset.u2, 'button');
    const icon = iconButton('cog', (e) => got.push(`icon ${e.type}`), {tooltip: 'Settings'});
    plain.click();
    icon.click();
  });
  const group = new ButtonGroup({toggle: 'single', items: [{id: 'a', label: 'A', onClick: (e, g) =>
    got.push(`group ${e.type} ${g === group}`)}, {id: 'b', label: 'B'}]});
  group.root.querySelector('[data-u2-item="a"]').click();
  const s = state(group.root);
  assert.deepEqual([s.value, s.current, s.selection], ['a', 'a', ['a']]);
  assert.deepEqual(got, ['save click', 'icon click', 'group click true']);
  group.dispose();
  scope.dispose();
});

ui('dropDownButton: its listeners and its open menu live in a scope of its own under the builder', () => {
  const scope = new Scope();
  const root = Scope.runWith(scope, () => dropDownButton('Export', (m) => m.item('CSV', () => {})));
  document.body.append(root);
  root.click();
  assert.notEqual(document.body.querySelector('.u2-menu'), null);
  scope.dispose();
  assert.equal(document.body.querySelector('.u2-menu'), null, 'the builder\'s end closes the menu');
  fire(root, 'keydown', {key: 'ArrowDown'});
  assert.equal(document.body.querySelector('.u2-menu'), null, 'and its listeners are gone');
});

// --- display ---------------------------------------------------------------------------------------

ui('ProgressBar and StatCard report what they show', () => {
  const bar = new ProgressBar({showPercent: true, description: 'Uploading'});
  bar.value.value = 0.25;
  assert.deepEqual([state(bar.root).value, state(bar.root).text, state(bar.root).busy, state(bar.root).label],
    [0.25, '25%', false, 'Uploading']);
  bar.indeterminate = true;
  assert.equal(state(bar.root).busy, true);
  bar.dispose();

  const revenue = signal(1200);
  const card = new StatCard({label: 'Revenue', value: revenue, delta: 0.12});
  const s = state(card.root);
  assert.deepEqual([s.label, s.value, s.text, s.extra.delta], ['Revenue', 1200, '1,200', 0.12]);
  card.dispose();
});

ui('AsyncView keeps the last answer through a reload, busy and dimmed, until the next one lands (P22)', async () => {
  const calls = [];
  const source = new AsyncSource(() => new Promise((resolve, reject) => calls.push(Object.assign(resolve, {reject}))),
    {debounceMs: 0});
  const view = new AsyncView(source, (items) => content(items.join(',')));
  view.refresh();
  await flush();
  calls[0](['a', 'b']);
  await flush();
  assert.equal(view.root.textContent, 'a,b');
  view.refresh();
  assert.equal(view.root.textContent, 'a,b', 'the reload keeps the content');
  assert.equal(view.root.getAttribute('aria-busy'), 'true');
  assert.equal(view.root.classList.contains('u2-async-stale'), true);
  await flush();
  calls[1](['c']);
  await flush();
  assert.equal(view.root.textContent, 'c');
  assert.equal(view.root.getAttribute('aria-busy'), 'false');
  assert.equal(view.root.classList.contains('u2-async-stale'), false);
  assert.equal(state(view.root).error, null);
  view.refresh();
  await flush();
  calls[2].reject(new Error('boom'));
  await flush();
  assert.deepEqual([state(view.root).error, state(view.root).busy], ['boom', false]);
  view.dispose();
  source.dispose();
});

ui('Breadcrumbs: an empty path keeps the strip\'s row height', () => {
  const crumbs = new Breadcrumbs({items: []});
  assert.equal(crumbs.root.children.length, 0);
  assert.equal(crumbs.getItems().count, 0);
  const sheet = readFileSync(new URL('../css/breadcrumbs.css', import.meta.url), 'utf8');
  const at = sheet.indexOf('.u2-breadcrumbs {');
  assert.match(sheet.slice(at, sheet.indexOf('}', at)), /min-height: var\(--dg-inline-control-height\);/);
  crumbs.dispose();
});
