/* u2/dg ViewHost.ribbon: a package's ribbon group over a platform TableView — after a rule at the ribbon's end,
   pinned while the ribbon overflows, the small icons folding into "…", the buttons always shown — plus the crumb
   back to the app (the shell renders it: `showParentCrumb`), Save left off (`showSaveButton`), a status line; all of
   it owned by the view. The double renders the ribbon as the stand does
   (probe 2026-10-27): an outer `.d4-ribbon` with the view ribbon (name area) and a ribbon of panels, each element
   wrapped in a `.d4-ribbon-item`. The shell closes every view a test showed when it ends. */

import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {layout} from './layout.mjs';
import {Scope} from '../src/core/scope.js';
import {Control} from '../src/core/component.js';
import {signal} from '../src/core/signals.js';
import {platform} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const {ViewHost} = await import('../src/dg/shell/view-host.js');
const grok = await import('datagrok-api/grok');
const ui = await import('datagrok-api/ui');

function el(tag, cls, children = [], text) {
  const e = document.createElement(tag);
  if (cls)
    e.className = cls;
  e.append(...children);
  if (text !== undefined)
    e.append(text);
  return e;
}

class FakeTableView {
  constructor(id, name = 'KHL-3') {
    const header = el('div', 'grok-ns-entity-name d4-ribbon-name', [el('span', '', [], name)]);
    const viewRibbon = el('div', 'd4-ribbon d4-view-ribbon', [el('div', 'd4-ribbon-panel'), header]);
    const panels = el('div', 'd4-ribbon');
    const outer = el('div', 'd4-ribbon', [viewRibbon, panels]);
    const container = el('div', 'grok-view-container', [outer]);
    this.dart = {id, name, parentView: null, showParentCrumb: false, showSaveButton: true, statusBarPanels: [],
      groups: [], container, outer, viewRibbon, header, panels};
    const save = el('button', 'ui-btn ui-btn-ok ui-btn-raised', [el('i', 'grok-icon fal fa-cloud-upload')], 'Save');
    this.setRibbonPanels([[el('i', 'grok-icon fal fa-plus'), save], [el('i', 'grok-icon fal fa-arrow-to-bottom')]]);
  }

  get name() { return this.dart.name; }
  set name(x) { this.dart.name = x; }
  get parentView() { return this.dart.parentView; }
  set parentView(x) { this.dart.parentView = x; }
  get showParentCrumb() { return this.dart.showParentCrumb; }
  set showParentCrumb(x) { this.dart.showParentCrumb = x; }
  get showSaveButton() { return this.dart.showSaveButton; }
  set showSaveButton(x) { this.dart.showSaveButton = x; }

  /** What the shell renders before the name of a current table view that shows its parent crumb; a press on it makes
   * the parent current, whose ribbon the shell renders at once — taking the crumb out of this one. */
  renderParentCrumb() {
    const crumb = el('span', 'd4-ribbon-name grok-table-view-parent', [], this.parentView.name);
    crumb.addEventListener('click', () => crumb.remove());
    this.dart.viewRibbon.insertBefore(crumb, this.dart.header);
    return crumb;
  }
  get statusBarPanels() { return this.dart.statusBarPanels; }
  set statusBarPanels(x) { this.dart.statusBarPanels = x; }
  get dataFrame() { return null; }

  getRibbonPanels() { return this.dart.groups.map((group) => group.slice()); }

  setRibbonPanels(panels, clear = true) {
    const item = (e) => e.classList.contains('d4-ribbon-item') ? e : el('div', 'd4-ribbon-item', [e]);
    const groups = panels.map((group) => group.map(item));
    this.dart.groups = clear ? groups : [...this.dart.groups, ...groups];
    this.dart.panels.replaceChildren(...this.dart.groups.map((group) => el('div', 'd4-ribbon-panel', group)));
  }
}

function chromed(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      for (const view of [...grok.shell.dart.views])
        removed(view);
      grok.shell.v = null;
      ui.hints.added.length = 0;
      ui.hints.removed.length = 0;
      platform.reset();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
    for (const stream of ['onViewRemoved', 'onCurrentViewChanged'])
      assert.equal(grok.events[stream].count, 0, `nothing left on ${stream}`);
  });
}

function shown(view) {
  document.body.append(view.dart.container);
  grok.shell.dart.views.push(view);
  return view;
}

/** What the shell does when a view closes: out of `grok.shell.views` first, then `onViewRemoved`. */
function removed(view) {
  grok.shell.dart.views.splice(grok.shell.dart.views.indexOf(view), 1);
  grok.events.onViewRemoved.fire({dart: view.dart});
}

const actions = (runs = []) => [
  {icon: 'sync', tooltip: 'Refresh Campaign', run: () => runs.push('refresh')},
  {icon: 'plus', tooltip: 'Add new row', run: () => runs.push('add')},
  {icon: 'code-merge', tooltip: 'Merge another table into the campaign',
    text: 'Merge another table into the campaign...', run: () => runs.push('merge')},
];

chromed('the group goes after the platform\'s panels: a rule, the icons by their tooltips, the buttons', () => {
  const view = shown(new FakeTableView('tv'));
  const runs = [];
  const chrome = ViewHost.of(view).ribbon({actions: actions(runs), name: 'hit-design-actions',
    buttons: [{text: 'Submit...', run: () => runs.push('submit'), cls: 'hit-design-submit-button'}]});
  const groups = view.getRibbonPanels();
  assert.equal(groups.length, 3, 'the platform\'s two, then the group');
  assert.ok(groups[2][0].contains(chrome.root));
  assert.ok(chrome.root.querySelector('.u2-ribbon-rule'));
  const refresh = chrome.root.querySelector('[aria-label="Refresh Campaign"]');
  fire(refresh, 'click');
  fire(chrome.root.querySelector('.hit-design-submit-button'), 'click');
  assert.deepEqual(runs, ['refresh', 'submit']);
  const state = chrome.getState();
  assert.equal(state.kind, 'ribbon-group');
  assert.equal(state.name, 'hit-design-actions');
  assert.deepEqual(state.items.list.map((i) => i.key),
    ['Refresh Campaign', 'Add new row', 'Merge another table into the campaign', 'Submit...']);
  assert.deepEqual(state.extra, {folded: false});
  chrome.dispose();
  assert.equal(view.getRibbonPanels().length, 2, 'gone before its view, it leaves the ribbon');
});

chromed('hideSave turns the platform\'s Save off and keeps the panels as they are; gone early, it puts Save back', () => {
  const view = shown(new FakeTableView('tv'));
  const chrome = ViewHost.of(view).ribbon({hideSave: true, actions: actions()});
  assert.equal(view.showSaveButton, false);
  assert.equal(view.getRibbonPanels().length, 3, 'the platform\'s panels untouched, then the group');
  assert.ok(view.dart.panels.querySelector('.fa-plus'));
  chrome.dispose();
  assert.equal(view.showSaveButton, true);
});

chromed('an action or a button offered by a signal comes and goes', () => {
  const view = shown(new FakeTableView('tv'));
  const canEdit = signal(false);
  const chrome = ViewHost.of(view).ribbon({
    actions: [...actions(), {icon: 'users', tooltip: 'Edit campaign permissions', run: () => {}, visible: canEdit}],
    buttons: [{text: 'Progress tracker', run: () => {}, visible: canEdit}, {text: 'Submit...', run: () => {}}]});
  const keys = () => chrome.getItems().list.map((i) => i.key);
  assert.deepEqual(keys(), ['Refresh Campaign', 'Add new row', 'Merge another table into the campaign', 'Submit...']);
  assert.equal(chrome.root.querySelector('[aria-label="Edit campaign permissions"]').hidden, true);
  canEdit.value = true;
  assert.ok(keys().includes('Edit campaign permissions'));
  assert.ok(keys().includes('Progress tracker'));
  chrome.dispose();
});

chromed('folded, the icons are one "…" menu — the labels, the icons — and the buttons stay', () => {
  const view = shown(new FakeTableView('tv'));
  const runs = [];
  const chrome = ViewHost.of(view).ribbon({actions: actions(runs), fold: 'always',
    buttons: [{text: 'Submit...', run: () => {}}]});
  chrome.refit();
  assert.equal(chrome.folded.value, true);
  const more = chrome.root.querySelector('[data-u2-part="more"]');
  assert.equal(more.hidden, false);
  assert.equal(more.getAttribute('aria-label'), 'More actions');
  assert.equal(chrome.root.querySelector('.u2-ribbon-icons').hidden, true);
  assert.equal(chrome.getItems().list[0].rendered, false, 'a folded icon is not rendered');
  assert.equal(chrome.getItems().list[3].rendered, true);
  fire(more, 'click');
  const items = [...document.querySelectorAll('.u2-menu-item .u2-menu-label')].map((l) => l.textContent);
  assert.deepEqual(items, ['Refresh Campaign', 'Add new row', 'Merge another table into the campaign...']);
  fire(document.querySelectorAll('.u2-menu-item')[2], 'click');
  assert.deepEqual(runs, ['merge']);
  chrome.dispose();
});

chromed('auto: folds when the ribbon overflows and unfolds when it fits; the group is pinned at the end', () => {
  const view = shown(new FakeTableView('tv'));
  const chrome = ViewHost.of(view).ribbon({actions: actions()});
  const outer = view.dart.outer;
  layout(outer, {clientWidth: 730, scrollWidth: 1035});
  chrome.refit();
  assert.equal(chrome.folded.value, true);
  const panel = chrome.root.closest('.d4-ribbon-panel');
  assert.ok(panel.classList.contains('u2-ribbon-pinned'));
  assert.ok(panel.classList.contains('u2-ribbon-overflow'), 'an edge over the platform icons that run under it');
  layout(outer, {scrollWidth: 730});
  chrome.refit();
  assert.equal(chrome.folded.value, false);
  assert.ok(!panel.classList.contains('u2-ribbon-overflow'));
  chrome.dispose();
});

chromed('panels the platform rebuilds keep the group pinned, before the next paint', async () => {
  const view = shown(new FakeTableView('tv'));
  const chrome = ViewHost.of(view).ribbon({actions: actions()});
  chrome.refit();
  const first = chrome.root.closest('.d4-ribbon-panel');
  view.dart.panels.replaceChildren(...view.dart.groups.map((group) => el('div', 'd4-ribbon-panel', group)));
  const rebuilt = chrome.root.closest('.d4-ribbon-panel');
  assert.ok(rebuilt !== first);
  await flush();
  assert.ok(rebuilt.classList.contains('u2-ribbon-pinned'));
  chrome.dispose();
});

chromed('a hint: the icon pulses and says it; folded, "…" carries the dot and its item the detail', () => {
  const view = shown(new FakeTableView('tv'));
  const hint = signal(null);
  const chrome = ViewHost.of(view).ribbon({actions: [{icon: 'sync', tooltip: 'Refresh Campaign', run: () => {},
    hint, hintDetail: 'New changes'}, ...actions().slice(1)]});
  const dot = chrome.root.querySelector('.u2-ribbon-hint-dot');
  const refresh = chrome.root.querySelector('[data-u2-item="Refresh Campaign"]');
  const more = chrome.root.querySelector('[data-u2-part="more"]');
  hint.value = 'New changes available. Click to refresh campaign.';
  assert.ok(ui.hints.added.at(-1).el === refresh, 'the platform hint pulse on the icon');
  assert.equal(refresh.getAttribute('aria-label'), 'New changes available. Click to refresh campaign.',
    'its name says the hint as its tooltip does');
  assert.equal(more.getAttribute('aria-label'), 'More actions — New changes');
  assert.equal(refresh.classList.contains('u2-ribbon-hinted'), true);
  assert.equal(dot.hidden, false);
  layout(view.dart.outer, {clientWidth: 400, scrollWidth: 900});
  chrome.refit();
  fire(more, 'click');
  assert.equal(document.querySelector('.u2-menu-item .u2-menu-detail')?.textContent, 'New changes');
  hint.value = null;
  assert.equal(dot.hidden, true);
  assert.equal(ui.hints.removed.length, 1, 'the pulse ends with its hint');
  assert.equal(refresh.getAttribute('aria-label'), 'Refresh Campaign');
  assert.equal(more.getAttribute('aria-label'), 'More actions');
  chrome.dispose();
});

chromed('an action that fails reports to the chrome\'s scope, not as an unhandled rejection', async () => {
  const view = shown(new FakeTableView('tv'));
  const reported = [];
  const chrome = ViewHost.of(view).ribbon({actions: [{icon: 'sync', tooltip: 'Refresh Campaign',
    run: async () => { throw new Error('Server is down'); }}],
  buttons: [{text: 'Submit...', run: async () => { throw new Error('No submit function'); }}]});
  chrome.scope.onError = (e) => reported.push(e.message);
  fire(chrome.root.querySelector('[data-u2-item="Refresh Campaign"]'), 'click');
  fire(chrome.root.querySelector('[data-u2-item="Submit..."]'), 'click');
  await flush();
  assert.deepEqual(reported, ['Server is down', 'No submit function']);
  chrome.dispose();
});

chromed('crumbs: the shell renders the parent before the name; its press runs onClick, nothing else does', async () => {
  const view = shown(new FakeTableView('tv'));
  const landing = {dart: {id: 'landing'}, name: 'Hit Design'};
  const clicks = [];
  const chrome = ViewHost.of(view).ribbon({breadcrumbs: {parent: landing, onClick: () => clicks.push('back')}});
  chrome.refit();
  assert.ok(view.parentView === landing, 'the platform\'s own parent link');
  assert.equal(view.showParentCrumb, true);
  assert.equal(view.dart.header.textContent, 'KHL-3', 'the shell\'s name area untouched');
  fire(view.renderParentCrumb(), 'click');
  fire(view.dart.header, 'click');
  await flush();
  assert.deepEqual(clicks, ['back']);
  chrome.dispose();
  assert.equal(view.showParentCrumb, false, 'gone early, the crumb goes with it');
  assert.equal(view.parentView, null, 'and the parent link is the one it found');
  view.parentView = landing;
  fire(view.renderParentCrumb(), 'click');
  await flush();
  assert.deepEqual(clicks, ['back'], 'nothing heard once disposed');
});

chromed('crumbs: a scripted press runs onClick after the shell made the parent current, as a real click does', async () => {
  const view = shown(new FakeTableView('tv'));
  const landing = {dart: {id: 'landing'}, name: 'Hit Design'};
  const order = [];
  const chrome = ViewHost.of(view).ribbon({breadcrumbs: {parent: landing, onClick: () => order.push('view closed')}});
  chrome.refit();
  const crumb = view.renderParentCrumb();
  // the shell's handler is a Dart `async` closure (simple_mode.dart): `shell.v = parent` runs a microtask later
  crumb.addEventListener('click', () => queueMicrotask(() => order.push('parent current')));
  fire(crumb, 'click');
  await flush();
  assert.deepEqual(order, ['parent current', 'view closed']);
  chrome.dispose();
});

chromed('gone early, the chrome puts back the Save, the parent and the crumb the view had', () => {
  const view = shown(new FakeTableView('tv'));
  const home = {dart: {id: 'home'}, name: 'Home'};
  Object.assign(view.dart, {parentView: home, showParentCrumb: true, showSaveButton: false});
  const chrome = ViewHost.of(view).ribbon({hideSave: true, breadcrumbs: {parent: {dart: {id: 'landing'}, name: 'Hit'}}});
  assert.equal(view.parentView.name, 'Hit');
  chrome.dispose();
  assert.ok(view.parentView === home);
  assert.equal(view.showParentCrumb, true);
  assert.equal(view.showSaveButton, false);
});

chromed('a status line, live; owned by the view and by whoever built it', () => {
  const view = shown(new FakeTableView('tv'));
  const status = signal('View only — changes are not saved');
  const session = new Control();
  const chrome = session.build(() => ViewHost.of(view).ribbon({status}));
  assert.equal(view.statusBarPanels.at(-1).textContent, 'View only — changes are not saved');
  status.value = 'Saving…';
  assert.equal(view.statusBarPanels.at(-1).textContent, 'Saving…');
  session.dispose();
  assert.equal(chrome.scope.isDisposed, true);
  assert.equal(view.statusBarPanels.length, 0);

  const second = ViewHost.of(view).ribbon({actions: actions()});
  removed(view);
  assert.equal(second.scope.isDisposed, true, 'the view gone, the chrome goes');
});

chromed('a session the view\'s removal disposes first takes the chrome along without writing to the closed view', () => {
  const view = shown(new FakeTableView('tv'));
  const parent = shown(new FakeTableView('landing', 'Hit Design'));
  const session = ViewHost.of(view).build(() => new Control());
  session.build(() => ViewHost.of(view).ribbon({hideSave: true, actions: actions(), status: signal('Ready'),
    breadcrumbs: {parent}}));
  const groups = view.dart.groups;
  const statusBar = view.statusBarPanels;
  removed(view);
  assert.equal(session.scope.isDisposed, true);
  assert.ok(view.dart.groups === groups, 'the closed view\'s ribbon is not rebuilt');
  assert.ok(view.statusBarPanels === statusBar);
  assert.equal(view.showSaveButton, false, 'nor its Save or its crumb put back');
  assert.equal(view.showParentCrumb, true);
});

