import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Control, component} from '../src/core/component.js';
import {signal} from '../src/core/signals.js';
import {Dialog} from '../src/components/containers/dialog.js';
import {DataFrame, Column, Stream, platform} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const {ViewHost} = await import('../src/dg/shell/view-host.js');
const grok = await import('datagrok-api/grok');
const DG = await import('datagrok-api/dg');

const made = new Set();

class FakeView {
  constructor(dart, extra = {}) {
    this.dart = dart;
    Object.assign(this, extra);
    made.add(dart);
  }
}

class FakeTableView extends FakeView {
  get dataFrame() { return this.dart.dataFrame; }
  get grid() { return this.dart.grid; }
  get dockManager() { return this.dart.dockManager; }
}

const wrap = (view) => view instanceof FakeTableView ? new FakeTableView(view.dart) : new FakeView(view.dart);
const removed = (view) => grok.events.onViewRemoved.fire({dart: view.dart});

/** dock_spawn's dock manager (dock_manager.dart, panel_dock_container.dart): its ✕ undocks and fires `onClosed`;
 * `close(el)` from code undocks and only hides the panel (`onPanelHidden`); a deselected tab is hidden and stays
 * in the layout. */
class FakeDockManager {
  constructor() {
    this.nodes = new Map();
    this.onClosed = new Stream();
    this.onPanelHidden = new Stream();
  }

  dock(el, type, ref, title, ratio) {
    const node = {el, type, ref, title, ratio};
    this.nodes.set(el, node);
    document.body.append(el);
    return node;
  }

  findNode(el) { return this.nodes.get(el); }

  close(el) {
    if (!this.nodes.delete(el))
      return;
    el.remove();
    this.onPanelHidden.fire(el);
  }

  pressClose(el) {
    this.nodes.delete(el);
    el.remove();
    this.onClosed.fire(el);
  }

  deselect(el) { this.onPanelHidden.fire(el); }
}

function hosted(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      for (const dart of made)
        removed({dart});
      made.clear();
      grok.shell.v = null;
      delete grok.shell.dockManager;
      platform.reset();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
    for (const stream of ['onViewRemoved', 'onCurrentViewChanged', 'onViewLayoutApplied'])
      assert.equal(grok.events[stream].count, 0, `nothing left on ${stream}`);
  });
}

hosted('one host per view, whichever wrapper the platform hands out; a new one once the view is gone', () => {
  const view = new FakeView({id: 'a'});
  const host = ViewHost.of(view);
  assert.ok(ViewHost.of(wrap(view)) === host, 'another wrapper of the same view');
  assert.ok(ViewHost.of(new FakeView({id: 'b'})) !== host, 'another view');
  assert.equal(host.view, view);
  removed(view);
  assert.equal(host.scope.isDisposed, true);
  const again = ViewHost.of(wrap(view));
  assert.ok(again !== host);
  assert.equal(again.scope.isDisposed, false);
});

hosted('the host goes with the view, and what was built under it; another builder around it does not own it', () => {
  const view = new FakeView({id: 'a'});
  const page = Scope.root('page');
  const host = page.build(() => ViewHost.of(view));
  const built = host.build(() => new Control());
  page.dispose();
  assert.equal(host.scope.isDisposed, false, 'the view, not the builder, ends the host');
  removed(new FakeView({id: 'b'}));
  assert.equal(built.scope.isDisposed, false, 'another view went');
  removed(view);
  assert.equal(built.scope.isDisposed, true);
  assert.equal(host.scope.isDisposed, true);
});

hosted('a created view: the content is built under the host; the platform\'s kill and the shell\'s removal end it', () => {
  let inner;
  let item;
  const status = signal('loading');
  const host = ViewHost.create({name: 'App', content: () => {
    inner = new Control();
    item = new Control();
    return {root: [inner, 'text'], ribbon: [[item]], status, path: '/apps/Pkg/App'};
  }});
  const view = host.view;
  assert.ok(view instanceof DG.View);
  assert.equal(view.name, 'App');
  assert.equal(view.path, '/apps/Pkg/App');
  assert.ok(ViewHost.of(view) === host, 'the view answers the host that made it');
  assert.ok(view.root.contains(inner.root));
  assert.equal(view.root.getAttribute('data-kill-on-close'), 'true');
  assert.equal(view.root.dataset.u2Name, 'App', 'the content root is the app\'s automation scope');
  assert.deepEqual(view.ribbonPanels, [[item.root]]);
  status.value = 'ready';
  assert.equal(view.statusBarPanels[0].textContent, 'ready');
  globalThis.grok_Widget_Kill(view.root);
  assert.equal(host.scope.isDisposed, true, 'the kill channel ends the host');
  assert.equal(inner.scope.isDisposed, true);
  assert.equal(item.scope.isDisposed, true, 'the chrome goes with it');

  const second = ViewHost.create({name: 'Other', content: () => new Control()});
  grok.events.onViewRemoved.fire(second.view);
  assert.equal(second.scope.isDisposed, true, 'the shell removing the view ends it too');
});

hosted('a created view\'s signal path is mirrored on change; the mirror dies with the host', () => {
  const path = signal('/apps/Pkg/App/a');
  const host = ViewHost.create({name: 'App', content: () => ({root: new Control(), path})});
  const view = host.view;
  assert.equal(view.path, '/apps/Pkg/App/a');
  path.value = '/apps/Pkg/App/b';
  assert.equal(view.path, '/apps/Pkg/App/b');
  globalThis.grok_Widget_Kill(view.root);
  path.value = '/apps/Pkg/App/c';
  assert.equal(view.path, '/apps/Pkg/App/b', 'the mirror died with the host');
});

hosted('a content that throws leaves nothing; one that reaches for the view before it exists is told why', () => {
  let built;
  assert.throws(() => ViewHost.create({name: 'App', content: () => {
    built = new Control();
    built.own(() => {});
    throw new Error('no data');
  }}), /no data/);
  assert.equal(built.scope.isDisposed, true);
  assert.throws(() => ViewHost.create({name: 'App', content: (host) => host.signals}),
    /view-not-ready: .*after create\(\) answers/);
});

hosted('signals: current, the grid through a layout, the table\'s row, columns and selection; owned by the host', () => {
  const frame = new DataFrame([new Column('Molecule'), new Column('Stage')], [{Molecule: 'C', Stage: 'Idea'},
    {Molecule: 'CC', Stage: 'Design'}]);
  const view = new FakeTableView({id: 'tv', dataFrame: frame, grid: {dart: {id: 'grid-1'}}});
  const host = ViewHost.of(view);
  const s = host.signals;
  assert.ok(ViewHost.of(wrap(view)).signals === s, 'made once');
  assert.equal(s.isCurrent.value, false);
  grok.shell.v = wrap(view);
  grok.events.onCurrentViewChanged.fire({});
  assert.equal(s.isCurrent.value, true);
  const grids = [];
  host.effect(() => grids.push(s.grid.value?.dart.id));
  view.dart.grid = {dart: {id: 'grid-2'}};
  grok.events.onViewLayoutApplied.fire({view: wrap(view)});
  assert.deepEqual(grids, ['grid-1', 'grid-2'], 'a layout swaps the grid');
  assert.equal(s.rowCount.value, 2);
  assert.deepEqual(s.columns.value, ['Molecule', 'Stage']);
  s.currentRowIdx.value = 1;
  assert.equal(frame.currentRowIdx, 1, 'a write moves the table\'s current row');
  removed(view);
  assert.equal(frame.liveSubscriptions(), 0, 'the view gone, the table is let go');
});

hosted('show: hidden while another view is current, closed with the view; a dialog closed by the user ends its pin', () => {
  const view = new FakeView({id: 'campaign'});
  const host = ViewHost.of(view);
  grok.shell.v = wrap(view);
  const dialog = host.show(new Dialog('Progress tracker').onOK(() => {}), {modal: true});
  assert.equal(dialog.root.hidden, false);
  grok.shell.v = new FakeView({id: 'landing'});
  grok.events.onCurrentViewChanged.fire({});
  assert.equal(dialog.root.hidden, true);
  assert.equal(dialog.isOpen.value, true, 'hidden, not closed');
  grok.shell.v = wrap(view);
  grok.events.onCurrentViewChanged.fire({});
  assert.equal(dialog.root.hidden, false);

  const tracker = new Dialog('Tracker', {keepAlive: true}).onOK(() => {});
  host.show(tracker);
  tracker.close();
  assert.equal(grok.events.onCurrentViewChanged.count, 1, 'the closed one follows nothing');
  host.show(tracker);
  host.show(tracker);
  assert.equal(grok.events.onCurrentViewChanged.count, 2, 'an open dialog is not pinned twice');

  removed(view);
  assert.equal(dialog.isOpen.value, false);
  assert.equal(dialog.scope.isDisposed, true);
  assert.equal(tracker.isOpen.value, false);
  tracker.dispose();
});

hosted('dock: gone on its ✕, on a close from code and with the view; a tab deselected keeps it; disposed, it undocks', () => {
  const shellDock = grok.shell.dockManager = new FakeDockManager();
  const view = new FakeView({id: 'app'});
  const host = ViewHost.of(view);
  const panel = () => host.build(() => new Control());

  const pressed = panel();
  host.dock(pressed, 'right', {title: 'Details', ratio: 0.3});
  const node = shellDock.findNode(pressed.root);
  assert.deepEqual([node.type, node.title, node.ratio], ['right', 'Details', 0.3]);
  assert.equal(pressed.root.getAttribute('data-kill-on-close'), 'true');
  shellDock.pressClose(pressed.root);
  assert.equal(pressed.scope.isDisposed, true, 'its ✕');

  const closed = panel();
  host.dock(closed, 'left');
  shellDock.close(closed.root);
  assert.equal(closed.scope.isDisposed, true, 'dockManager.close from code');

  const tabbed = panel();
  host.dock(tabbed, 'fill');
  shellDock.deselect(tabbed.root);
  assert.equal(tabbed.scope.isDisposed, false, 'a deselected tab is still in the layout');
  tabbed.dispose();
  assert.equal(shellDock.findNode(tabbed.root), undefined, 'disposed, its panel is closed');

  const kept = panel();
  host.dock(kept, 'down');
  removed(view);
  assert.equal(kept.scope.isDisposed, true, 'the view gone');
  assert.equal(shellDock.findNode(kept.root), undefined, 'and off the shell\'s dock');
  assert.equal(shellDock.onClosed.count + shellDock.onPanelHidden.count, 0, 'nothing left listening');
});

hosted('dock: a TableView docks into its own dock manager', () => {
  grok.shell.dockManager = new FakeDockManager();
  const own = new FakeDockManager();
  const view = new FakeTableView({id: 'tv', dataFrame: new DataFrame([new Column('a')], []), dockManager: own});
  const editor = component(() => 'Editor', {name: 'editor'});
  ViewHost.of(view).dock(editor, 'right');
  assert.ok(own.findNode(editor.root));
  assert.equal(grok.shell.dockManager.findNode(editor.root), undefined);
  own.pressClose(editor.root);
  assert.equal(editor.scope.isDisposed, true, 'a view\'s own panel has no kill: its ✕ ends it here');
  const details = component(() => 'Details', {name: 'details'});
  ViewHost.of(wrap(view)).dock(details, 'down');
  removed(view);
  assert.equal(details.scope.isDisposed, true, 'built on its own, docked: the view owns it');
});

hosted('task: what is built after an await is the view\'s; once the view goes the rest is skipped and aborted', async () => {
  const view = new FakeView({id: 'a'});
  const host = ViewHost.of(view);
  let open;
  let built;
  await host.task(async (t) => {
    await Promise.resolve();
    built = t.build(() => new Control());
  }).done;
  assert.ok(built instanceof Control);

  let late = 'not run';
  let aborted;
  const gate = new Promise((resolve) => open = resolve);
  const run = host.task(async (t) => {
    await gate;
    aborted = t.signal.aborted;
    late = t.build(() => new Control());
  });
  removed(view);
  assert.equal(built.scope.isDisposed, true, 'built after the await, owned by the view');
  open();
  await run.done;
  assert.equal(aborted, true);
  assert.equal(late, undefined, 'nothing built once the view is gone');
});
