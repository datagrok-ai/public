/* WO-V11 — the viewer node's Properties section is the platform property grid over the live look
   (plan-v25 VP-22..VP-24), here over the WO-V10 doubles: the grid is built through the
   feature-detected global with the registered props and the viewer's frame; a look write — what a
   grid editor's commit is — fires the viewer's property event and becomes one `set-prop` after the
   microtask; what never becomes a patch (a bound prop, a same value, an object, Run mode, a viewer
   mid-repoint); the grid re-pointed at the viewer a patch re-creates; the release with the panel. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {signal} from '../src/core/signals.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {SpecEditor} from '../src/spec/editor.js';
import {dfBindings} from '../src/sources/df-bindings.js';
import {DataFrame, Property, WidgetDescriptor, platform} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const {SpecNodeRef} = await import('../src/dg/designer/node-ref.js');
const {lookGrid} = await import('../src/dg/designer/look-grid.js');
const {SpecNodeHandler, disposePanel} = await import('../src/dg/designer/handler.js');
const {REPOINTING} = await import('../src/dg/viewers/viewer-control.js');
const {registerPlatformComponents} = await import('../src/dg/viewers/registrations.js');
const {shell} = await import('datagrok-api/grok');

const UPDATE = 'grok_PropertyGrid_Update';

/** Strings, a bool, a num and a list of maps: the spec types `string`, `bool`, `double` and the
 * `object` the document never carries from a grid edit — and the platform's own `table` row, which
 * the spec knows only as the bind-only frame. */
const descriptor = () => new WidgetDescriptor('Scatter plot', [
  new Property('table', 'string', {category: 'Data'}),
  new Property('xColumnName', 'string', {category: 'Data'}),
  new Property('yColumnName', 'string', {category: 'Data'}),
  new Property('showRegressionLine', 'bool', {category: 'Regression'}),
  new Property('markerMinSize', 'num', {category: 'Misc'}),
  new Property('formulaLines', 'list', {category: 'Misc', subType: 'map'}),
]);

/** The filter group: `columnNames` is the write-once alias the platform turns into `filters`. */
const filtersDescriptor = () => new WidgetDescriptor('Filters', [
  new Property('columnNames', 'list', {subType: 'string'}),
  new Property('filters', 'list', {subType: 'map'}),
]);

/** One spec-built scatter plot over `$.orders` with the registered platform components; the global
 * is wrapped to hand back the grid's handle, where the double records every update. */
function grid(name, body, {bind = {}, props = {xColumnName: 'total'}, tag = 'u2-viewer-scatter-plot'} = {}) {
  test(name, async () => {
    const live = Scope.liveCount;
    const warn = console.warn;
    console.warn = () => {};
    const warnings = [];
    const warning = shell.warning;
    shell.warning = (message) => warnings.push(message);
    WidgetDescriptor.registry = [descriptor(), filtersDescriptor()];
    const reg = new Registry();
    registerPlatformComponents(reg);
    const scope = new Scope();
    const df = new DataFrame([{name: 'city', type: 'string'}, {name: 'total', type: 'double'}],
      [{city: 'Kyiv', total: 1240}], 'orders');
    const plot = {tag, name: 'plot', bind: {table: '$.orders', ...bind}, props: {...props}};
    const spec = {$schema: 'dg-ui/1', root: {tag: 'u2-div-v', name: 'box', children: [plot]}};
    const instance = renderSpec(spec,
      new SpecContext({data: {orders: dfBindings(signal(df), scope), reagent: 'total'}}), reg, {designTime: true});
    const editor = new SpecEditor(instance);
    const update = globalThis[UPDATE];
    let handle;
    globalThis[UPDATE] = (dart, src, shown, table) => {
      handle = dart;
      update(dart, src, shown, table);
    };
    const panel = new Scope();
    const calls = [];
    const funnel = (patches) => {
      calls.push(patches);
      if (patches.length === 1)
        editor.apply(patches[0]);
      else
        editor.applyAll(patches);
    };
    try {
      await body({
        instance, editor, plot, df, panel, calls, warnings,
        viewer: () => instance.nodes().get(plot),
        open: () => lookGrid(new SpecNodeRef(instance, plot, editor), panel, funnel),
        updates: () => handle.updates,
        patches: () => calls.flat().map((p) => [p.op, p.name, p.value]),
        panelOf: () => new SpecNodeHandler().renderProperties(new SpecNodeRef(instance, plot, editor)),
      });
    } finally {
      disposePanel();
      panel.dispose();
      instance.dispose();
      scope.dispose();
      globalThis[UPDATE] = update;
      WidgetDescriptor.registry = [];
      platform.reset();
      shell.warning = warning;
      console.warn = warn;
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

grid('the grid is built through the global: the live look, the registered props, the viewer\'s frame',
  ({open, viewer, df, updates}) => {
    const look = open();
    assert.ok(look.root.classList.contains('d4-property-grid-widget') && look.root.classList.contains('u2-designer-look'));
    assert.equal(updates().length, 1);
    const [first] = updates();
    assert.equal(first.src, viewer().dart.look, 'the look, read through grok_Viewer_Get_Look');
    assert.deepEqual(first.props, viewer().getProperties().filter((p) => p.name !== 'table').map((p) => p.dart),
      'exactly the registered spec props ∩ getProperties(), as dart handles — look is no platform property, and ' +
      'the platform\'s table row is the bind-only frame (M4)');
    assert.equal(first.table, df.dart, 'the viewer\'s own frame: the designer has no table in the shell');
  });

grid('a look write is one set-prop, after the microtask', async ({open, viewer, plot, calls, patches}) => {
  open();
  viewer().props.xColumnName = 'city';
  assert.equal(calls.length, 0, 'coalesced: nothing reaches the editor inside the event');
  await flush();
  assert.deepEqual(patches(), [['set-prop', 'xColumnName', 'city']]);
  assert.equal(calls[0][0].node, plot);
  assert.equal(plot.props.xColumnName, 'city');
});

grid('a bound prop edited in the grid is live, never a patch', async ({open, viewer, patches}) => {
  open();
  viewer().props.yColumnName = 'city';
  await flush();
  assert.deepEqual(patches(), []);
}, {bind: {yColumnName: '$.reagent'}});

grid('a same value, an empty value for a prop the document lacks, an object: no patch',
  async ({open, viewer, patches, plot}) => {
    open();
    viewer().props.xColumnName = 'total';
    viewer().props.yColumnName = '';
    viewer().props.formulaLines = [{title: 'x = 1'}];
    await flush();
    assert.deepEqual(patches(), []);
    assert.deepEqual(plot.props, {xColumnName: 'total'});
  });

grid('null is "unset": the key is removed from the document', async ({open, viewer, patches, plot}) => {
  open();
  viewer().props.xColumnName = null;
  await flush();
  assert.deepEqual(patches(), [['set-prop', 'xColumnName', undefined]]);
  assert.equal(plot.props?.xColumnName, undefined);
});

grid('two names in one tick reach the funnel as one list', async ({open, viewer, calls, plot}) => {
  open();
  const before = viewer();
  before.props.xColumnName = 'city';
  before.props.showRegressionLine = true;
  await flush();
  assert.equal(calls.length, 1);
  assert.deepEqual(calls[0].map((p) => [p.name, p.value]), [['xColumnName', 'city'], ['showRegressionLine', true]]);
  assert.deepEqual(plot.props, {xColumnName: 'city', showRegressionLine: true});
  assert.notEqual(viewer(), before, 'one applyAll, one re-created viewer');
});

grid('through the panel: two names in one tick are one applyAll, one onDidApply, one undo entry',
  async ({panelOf, viewer, editor, plot}) => {
    const shown = panelOf();
    const kids = [...shown.children];
    const at = kids.findIndex((el) => el.classList.contains('u2-designer-look'));
    assert.ok(at > 0 && kids[at - 1].textContent === 'Properties', 'the first section');
    const all = [];
    const applyAll = editor.applyAll.bind(editor);
    editor.applyAll = (patches) => {
      all.push(patches);
      applyAll(patches);
    };
    let published = -1;
    const watch = new Scope().effect(() => {
      editor.onDidApply.value;
      published++;
    });
    try {
      viewer().props.xColumnName = 'city';
      viewer().props.showRegressionLine = true;
      await flush();
      assert.equal(all.length, 1);
      assert.equal(all[0].length, 2);
      assert.equal(published, 1);
      assert.deepEqual(plot.props, {xColumnName: 'city', showRegressionLine: true});
      editor.undo();
      assert.deepEqual(plot.props, {xColumnName: 'total'}, 'one undo entry for the whole commit');
      assert.equal(editor.canUndo.peek(), false);
    } finally {
      watch.dispose();
    }
  });

grid('through the panel: a refused member is dropped with a warning, the rest applied',
  async ({panelOf, viewer, warnings, plot}) => {
    panelOf();
    viewer().props.xColumnName = 'city';
    viewer().props.markerMinSize = 'big';
    await flush();
    assert.deepEqual(warnings, ['prop "markerMinSize" expects double']);
    assert.deepEqual(plot.props, {xColumnName: 'city'});
  });

grid('through the panel: what the apply throws is owned by the funnel — shell.error, no unhandled rejection',
  async ({panelOf, viewer, editor, plot}) => {
    panelOf();
    const errors = [];
    const error = shell.error;
    shell.error = (message) => errors.push(message);
    editor.applyAll = () => {
      throw new Error('u2 editor: boom');
    };
    try {
      viewer().props.xColumnName = 'city';
      await flush();
    } finally {
      shell.error = error;
    }
    assert.deepEqual(errors, ['u2 editor: boom']);
    assert.deepEqual(plot.props, {xColumnName: 'total'});
  });

grid('Run mode: the edit changes the viewer and nothing else (R-c)', async ({open, instance, viewer, patches, plot}) => {
  const look = open();
  instance.setDesignTime(false);
  look.refresh();
  viewer().props.xColumnName = 'city';
  await flush();
  assert.deepEqual(patches(), []);
  assert.equal(viewer().props.xColumnName, 'city');
  assert.equal(plot.props.xColumnName, 'total');
});

grid('a viewer mid-repoint echoes its look back: no patch', async ({open, viewer, patches}) => {
  open();
  REPOINTING.add(viewer());
  try {
    viewer().props.xColumnName = 'city';
    await flush();
  } finally {
    REPOINTING.delete(viewer());
  }
  assert.deepEqual(patches(), []);
});

grid('a patch re-creates the viewer: refresh re-points the grid at the new look and releases the old subscription',
  async ({open, viewer, updates, patches}) => {
    const look = open();
    const old = viewer();
    assert.equal(old.onPropertyValueChanged.count, 1);
    old.props.xColumnName = 'city';
    await flush();
    const fresh = viewer();
    assert.notEqual(fresh, old);
    assert.equal(updates().length, 1, 'not re-pointed until asked');
    look.refresh();
    assert.equal(updates().length, 2);
    assert.equal(updates()[1].src, fresh.dart.look, 'the NEW viewer\'s look');
    assert.equal(old.onPropertyValueChanged.count, 0, 'nothing lives on the old viewer');
    look.refresh();
    assert.equal(updates().length, 2, 'the same build: a no-op');
    fresh.props.showRegressionLine = true;
    await flush();
    assert.deepEqual(patches(), [['set-prop', 'xColumnName', 'city'], ['set-prop', 'showRegressionLine', true]]);
  });

grid('the panel\'s disposal empties the grid, which drops its look subscriptions', async ({open, panel, updates}) => {
  open();
  panel.dispose();
  const last = updates()[updates().length - 1];
  assert.deepEqual([last.src, last.props], [{}, []]);
});

grid('without the platform global there is a hint, no grid and no throw', ({open, updates}) => {
  const update = globalThis[UPDATE];
  delete globalThis[UPDATE];
  try {
    const look = open();
    assert.equal(look.root.className, 'u2-designer-hint');
    assert.match(look.root.textContent, /property grid is not available/);
    look.refresh();
  } finally {
    globalThis[UPDATE] = update;
  }
  assert.throws(updates, 'no grid was built');
});

/* UX WO (V-2.5) — B2: the look's own value written back (opening the platform color editor does
   that) is no edit; M2: focus stays on the edited row across the re-point, and the undo/redo chord
   works from a grid editor; M4/M5: no `table` row, no `columnNames` row on the filter group; the
   rows the document carries are marked. */

grid('an echo of the look\'s own value at attach is no patch; a change away from it is one (B2)',
  async ({open, viewer, patches, plot}) => {
    viewer().dart.look.markerMinSize = 5;
    open();
    viewer().props.markerMinSize = 5;
    await flush();
    assert.deepEqual(patches(), [], 'the document lacks it, and it is still not an edit');
    viewer().props.markerMinSize = 7;
    await flush();
    assert.deepEqual(patches(), [['set-prop', 'markerMinSize', 7]]);
    assert.deepEqual(plot.props, {xColumnName: 'total', markerMinSize: 7});
    viewer().props.markerMinSize = 7;
    await flush();
    assert.equal(patches().length, 1, 'the baseline follows the re-attach');
  });

grid('the platform color popup is removed on re-attach and on the panel\'s disposal (B2)',
  async ({open, viewer, panel}) => {
    const popup = () => {
      const host = document.createElement('div');
      host.className = 'property-grid-item-editor-color-picker-host';
      document.body.append(host);
      return host;
    };
    const stale = popup();
    const look = open();
    assert.equal(stale.isConnected, false, 'a leftover from an earlier panel goes when the grid opens');
    const before = popup();
    viewer().props.xColumnName = 'city';
    await flush();
    look.refresh();
    assert.equal(before.isConnected, false, 'gone with the re-point');
    const again = popup();
    panel.dispose();
    assert.equal(again.isConnected, false, 'gone with the panel');
  });

grid('refresh puts focus back on the edited row\'s editor, and marks the rows the document carries (M2, taste)',
  async ({open, viewer, plot}) => {
    const rebuild = globalThis[UPDATE];
    // the Dart grid rebuilds its rows on every update: the focused input is a new element after
    globalThis[UPDATE] = (dart, ...rest) => {
      rebuild(dart, ...rest);
      const rows = ['xColumnName', 'yColumnName'].map((name) => {
        const row = document.createElement('tr');
        row.dataset.propName = name;
        row.append(document.createElement('input'));
        return row;
      });
      dart.root.replaceChildren(...rows);
    };
    const look = open();
    document.body.append(look.root);
    const marked = () => [...look.root.querySelectorAll('tr')]
      .filter((r) => r.classList.contains('u2-designer-look-set')).map((r) => r.dataset.propName);
    assert.deepEqual(marked(), ['xColumnName']);
    const first = look.root.querySelector('tr[data-prop-name="yColumnName"] input');
    first.focus();
    viewer().props.yColumnName = 'city';
    await flush();
    look.refresh();
    const after = look.root.querySelector('tr[data-prop-name="yColumnName"] input');
    assert.notEqual(after, first, 'rebuilt');
    assert.equal(document.activeElement, after, 'and focused again');
    assert.deepEqual(marked(), ['xColumnName', 'yColumnName']);
    assert.deepEqual(plot.props, {xColumnName: 'total', yColumnName: 'city'});
  });

grid('Ctrl+Z / Ctrl+Y from a grid editor reach the designer\'s undo and redo; Delete does not (M2)',
  async ({open, viewer, editor, plot}) => {
    const look = open();
    const input = document.createElement('input');
    look.root.append(input);
    viewer().props.xColumnName = 'city';
    await flush();
    assert.equal(plot.props.xColumnName, 'city');
    fire(input, 'keydown', {key: 'Delete'});
    assert.equal(plot.props.xColumnName, 'city');
    fire(input, 'keydown', {key: 'z', ctrlKey: true});
    assert.equal(plot.props.xColumnName, 'total', 'undone');
    fire(input, 'keydown', {key: 'y', ctrlKey: true});
    assert.equal(plot.props.xColumnName, 'city', 'redone');
    fire(input, 'keydown', {key: 'z', ctrlKey: true, shiftKey: true});
    assert.equal(plot.props.xColumnName, 'city', 'Ctrl+Shift+Z is redo: nothing left to redo');
    assert.equal(editor.canRedo.peek(), false);
  });

grid('the filter group\'s grid hides `columnNames`, the alias the platform turns into `filters` (M5)',
  ({open, updates}) => {
    open();
    assert.deepEqual(updates()[0].props.map((p) => p.name), ['filters']);
  }, {tag: 'u2-viewer-filters', props: {}});

/* Final fixer (V-2.5) — U3: the Dart textbox editor removes itself on Enter BEFORE the look is
   written (editor_textbox.dart `_notifyFinishEditing` → `finishEditing` → `hideEditor` →
   `model.setValue`), so the browser has already moved focus to the body when the property event
   fires; the refresh after the patch puts focus back on the committed row — inside the grid root,
   whose keydown listener owns the undo chord. F4: the Dart grid is killed with the panel. U6: the
   filter group's panel surfaces the tag's designer actions as buttons. */

grid('a string commit drops focus to the body: refresh puts it on the committed row, where Ctrl+Z undoes (U3)',
  async ({open, viewer, plot, editor}) => {
    const rebuild = globalThis[UPDATE];
    // a string row shows a label, not an editor, once rebuilt
    globalThis[UPDATE] = (dart, ...rest) => {
      rebuild(dart, ...rest);
      const row = document.createElement('tr');
      row.dataset.propName = 'xColumnName';
      row.append(document.createElement('div'));
      dart.root.replaceChildren(row);
    };
    const look = open();
    document.body.append(look.root);
    const textbox = document.createElement('input');
    look.root.querySelector('tr').append(textbox);
    textbox.focus();
    textbox.blur();
    textbox.remove();
    assert.equal(document.activeElement, document.body, 'what hideEditor leaves behind');
    viewer().props.xColumnName = 'city';
    await flush();
    assert.equal(plot.props.xColumnName, 'city');
    look.refresh();
    const row = look.root.querySelector('tr[data-prop-name="xColumnName"]');
    assert.equal(document.activeElement, row, 'the committed row, made focusable');
    assert.equal(row.tabIndex, -1);
    fire(row, 'keydown', {key: 'z', ctrlKey: true});
    assert.equal(plot.props.xColumnName, 'total', 'undone from the grid, not the platform');
    assert.equal(editor.canUndo.peek(), false);

    // focus held elsewhere is nobody's to steal: a write that is not a grid commit leaves it
    const other = document.createElement('input');
    document.body.append(other);
    other.focus();
    viewer().props.xColumnName = 'city';
    await flush();
    look.refresh();
    assert.equal(document.activeElement, other);
  });

grid('the panel\'s disposal kills the Dart grid through the kill-walk global, once (F4)', ({open, panel, updates}) => {
  const look = open();
  assert.equal(platform.kills.includes(look.root), false);
  panel.dispose();
  assert.equal(platform.kills.filter((el) => el === look.root).length, 1);
  const last = updates()[updates().length - 1];
  assert.deepEqual([last.src, last.props], [{}, []], 'emptied first, then killed');
});

grid('the filter group\'s panel surfaces the tag\'s designer actions as buttons, run through the designer\'s verbs (U6)',
  ({instance, editor, plot}) => {
    const ran = [];
    const verbs = () => [{name: 'Delete', run: () => ran.push('delete')},
      {name: 'Add filter for column…', run: () => ran.push('add')}, {name: 'Remove filter…', run: () => ran.push('remove')}];
    const shown = new SpecNodeHandler().renderProperties(new SpecNodeRef(instance, plot, editor, verbs));
    const row = shown.querySelector('.u2-designer-verbs');
    assert.deepEqual([...row.querySelectorAll('button')].map((b) => b.textContent),
      ['Add filter for column…', 'Remove filter…'], 'the tag\'s own verbs, never the editing ones');
    const kids = [...shown.children];
    assert.ok(kids.indexOf(row) < kids.findIndex((el) => el.tagName === 'H3' && el.textContent === 'Properties'),
      'in the Node section, above the grid');
    fire(row.querySelector('button'), 'click');
    assert.deepEqual(ran, ['add']);
    assert.equal(new SpecNodeHandler().renderProperties(new SpecNodeRef(instance, plot, editor)).querySelector('.u2-designer-verbs'),
      null, 'a ref with no verbs to offer: no row');
  }, {tag: 'u2-viewer-filters', props: {}});

grid('a tag with no designer actions shows no verb row', ({instance, editor, plot}) => {
  const verbs = () => [{name: 'Delete', run: () => {}}, {name: 'Refresh', run: () => {}}];
  assert.equal(new SpecNodeHandler().renderProperties(new SpecNodeRef(instance, plot, editor, verbs))
    .querySelector('.u2-designer-verbs'), null);
});
