// What solutions reach for that u2's test doubles leave out; the TableView renders as table-view-chrome.test.js's.
// The page's `datagrok-api/dg` (lib/hooks.mjs) takes ViewBase, View and TableView from here.
import * as grok from 'datagrok-api/grok';
import {DataFrame, Stream, View as ViewDouble} from '../../tests/platform-doubles.mjs';

const tableViews = new WeakSet();
const isView = (x) => x instanceof ViewDouble || tableViews.has(x);

/** `DG.ViewBase`, `DG.View` and `DG.TableView` as the page has them: `instanceof` holds for the doubles' views and
 * for the TableViews {@link tableView} opens, as the platform's class hierarchy does. */
export class ViewBase {
  static [Symbol.hasInstance](x) { return isView(x); }
}

export class View extends ViewDouble {
  static [Symbol.hasInstance](x) { return isView(x); }
}

export class TableView extends View {
  static [Symbol.hasInstance](x) { return tableViews.has(x); }
}

/** A `DG.Row` as it was before js-api 2026-10: reading any name that is not a column, `then` included, throws. */
export function rowProxy(table, idx) {
  return new Proxy({table, idx}, {
    get(t, key) {
      if (key === 'table' || key === 'idx' || typeof key === 'symbol')
        return t[key];
      if (!table.columns.contains(key))
        throw new Error(`No such column: ${String(key)}`);
      return table.get(key, idx);
    },
  });
}

/** The js-api DataFrame events the doubles leave out (`data-frame.ts` `on…`): all present, firing where the page
 * can tell — a rename of the table, columns added or removed; the rest stay quiet. */
const MORE_EVENTS = ['onMouseOverRowChanged', 'onCurrentColChanged', 'onMouseOverColChanged', 'onCurrentCellChanged',
  'onMouseOverRowGroupChanged', 'onNameChanged', 'onMetadataChanged', 'onColumnSelectionChanged', 'onColumnsAdded',
  'onColumnsRemoved', 'onRowsFiltered', 'onRowsFiltering', 'onSemanticTypeDetecting', 'onSemanticTypeDetected'];
const DATA_EVENTS = ['onValuesChanged', 'onRowsAdded', 'onRowsRemoved', 'onColumnsChanged'];

const frames = DataFrame.prototype;
for (const event of MORE_EVENTS.filter((e) => !(e in frames)))
  Object.defineProperty(frames, event, {get() { return this.dart.events[event] ??= new Stream(); }});
Object.defineProperty(frames, 'onDataChanged', {get() {
  return {subscribe: (next) => {
    const subs = DATA_EVENTS.map((event) => this[event].subscribe(next));
    return {unsubscribe() {
      for (const sub of subs)
        sub.unsubscribe();
    }};
  }};
}});
const named = Object.getOwnPropertyDescriptor(frames, 'name');
Object.defineProperty(frames, 'name', {configurable: true, get: named.get, set(x) {
  named.set.call(this, x);
  this.onNameChanged.fire(x);
}});

/** `rows.addNew(values, notify)` answers a {@link rowProxy}. */
Object.defineProperties(frames, {
  col: {value(name) { return this.columns.byName(name); }},
  rows: {get() {
    return this.dart.rowList ??= {addNew: (values = null, notify = true) => {
      const row = {};
      for (const [i, column] of this.columns.toList().entries())
        row[column.name] = values?.[i] ?? null;
      this.dart.rows.push(row);
      if (notify)
        this.onRowsAdded.fire();
      return rowProxy(this, this.rowCount - 1);
    }};
  }},
});

const columnLists = Object.getPrototypeOf(new DataFrame().columns);
const {add, remove} = columnLists;
Object.assign(columnLists, {
  bySemType(semType) { return this.toList().find((c) => c.semType === semType) ?? null; },
  bySemTypeAll(semType) { return this.toList().filter((c) => c.semType === semType); },
  add(column) {
    const added = add.call(this, column);
    this.dart.frame.onColumnsAdded.fire({args: {columns: [added]}});
    return added;
  },
  remove(name) {
    const gone = this.byName(name);
    remove.call(this, name);
    if (gone !== null)
      this.dart.frame.onColumnsRemoved.fire({args: {columns: [gone]}});
  },
});

/** A table: the doubles' frame, which carries everything above. */
export const frame = (columns, rows = [], name = 'table') => new DataFrame(columns, rows, name);

/** The session a user works in: tables open (`demog`, one with no name, `refs`), and another view current. */
export function busySession() {
  grok.shell.dart.tables = [];
  for (const name of ['demog', null, 'refs']) {
    const table = grok.shell.addTable(frame([{name: 'smiles', type: 'string', semType: 'Molecule'}],
      [{smiles: 'CCO'}], name ?? 'unnamed'));
    if (name === null)
      table.name = undefined;
  }
  const browse = ViewDouble.fromRoot(document.createElement('div'));
  browse.name = 'Browse';
  grok.shell.dart.views = [browse];
  grok.shell.v = browse;
  grok.events.onCurrentViewChanged.fire(browse);
}

function el(tag, cls, children = [], text) {
  const e = document.createElement(tag);
  e.className = cls;
  e.append(...children);
  if (text !== undefined)
    e.append(text);
  return e;
}

/** A platform TableView over `table`, shown and current; `close()` is the shell's close. */
export function tableView(table, name = 'KHL-3') {
  const header = el('div', 'grok-ns-entity-name d4-ribbon-name', [el('span', '', [], name)]);
  const viewRibbon = el('div', 'd4-ribbon d4-view-ribbon', [el('div', 'd4-ribbon-panel'), header]);
  const panels = el('div', 'd4-ribbon');
  const outer = el('div', 'd4-ribbon', [viewRibbon, panels]);
  const container = el('div', 'grok-view-container', [outer]);
  const dart = {id: name, name, parentView: null, showParentCrumb: false, showSaveButton: true, statusBarPanels: [],
    groups: [], container, panels, table};
  const view = {
    dart,
    get name() { return dart.name; },
    set name(x) { dart.name = x; },
    get dataFrame() { return dart.table; },
    get parentView() { return dart.parentView; },
    set parentView(x) { dart.parentView = x; },
    get showParentCrumb() { return dart.showParentCrumb; },
    set showParentCrumb(x) { dart.showParentCrumb = x; },
    get showSaveButton() { return dart.showSaveButton; },
    set showSaveButton(x) { dart.showSaveButton = x; },
    get statusBarPanels() { return dart.statusBarPanels; },
    set statusBarPanels(x) { dart.statusBarPanels = x; },
    get root() { return container; },
    getRibbonPanels() { return dart.groups.map((group) => group.slice()); },
    setRibbonPanels(groups, clear = true) {
      const item = (e) => e.classList.contains('d4-ribbon-item') ? e : el('div', 'd4-ribbon-item', [e]);
      dart.groups = [...clear ? [] : dart.groups, ...groups.map((group) => group.map(item))];
      panels.replaceChildren(...dart.groups.map((group) => el('div', 'd4-ribbon-panel', group)));
    },
    close() {
      const at = grok.shell.dart.views.indexOf(view);
      if (at >= 0)
        grok.shell.dart.views.splice(at, 1);
      if (grok.shell.v === view)
        grok.shell.v = null;
      container.remove();
      grok.events.onViewRemoved.fire({dart});
    },
  };
  tableViews.add(view);
  view.setRibbonPanels([[el('i', 'grok-icon fal fa-filter')], [el('i', 'grok-icon fal fa-arrow-to-bottom')]]);
  document.body.append(container);
  grok.shell.dart.views.push(view);
  grok.events.onViewAdded.fire(view);
  grok.shell.v = view;
  grok.events.onCurrentViewChanged.fire(view);
  return view;
}
