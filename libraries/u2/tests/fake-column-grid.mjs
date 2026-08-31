/* The js-api ColumnGrid surface the u2 platform backends drive (grid.ts:1454-1535) — a test-local
   stand-in for a js-api proxy, not a platform entity double, so nothing under test spreads it.
   `dfColumns` is a real frame double: a `currentRowIdx` write fires `onCurrentRowChanged`, the
   seam the single-pick backend answers picks from; the selector flavor buffers checks the way the
   Dart grid's checkColumn does (`column_grid.dart:693-728`), read back through
   `getCheckedColumnNames`. Installed through the dg stub's live binding (`_installColumnGrid`). */
import {DataFrame} from './platform-doubles.mjs';

export class FakeColumnGrid {
  static calls = [];

  static popup(table, options = {}) {
    const cg = new FakeColumnGrid(table, options);
    FakeColumnGrid.calls.push(cg);
    return cg;
  }

  /** The checkbox mode (`ColumnGrid.columnSelector`, grid.ts:1467-1475): checks over the filtered
   * columns. The constructor's `isChecked` is accepted but UNUSED by the Dart side
   * (`column_grid.dart:420-425` never reads it) — modeled faithfully, so seeding must go through
   * `setSelectedColumns`; only `checkAll` seeds at construction. */
  static columnSelector(table, options = {}) {
    const cg = new FakeColumnGrid(table, options);
    if (options.checkAll)
      cg.checkAll(true);
    FakeColumnGrid.calls.push(cg);
    return cg;
  }

  constructor(table, options) {
    this.table = table;
    this.options = options;
    this.columns = table.columns.toList().filter((c) => options.filter === undefined || options.filter(c));
    if (options.addEmpty)
      this.columns.unshift(null);
    this.dfColumns = new DataFrame([{name: 'name', type: 'string'}],
      this.columns.map((c) => ({name: c?.name ?? ''})), 'columns');
    this.dfColumns.dart.currentRowIdx = -1;
    const permissive = {length: this.columns.length, trueCount: this.columns.length,
      get: () => true, set: () => {}, setAll: () => {}};
    this.dfColumns.dart.selection = permissive;
    this.dfColumns.dart.filter = permissive;
    this.root = document.createElement('div');
    this.root.className = 'd4-root d4-column-grid';
    this.search = document.createElement('input');
    // the Dart search box swallows Esc at bubble even when empty (text_input.dart:148-151)
    this.search.addEventListener('keydown', (e) => {
      if (e.key === 'Escape') {
        this.search.value = '';
        e.stopImmediatePropagation();
      }
    });
    this.root.append(this.search);
    this.grid = {autoSize: () => {}, props: {}, root: document.createElement('div')};
    this.showSearch = false;
    this.closed = 0;
    this._checked = new Set();
  }

  get currentColumn() { return this.columns[this.dfColumns.currentRowIdx] ?? null; }
  get mouseOverColumn() { return this.columns[this.dfColumns.mouseOverRowIdx ?? -1] ?? null; }

  getRow(column) { return this.columns.indexOf(column); }

  /** The grid's own verdict (column_grid.dart:98-106): the search text against name OR
   * friendlyName, ANDed with the grid's filter. */
  passesFilter(column, columnName) {
    const search = this.search.value.toLowerCase();
    return (column != null || columnName != null)
      && ((columnName ?? column.name).toLowerCase().includes(search)
        || (column?.friendlyName?.toLowerCase()?.includes(search) ?? false))
      && (this.options.filter === undefined || column == null || this.options.filter(column));
  }

  /** Drives a pick the way the Dart grid reports one. The empty name is the addEmpty row. */
  pick(name) {
    this.dfColumns.currentRowIdx = this.columns.findIndex((c) => (c?.name ?? '') === name);
  }

  // the checked-set surface (grid.ts:1524-1535), read back in grid order
  checkAll(flag = false) {
    this._checked = new Set(flag ? this.columns.filter((c) => c != null).map((c) => c.name) : []);
  }

  setSelectedColumns(columnIds) {
    this._checked = new Set(this.columns
      .filter((c) => c != null && columnIds.includes(c.name)).map((c) => c.name));
  }

  getCheckedColumnNames() {
    return this.columns.filter((c) => c != null && this._checked.has(c.name)).map((c) => c.name);
  }

  getCheckedColumns() {
    return this.columns.filter((c) => c != null && this._checked.has(c.name));
  }

  /** Drives one user checkbox click. */
  toggle(name) {
    if (!this._checked.delete(name))
      this._checked.add(name);
  }

  close() { this.closed++; }
  detach() {}
}
