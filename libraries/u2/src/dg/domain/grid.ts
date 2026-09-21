/* `domains.grid` — the platform grid over a source's frame with the source's writer attached
   (GOAL "Grid, subclass or hook": `Grid.attachEditor`, no subclass): in-cell edits go to the
   js-api editor the source already holds, its state paints the cells, the table's handler
   decorates the columns. Selection and the current row are the frame's, so the source sees them
   without a wire; Save and Discard are the session's, not the grid's. */
import * as DG from 'datagrok-api/dg';
import {Control} from '../../core/component.js';
import {signal} from '../../core/signals.js';
import type {ReadonlySignal} from '../../core/signals.js';
import {text} from '../../core/text.js';
import {Dates} from '../../core/dates.js';
import {span, timestamp} from '../../core/elements.js';
import type {IProperty} from '../../core/property-like.js';
import type {ObjectRenderer} from '../../core/object-renderer.js';
import {DataTable} from '../../components/collections/data-table.js';
import type {CellStateLike, DataTableColumn} from '../../components/collections/data-table.js';
import type {DomainSource} from '../../sources/domain-source.js';
import {Rows} from '../../sources/rows-like.js';
import type {RowView} from '../../sources/rows-like.js';
import {viewers} from '../viewers/viewers.js';
import type {Bindable} from '../viewers/viewer-control.js';
import {SYSTEM_COLUMNS} from './backend.js';
import {EditorEditState} from './editor-state.js';
import {DomainSelection} from './selection.js';
import {DomainTable} from './index.js';
import {DomainForm} from './form.js';

export interface DomainGridOptions {
  /** The grid's look, each option a value or a signal the grid follows. */
  look?: Bindable<DG.IGridSettings>;
  /** Columns the grid never shows, on top of the system ones — a child collection's FK to the
   * parent every row repeats. */
  hiddenColumns?: string[];
}

/** What a datetime column without a format of its own is drawn with: a `DateTime.toString()`
 * ("2026-02-03 00:00:00.000Z") is not a date. A column whose values all fall on midnight is a
 * date, the rest keep the time. */
const DATE_FORMAT = 'MMM d, yyyy';
const DATE_TIME_FORMAT = 'MMM d, yyyy HH:mm';

export class DomainGrid extends Control {
  readonly grid: DG.Grid;
  private _refused: {unsubscribe(): void} | undefined;

  constructor(readonly source: DomainSource, private readonly _options: DomainGridOptions = {}) {
    super();
    this.root.classList.add('u2-domain-grid');
    this.root.dataset.u2 = 'domain-grid';
    this.grid = this.runInScope(() =>
      viewers.grid(source.df as unknown as ReadonlySignal<DG.DataFrame | undefined>, _options.look));
    this.root.append(this.grid.root);
    // the editor follows the source's writer: attached over each frame the source loads, released
    // when the source lets the frame go
    this.effect(() => {
      const edit = source.edit.value;
      this._refused?.unsubscribe();
      this._refused = undefined;
      if (edit instanceof EditorEditState) {
        this._decorate(edit.df);
        this.grid.attachEditor(edit.editor);
        // the grid balloons a refused cell edit itself; the source carries it to the status bar
        this._refused = edit.editor.onRefused.subscribe((r) => source.refuse(r.message));
      } else
        this.grid.detachEditor();
    });
    // the grid columns are the viewer's, and it repoints after the effect above — and again
    // whenever the editor rebuilds the frame, so decoration runs here too
    const repointed = this.grid.onDataFrameChanged.subscribe(() => {
      const edit = source.edit.peek();
      if (edit instanceof EditorEditState)
        this._decorate(edit.df);
      this._hideSystem();
    });
    this.own(() => {
      repointed.unsubscribe();
      this._refused?.unsubscribe();
    });
  }

  /** The platform's decoration over the frame the grid HOLDS — the ref-cell renderers that draw a
   * `<schema>.<table>` cell as the target row's name are resolved against the grid's own frame,
   * and the grid re-resolves them on its own only at its first render (`grid_core` render latch):
   * a frame decorated before the viewer repoints to it keeps drawing raw ids. */
  private _decorate(df: DG.DataFrame): void {
    if (this.grid.dataFrame?.dart !== df.dart)
      return;
    DG.DomainObjectHandler.decorateGrid(this.grid, this.source.table, df);
    this._caption(df);
    this._dateFormat(df);
  }

  /** The schema's captions where `decorateGrid` found no per-table meta to take them from —
   * written as the column's friendly name, not the grid column's, so the editor still matches
   * its writable columns by name. */
  private _caption(df: DG.DataFrame): void {
    for (const p of this.source.schema.properties) {
      const column = df.columns.byName(p.name);
      if (column !== null && p.friendlyName !== undefined && p.friendlyName !== '')
        column.setTag(DG.TAGS.FRIENDLY_NAME, p.friendlyName);
    }
  }

  /** The columns a form keeps in its footer, and whatever the caller hid, are not the grid's to show. */
  private _hideSystem(): void {
    for (const name of [...SYSTEM_COLUMNS.map(([column]) => column), ...(this._options.hiddenColumns ?? [])]) {
      const column = this.grid.columns.byName(name);
      if (column !== null)
        column.visible = false;
    }
  }

  /** A datetime column with no format of its own: a date where every value it holds falls on
   * midnight, date and time otherwise. */
  private _dateFormat(df: DG.DataFrame): void {
    for (const p of this.source.schema.properties) {
      if ((p.propertyType ?? p.type) !== 'datetime' || (p.format ?? '') !== '')
        continue;
      const column = df.columns.byName(p.name!);
      if (column === null || (column.getTag(DG.TAGS.FORMAT) ?? '') !== '')
        continue;
      column.setTag(DG.TAGS.FORMAT, DomainGrid.isDateOnly(column) ? DATE_FORMAT : DATE_TIME_FORMAT);
    }
  }

  /** Whether every value the column holds falls on midnight UTC — the shape of a date column
   * over a `datetime` type, which is the only one the schema has. */
  static isDateOnly(column: DG.Column): boolean {
    const raw = column.getRawData();
    for (let i = 0; i < raw.length; i++) {
      if (!Dates.isDateOnlyMicros(raw[i]))
        return false;
    }
    return true;
  }
}

/** What `domains.dataTable` shows and hides — the grid's options in the HTML table's terms. */
export interface DomainDataTableOptions {
  /** The columns to show, in this order; the schema's own selection by default (the name column
   * first, the system and service ones out). */
  columns?: string[];
  /** Columns the table never shows, on top of those — a child collection's FK to the parent. */
  hiddenColumns?: string[];
  /** Row and header height in pixels. */
  rowHeight?: number;
  /** Enter on the selected row, or a double-click on it; by default the row is handed to the
   * form paired through the source, as a list's Enter hands it over. */
  onActivate?: (row: RowView, index: number) => void;
}

/** `domains.dataTable` — the source's rows as the HTML {@link DataTable}: virtualized rows of
 * pooled cells, the columns derived from the table's schema the way the grid's decoration
 * derives them, and the source's writer painting the cells (pending amber, refused red, a cell
 * the caller may not write muted). Read-only by design (ruling 5: the Dart grid stays THE
 * editor) — what it is for is a wide collection beside a tree or a form, where a canvas viewer
 * is more than the page needs. */
export class DomainDataTable extends Control {
  /** The table, once the source's schema is known — it is what the columns are derived from. */
  readonly table: ReadonlySignal<DataTable<RowView> | null>;

  private readonly _table = signal<DataTable<RowView> | null>(null);
  /** `<semType>|<id>` → the caption a lookup answered for a ref the query did not project. */
  private readonly _looked = new Map<string, string>();
  private _tick: {unsubscribe(): void} | undefined;
  private _plain: ObjectRenderer<RowView> | undefined;

  constructor(readonly source: DomainSource, private readonly _options: DomainDataTableOptions = {}) {
    super();
    this.table = this._table;
    this.root.classList.add('u2-domain-data-table');
    this.root.dataset.u2 = 'domain-data-table';
    this.effect(() => {
      source.state.value;
      if (this._table.peek() === null && source.schema.properties.length > 0)
        this._build();
    });
    // the writer's verdicts are not the items: a pending edit repaints the window in place
    this.effect(() => {
      const edit = source.edit.value;
      const table = this._table.value;
      this._tick?.unsubscribe();
      this._tick = edit === undefined || table === null ? undefined :
        edit.onChanged.subscribe(() => table.refresh());
    });
    this.own(() => this._tick?.unsubscribe());
  }

  /** The columns the table shows: what the caller asked for, else every property the caller may
   * see that is neither a system nor a service column, the name column first. */
  columns(): string[] {
    const hidden = new Set([...SYSTEM_COLUMNS.map(([column]) => column), ...(this._options.hiddenColumns ?? [])]);
    const access = this.source.access.peek();
    const shown = (name: string) => !hidden.has(name) && !Rows.isService(name) && access.field(name) !== 'hidden';
    if (this._options.columns !== undefined)
      return this._options.columns.filter(shown);
    const names = this.source.schema.properties.map((p) => p.name!).filter(shown);
    const name = this.source.schema.info.nameColumn;
    return name === null || !names.includes(name) ? names : [name, ...names.filter((n) => n !== name)];
  }

  private _build(): void {
    const source = this.source;
    const byName = new Map(source.schema.properties.map((p): [string, IProperty] => [p.name!, p]));
    const nameColumn = source.schema.info.nameColumn;
    const columns = this.columns().map((name): DataTableColumn<RowView> => ({
      name,
      header: byName.get(name)?.friendlyName || name,
      width: DomainDataTable.width(byName.get(name), name === nameColumn),
      ...(DomainDataTable.isNumeric(byName.get(name)) ? {align: 'right' as const} : {}),
      render: (row) => this._cell(row, name, byName.get(name)),
    }));
    const table = this.runInScope(() => new DataTable<RowView>({columns, rowHeight: this._options.rowHeight,
      items: source.rows, cellState: this._cellState(),
      onActivate: this._options.onActivate ??
        (() => source.activate.value = source.activate.peek() + 1)}));
    DomainSelection.bind(this, source, table);
    this.root.replaceChildren(table.root);
    this._table.value = table;
  }

  /** The source's writer as the table reads a cell: keyed by ROW KEY, which is what the edit
   * state is keyed by too — the frame index never leaves it. */
  private _cellState(): CellStateLike {
    const source = this.source;
    return {
      isChanged: (key, column) => source.edit.peek()?.isChanged(key, column) ?? false,
      canEdit: (key, column) => {
        const row = source.rows.byKey(key);
        return row !== undefined && source.access.peek().row(row).field(column) === 'editable';
      },
      errorOf: (key, column) => {
        const message = source.edit.peek()?.errorOf(key, column) ?? null;
        return message === null ? null : {message, kind: 'error'};
      },
    };
  }

  /** A CSS grid track per column type: the name column takes the most room, a flag or a number
   * the least. Without this every column of a wide table gets an equal share and all of them
   * truncate together. */
  static width(prop: IProperty | undefined, isName: boolean): string {
    if (isName)
      return 'minmax(120px, 1.5fr)';
    switch (prop?.propertyType ?? prop?.type) {
      case 'bool': case 'int': case 'bigint': case 'double': case 'float': case 'num':
        return 'minmax(48px, 0.5fr)';
      case 'datetime':
        return 'minmax(96px, 0.8fr)';
      default:
        // a reference reads as the target's name, not as the id behind it: it grows like the name
        // column, but keeps a smaller floor — a wide table of refs must still fit a pane
        return prop !== undefined && DomainTable.isReference(prop) ? 'minmax(96px, 1.5fr)' :
          'minmax(72px, 1fr)';
    }
  }

  /** A column whose values are read right-aligned, the way every table of numbers is read. */
  static isNumeric(prop: IProperty | undefined): boolean {
    switch (prop?.propertyType ?? prop?.type) {
      case 'int': case 'bigint': case 'double': case 'float': case 'num': return true;
      default: return false;
    }
  }

  /** What a cell shows: the name column as the table's renderer captions the row (a draft as
   * "New <singular>"), a reference as the row it points at, a datetime through the same
   * `timestamp` every other u2 surface uses (a value at midnight is a date, the rest keep the
   * time), a number through the platform's formatter where the schema declares a format — and a
   * raw uuid or a `1.7999999999999998` nowhere. */
  private _cell(row: RowView, column: string, prop: IProperty | undefined): HTMLElement | string {
    if (column === this.source.schema.info.nameColumn) {
      const renderer = DomainTable.of(this.source)?.renderer ??
        (this._plain ??= DomainTable.schemaRenderer(() => this.source.schema));
      return renderer.caption(row);
    }
    const value = row[column];
    if (prop === undefined || value === null || value === undefined || value === '')
      return text(value);
    if (DomainTable.isReference(prop))
      return DomainDataTable.caption(prop, row, column, String(value), this._looked);
    switch (prop.propertyType ?? prop.type) {
      case 'datetime':
        // a domain `datetime` carrying a date is stamped at midnight UTC, and reading it locally
        // would move it a day back west of Greenwich
        return timestamp(value as string, undefined, {utcDates: true});
      case 'double': case 'float': case 'num':
        return DomainDataTable.number(value as number, prop.format);
      default:
        return text(value);
    }
  }

  /** The referenced row's caption: the one the query projected for the column
   * (`~caption_<column>`), which is the target's display name as the server computed it — null
   * where the caller may not see the target, and that is an empty cell, not a miss. Where the
   * frame carries none — a draft, a `User`, a `Group` — the one lookup a form's readonly
   * reference does (`DomainForm.captionOf`), patched into the cell it was drawn into; never the
   * uuid, which read as the value of the column. `cache` holds what that lookup answered, so a
   * pooled cell redrawn on every scroll asks once per id instead of once per paint. */
  static caption(prop: IProperty, row: RowView, column: string, id: string,
    cache?: Map<string, string>): HTMLElement | string {
    const projected = row[Rows.caption(column)];
    if (projected !== undefined)
      return typeof projected === 'string' ? projected : '';
    const key = `${prop.semType}|${id}`;
    const known = cache?.get(key);
    if (known !== undefined)
      return known;
    const el = span('…');
    void DomainForm.captionOf(prop, id).then((caption) => {
      el.textContent = caption || id;
      if (caption)
        cache?.set(key, caption);
    }, () => el.textContent = id);
    return el;
  }

  /** A number as the schema wants it read: the platform's formatter under a declared `format`,
   * else the value with the float noise rounded off (`1.7999999999999998` is `1.8`). */
  static number(value: number, format?: string | null): string {
    if (typeof value !== 'number' || !Number.isFinite(value))
      return text(value);
    return format ? DG.format(value, format) : String(Number(value.toFixed(4)));
  }
}
