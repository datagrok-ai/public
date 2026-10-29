import * as DG from 'datagrok-api/dg';
import {Control, signal, text, textWidth, Dates, span, timestamp, DataTable, Rows} from '@datagrok-libraries/u2';
import type {ReadonlySignal, IProperty, ObjectRenderer, CellStateLike, DataTableColumn,
  RowView} from '@datagrok-libraries/u2';
import type {DomainSource} from '../core/domain-source.js';
import {viewers} from '@datagrok-libraries/u2/dg';
import {collectionSelection, Editors} from '@datagrok-libraries/u2/dg';
import type {Bindable} from '@datagrok-libraries/u2/dg';
import {SYSTEM_COLUMNS} from './backend.js';
import {EditorEditState} from './editor-state.js';
import {DomainTable} from './facade.js';
import type {DomainAction} from './facade.js';
import {DomainForm} from './form.js';
import {pageOnScroll} from './list.js';

export interface DomainGridOptions {
  /** The grid's look, each option a value or a signal the grid follows. */
  look?: Bindable<DG.IGridSettings>;
  /** Columns the grid never shows, on top of the system ones — a child collection's FK to the
   * parent every row repeats. */
  hiddenColumns?: string[];
  /** The columns to show, in this order; every column but the system and hidden ones by default. */
  columns?: string[];
  /** A click on a cell of `column` hands its saved row to `run`. */
  open?: {column: string, run: (row: RowView) => void};
}

/** What a datetime column without a format of its own is drawn with: a `DateTime.toString()`
 * ("2026-02-03 00:00:00.000Z") is not a date. A column whose values all fall on midnight is a
 * date, the rest keep the time. */
const DATE_FORMAT = 'MMM d, yyyy';
const DATE_TIME_FORMAT = 'MMM d, yyyy HH:mm';
/** Two lines, so a unit gets its own: "Dosed Concentration (uM)". */
const HEADER_HEIGHT = 44;
/** Floors by type, so an empty frame still gets readable columns, and a cap. */
const COLUMN_WIDTH = {flag: 48, number: 64, date: 96, text: 96, name: 160, reference: 160, max: 320};
const SCROLLBAR = 12;
const HEADER_MARGIN = 8;

/** The platform grid over a source's frame with the source's writer attached (`Grid.attachEditor`, no subclass).
 * Selection and the current row are the frame's; Save and Discard are the session's, not the grid's. */
export class DomainGrid extends Control {
  readonly grid: DG.Grid;
  private _refused: {unsubscribe(): void} | undefined;
  private _rows: {unsubscribe(): void}[] = [];
  // fitted once per frame, so a column the user resized keeps its width
  private _fittedTo: unknown;

  constructor(readonly source: DomainSource, private readonly _options: DomainGridOptions = {}) {
    super();
    this.root.classList.add('u2-domain-grid');
    this.root.dataset.u2 = 'domain-grid';
    this.grid = this.build(() =>
      viewers.grid(source.df as unknown as ReadonlySignal<DG.DataFrame | undefined>, _options.look));
    this.root.append(this.grid.root);
    if (_options.look?.colHeaderHeight === undefined)
      this.grid.setOptions({colHeaderHeight: HEADER_HEIGHT});
    this.effect(() => this.grid.setOptions({allowEdit: !source.readOnly.value}));
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
      this._singleLine();
      this._fitColumns();
      this._followRows();
    });
    const open = _options.open;
    const clicked = this.grid.onCellClick.subscribe((cell) => {
      if (open === undefined || !cell.isTableCell || cell.tableColumn?.name !== open.column)
        return;
      const row = source.rows.byKey(String(cell.grid.dataFrame.get('id', cell.tableRowIndex!)));
      if (row !== undefined && !Rows.isDraft(row))
        open.run(row);
    });
    this.own(() => {
      repointed.unsubscribe();
      clicked.unsubscribe();
      this._refused?.unsubscribe();
      for (const s of this._rows)
        s.unsubscribe();
    });
  }

  // only over the frame the grid holds: it re-resolves the ref-cell renderers on its own only at its first
  // render (`grid_core` render latch), so a frame decorated before the viewer repoints keeps drawing raw ids
  private _decorate(df: DG.DataFrame): void {
    if (this.grid.dataFrame?.dart !== df.dart)
      return;
    DG.DomainObjectHandler.decorateGrid(this.grid, this.source.table, df);
    DomainGrid.sortByCaption(df);
    this._caption(df);
    this._formats(df);
  }

  private static sortByCaption(df: DG.DataFrame): void {
    for (const column of df.columns.toList()) {
      const captions = df.columns.byName(Rows.caption(column.name));
      if (captions === null)
        continue;
      const names = new Map<unknown, string>();
      for (let i = 0; i < df.rowCount; i++) {
        const name = captions.get(i);
        if (name !== null && name !== '')
          names.set(column.get(i), String(name));
      }
      const text = (id: unknown) => names.get(id) ?? String(id ?? '');
      column.valueComparer = (a, b) => text(a).localeCompare(text(b));
    }
  }

  // a wrapped text is clipped mid-line by the row height; the tooltip shows it whole
  private _singleLine(): void {
    for (const column of this._columns()) {
      if (column.column === null)
        continue;
      column.contentCellStyle.textWrap = 'ellipsis';
      column.headerCellStyle.vertAlign = 'bottom';
    }
  }

  // what a child pane caps the grid at, so a short table has no empty rules
  private _fitHeight(): void {
    const rows = Math.max(1, this.grid.dataFrame?.rowCount ?? 0);
    const look = this.grid.props;
    this.root.style.setProperty('--u2-domain-grid-height',
      `${look.colHeaderHeight + rows * look.rowHeight + SCROLLBAR}px`);
  }

  private _followRows(): void {
    for (const s of this._rows)
      s.unsubscribe();
    const df = this.grid.dataFrame;
    this._rows = !df ? [] : [df.onRowsAdded, df.onRowsRemoved].map((e) => e.subscribe(() => this._fitHeight()));
    this._fitHeight();
  }

  // the platform sizes a column by the first frame it saw, for a child often an empty one
  private _fitColumns(): void {
    const df = this.grid.dataFrame;
    if (df === null || df === undefined || this._fittedTo === df.dart)
      return;
    this._fittedTo = df.dart;
    const byName = new Map(this.source.schema.properties.map((p): [string, IProperty] => [p.name!, p]));
    const nameColumn = this.source.schema.info.nameColumn;
    const font = this.grid.props.colHeaderFont;
    for (const c of this._columns()) {
      if (c.column === null || !c.visible)
        continue;
      const floor = DomainGrid._floorOf(byName.get(c.column.name), c.column.name === nameColumn);
      const caption = DomainGrid._captionWidth(c.column.getTag(DG.TAGS.FRIENDLY_NAME) || c.name, font);
      c.width = Math.round(Math.min(COLUMN_WIDTH.max, Math.max(floor, c.getDataWidth(), caption + HEADER_MARGIN)));
    }
  }

  private static _captionWidth(caption: string, font: string): number {
    const words = caption.split(/\s+/).filter((w) => w !== '');
    let best = textWidth(words.join(' '), font);
    for (let at = 1; at < words.length; at++) {
      best = Math.min(best, Math.max(textWidth(words.slice(0, at).join(' '), font),
        textWidth(words.slice(at).join(' '), font)));
    }
    return best;
  }

  private static _floorOf(prop: IProperty | undefined, isName: boolean): number {
    const type = prop?.propertyType ?? prop?.type;
    return isName ? COLUMN_WIDTH.name : prop !== undefined && DomainTable.isReference(prop) ? COLUMN_WIDTH.reference :
      type === 'bool' ? COLUMN_WIDTH.flag : DomainDataTable.isNumeric(prop) ? COLUMN_WIDTH.number :
        type === 'datetime' ? COLUMN_WIDTH.date : COLUMN_WIDTH.text;
  }

  private _columns(): DG.GridColumn[] {
    const columns = this.grid.columns;
    return Array.from({length: columns.length}, (_, i) => columns.byIndex(i))
      .filter((c): c is DG.GridColumn => c !== null);
  }

  // the column's friendly name, not the grid column's, so the editor still matches its writable columns by name
  private _caption(df: DG.DataFrame): void {
    for (const p of this.source.schema.properties) {
      const column = df.columns.byName(p.name);
      if (column !== null && p.friendlyName !== undefined && p.friendlyName !== '')
        column.setTag(DG.TAGS.FRIENDLY_NAME, p.friendlyName);
    }
  }

  private _hideSystem(): void {
    const shown = this._options.columns;
    if (shown !== undefined) {
      this.grid.columns.setVisible(shown);
      this.grid.columns.setOrder(shown);
      return;
    }
    for (const name of [...SYSTEM_COLUMNS.map(([column]) => column), ...(this._options.hiddenColumns ?? [])]) {
      const column = this.grid.columns.byName(name);
      if (column !== null)
        column.visible = false;
    }
  }

  private _formats(df: DG.DataFrame): void {
    for (const p of this.source.schema.properties) {
      const column = df.columns.byName(p.name!);
      if (column === null || (column.getTag(DG.TAGS.FORMAT) ?? '') !== '')
        continue;
      if ((p.format ?? '') !== '')
        column.setTag(DG.TAGS.FORMAT, p.format!);
      else if ((p.propertyType ?? p.type) === 'datetime')
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
  /** Row height in pixels; the header keeps its own, at most 28. */
  rowHeight?: number;
  /** Enter on the selected row, or a double-click on it; by default the row is handed to the
   * form paired through the source, as a list's Enter hands it over. */
  onActivate?: (row: RowView, index: number) => void;
  /** Context menu actions on top of Open and the table's registry. */
  actions?: DomainAction[];
  /** Columns after the table's own — values that are not the row's, such as aggregations of its children. */
  extraColumns?: DataTableColumn<RowView>[];
  /** Per column, how a cell is drawn where the app knows better than the type (a structure as a picture). */
  cells?: Record<string, (row: RowView) => HTMLElement | string>;
}

/** `domains.dataTable` — the source's rows as the HTML {@link DataTable}: virtualized rows of
 * pooled cells, the columns derived from the table's schema the way the grid's decoration
 * derives them, the headers sorting through the source, and the source's writer painting the
 * cells (pending amber, refused red). Read-only by design (the Dart grid stays THE editor) —
 * what it is for is a wide collection beside a tree or a form, where a canvas viewer is more than the page needs. */
export class DomainDataTable extends Control {
  static readonly PICTURE_ROW_HEIGHT = 60;
  /** The table, once the source's schema is known — it is what the columns are derived from. */
  readonly table: ReadonlySignal<DataTable<RowView> | null>;

  private readonly _table = signal<DataTable<RowView> | null>(null);
  // `<semType>|<id>` → the caption a lookup answered for a ref the query did not project
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
      render: (row) => this._options.cells?.[name]?.(row) ?? this._cell(row, name, byName.get(name)),
    })).concat(this._options.extraColumns ?? []);
    const pictured = this.columns().some((name) => byName.has(name) && Editors.displays(byName.get(name)!));
    const rowHeight = this._options.rowHeight ?? (pictured ? DomainDataTable.PICTURE_ROW_HEIGHT : undefined);
    const table = this.build(() => new DataTable<RowView>({columns, rowHeight,
      items: source.rows, sort: source.sort, cellState: this._cellState(), fitColumns: true,
      contextActions: (row) => DomainTable.of(source)?.rowActions(row, this._options.actions) ?? [],
      onActivate: this._options.onActivate ??
        (() => source.activate.value = source.activate.peek() + 1)}));
    collectionSelection(source, table, {owner: this});
    pageOnScroll(this, table, source);
    this.root.replaceChildren(table.root);
    this._table.value = table;
  }

  // keyed by row key, as the edit state is: the frame index never leaves it
  private _cellState(): CellStateLike {
    const source = this.source;
    return {
      isChanged: (key, column) => source.edit.peek()?.isChanged(key, column) ?? false,
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
    if (prop !== undefined && Editors.displays(prop))
      return `${Math.round((DomainDataTable.PICTURE_ROW_HEIGHT - 4) * 1.5) + 8}px`;
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

  // never a raw uuid or a `1.7999999999999998`: the renderer's caption, the referenced row, `timestamp`, the format
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
    if (Editors.displays(prop)) {
      const height = (this._options.rowHeight ?? DomainDataTable.PICTURE_ROW_HEIGHT) - 4;
      return Editors.display(prop, value, {width: Math.round(height * 1.5), height})!;
    }
    switch (prop.propertyType ?? prop.type) {
      case 'datetime':
        // a domain `datetime` carrying a date is stamped at midnight UTC, and reading it locally
        // would move it a day back west of Greenwich
        return timestamp(value as string, undefined, {utcDates: true});
      case 'double': case 'float': case 'num':
        return DomainDataTable.number(value as number, prop.format);
      case 'int': case 'bigint':
        return prop.format ? DomainDataTable.number(value as number, prop.format) : text(value);
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
