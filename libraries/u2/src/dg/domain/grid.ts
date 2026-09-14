/* `domains.grid` — the platform grid over a source's frame with the source's writer attached
   (GOAL "Grid, subclass or hook": `Grid.attachEditor`, no subclass): in-cell edits go to the
   js-api editor the source already holds, its state paints the cells, the table's handler
   decorates the columns. Selection and the current row are the frame's, so the source sees them
   without a wire; Save and Discard are the session's, not the grid's. */
import * as DG from 'datagrok-api/dg';
import {Control} from '../../core/component.js';
import type {ReadonlySignal} from '../../core/signals.js';
import {Filters} from '../../core/filter/index.js';
import type {DomainSource} from '../../sources/domain-source.js';
import {viewers} from '../viewers/viewers.js';
import type {Bindable} from '../viewers/viewer-control.js';
import {SYSTEM_COLUMNS} from './backend.js';
import {EditorEditState} from './editor-state.js';

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
const MICROSECONDS_PER_DAY = 86400000000;

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
      const value = raw[i];
      if (value !== Filters.FLOAT_NULL && Number.isFinite(value) && value % MICROSECONDS_PER_DAY !== 0)
        return false;
    }
    return true;
  }
}
