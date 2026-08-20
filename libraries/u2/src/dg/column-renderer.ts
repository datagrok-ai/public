/* Column names on u2 surfaces: type glyph, name, semantic-type hint — the affordances Dart's
   ColumnComboBox carries (`column_input.dart:13-17`), packaged as an `ObjectRenderer` so any u2
   control that takes one (TypeAhead rows, list rows, chips) shows a column the same way.
   Renders NAMES, not `DG.Column`s, because that is what the pickers' value is (pickers.ts). */
import * as DG from 'datagrok-api/dg';
import {ObjectRenderer} from '../core/object-renderer.js';
import {divH, span} from '../core/elements.js';
import {icon} from '../components/display/icon.js';

// the platform has no glyph table to port (its column selector shows none); these are the FA
// names the app uses for the same notions elsewhere
const GLYPHS: {[type: string]: string} = {
  [DG.COLUMN_TYPE.STRING]: 'font',
  [DG.COLUMN_TYPE.INT]: 'hashtag',
  [DG.COLUMN_TYPE.FLOAT]: 'hashtag',
  [DG.COLUMN_TYPE.BIG_INT]: 'hashtag',
  [DG.COLUMN_TYPE.QNUM]: 'hashtag',
  [DG.COLUMN_TYPE.BOOL]: 'check-square',
  [DG.COLUMN_TYPE.DATE_TIME]: 'calendar-alt',
  [DG.COLUMN_TYPE.DATA_FRAME]: 'table',
};
const UNKNOWN_GLYPH = 'question-circle';

/** Presents a column of one table by name. A name the table no longer has still renders — as
 * the bare name with the unknown glyph — so a stale value never blanks the control. */
export class ColumnRenderer implements ObjectRenderer<string> {
  private readonly _table: DG.DataFrame;

  constructor(table: DG.DataFrame) {
    this._table = table;
  }

  column(name: string): DG.Column | null {
    return this._table.columns.contains(name) ? this._table.columns.byName(name) : null;
  }

  caption(name: string): string {
    return name;
  }

  icon(name: string): HTMLElement {
    const column = this.column(name);
    return icon(GLYPHS[column?.type ?? ''] ?? UNKNOWN_GLYPH, {cls: 'u2-column-glyph'});
  }

  listItem(name: string): HTMLElement {
    const semType = this.column(name)?.semType;
    const parts = [this.icon(name), span(name, 'u2-column-name')];
    if (semType)
      parts.push(span(semType, 'u2-column-semtype'));
    return divH(parts, 'u2-column-row');
  }

  markup(name: string): HTMLElement {
    return divH([this.icon(name), span(name, 'u2-column-name')], 'u2-column-row');
  }

  tooltip(name: string): string {
    const column = this.column(name);
    if (!column)
      return name;
    return column.semType ? `${name} · ${column.type} · ${column.semType}` : `${name} · ${column.type}`;
  }
}

export function columnRenderer(table: DG.DataFrame): ColumnRenderer {
  return new ColumnRenderer(table);
}
