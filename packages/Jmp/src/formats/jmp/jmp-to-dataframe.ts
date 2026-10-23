import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {JmpColumn, JmpTable, JMP_EPOCH_OFFSET_S, ROWSTATE_MARKERS} from './jmp-types';
import {JmpColumnRef, JmpScriptConversion, convertScripts} from './jsl-to-viewers';
import {isSourceScript} from './jmp-scripts';

/** The columns that come across: all but those the reader couldn't decode, which would hold nothing. */
export function dataColumns(table: JmpTable): JmpColumn[] {
  return table.columns.filter((c) => !c.decodeError);
}

/** Columns the reader could not decode. */
export function skippedColumns(table: JmpTable): {name: string, reason: string}[] {
  return table.columns
    .filter((c) => c.decodeError)
    .map((c) => ({name: c.info.name, reason: c.decodeError!}));
}

/** The saved analyses, converted to viewers where we can; the "Source" import script isn't one. */
export function savedAnalyses(table: JmpTable): JmpScriptConversion[] {
  const columns: JmpColumnRef[] = dataColumns(table).filter((c) => c.kind !== 'rowstate')
    .map((c) => ({name: c.info.name, numeric: c.kind !== 'string'}));
  return convertScripts(table.scripts.filter((s) => !isSourceScript(s)), columns);
}

/** JMP marker index to the nearest Datagrok marker; markers with no counterpart are left to the default. */
const MARKERS: Record<number, DG.MarkerCodingType> = {
  0: 'circle', 1: 'cross border', 2: 'cross x border', 3: 'square border', 4: 'diamond border',
  5: 'triangle top border', 8: 'circle border', 11: 'asterisk', 12: 'circle', 15: 'square', 16: 'diamond',
  17: 'triangle bottom border', 18: 'triangle left border', 19: 'triangle right border', 20: 'triangle top',
  21: 'triangle bottom', 22: 'triangle left', 23: 'triangle right',
};

/** Converts a parsed JMP table to a DataFrame; see JMP_PLUGIN.md §5 for the mapping. */
export function jmpTableToDataFrame(table: JmpTable, name?: string): DG.DataFrame {
  const df = DG.DataFrame.create(table.info.nrows, name);
  for (const c of dataColumns(table)) {
    if (c.kind === 'rowstate')
      addRowStateColumns(df, c);
    else
      df.columns.add(toColumn(c, df.columns.getUnusedName(c.info.name)));
  }
  df.setTag('jmp.version', table.info.jmpVersion);
  df.setTag('jmp.buildString', table.info.buildString);
  df.setTag('jmp.saveTime', table.info.savetime.toISOString());
  df.setTag('jmp.charset', table.info.charset);
  if (table.scripts.length > 0)
    df.setTag('jmp.scripts', JSON.stringify(table.scripts.map((s) => ({name: s.name, text: s.text}))));
  return df;
}

/** A row-state column becomes two: "<name> marker" with the JMP glyph, drawn with the nearest Datagrok
 * marker, and "<name> color" with the colour as #RRGGBB, colour-coded with itself. */
function addRowStateColumns(df: DG.DataFrame, c: JmpColumn): void {
  const states = c.rowStates ?? [];
  const glyph = (m: number): string => m <= 0x1f ? ROWSTATE_MARKERS[m] : String.fromCodePoint(m);
  const marker = df.columns.addNewString(df.columns.getUnusedName(`${c.info.name} marker`));
  const color = df.columns.addNewString(df.columns.getUnusedName(`${c.info.name} color`));
  marker.init((i) => glyph(states[i].marker));
  color.init((i) => states[i].color);
  for (const m of new Set(states.map((s) => s.marker))) {
    if (MARKERS[m])
      marker.meta.markers.assign(glyph(m), MARKERS[m]);
  }
  color.meta.colors.setCategorical(Object.fromEntries([...new Set(states.map((s) => s.color))].map((h) => [h, h])));
  for (const col of [marker, color])
    col.setTag('jmp.rowState', c.info.name);
}

function toColumn(c: JmpColumn, name: string): DG.Column {
  let col: DG.Column;
  switch (c.kind) {
  case 'string':
    col = DG.Column.fromStrings(name, c.values.map((v) => v === null ? '' : v as string));
    break;
  case 'int':
    col = DG.Column.fromInt32Array(name, Int32Array.from(c.values, (v) => v === null ? DG.INT_NULL : v as number));
    break;
  case 'bigint':
    col = c.values.every((v) => v === null || (v as number) > DG.INT_NULL && (v as number) <= 2147483647) ?
      DG.Column.fromInt32Array(name, Int32Array.from(c.values, (v) => v === null ? DG.INT_NULL : v as number)) :
      DG.Column.fromList(DG.COLUMN_TYPE.BIG_INT, name, c.values.map((v) => v === null ? null : BigInt(v)));
    break;
  case 'date':
  case 'datetime':
    col = DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, name,
      c.values.map((v) => v === null ? null : ((v as number) - JMP_EPOCH_OFFSET_S) * 1000));
    if (c.kind === 'date')
      col.meta.format = 'yyyy-MM-dd';
    break;
  default:
    // float, currency, latitude, longitude, duration
    col = DG.Column.fromFloat64Array(name,
      Float64Array.from(c.values, (v) => v === null ? DG.FLOAT_NULL : v as number));
  }
  if (c.kind === 'duration')
    col.meta.units = 's';
  else if (c.kind === 'latitude')
    col.semType = DG.SEMTYPE.LATITUDE;
  else if (c.kind === 'longitude')
    col.semType = DG.SEMTYPE.LONGITUDE;
  else if (c.kind === 'currency')
    col.meta.format = '#,##0.00';
  col.setTag('jmp.format', c.descriptor.map((b) => b.toString(16).padStart(2, '0')).join(' '));
  col.setTag('jmp.displayWidth', `${c.info.width}`);
  return col;
}

/** Adds the viewers converted from saved analyses to a table view. Returns the number added; a viewer the
 * platform refuses is skipped with a warning, not an error, since the table itself opened fine. */
export function addSavedAnalyses(tv: DG.TableView, analyses: JmpScriptConversion[]): number {
  let added = 0;
  for (const a of analyses) {
    for (const v of a.viewers) {
      try {
        tv.addViewer(v.type, v.options);
        added++;
      } catch (e) {
        grok.shell.warning(`JMP: couldn't add the ${v.type} for "${a.name}": ${e instanceof Error ? e.message : e}`);
      }
    }
  }
  return added;
}
