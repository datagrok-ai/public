import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {JMP_EPOCH_OFFSET_S, JMP_INT_MISSING, JmpColumn, JmpData, JmpTable, hexColor, rowStateGlyph} from './jmp-types';
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

/** Columns hidden in JMP. The hidden-column index is taken to hold column positions; one fixture
 * supports that and none contradicts it (JMP_PLUGIN.md §4.4). */
function hiddenColumns(table: JmpTable): Set<JmpColumn> {
  return new Set(table.info.hiddenIndex.filter((i) => i < table.columns.length).map((i) => table.columns[i]));
}

/** Hides, in a table view's grid, the columns JMP hides. They stay in the table, as in JMP. */
export function hideJmpHiddenColumns(tv: DG.TableView): void {
  for (const col of tv.dataFrame.columns.toList()) {
    if (col.getTag('jmp.hidden') === 'true') {
      const gc = tv.grid.columns.byName(col.name);
      if (gc)
        gc.visible = false;
    }
  }
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
  const hidden = hiddenColumns(table);
  for (const c of dataColumns(table)) {
    if (c.kind === 'rowstate')
      addRowStateColumns(df, c);
    else {
      const col = df.columns.add(toColumn(c, df.columns.getUnusedName(c.info.name)));
      if (hidden.has(c))
        col.setTag('jmp.hidden', 'true');
    }
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
  if (c.data.type !== 'rowstate')
    return;
  const {markers, colors} = c.data;
  const marker = df.columns.add(categoryColumn(df.columns.getUnusedName(`${c.info.name} marker`), markers,
    rowStateGlyph));
  const color = df.columns.add(categoryColumn(df.columns.getUnusedName(`${c.info.name} color`), colors, hexColor));
  for (const m of new Set(markers)) {
    if (MARKERS[m])
      marker.meta.markers.assign(rowStateGlyph(m), MARKERS[m]);
  }
  color.meta.colors.setCategorical(Object.fromEntries(color.categories.map((h) => [h, h])));
  for (const col of [marker, color])
    col.setTag('jmp.rowState', c.info.name);
}

/** A string column of the distinct `codes`, each shown as `label(code)`. */
function categoryColumn(name: string, codes: Uint32Array, label: (code: number) => string): DG.Column {
  const index = new Map<number, number>();
  const categories: string[] = [];
  const indexes = new Int32Array(codes.length);
  for (let i = 0; i < codes.length; i++) {
    let k = index.get(codes[i]);
    if (k === undefined) {
      k = categories.length;
      categories.push(label(codes[i]));
      index.set(codes[i], k);
    }
    indexes[i] = k;
  }
  return DG.Column.fromIndexes(name, categories, indexes);
}

/** The column's numbers as doubles, missing as NaN. */
function asFloat64(d: JmpData): Float64Array {
  switch (d.type) {
  case 'float': return d.values;
  case 'int': return Float64Array.from(d.values, (v) => v === JMP_INT_MISSING ? NaN : v);
  case 'bigint': return Float64Array.from(d.values, (v) => Number(v));
  default: return new Float64Array(0);
  }
}

/** A copy of `values` for a Datagrok column's storage, with `map` applied and NaN made FLOAT_NULL. */
function toStorage(values: Float64Array, map?: (v: number) => number): Float64Array {
  const out = new Float64Array(values.length);
  for (let i = 0; i < values.length; i++) {
    const v = values[i];
    out[i] = v !== v ? DG.FLOAT_NULL : map ? map(v) : v;
  }
  return out;
}

function toColumn(c: JmpColumn, name: string): DG.Column {
  const d = c.data;
  let col: DG.Column;
  if (d.type === 'string')
    // fromIndexes, not fromStrings: text stays text even when it reads as numbers
    col = DG.Column.fromIndexes(name, d.categories, d.indexes.slice());
  else if (d.type === 'int')
    // JMP_INT_MISSING is DG.INT_NULL, so the array is the column's storage as it is
    col = DG.Column.fromInt32Array(name, d.values.slice());
  else if (d.type === 'bigint') {
    const fits = d.values.every((v) => v > BigInt(DG.INT_NULL) && v <= BigInt(2147483647));
    col = fits ? DG.Column.fromInt32Array(name, Int32Array.from(d.values, (v) => Number(v))) :
      DG.Column.fromBigInt64Array(name, d.values);
  } else if (c.kind === 'date' || c.kind === 'datetime') {
    // Datagrok keeps date-times as microseconds since 1970
    const micros = toStorage(asFloat64(d), (s) => (s - JMP_EPOCH_OFFSET_S) * 1e6);
    col = DG.Column.dateTime(name, micros.length);
    col.setRawData(micros, false);
    if (c.kind === 'date')
      col.meta.format = 'yyyy-MM-dd';
  } else
    // float, currency, latitude, longitude, duration
    col = DG.Column.fromFloat64Array(name, toStorage(asFloat64(d)));
  if (c.kind === 'duration')
    col.meta.units = 's';
  else if (c.kind === 'latitude')
    col.semType = DG.SEMTYPE.LATITUDE;
  else if (c.kind === 'longitude')
    col.semType = DG.SEMTYPE.LONGITUDE;
  else if (c.kind === 'currency') {
    col.meta.format = '#,##0.00';
    if (c.currency) {
      col.meta.units = c.currency;
      col.setTag('jmp.currency', c.currency);
    }
  }
  if (c.fixedDecimals !== undefined)
    col.meta.format = c.fixedDecimals === 0 ? '0' : `0.${'0'.repeat(c.fixedDecimals)}`;
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
