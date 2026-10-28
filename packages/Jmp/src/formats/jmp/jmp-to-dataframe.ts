import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {
  JMP_EPOCH_OFFSET_S, JMP_INT_MISSING, JmpColumn, JmpData, JmpTable, hexColor, jmpValue, rowStateGlyph,
} from './jmp-types';
import {JmpValueLabels} from './jmp-properties';
import {JslNode, isCall, jslKey} from './jsl';
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

/** Columns hidden in JMP: the hidden-column index holds column positions (JMP_PLUGIN.md §4.4). */
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
      const name = df.columns.getUnusedName(c.info.name);
      const labels = c.properties?.valueLabels;
      const col = df.columns.add(labels ? labelledColumn(c, name, labels, table.info.nrows) : toColumn(c, name));
      setColumnTags(col, c);
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
  return col;
}

/** A column with value labels shows the labels, as JMP does: a string column of each value's label, or of
 * the value itself where it has none. The values and their labels are kept in the `jmp.valueLabels` tag. */
function labelledColumn(c: JmpColumn, name: string, labels: JmpValueLabels, nrows: number): DG.Column {
  const label = new Map<number | string, string>(labels.keys.map((k, i) => [k, labels.labels[i]]));
  const index = new Map<string, number>();
  const categories: string[] = [];
  const indexes = new Int32Array(nrows);
  for (let i = 0; i < nrows; i++) {
    const v = jmpValue(c, i);
    const text = v === null ? '' : label.get(v) ?? String(v);
    let k = index.get(text);
    if (k === undefined) {
      k = categories.length;
      categories.push(text);
      index.set(text, k);
    }
    indexes[i] = k;
  }
  const col = DG.Column.fromIndexes(name, categories, indexes);
  col.setTag('jmp.valueLabels', JSON.stringify(labels.keys.map((k, i) => [k, labels.labels[i]])));
  return col;
}

/** Tags every imported column gets, and what its JMP column properties map to (JMP_PLUGIN.md §5). */
function setColumnTags(col: DG.Column, c: JmpColumn): void {
  col.setTag('jmp.format', c.descriptor.map((b) => b.toString(16).padStart(2, '0')).join(' '));
  col.setTag('jmp.displayWidth', `${c.info.width}`);
  const props = c.properties?.properties ?? [];
  if (props.length === 0)
    return;
  col.setTag('jmp.properties', JSON.stringify(props.map((p) => ({name: p.name, text: p.text}))));
  const labels = c.properties?.valueLabels;
  const shown = (v: string | number): string => {
    const i = labels ? labels.keys.indexOf(v) : -1;
    return i >= 0 ? labels!.labels[i] : String(v);
  };
  for (const p of props) {
    const node = p.node;
    if (!node)
      continue;
    const text = node.args[0]?.type === 'string' ? node.args[0].value : undefined;
    if (isCall(node, 'Notes') && text !== undefined)
      col.setTag(DG.TAGS.DESCRIPTION, text);
    else if (isCall(node, 'Units') && text !== undefined)
      col.meta.units = text;
    else if (isCall(node, 'Formula'))
      col.setTag('jmp.formula', p.text.slice(p.text.indexOf('(') + 1, p.text.lastIndexOf(')')).trim());
    else if (isCall(node, 'Spec Limits')) {
      const {limits, show} = specLimits(node);
      if (Object.keys(limits).length > 0) {
        col.setTag('jmp.specLimits', JSON.stringify(limits));
        addSpecLimitLines(col, limits, show);
      }
    } else if (col.type === DG.TYPE.STRING) {
      const order = categoryOrder(node);
      if (order)
        col.setCategoryOrder(order.map(shown));
    }
  }
}

/** The order of `Value Order({Custom Order({...})})`, or of `"Value Ordering"({...})`, which JMP writes
 * alongside it. */
function categoryOrder(node: Extract<JslNode, {type: 'call'}>): (string | number)[] | undefined {
  let list: JslNode | undefined;
  const key = jslKey(node.name);
  if (key === 'valueorder') {
    const items = node.args[0]?.type === 'list' ? node.args[0].items : node.args;
    list = items.find((a): a is Extract<JslNode, {type: 'call'}> => isCall(a, 'Custom Order'))?.args[0];
  } else if (key === 'valueordering')
    list = node.args[0];
  if (list?.type !== 'list')
    return undefined;
  const values = list.items.map((a) => a.type === 'string' || a.type === 'number' ? a.value : undefined);
  return values.every((v) => v !== undefined) ? values as (string | number)[] : undefined;
}

type SpecLimits = {lsl?: number, usl?: number, target?: number};

/** `Spec Limits({LSL(1), USL(5), Target(3), Show Limits(0)})` as {lsl, usl, target}, missing ones left out,
 * and whether JMP draws them on its graphs. */
function specLimits(node: Extract<JslNode, {type: 'call'}>): {limits: SpecLimits, show: boolean} {
  const items = node.args[0]?.type === 'list' ? node.args[0].items : node.args;
  const arg = (name: string) => items.find((a): a is Extract<JslNode, {type: 'call'}> => isCall(a, name))?.args[0];
  // A negative limit is a minus applied to a number, not a number literal.
  const num = (v?: JslNode): number | undefined => v?.type === 'number' ? v.value :
    v?.type === 'unary' && v.op === '-' && v.operand.type === 'number' ? -v.operand.value : undefined;
  const limits: SpecLimits = {};
  for (const [key, name] of [['lsl', 'LSL'], ['usl', 'USL'], ['target', 'Target']] as const) {
    const v = num(arg(name));
    if (v !== undefined && !Number.isNaN(v))
      limits[key] = v;
  }
  const show = arg('Show Limits');
  return {limits, show: show?.type === 'number' && show.value === 1};
}

/** Spec limits as table formula lines, so every chart with the column on its value axis draws them; hidden,
 * but still listed, when JMP doesn't show them. */
function addSpecLimitLines(col: DG.Column, limits: SpecLimits, show: boolean): void {
  const lines = [['lsl', 'LSL', '#d62728'], ['usl', 'USL', '#d62728'], ['target', 'Target', '#2ca02c']] as const;
  for (const [key, title, color] of lines) {
    if (limits[key] !== undefined) {
      col.dataFrame.meta.formulaLines.addLine({title: `${col.name} ${title}`,
        formula: `\${${col.name}} = ${limits[key]}`, color, width: 1, style: 'dashed', visible: show});
    }
  }
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
