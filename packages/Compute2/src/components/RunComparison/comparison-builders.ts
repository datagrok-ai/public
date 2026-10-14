// Builds the comparison chart dataframes from already-extracted entries and matched targets.
// DG-dependent, no server calls.

import * as DG from 'datagrok-api/dg';
import {
  ComparisonEntry, ComparisonTarget, ScalarTarget, ColumnTarget, RUN_COLUMN, isNumericType, isTimeIndexType,
  TimeUnit, TIME_UNIT_MS, AxisConfigMap, toTargetUnits,
} from './types';
import {resolveDisplayUnit, pickAutoUnit, targetAxisMode, targetContexts} from './selection';
import {normalizeUnits} from './matching';

// appends (2), (3)... while the name is taken; "claim" — mutates the set,
// registering the result so successive calls stay collision-free
const claimUniqueName = (base: string, taken: Set<string>) => {
  let name = base;
  for (let n = 2; taken.has(name); n++)
    name = `${base} (${n})`;
  taken.add(name);
  return name;
};

// dataframes may carry bigint (int64) columns; those values cannot enter a float column as-is
const toFloatValues = (list: any[]) => list.map((v) => typeof v === 'bigint' ? Number(v) : v);

export interface ScalarComparisonResult {
  chartDf: DG.DataFrame;
  // actual value column label: the display name deduped against the built-in columns
  valueColumnName: string;
}

/** One row per run: the scalar value with its source path. */
export function buildScalarComparison(
  target: ScalarTarget,
  entries: ComparisonEntry[],
): ScalarComparisonResult {
  const ordered = entries.filter((entry) => target.bindings.some((b) => b.entryId === entry.id));
  const bindings = ordered.map((entry) => target.bindings.find((b) => b.entryId === entry.id)!);
  const valueColumnName = claimUniqueName(withUnits(target.displayName, target.units), new Set([RUN_COLUMN, 'Path']));

  // fromList with an explicit string type: fromStrings infers the type from values,
  // so numeric-looking run names would produce a numeric column and break chart legends
  const chartDf = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, RUN_COLUMN, ordered.map((entry) => entry.name)),
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Path', bindings.map((b) => b.friendlyPath ?? b.path)),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, valueColumnName,
      toFloatValues(bindings.map((b) => toTargetUnits(b.value, b.factor)))),
  ]);
  chartDf.name = `Comparison: ${target.displayName}`;
  return {chartDf, valueColumnName};
}

// units go into the column name, so they survive any export (metadata stays in Datagrok)
const withUnits = (name: string, units?: string) => units ? `${name} (${units})` : name;

// names repeated in the selection get their step, as in the compare list
function dedupeLabels(targets: ComparisonTarget[], reserved: string[] = []): Map<string, string> {
  const contexts = targetContexts(targets);
  const taken = new Set(reserved);
  const labels = new Map<string, string>();
  for (const target of targets) {
    const context = contexts.get(target.key);
    const name = context ? `${context} · ${target.displayName}` : target.displayName;
    labels.set(target.key, claimUniqueName(withUnits(name, target.units), taken));
  }
  return labels;
}

export interface MultiScalarComparisonResult {
  chartDf: DG.DataFrame;
  valueColumnNames: string[];
}

/** One row per run, one value column per target — wide shape for radar/PC plot axes. */
export function buildMultiScalarComparison(
  targets: ScalarTarget[],
  entries: ComparisonEntry[],
): MultiScalarComparisonResult | null {
  const labels = dedupeLabels(targets, [RUN_COLUMN]);
  const participating = entries.filter((entry) =>
    targets.some((target) => target.bindings.some((b) => b.entryId === entry.id)));
  if (participating.length < 2)
    return null;

  const valueColumnNames = targets.map((target) => labels.get(target.key)!);
  const chartDf = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, RUN_COLUMN, participating.map((entry) => entry.name)),
    ...targets.map((target) => DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, labels.get(target.key)!,
      toFloatValues(participating.map((entry) => {
        const binding = target.bindings.find((b) => b.entryId === entry.id);
        return binding ? toTargetUnits(binding.value, binding.factor) : null;
      })))),
  ]);
  chartDf.name = 'Comparison: multiple values';
  return {chartDf, valueColumnNames};
}

export interface ColumnComparisonResult {
  chartDf: DG.DataFrame;
  // chart y columns, one per target: display names deduped against the run/index/split
  // columns; more than one marks a multi-value result; absent from the melted frame
  valueColumnNames: string[];
  // the first of valueColumnNames, charted by single-value results
  valueColumnName: string;
  indexColumnName: string;
  // display name of the inner split column, when set on any participating table
  splitColumnName?: string;
  // set when a run charts picks from several tables
  stepColumnName?: string;
  isKeyIndex: boolean;
  // the key index comes from mixing datetime and numeric indexes, which only the elapsed axis can chart
  mixedTimeIndex?: boolean;
  // every target has the same explicit units, so the values can share one y scale
  sameUnits: boolean;
  // every participating table is in independent-points mode — chart rows as points, not lines
  isScatter: boolean;
  // multi-value scatter melts values into one column keyed by target name
  melted?: {seriesColumnName: string, valueColumnName: string};
  // display unit of the elapsed axis; set only when the target charted in time-series mode
  timeSeriesUnit?: TimeUnit;
}


interface RowBlock {
  entry: ComparisonEntry;
  binding: ColumnTarget['bindings'][number];
}

// one block of rows per (run, table) any target picks from, in selection order; enabled
// bindings only, matching guarantees at most one per target and run
function getBlocks(targets: ColumnTarget[], entries: ComparisonEntry[]): RowBlock[] {
  return entries.flatMap((entry) => {
    const blocks = new Map<string, RowBlock>();
    for (const target of targets) {
      const binding = target.bindings.find((b) => b.entryId === entry.id && entry.dataFrames.has(b.tablePath));
      if (binding && !blocks.has(binding.tablePath))
        blocks.set(binding.tablePath, {entry, binding});
    }
    return [...blocks.values()];
  });
}

const allInMode = (blocks: RowBlock[], axisModes: AxisConfigMap, mode: 'timeseries' | 'points') =>
  targetAxisMode({bindings: blocks.map(({binding}) => binding)}, axisModes, mode) === 'full';

type IndexAxis =
  | {kind: 'datetime' | 'float'}
  | {kind: 'key', mixedTime: boolean}
  | {kind: 'elapsed', unit: TimeUnit | 'auto'};

// time-series mode requires every charted table to be configured (require-all) and
// every index column to actually be numeric or datetime; otherwise the axis follows
// the index column types
function getIndexAxis(blocks: RowBlock[], axisModes?: AxisConfigMap): IndexAxis {
  const indexTypes = blocks.map(({entry, binding}) =>
    entry.dataFrames.get(binding.tablePath)?.getCol(binding.indexColumnName).type);
  const isDatetimeOnly = indexTypes.every((type) => type === DG.TYPE.DATE_TIME);
  const timeChartable = indexTypes.every(isTimeIndexType);
  const configured = axisModes != null && allInMode(blocks, axisModes, 'timeseries');
  if (configured && timeChartable) {
    return {
      kind: 'elapsed',
      unit: resolveDisplayUnit(blocks.map(({binding}) => binding), axisModes!),
    };
  }
  if (isDatetimeOnly)
    return {kind: 'datetime'};
  if (indexTypes.some((type) => type != null && !isNumericType(type)))
    return {kind: 'key', mixedTime: timeChartable};
  return {kind: 'float'};
}

function getSplitName(blocks: RowBlock[]): string | undefined {
  const rawName = blocks.map(({binding}) => binding.splitColumnName).find((name) => !!name);
  return rawName === RUN_COLUMN ? `${rawName} (split)` : rawName;
}

// index units every charted table agrees on (up to spelling), in the first one's spelling
function sharedIndexUnits(blocks: RowBlock[]): string | undefined {
  const units = blocks.map(({entry, binding}) => entry.nodes.tables
    .find((table) => table.path === binding.tablePath)?.columns
    .find((col) => col.name === binding.indexColumnName)?.units);
  return units[0] && units.every((item) => normalizeUnits(item) === normalizeUnits(units[0])) ? units[0] : undefined;
}

const makeIndexColumn = (name: string, list: any[], axis: IndexAxis) => {
  if (axis.kind === 'key')
    return DG.Column.fromList(DG.COLUMN_TYPE.STRING, name, list);
  if (axis.kind === 'datetime')
    return DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, name, list);
  if (axis.kind === 'elapsed') {
    // explicit 64-bit storage: fromList(FLOAT) may store 32-bit, which quantizes
    // epoch-scale values (unaligned mode) by ~±100 s
    const data = new Float64Array(list.length);
    for (let i = 0; i < list.length; i++)
      data[i] = list[i] == null ? DG.FLOAT_NULL : list[i];
    return DG.Column.fromFloat64Array(name, data);
  }
  return DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, toFloatValues(list));
};

/** Column values as primitives: nulls for missing cells, datetimes as epoch ms. */
export function columnPrimitives(col: DG.Column): unknown[] {
  const values = new Array(col.length);
  for (let i = 0; i < col.length; i++) {
    const value = col.isNone(i) ? null : col.get(i);
    values[i] = value instanceof Object ? value.valueOf() : value;
  }
  return values;
}

// per-table index values in milliseconds: datetime via epoch ms, numeric scaled by its
// declared units; nulls preserved
function indexMsValues(
  {entry, binding}: RowBlock,
  configs: AxisConfigMap,
): (number | null)[] {
  const col = entry.dataFrames.get(binding.tablePath)!.getCol(binding.indexColumnName);
  const units = configs.get(binding.entryId)?.get(binding.tablePath)?.units;
  const scale = col.type === DG.TYPE.DATE_TIME ? 1 : TIME_UNIT_MS[units ?? 's'];
  // Number(): bigint values cannot mix with number arithmetic
  return columnPrimitives(col).map((value) => value == null ? null : Number(value) * scale);
}

/**
 * Concatenates raw rows of all (run, table) blocks into one long dataframe, one value column
 * per target; no row matching, the chart splits by run, by step when a run charts several
 * tables, and by the inner split column, if set.
 */
export function buildColumnComparison(
  targets: ColumnTarget[],
  entries: ComparisonEntry[],
  axisModes?: AxisConfigMap,
): ColumnComparisonResult | null {
  const blocks = getBlocks(targets, entries);
  const runCount = new Set(blocks.map(({entry}) => entry)).size;
  if (runCount < 2)
    return null;
  const splitColumnName = getSplitName(blocks);
  const builtIn = [RUN_COLUMN, blocks[0].binding.indexColumnName, ...splitColumnName ? [splitColumnName] : []];
  const stepColumnName = blocks.length > runCount ? claimUniqueName('Step', new Set(builtIn)) : undefined;
  // value labels must not shadow the built-in run/index/split/step columns
  const labels = dedupeLabels(targets, [...builtIn, ...stepColumnName ? [stepColumnName] : []]);
  const axis = getIndexAxis(blocks, axisModes);
  const isKeyIndex = axis.kind === 'key';
  // same require-all shape as the elapsed axis; a mixed numeric/datetime index falls
  // back to the key axis above, which a scatterplot cannot chart
  const isScatter = axisModes != null && !isKeyIndex && allInMode(blocks, axisModes, 'points');

  // elapsed axis: convert every block's index to ms, subtract its own start,
  // and express the result in one shared display unit
  let displayUnit: TimeUnit | undefined;
  let elapsedByBlock: Map<RowBlock, (number | null)[]> | undefined;
  if (axis.kind === 'elapsed') {
    const msByBlock = blocks.map((block) => ({block, ms: indexMsValues(block, axisModes!)}));
    let maxElapsedMs = 0;
    const shifted = msByBlock.map(({block, ms}) => {
      // per-table start: min of non-null values (rows may be unsorted); all-null -> 0
      const t0 = ms.reduce<number | null>((acc, value) =>
        value == null ? acc : (acc == null ? value : Math.min(acc, value)), null) ?? 0;
      const elapsed = ms.map((value) => value == null ? null : value - t0);
      for (const value of elapsed) {
        if (value != null)
          maxElapsedMs = Math.max(maxElapsedMs, value);
      }
      return {block, elapsed};
    });
    displayUnit = axis.unit === 'auto' ? pickAutoUnit(maxElapsedMs) : axis.unit;
    elapsedByBlock = new Map(shifted.map(({block, elapsed}) =>
      [block, elapsed.map((value) => value == null ? null : value / TIME_UNIT_MS[displayUnit!])]));
  }

  // the axis is labeled with its unit (the display unit on the elapsed axis, the runs' shared
  // index units on a float axis); keep it clear of the other column labels
  let indexColumnName = blocks[0].binding.indexColumnName;
  const axisUnits = axis.kind === 'elapsed' ? displayUnit :
    axis.kind === 'float' ? sharedIndexUnits(blocks) : undefined;
  if (axisUnits) {
    const taken = new Set([RUN_COLUMN, ...splitColumnName ? [splitColumnName] : [],
      ...stepColumnName ? [stepColumnName] : [], ...targets.map((target) => labels.get(target.key)!)]);
    indexColumnName = claimUniqueName(`${indexColumnName} (${axisUnits})`, taken);
  }

  const longIndex: any[] = [];
  const longSplits: string[] = [];
  const longRuns: string[] = [];
  const longSteps: string[] = [];
  const longValues = new Map<string, (number | null)[]>(targets.map((target) => [labels.get(target.key)!, []]));
  for (const block of blocks) {
    const {entry, binding} = block;
    const df = entry.dataFrames.get(binding.tablePath)!;
    const step = binding.tableFriendlyPath ?? binding.tablePath;
    const index = elapsedByBlock?.get(block) ?? df.getCol(binding.indexColumnName).toList();
    const splits = binding.splitColumnName ?
      df.col(binding.splitColumnName)?.toList() : undefined;
    for (let i = 0; i < index.length; i++) {
      longIndex.push(isKeyIndex && index[i] != null ? `${index[i]}` : index[i]);
      longRuns.push(entry.name);
      if (stepColumnName)
        longSteps.push(step);
      if (splitColumnName)
        longSplits.push(splits?.[i] == null ? '' : `${splits[i]}`);
    }
    // a pick from another table charts in that table's block — padded here like a missing one
    for (const target of targets) {
      const targetBinding = target.bindings.find((b) => b.entryId === entry.id);
      const usable = targetBinding != null && targetBinding.tablePath === binding.tablePath;
      const values = usable ?
        df.getCol(targetBinding.columnName).toList().map((v) => toTargetUnits(v, targetBinding.factor)) :
        index.map(() => null);
      longValues.get(labels.get(target.key)!)!.push(...values);
    }
  }

  const valueColumnNames = targets.map((target) => labels.get(target.key)!);

  // multi-value scatter: one point cloud per target on a shared y axis, keyed by
  // the target name — the scatterplot charts a single y column
  let melted: {seriesColumnName: string, valueColumnName: string} | undefined;
  let chartDf: DG.DataFrame;
  if (isScatter && targets.length > 1) {
    const taken = new Set([indexColumnName, RUN_COLUMN, ...splitColumnName ? [splitColumnName] : [],
      ...stepColumnName ? [stepColumnName] : []]);
    melted = {seriesColumnName: claimUniqueName('Data', taken), valueColumnName: claimUniqueName('Value', taken)};
    const meltedIndex: any[] = [];
    const meltedSplits: string[] = [];
    const meltedRuns: string[] = [];
    const meltedSteps: string[] = [];
    const meltedSeries: string[] = [];
    const meltedValues: (number | null)[] = [];
    for (const label of valueColumnNames) {
      meltedIndex.push(...longIndex);
      meltedRuns.push(...longRuns);
      meltedSteps.push(...longSteps);
      if (splitColumnName)
        meltedSplits.push(...longSplits);
      for (const value of longValues.get(label)!) {
        meltedSeries.push(label);
        meltedValues.push(value);
      }
    }
    chartDf = DG.DataFrame.fromColumns([
      makeIndexColumn(indexColumnName, meltedIndex, axis),
      ...splitColumnName ? [DG.Column.fromList(DG.COLUMN_TYPE.STRING, splitColumnName, meltedSplits)] : [],
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, RUN_COLUMN, meltedRuns),
      ...stepColumnName ? [DG.Column.fromList(DG.COLUMN_TYPE.STRING, stepColumnName, meltedSteps)] : [],
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, melted.seriesColumnName, meltedSeries),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, melted.valueColumnName, toFloatValues(meltedValues)),
    ]);
  } else {
    chartDf = DG.DataFrame.fromColumns([
      makeIndexColumn(indexColumnName, longIndex, axis),
      ...splitColumnName ? [DG.Column.fromList(DG.COLUMN_TYPE.STRING, splitColumnName, longSplits)] : [],
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, RUN_COLUMN, longRuns),
      ...stepColumnName ? [DG.Column.fromList(DG.COLUMN_TYPE.STRING, stepColumnName, longSteps)] : [],
      ...valueColumnNames.map((label) =>
        DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, label, toFloatValues(longValues.get(label)!))),
    ]);
  }
  chartDf.name = targets.length === 1 ? `Comparison: ${targets[0].displayName}` : 'Comparison: multiple values';
  return {
    chartDf, indexColumnName, splitColumnName, stepColumnName, isKeyIndex, isScatter, valueColumnNames,
    sameUnits: targets.every((target) =>
      !!target.units && normalizeUnits(target.units) === normalizeUnits(targets[0].units)),
    valueColumnName: valueColumnNames[0],
    ...axis.kind === 'key' && axis.mixedTime ? {mixedTimeIndex: true} : {},
    ...melted ? {melted} : {},
    ...displayUnit ? {timeSeriesUnit: displayUnit} : {},
  };
}

export type Comparison =
  | {kind: 'scalar'} & ScalarComparisonResult
  | {kind: 'multi-scalar'} & MultiScalarComparisonResult
  | {kind: 'column'} & ColumnComparisonResult;

/** Chart data for the selected targets (anchor first); targets of the anchor's kind only. */
export function buildComparison(
  selected: ComparisonTarget[],
  entries: ComparisonEntry[],
  axisModes?: AxisConfigMap,
): Comparison | null {
  if (selected.length === 0 || entries.length < 2)
    return null;
  const scalars = selected.filter((target): target is ScalarTarget => target.kind === 'scalar');
  if (selected[0].kind === 'scalar') {
    if (scalars.length === 1)
      return {kind: 'scalar', ...buildScalarComparison(scalars[0], entries)};
    const result = buildMultiScalarComparison(scalars, entries);
    return result && {kind: 'multi-scalar', ...result};
  }
  const columns = selected.filter((target): target is ColumnTarget => target.kind === 'column');
  const result = buildColumnComparison(columns, entries, axisModes);
  return result && {kind: 'column', ...result};
}
