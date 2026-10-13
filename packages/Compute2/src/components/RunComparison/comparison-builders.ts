// Builds the comparison chart dataframes from already-extracted entries and matched targets.
// DG-dependent, no server calls.

import * as DG from 'datagrok-api/dg';
import {
  ComparisonEntry, ComparisonTarget, ScalarTarget, ColumnTarget, RUN_COLUMN, isNumericType, isTimeIndexType,
  TimeUnit, TIME_UNIT_MS, AxisConfigMap, toTargetUnits,
} from './types';
import {resolveDisplayUnit, pickAutoUnit, targetAxisMode} from './selection';

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

function dedupeLabels(
  targets: {key: string, displayName: string, units?: string}[], reserved: string[] = [],
): Map<string, string> {
  const taken = new Set(reserved);
  const labels = new Map<string, string>();
  for (const target of targets)
    labels.set(target.key, claimUniqueName(withUnits(target.displayName, target.units), taken));
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
  isKeyIndex: boolean;
  // every participating table is in independent-points mode — chart rows as points, not lines
  isScatter: boolean;
  // multi-value scatter melts values into one column keyed by target name
  melted?: {seriesColumnName: string, valueColumnName: string};
  // display unit of the elapsed axis; set only when the target charted in time-series mode
  timeSeriesUnit?: TimeUnit;
}


interface ParticipatingRun {
  entry: ComparisonEntry;
  binding: ColumnTarget['bindings'][number];
}

// every run any target charts; its rows come from the first target (in selection order)
// with a pick there — enabled bindings only, matching guarantees at most one per run
function getParticipating(targets: ColumnTarget[], entries: ComparisonEntry[]): ParticipatingRun[] {
  return entries.flatMap((entry) => {
    for (const target of targets) {
      const binding = target.bindings.find((b) => b.entryId === entry.id && entry.dataFrames.has(b.tablePath));
      if (binding)
        return [{entry, binding}];
    }
    return [];
  });
}

const allInMode = (participating: ParticipatingRun[], axisModes: AxisConfigMap, mode: 'timeseries' | 'points') =>
  targetAxisMode({bindings: participating.map(({binding}) => binding)}, axisModes, mode) === 'full';

type IndexAxis =
  | {kind: 'key' | 'datetime' | 'float'}
  | {kind: 'elapsed', unit: TimeUnit | 'auto'};

// time-series mode requires every participating table to be configured (require-all) and
// every index column to actually be numeric or datetime; otherwise the axis follows
// the index column types
function getIndexAxis(participating: ParticipatingRun[], axisModes?: AxisConfigMap): IndexAxis {
  const indexTypes = participating.map(({entry, binding}) =>
    entry.dataFrames.get(binding.tablePath)?.getCol(binding.indexColumnName).type);
  const isDatetimeOnly = indexTypes.every((type) => type === DG.TYPE.DATE_TIME);
  const timeChartable = indexTypes.every(isTimeIndexType);
  const configured = axisModes != null && allInMode(participating, axisModes, 'timeseries');
  if (configured && timeChartable) {
    return {
      kind: 'elapsed',
      unit: resolveDisplayUnit(participating.map(({binding}) => binding), axisModes!),
    };
  }
  if (isDatetimeOnly)
    return {kind: 'datetime'};
  if (indexTypes.some((type) => type != null && !isNumericType(type)))
    return {kind: 'key'};
  return {kind: 'float'};
}

function getSplitName(participating: ParticipatingRun[]): string | undefined {
  const rawName = participating.map(({binding}) => binding.splitColumnName).find((name) => !!name);
  return rawName === RUN_COLUMN ? `${rawName} (split)` : rawName;
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

// per-run index values in milliseconds: datetime via epoch ms, numeric scaled by its
// declared units; nulls preserved
function indexMsValues(
  {entry, binding}: ParticipatingRun,
  configs: AxisConfigMap,
): (number | null)[] {
  const col = entry.dataFrames.get(binding.tablePath)!.getCol(binding.indexColumnName);
  const units = configs.get(binding.entryId)?.get(binding.tablePath)?.units;
  const scale = col.type === DG.TYPE.DATE_TIME ? 1 : TIME_UNIT_MS[units ?? 's'];
  // Number(): bigint values cannot mix with number arithmetic
  return columnPrimitives(col).map((value) => value == null ? null : Number(value) * scale);
}

/**
 * Concatenates raw rows of all participating runs into one long dataframe, one value column
 * per target; no row matching, the chart splits by run (and by the inner split column, if set).
 */
export function buildColumnComparison(
  targets: ColumnTarget[],
  entries: ComparisonEntry[],
  axisModes?: AxisConfigMap,
): ColumnComparisonResult | null {
  const participating = getParticipating(targets, entries);
  if (participating.length < 2)
    return null;
  const splitColumnName = getSplitName(participating);
  // value labels must not shadow the built-in run/index/split columns
  const labels = dedupeLabels(targets, [
    RUN_COLUMN,
    participating[0].binding.indexColumnName,
    ...splitColumnName ? [splitColumnName] : [],
  ]);
  const axis = getIndexAxis(participating, axisModes);
  const isKeyIndex = axis.kind === 'key';
  // same require-all shape as the elapsed axis; a mixed numeric/datetime index falls
  // back to the key axis above, which a scatterplot cannot chart
  const isScatter = axisModes != null && !isKeyIndex && allInMode(participating, axisModes, 'points');

  // elapsed axis: convert every run's index to ms, subtract its own start,
  // and express the result in one shared display unit
  let displayUnit: TimeUnit | undefined;
  let elapsedByRun: Map<ParticipatingRun, (number | null)[]> | undefined;
  if (axis.kind === 'elapsed') {
    const msByRun = participating.map((run) => ({run, ms: indexMsValues(run, axisModes!)}));
    let maxElapsedMs = 0;
    const shifted = msByRun.map(({run, ms}) => {
      // per-table start: min of non-null values (rows may be unsorted); all-null -> 0
      const t0 = ms.reduce<number | null>((acc, value) =>
        value == null ? acc : (acc == null ? value : Math.min(acc, value)), null) ?? 0;
      const elapsed = ms.map((value) => value == null ? null : value - t0);
      for (const value of elapsed) {
        if (value != null)
          maxElapsedMs = Math.max(maxElapsedMs, value);
      }
      return {run, elapsed};
    });
    displayUnit = axis.unit === 'auto' ? pickAutoUnit(maxElapsedMs) : axis.unit;
    elapsedByRun = new Map(shifted.map(({run, elapsed}) =>
      [run, elapsed.map((value) => value == null ? null : value / TIME_UNIT_MS[displayUnit!])]));
  }

  // the axis is labeled with its unit (the display unit on the elapsed axis, the first run's
  // index units on a float axis); keep it clear of the other column labels
  let indexColumnName = participating[0].binding.indexColumnName;
  const axisUnits = axis.kind === 'elapsed' ? displayUnit : axis.kind === 'float' ? targets[0].indexUnits : undefined;
  if (axisUnits) {
    const taken = new Set([RUN_COLUMN, ...splitColumnName ? [splitColumnName] : [],
      ...targets.map((target) => labels.get(target.key)!)]);
    indexColumnName = claimUniqueName(`${indexColumnName} (${axisUnits})`, taken);
  }

  const longIndex: any[] = [];
  const longSplits: string[] = [];
  const longRuns: string[] = [];
  const longValues = new Map<string, (number | null)[]>(targets.map((target) => [labels.get(target.key)!, []]));
  for (const run of participating) {
    const {entry, binding} = run;
    const df = entry.dataFrames.get(binding.tablePath)!;
    // the elapsed axis converts through the time units picker instead of unit mappings
    const index = elapsedByRun?.get(run) ?? (axis.kind === 'float' ?
      df.getCol(binding.indexColumnName).toList().map((v) => toTargetUnits(v, binding.indexFactor)) :
      df.getCol(binding.indexColumnName).toList());
    const splits = binding.splitColumnName ?
      df.col(binding.splitColumnName)?.toList() : undefined;
    for (let i = 0; i < index.length; i++) {
      longIndex.push(isKeyIndex ? `${index[i]}` : index[i]);
      longRuns.push(entry.name);
      if (splitColumnName)
        longSplits.push(splits?.[i] == null ? '' : `${splits[i]}`);
    }
    // a pick from another table cannot share this block's rows — padded like a missing one
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
    const taken = new Set([indexColumnName, RUN_COLUMN, ...splitColumnName ? [splitColumnName] : []]);
    melted = {seriesColumnName: claimUniqueName('Data', taken), valueColumnName: claimUniqueName('Value', taken)};
    const meltedIndex: any[] = [];
    const meltedSplits: string[] = [];
    const meltedRuns: string[] = [];
    const meltedSeries: string[] = [];
    const meltedValues: (number | null)[] = [];
    for (const label of valueColumnNames) {
      meltedIndex.push(...longIndex);
      meltedRuns.push(...longRuns);
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
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, melted.seriesColumnName, meltedSeries),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, melted.valueColumnName, toFloatValues(meltedValues)),
    ]);
  } else {
    chartDf = DG.DataFrame.fromColumns([
      makeIndexColumn(indexColumnName, longIndex, axis),
      ...splitColumnName ? [DG.Column.fromList(DG.COLUMN_TYPE.STRING, splitColumnName, longSplits)] : [],
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, RUN_COLUMN, longRuns),
      ...valueColumnNames.map((label) =>
        DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, label, toFloatValues(longValues.get(label)!))),
    ]);
  }
  chartDf.name = targets.length === 1 ? `Comparison: ${targets[0].displayName}` : 'Comparison: multiple values';
  return {
    chartDf, indexColumnName, splitColumnName, isKeyIndex, isScatter, valueColumnNames,
    valueColumnName: valueColumnNames[0],
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
