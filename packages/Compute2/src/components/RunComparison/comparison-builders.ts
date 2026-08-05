// Builds the comparison chart dataframes from already-extracted entries and matched targets.
// DG-dependent, no server calls.

import * as DG from 'datagrok-api/dg';
import {
  ComparisonEntry, ScalarTarget, ColumnTarget, RUN_COLUMN, isNumericType,
  TimeUnit, TIME_UNIT_MS, AxisConfigMap,
} from './types';
import {resolveDisplayUnit, pickAutoUnit} from './selection';

export interface ScalarComparisonResult {
  chartDf: DG.DataFrame;
}

/** One row per run: the scalar value with its source path. */
export function buildScalarComparison(
  target: ScalarTarget,
  entries: ComparisonEntry[],
): ScalarComparisonResult {
  const ordered = entries.filter((entry) => target.bindings.some((b) => b.entryId === entry.id));
  const bindings = ordered.map((entry) => target.bindings.find((b) => b.entryId === entry.id)!);

  // fromList with an explicit string type: fromStrings infers the type from values,
  // so numeric-looking run names would produce a numeric column and break chart legends
  const chartDf = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, RUN_COLUMN, ordered.map((entry) => entry.name)),
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Path', bindings.map((b) => b.friendlyPath ?? b.path)),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, target.displayName, bindings.map((b) => b.value)),
  ]);
  chartDf.name = `Comparison: ${target.displayName}`;
  return {chartDf};
}

function dedupeLabels(targets: {key: string, displayName: string}[]): Map<string, string> {
  const labels = new Map<string, string>();
  const seenNames = new Map<string, number>();
  for (const target of targets) {
    const count = (seenNames.get(target.displayName) ?? 0) + 1;
    seenNames.set(target.displayName, count);
    labels.set(target.key, count > 1 ? `${target.displayName} (${count})` : target.displayName);
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
  const labels = dedupeLabels(targets);
  const participating = entries.filter((entry) =>
    targets.some((target) => target.bindings.some((b) => b.entryId === entry.id)));
  if (participating.length < 2)
    return null;

  const valueColumnNames = targets.map((target) => labels.get(target.key)!);
  const chartDf = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, RUN_COLUMN, participating.map((entry) => entry.name)),
    ...targets.map((target) => DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, labels.get(target.key)!,
      participating.map((entry) =>
        target.bindings.find((b) => b.entryId === entry.id)?.value ?? null))),
  ]);
  chartDf.name = 'Comparison: multiple values';
  return {chartDf, valueColumnNames};
}

export interface ColumnComparisonResult {
  chartDf: DG.DataFrame;
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

export interface MultiColumnComparisonResult extends ColumnComparisonResult {
  // chart y-columns, one per selected target; not present in the melted frame
  valueColumnNames: string[];
}


interface ParticipatingRun {
  entry: ComparisonEntry;
  binding: ColumnTarget['bindings'][number];
}

// enabled bindings only; matching guarantees at most one per run
function getParticipating(target: ColumnTarget, entries: ComparisonEntry[]): ParticipatingRun[] {
  return entries
    .map((entry) => ({entry, binding: target.bindings.find((b) =>
      b.entryId === entry.id && entry.dataFrames.has(b.tablePath))}))
    .filter((item) => !!item.binding) as ParticipatingRun[];
}

type IndexAxis =
  | {kind: 'key' | 'datetime' | 'float'}
  | {kind: 'elapsed', unit: TimeUnit | 'auto'};

// time-series mode requires every participating table to be configured (require-all) and
// every index column to actually be numeric or datetime; otherwise the legacy rules apply
function getIndexAxis(participating: ParticipatingRun[], axisModes?: AxisConfigMap): IndexAxis {
  const indexTypes = participating.map(({entry, binding}) =>
    entry.dataFrames.get(binding.tablePath)?.getCol(binding.indexColumnName).type);
  const isDatetimeOnly = indexTypes.every((type) => type === DG.TYPE.DATE_TIME);
  const timeChartable = indexTypes.every((type) =>
    type != null && (isNumericType(type) || type === DG.TYPE.DATE_TIME));
  const configured = axisModes != null && participating.every(({binding}) =>
    axisModes.get(binding.entryId)?.get(binding.tablePath)?.mode === 'timeseries');
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
  return DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, list);
};

// per-run index values in milliseconds: datetime via epoch ms, numeric scaled by its
// declared units; nulls preserved
function indexMsValues(
  {entry, binding}: ParticipatingRun,
  configs: AxisConfigMap,
): (number | null)[] {
  const col = entry.dataFrames.get(binding.tablePath)!.getCol(binding.indexColumnName);
  const isDatetime = col.type === DG.TYPE.DATE_TIME;
  const units = configs.get(binding.entryId)?.get(binding.tablePath)?.units;
  const scale = TIME_UNIT_MS[units ?? 's'];
  const values: (number | null)[] = new Array(col.length);
  for (let i = 0; i < col.length; i++) {
    if (col.isNone(i)) {
      values[i] = null;
      continue;
    }
    const value = col.get(i);
    // Number(): bigint columns return BigInt values, which cannot mix with number arithmetic
    values[i] = isDatetime ? (value as {valueOf(): number}).valueOf() : Number(value) * scale;
  }
  return values;
}

/**
 * Concatenates raw rows of all participating runs into one long dataframe;
 * no row matching, the chart splits by run (and by the inner split column, if set).
 */
export function buildColumnComparison(
  target: ColumnTarget,
  entries: ComparisonEntry[],
  axisModes?: AxisConfigMap,
): ColumnComparisonResult | null {
  const multi = buildMultiColumnComparison([target], entries, axisModes);
  if (!multi)
    return null;
  // the presence of valueColumnNames is what marks a multi-value result downstream
  const {valueColumnNames: _valueColumnNames, ...single} = multi;
  single.chartDf.name = `Comparison: ${target.displayName}`;
  return single;
}

/**
 * Same as buildColumnComparison for several targets sharing the bindings signature:
 * one value column per target over the shared concatenated rows.
 */
export function buildMultiColumnComparison(
  targets: ColumnTarget[],
  entries: ComparisonEntry[],
  axisModes?: AxisConfigMap,
): MultiColumnComparisonResult | null {
  const labels = dedupeLabels(targets);
  const participating = getParticipating(targets[0], entries);
  if (participating.length < 2)
    return null;
  const axis = getIndexAxis(participating, axisModes);
  const isKeyIndex = axis.kind === 'key';
  // same require-all shape as the elapsed axis; a mixed numeric/datetime index falls
  // back to the key axis above, which a scatterplot cannot chart
  const isScatter = axisModes != null && !isKeyIndex &&
    participating.every(({binding}) =>
      axisModes.get(binding.entryId)?.get(binding.tablePath)?.mode === 'points');
  const splitColumnName = getSplitName(participating);

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

  // the elapsed axis is labeled with its unit; keep it clear of the value column labels
  let indexColumnName = participating[0].binding.indexColumnName;
  if (axis.kind === 'elapsed') {
    const taken = new Set(targets.map((target) => labels.get(target.key)!));
    const base = `${indexColumnName} (${displayUnit})`;
    indexColumnName = base;
    for (let n = 2; taken.has(indexColumnName); n++)
      indexColumnName = `${base} (${n})`;
  }

  const longIndex: any[] = [];
  const longSplits: string[] = [];
  const longRuns: string[] = [];
  const longValues = new Map<string, (number | null)[]>(targets.map((target) => [labels.get(target.key)!, []]));
  for (const run of participating) {
    const {entry, binding} = run;
    const df = entry.dataFrames.get(binding.tablePath)!;
    const index = elapsedByRun?.get(run) ?? df.getCol(binding.indexColumnName).toList();
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
      const values = usable ? df.getCol(targetBinding.columnName).toList() : index.map(() => null);
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
    const unique = (base: string) => {
      let name = base;
      for (let n = 2; taken.has(name); n++)
        name = `${base} (${n})`;
      taken.add(name);
      return name;
    };
    melted = {seriesColumnName: unique('Data'), valueColumnName: unique('Value')};
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
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, melted.valueColumnName, meltedValues),
    ]);
  } else {
    chartDf = DG.DataFrame.fromColumns([
      makeIndexColumn(indexColumnName, longIndex, axis),
      ...splitColumnName ? [DG.Column.fromList(DG.COLUMN_TYPE.STRING, splitColumnName, longSplits)] : [],
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, RUN_COLUMN, longRuns),
      ...valueColumnNames.map((label) =>
        DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, label, longValues.get(label)!)),
    ]);
  }
  chartDf.name = 'Comparison: multiple values';
  return {
    chartDf, indexColumnName, splitColumnName, isKeyIndex, isScatter, valueColumnNames,
    ...melted ? {melted} : {},
    ...displayUnit ? {timeSeriesUnit: displayUnit} : {},
  };
}
