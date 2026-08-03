// Builds the comparison chart dataframes from already-extracted entries and matched targets.
// DG-dependent, no server calls.

import * as DG from 'datagrok-api/dg';
import {ComparisonEntry, ScalarTarget, ColumnTarget, RUN_COLUMN, isNumericType} from './types';

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

export interface ColumnComparisonResult {
  chartDf: DG.DataFrame;
  indexColumnName: string;
  // display name of the inner split column, when set on any participating table
  splitColumnName?: string;
  isKeyIndex: boolean;
}

export interface MultiColumnComparisonResult extends ColumnComparisonResult {
  // chart y-columns, one per selected target
  valueColumnNames: string[];
}

interface ParticipatingRun {
  entry: ComparisonEntry;
  binding: ColumnTarget['bindings'][number];
}

function getParticipating(target: ColumnTarget, entries: ComparisonEntry[]): ParticipatingRun[] {
  return entries
    .map((entry) => ({entry, binding: target.bindings.find((b) => b.entryId === entry.id)}))
    .filter((item) => !!item.binding) as ParticipatingRun[];
}

function getIndexKind(participating: ParticipatingRun[]) {
  const indexTypes = participating.map(({entry, binding}) =>
    entry.dataFrames.get(binding.tablePath)?.getCol(binding.indexColumnName).type);
  const isDatetimeIndex = indexTypes.every((type) => type === DG.TYPE.DATE_TIME);
  const isKeyIndex = !isDatetimeIndex && indexTypes.some((type) => type != null && !isNumericType(type));
  return {isDatetimeIndex, isKeyIndex};
}

function getSplitName(participating: ParticipatingRun[]): string | undefined {
  const rawName = participating.map(({binding}) => binding.splitColumnName).find((name) => !!name);
  return rawName === RUN_COLUMN ? `${rawName} (split)` : rawName;
}

const makeIndexColumn = (name: string, list: any[], isKeyIndex: boolean, isDatetimeIndex: boolean) => {
  if (isKeyIndex)
    return DG.Column.fromList(DG.COLUMN_TYPE.STRING, name, list);
  if (isDatetimeIndex)
    return DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, name, list);
  return DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, list);
};

/**
 * Concatenates raw rows of all participating runs into one long dataframe;
 * no row matching, the chart splits by run (and by the inner split column, if set).
 */
export function buildColumnComparison(
  target: ColumnTarget,
  entries: ComparisonEntry[],
): ColumnComparisonResult | null {
  const multi = buildMultiColumnComparison([target], entries);
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
): MultiColumnComparisonResult | null {
  const labels = new Map<string, string>();
  const seenNames = new Map<string, number>();
  for (const target of targets) {
    const count = (seenNames.get(target.displayName) ?? 0) + 1;
    seenNames.set(target.displayName, count);
    labels.set(target.key, count > 1 ? `${target.displayName} (${count})` : target.displayName);
  }

  const participating = getParticipating(targets[0], entries);
  if (participating.length < 2)
    return null;
  const {isDatetimeIndex, isKeyIndex} = getIndexKind(participating);
  const indexColumnName = participating[0].binding.indexColumnName;
  const splitColumnName = getSplitName(participating);

  const longIndex: any[] = [];
  const longSplits: string[] = [];
  const longRuns: string[] = [];
  const longValues = new Map<string, (number | null)[]>(targets.map((target) => [labels.get(target.key)!, []]));
  for (const {entry, binding} of participating) {
    const df = entry.dataFrames.get(binding.tablePath)!;
    const index = df.getCol(binding.indexColumnName).toList();
    const splits = binding.splitColumnName ?
      df.col(binding.splitColumnName)?.toList() : undefined;
    for (let i = 0; i < index.length; i++) {
      longIndex.push(isKeyIndex ? `${index[i]}` : index[i]);
      longRuns.push(entry.name);
      if (splitColumnName)
        longSplits.push(splits?.[i] == null ? '' : `${splits[i]}`);
    }
    for (const target of targets) {
      const targetBinding = target.bindings.find((b) => b.entryId === entry.id);
      const values = targetBinding ? df.getCol(targetBinding.columnName).toList() : index.map(() => null);
      longValues.get(labels.get(target.key)!)!.push(...values);
    }
  }

  const valueColumnNames = targets.map((target) => labels.get(target.key)!);
  const chartDf = DG.DataFrame.fromColumns([
    makeIndexColumn(indexColumnName, longIndex, isKeyIndex, isDatetimeIndex),
    ...splitColumnName ? [DG.Column.fromList(DG.COLUMN_TYPE.STRING, splitColumnName, longSplits)] : [],
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, RUN_COLUMN, longRuns),
    ...valueColumnNames.map((label) =>
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, label, longValues.get(label)!)),
  ]);
  chartDf.name = 'Comparison: multiple values';
  return {chartDf, indexColumnName, splitColumnName, isKeyIndex, valueColumnNames};
}
