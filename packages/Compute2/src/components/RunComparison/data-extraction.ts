import * as DG from 'datagrok-api/dg';
import {historyUtils} from '@datagrok-libraries/compute-utils';
import {deserialize} from '@datagrok-libraries/utils/src/json-serialization';
import {
  ComparisonEntryNodes, ScalarNodeInfo, TableNodeInfo, ColumnTarget, ScalarTarget,
  NumericSeries, KeyedSeries, ExclusionReason,
  isNumericType, alignSeriesByIndex, alignSeriesByKey, computeDelta, computeDeltaPct,
} from './comparison-core';

// mirrors CONFIG_PATH in reactive-tree-driver/src/runtime/funccall-utils.ts (not exported there)
const PIPELINE_CONFIG_OPTION = 'PIPELINE_CONFIG';

export const RUN_COLUMN = 'Run';

export type EntrySourceKind = 'workflow' | 'function' | 'raw';

export interface ComparisonEntry {
  id: string;
  name: string;
  sourceKind: EntrySourceKind;
  modelName: string;
  nodes: ComparisonEntryNodes;
  dataFrames: Map<string, DG.DataFrame>;
}

const SCALAR_PROPERTY_TYPES = new Set<string>([DG.TYPE.INT, DG.TYPE.FLOAT, DG.TYPE.BIG_INT]);

function getRunDisplayName(call: DG.FuncCall): string {
  const title = call.options?.['title'];
  if (title)
    return title;
  const funcName = call.func?.friendlyName ?? call.func?.name ?? 'Run';
  const started = (call.started as any)?.format?.('MMM D, HH:mm');
  return started ? `${funcName} — ${started}` : funcName;
}

function extractCallNodes(
  call: DG.FuncCall,
  pathPrefix: string,
  scalars: ScalarNodeInfo[],
  tables: TableNodeInfo[],
  dataFrames: Map<string, DG.DataFrame>,
) {
  const io = [
    ...[...call.inputParams.values()].map((p) => ({param: p, value: call.inputs[p.property.name]})),
    ...[...call.outputParams.values()].map((p) => ({param: p, value: call.outputs[p.property.name]})),
  ];
  for (const {param, value} of io) {
    const prop = param.property;
    const path = pathPrefix ? `${pathPrefix}/${prop.name}` : prop.name;
    const name = prop.caption ?? prop.name;
    if (SCALAR_PROPERTY_TYPES.has(prop.propertyType)) {
      scalars.push({
        path,
        name,
        valueType: prop.propertyType,
        units: prop.options?.['units'] || undefined,
        value: (value == null || isNaN(Number(value))) ? null : Number(value),
      });
    } else if (prop.propertyType === DG.TYPE.DATA_FRAME && value != null) {
      const df = value as DG.DataFrame;
      tables.push({
        path,
        name,
        columns: [...df.columns].map((col) => ({
          name: col.name,
          type: col.type,
          units: col.meta?.units || undefined,
        })),
        rowCount: df.rowCount,
      });
      dataFrames.set(path, df);
    }
  }
}

function collectWorkflowSteps(state: any, prefix: string[], acc: {path: string, funcCallId: string}[]) {
  if (state == null || typeof state !== 'object')
    return;
  if (state.type === 'funccall') {
    if (state.funcCallId)
      acc.push({path: [...prefix, state.configId].join('/'), funcCallId: state.funcCallId});
    return;
  }
  if (Array.isArray(state.steps)) {
    for (const step of state.steps)
      collectWorkflowSteps(step, state.configId ? [...prefix, state.configId] : prefix, acc);
  }
}

/** Builds a comparison entry from a loaded run: a workflow meta-call is flattened into its steps. */
export async function entryFromFuncCall(call: DG.FuncCall): Promise<ComparisonEntry> {
  const scalars: ScalarNodeInfo[] = [];
  const tables: TableNodeInfo[] = [];
  const dataFrames = new Map<string, DG.DataFrame>();
  const serializedConfig = call.options?.[PIPELINE_CONFIG_OPTION];
  const modelName = call.func?.friendlyName ?? call.func?.name ?? '';

  if (serializedConfig) {
    const state = deserialize(serializedConfig);
    const steps: {path: string, funcCallId: string}[] = [];
    collectWorkflowSteps(state, [], steps);
    for (const step of steps) {
      try {
        const stepCall = await historyUtils.loadRun(step.funcCallId);
        extractCallNodes(stepCall, step.path, scalars, tables, dataFrames);
      } catch (e) {
        console.warn(`Run comparison: failed to load step run ${step.funcCallId}`, e);
      }
    }
    return {
      id: call.id,
      name: getRunDisplayName(call),
      sourceKind: 'workflow',
      modelName,
      nodes: {entryId: call.id, entryName: getRunDisplayName(call), scalars, tables},
      dataFrames,
    };
  }

  extractCallNodes(call, '', scalars, tables, dataFrames);
  return {
    id: call.id,
    name: getRunDisplayName(call),
    sourceKind: 'function',
    modelName,
    nodes: {entryId: call.id, entryName: getRunDisplayName(call), scalars, tables},
    dataFrames,
  };
}

/** Builds a raw-data comparison entry from an open workspace table. */
export function entryFromDataFrame(df: DG.DataFrame): ComparisonEntry {
  const id = `table:${df.name}:${df.rowCount}`;
  const table: TableNodeInfo = {
    path: df.name,
    name: df.name,
    columns: [...df.columns].map((col) => ({
      name: col.name,
      type: col.type,
      units: col.meta?.units || undefined,
    })),
    rowCount: df.rowCount,
  };
  return {
    id,
    name: df.name,
    sourceKind: 'raw',
    modelName: '',
    nodes: {entryId: id, entryName: df.name, scalars: [], tables: [table]},
    dataFrames: new Map([[df.name, df]]),
  };
}

export type ComparisonMode = 'values' | 'delta' | 'deltaPct';

export interface ExcludedEntry {
  entryId: string;
  reason: ExclusionReason;
}

export interface ScalarComparisonResult {
  gridDf: DG.DataFrame;
  chartDf: DG.DataFrame;
  excluded: ExcludedEntry[];
}

/** One row per run: value plus Δ/Δ% against the baseline run. */
export function buildScalarComparison(
  target: ScalarTarget,
  entries: ComparisonEntry[],
  baselineEntryId: string,
): ScalarComparisonResult {
  const ordered = entries.filter((entry) => target.bindings.some((b) => b.entryId === entry.id));
  const bindings = ordered.map((entry) => target.bindings.find((b) => b.entryId === entry.id)!);
  const values = bindings.map((b) => b.value);
  const baselineIdx = Math.max(0, ordered.findIndex((entry) => entry.id === baselineEntryId));
  const baseline = values.map(() => values[baselineIdx]);
  const delta = computeDelta(values, baseline);
  const deltaPct = computeDeltaPct(values, baseline);

  // fromList with an explicit string type: fromStrings infers the type from values,
  // so numeric-looking run names would produce a numeric column and break chart legends
  const gridDf = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, RUN_COLUMN, ordered.map((entry) => entry.name)),
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Path', bindings.map((b) => b.path)),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, target.displayName, values),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Δ', delta),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Δ%', deltaPct),
  ]);
  gridDf.name = `Comparison: ${target.displayName}`;
  return {gridDf, chartDf: gridDf, excluded: []};
}

export interface ColumnComparisonResult {
  gridDf: DG.DataFrame;
  chartDf: DG.DataFrame;
  indexColumnName: string;
  isKeyIndex: boolean;
  excluded: ExcludedEntry[];
}

/**
 * Aligns the target column across runs on their user-defined index columns
 * (intersection join) and builds the wide grid df and the long chart df.
 */
export function buildColumnComparison(
  target: ColumnTarget,
  entries: ComparisonEntry[],
  baselineEntryId: string,
  mode: ComparisonMode,
): ColumnComparisonResult | null {
  const participating = entries
    .map((entry) => ({entry, binding: target.bindings.find((b) => b.entryId === entry.id)}))
    .filter((item) => !!item.binding) as {entry: ComparisonEntry, binding: ColumnTarget['bindings'][number]}[];
  if (participating.length < 2)
    return null;

  const isKeyIndex = participating.some(({entry, binding}) => {
    const df = entry.dataFrames.get(binding.tablePath);
    return df && !isNumericType(df.getCol(binding.indexColumnName).type);
  });

  const baselineIdx = Math.max(0, participating.findIndex(({entry}) => entry.id === baselineEntryId));
  const indexColumnName = participating[baselineIdx].binding.indexColumnName;
  const excluded: ExcludedEntry[] = [];

  let index: (number | string)[];
  let values: (number | null)[][];
  if (isKeyIndex) {
    const seriesList: KeyedSeries[] = participating.map(({entry, binding}) => {
      const df = entry.dataFrames.get(binding.tablePath)!;
      return {
        keys: df.getCol(binding.indexColumnName).toList().map((v: any) => `${v}`),
        values: df.getCol(binding.columnName).toList(),
      };
    });
    const aligned = alignSeriesByKey(seriesList, 'intersection');
    index = aligned.keys;
    values = aligned.values;
  } else {
    const seriesList: NumericSeries[] = participating.map(({entry, binding}) => {
      const df = entry.dataFrames.get(binding.tablePath)!;
      return {
        index: df.getCol(binding.indexColumnName).toList(),
        values: df.getCol(binding.columnName).toList(),
      };
    });
    const aligned = alignSeriesByIndex(seriesList, 'intersection');
    index = aligned.index;
    values = aligned.values;
  }

  if (index.length === 0) {
    for (const {entry} of participating) {
      if (entry.id !== participating[baselineIdx].entry.id)
        excluded.push({entryId: entry.id, reason: isKeyIndex ? 'no matching rows' : 'index grids differ'});
    }
    return {
      gridDf: DG.DataFrame.create(0),
      chartDf: DG.DataFrame.create(0),
      indexColumnName,
      isKeyIndex,
      excluded,
    };
  }

  const baselineValues = values[baselineIdx];
  const displayValues = values.map((series) => {
    if (mode === 'delta')
      return computeDelta(series, baselineValues);
    if (mode === 'deltaPct')
      return computeDeltaPct(series, baselineValues);
    return series;
  });

  const indexColumn = isKeyIndex ?
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, indexColumnName, index) :
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, indexColumnName, index);

  const gridColumns: DG.Column[] = [indexColumn];
  participating.forEach(({entry}, i) => {
    gridColumns.push(DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, entry.name, values[i]));
  });
  participating.forEach(({entry}, i) => {
    if (i === baselineIdx)
      return;
    gridColumns.push(DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, `Δ ${entry.name}`, computeDelta(values[i], baselineValues)));
    gridColumns.push(DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, `Δ% ${entry.name}`, computeDeltaPct(values[i], baselineValues)));
  });
  const gridDf = DG.DataFrame.fromColumns(gridColumns);
  gridDf.name = `Comparison: ${target.displayName}`;

  const longIndex: (number | string)[] = [];
  const longValues: (number | null)[] = [];
  const longRuns: string[] = [];
  participating.forEach(({entry}, i) => {
    if (mode !== 'values' && i === baselineIdx)
      return;
    for (let j = 0; j < index.length; j++) {
      longIndex.push(index[j]);
      longValues.push(displayValues[i][j]);
      longRuns.push(entry.name);
    }
  });
  const chartDf = DG.DataFrame.fromColumns([
    isKeyIndex ?
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, indexColumnName, longIndex) :
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, indexColumnName, longIndex),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, target.displayName, longValues),
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, RUN_COLUMN, longRuns),
  ]);
  chartDf.name = gridDf.name;

  return {gridDf, chartDf, indexColumnName, isKeyIndex, excluded};
}
