import * as DG from 'datagrok-api/dg';
import {historyUtils} from '@datagrok-libraries/compute-utils';
import {getRunTitle} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {deserialize} from '@datagrok-libraries/utils/src/json-serialization';
import {
  ComparisonEntryNodes, ScalarNodeInfo, TableNodeInfo, ColumnTarget, ScalarTarget,
  isNumericType,
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

const columnInfos = (df: DG.DataFrame) => [...df.columns].map((col) => ({
  name: col.name,
  type: col.type,
  units: col.meta?.units || undefined,
}));

function extractCallNodes(
  call: DG.FuncCall,
  pathPrefix: string,
  friendlyPrefix: string,
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
    const friendlyPath = friendlyPrefix ? `${friendlyPrefix} · ${name}` : name;
    if (SCALAR_PROPERTY_TYPES.has(prop.propertyType)) {
      scalars.push({
        path,
        name,
        friendlyPath,
        valueType: prop.propertyType,
        units: prop.options?.['units'] || undefined,
        value: (value == null || isNaN(Number(value))) ? null : Number(value),
      });
    } else if (prop.propertyType === DG.TYPE.DATA_FRAME && value != null) {
      const df = value as DG.DataFrame;
      tables.push({
        path,
        name,
        friendlyPath,
        nqName: call.func?.nqName,
        columns: columnInfos(df),
        rowCount: df.rowCount,
      });
      dataFrames.set(path, df);
    }
  }
}

interface WorkflowStep {
  path: string;
  friendlyPath: string;
  funcCallId: string;
}

function collectWorkflowSteps(
  state: any,
  prefix: string[],
  friendlyPrefix: string[],
  acc: WorkflowStep[],
  isRoot = true,
) {
  if (state == null || typeof state !== 'object')
    return;
  const friendlyName = state.friendlyName ?? state.configId;
  if (state.type === 'funccall') {
    if (state.funcCallId) {
      acc.push({
        path: [...prefix, state.configId].join('/'),
        friendlyPath: [...friendlyPrefix, friendlyName].join(' · '),
        funcCallId: state.funcCallId,
      });
    }
    return;
  }
  if (Array.isArray(state.steps)) {
    // the root pipeline name is shared by every step, so it is left out of friendly paths
    const nestedFriendly = isRoot ? friendlyPrefix : [...friendlyPrefix, friendlyName];
    for (const step of state.steps) {
      collectWorkflowSteps(
        step, state.configId ? [...prefix, state.configId] : prefix, nestedFriendly, acc, false);
    }
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
    const steps: WorkflowStep[] = [];
    collectWorkflowSteps(state, [], [], steps);
    for (const step of steps) {
      try {
        const stepCall = await historyUtils.loadRun(step.funcCallId);
        extractCallNodes(stepCall, step.path, step.friendlyPath, scalars, tables, dataFrames);
      } catch (e) {
        console.warn(`Run comparison: failed to load step run ${step.funcCallId}`, e);
      }
    }
  } else {
    extractCallNodes(call, '', '', scalars, tables, dataFrames);
  }

  const name = getRunTitle(call);
  return {
    id: call.id,
    name,
    sourceKind: serializedConfig ? 'workflow' : 'function',
    modelName,
    nodes: {entryId: call.id, entryName: name, scalars, tables},
    dataFrames,
  };
}

/** Builds a raw-data comparison entry from an open workspace table. */
export function entryFromDataFrame(df: DG.DataFrame): ComparisonEntry {
  const id = `table:${df.name}:${df.rowCount}`;
  const table: TableNodeInfo = {
    path: df.name,
    name: df.name,
    friendlyPath: df.name,
    columns: columnInfos(df),
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
