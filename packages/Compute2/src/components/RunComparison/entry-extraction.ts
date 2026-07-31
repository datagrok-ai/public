// Loads runs from the platform and flattens them into ComparisonEntry objects:
// the only file in the comparison pipeline with server calls.

import * as DG from 'datagrok-api/dg';
import {historyUtils} from '@datagrok-libraries/compute-utils';
import {getRunTitle} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {deserialize} from '@datagrok-libraries/utils/src/json-serialization';
import {CONFIG_PATH} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {
  PipelineSerializedState, isFuncCallSerializedState,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {buildTraverseD} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/graph-traverse-utils';
import {ComparisonEntry, ScalarNodeInfo, TableNodeInfo} from './types';

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
        defaultIndexColumn: prop.options?.['comparisonIndex'] || undefined,
        defaultSplitColumn: prop.options?.['comparisonSplit'] || undefined,
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

interface StepPath {
  ids: string[];
  names: string[];
}

const traverseSteps = buildTraverseD<PipelineSerializedState, StepPath, boolean>(
  {ids: [], names: []},
  (state, path, isRoot) => {
    if (isFuncCallSerializedState(state))
      return [];
    // the root pipeline name is shared by every step, so it is left out of friendly paths
    const names = isRoot ? path.names : [...path.names, state.friendlyName ?? state.configId];
    return state.steps.map((step) =>
      [step, {ids: [...path.ids, state.configId], names}, false] as const);
  }, true);

function collectWorkflowSteps(state: PipelineSerializedState): WorkflowStep[] {
  return traverseSteps(state, (acc: WorkflowStep[], item, path) => {
    if (isFuncCallSerializedState(item) && item.funcCallId) {
      acc.push({
        path: [...path.ids, item.configId].join('/'),
        friendlyPath: [...path.names, item.friendlyName ?? item.configId].join(' · '),
        funcCallId: item.funcCallId,
      });
    }
    return acc;
  }, []);
}

/** Builds a comparison entry from a loaded run: a workflow meta-call is flattened into its steps. */
export async function entryFromFuncCall(call: DG.FuncCall): Promise<ComparisonEntry> {
  const scalars: ScalarNodeInfo[] = [];
  const tables: TableNodeInfo[] = [];
  const dataFrames = new Map<string, DG.DataFrame>();
  const serializedConfig = call.options?.[CONFIG_PATH];
  const modelName = call.func?.friendlyName ?? call.func?.name ?? '';

  if (serializedConfig) {
    let steps: WorkflowStep[] = [];
    try {
      steps = collectWorkflowSteps(deserialize(serializedConfig));
    } catch (e) {
      // legacy or malformed serialized configs degrade to an entry with no steps
      console.warn('Run comparison: failed to parse workflow config', e);
    }
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
