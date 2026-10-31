// Loads runs from the platform, or snapshots live state, and flattens them into ComparisonEntry
// objects: the only file in the comparison pipeline with server calls.

import * as DG from 'datagrok-api/dg';
import {historyUtils} from '@datagrok-libraries/compute-utils/history-utils/src/history-utils';
import {getRunTitle} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {tryParseJson} from '@datagrok-libraries/utils/src/string-utils';
import {readOutputOutdated} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {
  StateSnapshot, loadStateSnapshot, saveStateSnapshot, snapshotCall,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/state-snapshot';
import {
  PipelineStateRec, StepFunCallStateBase, isFuncCallState,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {buildTraverseD} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/graph-traverse-utils';
import type {ItemMetadata} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/view/ViewCommunication';
import {STEP_HISTORY_OPTION} from '../RunStorage/run-storage';
import {
  ComparisonEntry, EntrySnapshot, EntrySourceKind, ScalarNodeInfo, TableNodeInfo, AxisMode, TimeUnit, TIME_UNITS,
} from './types';

const SCALAR_PROPERTY_TYPES = new Set<string>([DG.TYPE.INT, DG.TYPE.FLOAT, DG.TYPE.BIG_INT]);

const AXIS_MODES = new Set<AxisMode>(['series', 'timeseries', 'points']);

export interface ComparisonDefaults {
  index?: string;
  split?: string;
  mode?: AxisMode;
  units?: TimeUnit;
}

/** Index defaults of a dataframe IO: the {comparison: {...}} JSON annotation. */
export function parseComparisonDefaults(options?: Record<string, any>): ComparisonDefaults {
  const raw = options?.['comparison'];
  // 'null' parses to null, property access below must survive it
  const parsed: any = (raw ? tryParseJson(raw, {warning: 'Run comparison: comparison annotation'}) : null) ?? {};
  return {
    index: typeof parsed.index === 'string' ? parsed.index : undefined,
    split: typeof parsed.split === 'string' ? parsed.split : undefined,
    ...AXIS_MODES.has(parsed.mode) ? {mode: parsed.mode as AxisMode} : {},
    ...TIME_UNITS.includes(parsed.units) ? {units: parsed.units as TimeUnit} : {},
  };
}

const columnInfos = (df: DG.DataFrame) => [...df.columns].map((col) => ({
  name: col.name,
  type: col.type,
  units: col.meta?.units || undefined,
}));

interface CallNodes {
  scalars: ScalarNodeInfo[];
  tables: TableNodeInfo[];
  dataFrames: Map<string, DG.DataFrame>;
}

function extractCallNodes(
  call: DG.FuncCall,
  pathPrefix: string,
  friendlyPrefix: string,
  isOutputOutdated: boolean,
): CallNodes {
  const scalars: ScalarNodeInfo[] = [];
  const tables: TableNodeInfo[] = [];
  const dataFrames = new Map<string, DG.DataFrame>();
  const io = [
    ...[...call.inputParams.values()].map((p) => ({param: p, value: call.inputs[p.property.name], outdated: false})),
    ...[...call.outputParams.values()].map((p) =>
      ({param: p, value: call.outputs[p.property.name], outdated: isOutputOutdated})),
  ].filter(({param, value}) => SCALAR_PROPERTY_TYPES.has(param.property.propertyType) ||
    (param.property.propertyType === DG.TYPE.DATA_FRAME && value != null));
  // captions are free text, so one call's input and output can share a display name;
  // property names are unique per call and disambiguate everything downstream
  const nameCounts = new Map<string, number>();
  for (const {param} of io) {
    const caption = param.property.caption ?? param.property.name;
    nameCounts.set(caption, (nameCounts.get(caption) ?? 0) + 1);
  }
  for (const {param, value, outdated} of io) {
    const prop = param.property;
    const path = pathPrefix ? `${pathPrefix}/${prop.name}` : prop.name;
    const caption = prop.caption ?? prop.name;
    const name = nameCounts.get(caption)! > 1 ? `${caption} (${prop.name})` : caption;
    const friendlyPath = friendlyPrefix ? `${friendlyPrefix} · ${name}` : name;
    if (SCALAR_PROPERTY_TYPES.has(prop.propertyType)) {
      scalars.push({
        path,
        name,
        friendlyPath,
        valueType: prop.propertyType,
        units: prop.options?.['units'] || undefined,
        value: (value == null || isNaN(Number(value))) ? null : Number(value),
        ...outdated ? {outdated} : {},
      });
    } else if (prop.propertyType === DG.TYPE.DATA_FRAME && value != null) {
      const df = value as DG.DataFrame;
      const defaults = parseComparisonDefaults(prop.options);
      tables.push({
        path,
        name,
        friendlyPath,
        nqName: call.func?.nqName,
        defaultIndexColumn: defaults.index,
        defaultSplitColumn: defaults.split,
        defaultAxisMode: defaults.mode,
        defaultTimeUnits: defaults.units,
        columns: columnInfos(df),
        rowCount: df.rowCount,
        ...outdated ? {outdated} : {},
      });
      dataFrames.set(path, df);
    }
  }
  return {scalars, tables, dataFrames};
}

interface WorkflowStep<S> {
  path: string;
  friendlyPath: string;
  step: S;
}

// a node's own segment included; repeated steps (dynamic items) are told apart by position
interface StepPath {
  ids: string[];
  names: string[];
}

/** Leaf steps of a saved or a live workflow state, with paths shared by both. */
export function collectWorkflowSteps<S extends StepFunCallStateBase, T>(
  state: PipelineStateRec<S, T>,
): WorkflowStep<S>[] {
  const traverse = buildTraverseD<PipelineStateRec<S, T>, StepPath, boolean>(
    {ids: [], names: []},
    (item, path, isRoot) => {
      if (isFuncCallState(item))
        return [];
      // the root pipeline name is shared by every step, so it is left out of friendly paths
      const base = isRoot ? {ids: [item.configId], names: []} : path;
      const steps: PipelineStateRec<S, T>[] = item.steps;
      return steps.map((step, i) => {
        const name = step.friendlyName ?? step.configId;
        const repeats = steps.filter((other) => other.configId === step.configId).length > 1;
        const position = steps.slice(0, i + 1).filter((other) => other.configId === step.configId).length;
        return [step, {
          ids: [...base.ids, repeats ? `${step.configId}[${position}]` : step.configId],
          names: [...base.names, repeats ? `${name} ${position}` : name],
        }, false] as const;
      });
    }, true);
  return traverse(state, (acc: WorkflowStep<S>[], item, path) => {
    if (isFuncCallState(item)) {
      acc.push({
        path: path.ids.length ? path.ids.join('/') : item.configId,
        friendlyPath: path.names.length ? path.names.join(' · ') : item.friendlyName ?? item.configId,
        step: item,
      });
    }
    return acc;
  }, []);
}

const mergeNodes = (parts: CallNodes[]): CallNodes => ({
  scalars: parts.flatMap((part) => part.scalars),
  tables: parts.flatMap((part) => part.tables),
  dataFrames: new Map(parts.flatMap((part) => [...part.dataFrames])),
});

const callModelName = (call: DG.FuncCall) => call.func?.friendlyName ?? call.func?.name ?? '';

function makeEntry(
  id: string, name: string, sourceKind: EntrySourceKind, modelName: string, nodes: CallNodes,
  lists: {failedSteps?: string[], outdatedSteps?: string[]} = {}, snapshot?: EntrySnapshot,
): ComparisonEntry {
  return {
    id,
    name,
    sourceKind,
    modelName,
    nodes: {entryId: id, entryName: name, sourceKind, scalars: nodes.scalars, tables: nodes.tables},
    dataFrames: nodes.dataFrames,
    ...lists.failedSteps?.length ? {failedSteps: lists.failedSteps} : {},
    ...lists.outdatedSteps?.length ? {outdatedSteps: lists.outdatedSteps} : {},
    ...snapshot ? {snapshot} : {},
  };
}

function callEntry(id: string, name: string, call: DG.FuncCall, isOutputOutdated: boolean, snapshot?: EntrySnapshot) {
  const modelName = callModelName(call);
  return makeEntry(id, name, 'function', modelName, extractCallNodes(call, '', '', isOutputOutdated),
    {outdatedSteps: isOutputOutdated ? [modelName] : []}, snapshot);
}

// steps without a loaded call are left out, the ones that failed to load are listed
function workflowEntry(id: string, name: string, modelName: string, state: StateSnapshot, snapshot?: EntrySnapshot) {
  const steps = collectWorkflowSteps(state);
  const loaded = steps.filter(({step}) => step.funcCall);
  const failed = steps.filter(({step}) => step.loadError);
  for (const {step} of failed)
    console.warn(`Run comparison: failed to load step run ${step.funcCallId}`, step.loadError);
  const nodes = mergeNodes(loaded.map(({path, friendlyPath, step}) =>
    extractCallNodes(step.funcCall!, path, friendlyPath, step.isOutputOutdated)));
  const outdatedSteps = loaded.filter(({step}) => step.isOutputOutdated).map(({friendlyPath}) => friendlyPath);
  return makeEntry(id, name, 'workflow', modelName, nodes,
    {failedSteps: failed.map(({friendlyPath}) => friendlyPath), outdatedSteps}, snapshot);
}

/** Builds a comparison entry from a loaded run: a workflow meta-call is flattened into its steps. */
export async function entryFromFuncCall(call: DG.FuncCall): Promise<ComparisonEntry> {
  const name = getRunTitle(call);
  const modelName = callModelName(call);
  let state: StateSnapshot | undefined;
  let isWorkflow = true;
  try {
    state = await loadStateSnapshot(call);
    isWorkflow = !!state;
  } catch (e) {
    // legacy or malformed serialized configs degrade to an entry with no steps
    console.warn('Run comparison: failed to parse workflow config', e);
  }
  if (!isWorkflow)
    return callEntry(call.id, name, call, readOutputOutdated(call));
  return state ?
    workflowEntry(call.id, name, modelName, state) :
    makeEntry(call.id, name, 'workflow', modelName, mergeNodes([]));
}

/** Loads runs by id; runs that fail to load are reported by id instead of thrown. */
export async function loadRunEntries(ids: string[]) {
  const results = await Promise.allSettled(ids.map(async (id) => entryFromFuncCall(await historyUtils.loadRun(id))));
  const entries: ComparisonEntry[] = [];
  const failed: string[] = [];
  results.forEach((result, i) => {
    if (result.status === 'fulfilled') {
      entries.push(result.value);
      return;
    }
    console.error(`Run comparison: failed to load run ${ids[i]}`, result.reason);
    failed.push(ids[i]);
  });
  return {entries, failed};
}

let snapshotCounter = 0;

/** An unsaved entry of a snapshot, with the same paths as its saved run: an edited step takes the paths of the plain
 * run it is based on. A step saves to the step history. */
export function entryFromSnapshot(state: StateSnapshot, name: string, {
  metadata = {}, stepHistory = false, saveable = true,
}: {metadata?: ItemMetadata, stepHistory?: boolean, saveable?: boolean} = {}): ComparisonEntry {
  const id = `snapshot:${++snapshotCounter}`;
  const save = async (edited: ItemMetadata) => {
    if (stepHistory && isFuncCallState(state))
      state.funcCall!.options[STEP_HISTORY_OPTION] = 'true';
    return saveStateSnapshot(state, edited);
  };
  if (isFuncCallState(state)) {
    return state.funcCall ?
      callEntry(id, name, state.funcCall, state.isOutputOutdated,
        {metadata, save: saveable ? save : undefined, model: state.funcCall.func}) :
      makeEntry(id, name, 'function', '', mergeNodes([]), {}, {metadata});
  }
  return workflowEntry(id, name, state.friendlyName ?? state.nqName ?? '', state,
    {metadata, save: state.nqName ? save : undefined, model: state.nqName ? DG.Func.byName(state.nqName) : undefined});
}

/** Snapshot of a standalone call; like the RFV's own save, outdated or failed results can't be saved. */
export const entryFromLiveCall = (
  call: DG.FuncCall, {isOutputOutdated, runError}: {isOutputOutdated: boolean, runError?: string}, name: string,
  metadata?: ItemMetadata,
) => entryFromSnapshot(snapshotCall(call, {isOutputOutdated, runError}), name,
  {metadata, saveable: !isOutputOutdated && !runError});

// same underlying table -> same entry id (so re-adding dedupes); keyed by the Dart
// handle because names and row counts alone can collide across distinct tables
const rawEntryIds = new WeakMap<object, string>();
let rawEntryCounter = 0;

/** Builds a raw-data comparison entry from an open workspace table. */
export function entryFromDataFrame(df: DG.DataFrame): ComparisonEntry {
  let id = rawEntryIds.get(df.dart);
  if (!id) {
    id = `table:${df.name}:${++rawEntryCounter}`;
    rawEntryIds.set(df.dart, id);
  }
  const table: TableNodeInfo = {
    path: df.name,
    name: df.name,
    friendlyPath: df.name,
    columns: columnInfos(df),
    rowCount: df.rowCount,
  };
  return makeEntry(id, df.name, 'raw', '', {scalars: [], tables: [table], dataFrames: new Map([[df.name, df]])});
}
