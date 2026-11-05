import * as DG from 'datagrok-api/dg';
import {historyUtils} from '../../../history-utils';
import {deepCopy} from '../../../shared-utils/utils';
import {
  PipelineSerializedState, PipelineState, PipelineStateRec, StepFunCallSerializedState,
  isFuncCallSerializedState, isFuncCallState,
} from '../config/PipelineInstance';
import {ItemMetadata} from '../view/ViewCommunication';
import {RestrictionState} from './FuncCallInstancesBridge';
import {ConsistencyInfo, FuncCallStateInfo} from './StateTreeNodes';
import {
  RUN_ERROR_PATH, readInstanceState, readOutputOutdated, saveCallState, saveInstanceState,
} from './funccall-utils';
import {copyDataFrame} from '../utils';

/** A saved step state with the call it holds; restrictions are only taken from live state. */
export type StepSnapshot = StepFunCallSerializedState & {
  funcCall?: DG.FuncCall;
  isOutputOutdated: boolean;
  runError?: string;
  restrictions?: Record<string, RestrictionState | undefined>;
  loadError?: unknown;
};

/** A workflow or step state that holds its calls: a copy of live state, or a saved run with its step runs. */
export type StateSnapshot = PipelineStateRec<StepSnapshot, {}>;

const copyValue = (value: any) => value instanceof DG.DataFrame ? copyDataFrame(value) : value;

/** A step snapshot of a call; the call is copied, so later edits don't reach it. */
export function snapshotCall(
  funcCall: DG.FuncCall,
  status: {isOutputOutdated: boolean, runError?: string, restrictions?: Record<string, RestrictionState | undefined>},
  step: Partial<StepFunCallSerializedState> = {},
): StepSnapshot {
  return {
    type: 'funccall', uuid: '', configId: funcCall.func.name, isReadonly: false, ...step,
    nqName: funcCall.func.nqName, funcCallId: funcCall.id, funcCall: deepCopy(funcCall), ...status,
  };
}

/** A copy of live view state, with the run status and restrictions the driver publishes for its steps. */
export function snapshotState(
  state: PipelineState,
  calls: Record<string, FuncCallStateInfo | undefined>,
  consistency: Record<string, Record<string, ConsistencyInfo> | undefined>,
): StateSnapshot {
  if (isFuncCallState(state)) {
    const {type, uuid, configId, friendlyName, isReadonly, sourceRunId, funcCall} = state;
    const step = {type, uuid, configId, friendlyName, isReadonly, sourceRunId};
    if (!funcCall)
      return {...step, nqName: '', isOutputOutdated: true};
    const restrictions = Object.fromEntries(Object.entries(consistency[uuid] ?? {}).map(
      ([name, {restriction, assignedValue}]) => [name, {type: restriction, assignedValue: copyValue(assignedValue)}]));
    const {isOutputOutdated = false, runError} = calls[uuid] ?? {};
    return snapshotCall(funcCall, {isOutputOutdated, runError, restrictions}, step);
  }
  const {uuid, configId, isReadonly, friendlyName, description, version, nqName, sourceRunId} = state;
  const base = {uuid, configId, isReadonly, friendlyName, description, version, nqName, sourceRunId};
  const steps = state.steps.map((step: PipelineState) => snapshotState(step, calls, consistency));
  return (state.type === 'static' ?
    {...base, type: state.type, steps} :
    {...base, type: state.type, stepTypes: state.stepTypes, steps}) as StateSnapshot;
}

async function loadSnapshotRec(state: PipelineSerializedState): Promise<StateSnapshot> {
  if (isFuncCallSerializedState(state)) {
    if (!state.funcCallId)
      return {...state, isOutputOutdated: true};
    try {
      const funcCall = await historyUtils.loadRun(state.funcCallId);
      const runError = funcCall.options[RUN_ERROR_PATH];
      return {...state, funcCall, isOutputOutdated: readOutputOutdated(funcCall), runError};
    } catch (loadError) {
      return {...state, isOutputOutdated: true, loadError};
    }
  }
  const steps = await Promise.all(state.steps.map((step) => loadSnapshotRec(step)));
  return {...state, steps} as StateSnapshot;
}

/** The saved workflow of a meta-call with its step runs, undefined for a plain run; a step run that fails to load
 * keeps the error in `loadError`. */
export async function loadStateSnapshot(metaCall: DG.FuncCall): Promise<StateSnapshot | undefined> {
  const state = readInstanceState(metaCall);
  return state ? loadSnapshotRec(state) : undefined;
}

async function saveStep(step: StepSnapshot, metadata?: ItemMetadata, nested = false):
  Promise<StepFunCallSerializedState> {
  const {funcCall, isOutputOutdated, runError, restrictions, loadError: _loadError, ...state} = step;
  if (!funcCall)
    throw new Error(`Step ${step.configId} has no call to save`);
  const saved = await saveCallState(funcCall, restrictions ?? {}, isOutputOutdated, runError, metadata, nested);
  return {...state, uuid: '', funcCallId: saved.id};
}

async function saveSnapshotRec(state: StateSnapshot): Promise<PipelineSerializedState> {
  if (isFuncCallState(state))
    return saveStep(state, undefined, true);
  const steps: PipelineSerializedState[] = [];
  for (const step of state.steps)
    steps.push(await saveSnapshotRec(step));
  return {...state, uuid: '', steps} as PipelineSerializedState;
}

/** Saves a snapshot like the driver saves a tree: every step call as a new run, then the workflow meta-call.
 * A step snapshot saves its call only. Resolves the saved meta-call or step call. */
export async function saveStateSnapshot(snapshot: StateSnapshot, metadata: ItemMetadata = {}): Promise<DG.FuncCall> {
  if (isFuncCallState(snapshot)) {
    await saveStep(snapshot, metadata);
    return snapshot.funcCall!;
  }
  if (!snapshot.nqName)
    throw new Error(`Attempting to save pipeline with no nqName`);
  return saveInstanceState(snapshot.nqName, await saveSnapshotRec(snapshot), metadata, snapshot.version);
}
