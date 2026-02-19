import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {
  isFuncCallState,
  isParallelPipelineState,
  isSequentialPipelineState,
  isStaticPipelineState,
  PipelineInstanceRuntimeData,
  PipelineState,
  PipelineStateParallel,
  PipelineStateSequential,
  StepFunCallState,
  ViewAction,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {zipSync, Zippable} from 'fflate';
import {dfToViewerMapping, getStartedOrNull, replaceForWindowsPath, richFunctionViewReport, ValidationResult} from '@datagrok-libraries/compute-utils';
import {ConsistencyInfo, FuncCallStateInfo, MetaCallInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import type Dayjs from 'dayjs';
import {ExportCbInput} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';

export type NodeWithPath = {
  state: PipelineState,
  pathSegments: number[],
}

function hasSteps(state: PipelineState) {
  return isParallelPipelineState(state) ||
  isSequentialPipelineState(state) ||
  isStaticPipelineState(state);
}

function _findTreeNode(
  steps: PipelineState[],
  pred: (state: PipelineState) => boolean,
  pathSegments: number[] = [],
): NodeWithPath | undefined {
  for (const [idx, stepState] of steps.entries()) {
    const currentPathSegments = [...pathSegments, idx];
    if (pred(stepState))
      return {state: stepState, pathSegments: currentPathSegments};

    if (hasSteps(stepState)) {
      const t = _findTreeNode(stepState.steps, pred, currentPathSegments);
      if (t)
        return t;
    }
  }
};

export function findNodeWithPathByUuid(uuid: string, state: PipelineState): NodeWithPath | undefined {
  return _findTreeNode([state], (state: PipelineState) => state.uuid === uuid);
};

export function findTreeNodeByPath(pathSegments: number[], state: PipelineState): NodeWithPath | undefined {
  const node = pathSegments.slice(1).reduce((acc, segment) => {
    if (acc && hasSteps(acc)) {
      acc = acc.steps[segment];

      return acc;
    }
    return undefined;
  }, state as PipelineState | undefined);

  return node ? {
    state: node,
    pathSegments,
  }: undefined;
};

export function findTreeNodeParrent(uuid: string, state: PipelineState): PipelineState | undefined {
  const notVisitedStates = [state];

  while (notVisitedStates.length > 0) {
    const currentState = notVisitedStates.pop()!;

    if (
      isParallelPipelineState(currentState) ||
        isSequentialPipelineState(currentState) ||
        isStaticPipelineState(currentState)
    ) {
      for (const item of currentState.steps) {
        if (item.uuid == uuid)
          return currentState;
      }
      notVisitedStates.push(...currentState.steps);
    }
  }
};

const suitableForNavStep = (state: PipelineState) => {
  return (state.type === 'funccall' || state.steps.length === 0 || state.forceNavigate) && !state.isReadonly;
};

export function findPrevStep(uuid: string, state: PipelineState): NodeWithPath | undefined {
  let prevUuid = '';
  const pred = (state: PipelineState) => {
    if (!suitableForNavStep(state))
      return false;
    if (state.uuid === uuid)
      return true;
    prevUuid = state.uuid;
    return false;
  };
  return _findTreeNode([state], pred) ? findNodeWithPathByUuid(prevUuid, state) : undefined;
}

export function findNextStep(uuid: string, state: PipelineState): NodeWithPath | undefined {
  let prevUuid = '';
  const pred = (state: PipelineState) => {
    if (!suitableForNavStep(state))
      return false;
    if (prevUuid === uuid)
      return true;
    prevUuid = state.uuid;
    return false;
  };
  return _findTreeNode([state], pred);
}

export function findNextSubStep(state: PipelineState): NodeWithPath | undefined {
  return _findTreeNode([state], suitableForNavStep);
}

export type PipelineWithAdd = PipelineStateSequential<StepFunCallState, PipelineInstanceRuntimeData> |
PipelineStateParallel<StepFunCallState, PipelineInstanceRuntimeData>;

export const hasRunnableSteps = (data: PipelineState) =>
  (isParallelPipelineState(data) || isSequentialPipelineState(data)) && !data.isReadonly && data.steps.length > 0;

export const hasAddControls = (data: PipelineState): data is PipelineWithAdd =>
  (isParallelPipelineState(data) || isSequentialPipelineState(data)) && !data.isReadonly &&
    data.stepTypes.filter((item) => !item.disableUIAdding).length > 0;

export const couldBeSaved = (data: PipelineState) => !isFuncCallState(data) && !!data.nqName && !data.disableHistory;

export const hasSubtreeFixableInconsistencies = (
  data: PipelineState,
  callStates: Record<string, FuncCallStateInfo | undefined>,
  consistencyStates: Record<string, Record<string, ConsistencyInfo> | undefined>,
) => {
  return _findTreeNode(
    [data],
    (state: PipelineState) => isFuncCallState(state) ?
      (!state.isReadonly && hasInconsistencies(consistencyStates[state.uuid]) && !callStates[state.uuid]?.pendingDependencies?.length) :
      false,
  );
};

export function getRelevantGlobalActions(data: PipelineState, currentStepUuid: string) : ViewAction[] {
  const nodePaths = _findTreeNode([data], (state) => state.uuid === currentStepUuid);
  const segments = nodePaths?.pathSegments ?? [];
  const states: PipelineState[] = [data];
  for (let idx = 1, currentState = data; idx < segments.length; idx++) {
    if (isFuncCallState(currentState))
      break;
    currentState = currentState.steps[segments[idx]];
    states.push(currentState);
  }
  const globalActions = states.flatMap(state => state?.actions?.filter(action => action.position === 'globalmenu') ?? []);
  return globalActions;
}

export const hasInconsistencies = (consistencyStates?: Record<string, ConsistencyInfo>) => {
  const firstInconsistency = Object.values(consistencyStates || {}).find(
    (val) => val.inconsistent && (val.restriction === 'disabled' || val.restriction === 'restricted'));
  return !!firstInconsistency;
};

export async function reportTree(
  {
    startDownload,
    treeState,
    meta = {},
    callInfoStates,
    validationStates,
    consistencyStates,
    descriptions,
    hasNotSavedEdits,
    cb,
  }: {
    startDownload: boolean;
    treeState: PipelineState;
    meta?: MetaCallInfo;
    callInfoStates?: Record<string, FuncCallStateInfo | undefined>,
    validationStates?: Record<string, Record<string, ValidationResult> | undefined>,
    consistencyStates?: Record<string, Record<string, ConsistencyInfo> | undefined>,
    descriptions?: Record<string, Record<string, string | string[]> | undefined>,
    hasNotSavedEdits?: boolean;
    cb?: (input: ExportCbInput) => Promise<void>,
  }) {
  const zipConfig: Zippable = {};

  const q = [{ state: treeState, idx: 0, path: [] as string[] }];

  while (q.length > 0) {
    const { state, idx, path } = q.shift()!;
    if (isFuncCallState(state)) {
      const funcCall = state.funcCall;
      const callInfo = callInfoStates?.[state.uuid];
      if (!funcCall || !callInfo)
        continue;

      const validation = validationStates?.[state.uuid];
      const consistency = consistencyStates?.[state.uuid];
      const description = descriptions?.[state.uuid];
      const {isOutputOutdated, runError} = callInfo;

      const [blob, wb] = await richFunctionViewReport(
        'Excel',
        funcCall.func,
        funcCall,
        dfToViewerMapping(funcCall),
        validation,
        consistency,
      );

      const rawFileName = getExportName(state, isOutputOutdated, description?.title as string, getStartedOrNull(funcCall), runError);
      const fileName = `${String(idx + 1).padStart(3, '0')}_${replaceForWindowsPath(rawFileName)}.xlsx`;
      const configKey = [...path, fileName].join('/')
      zipConfig[configKey] = [new Uint8Array(await blob.arrayBuffer()), { level: 0 }];
      if (cb) {
        await cb({
          fc: funcCall,
          wb,
          archive: zipConfig,
          path,
          fileName,
          isOutputOutdated,
          runError,
          validation,
          consistency,
          description
        });
      }
    } else {
      const dirName = `${String(idx + 1).padStart(3, '0')}_${replaceForWindowsPath(state.friendlyName ?? state.nqName ?? '')}`;
      const nPath = state === treeState ? [] : [...path, dirName]
      for (const [idx, stepState] of state.steps.entries()) {
        q.push({ state: stepState, idx, path: nPath })
      }
    }
  }

  const rawFileName = getExportName(treeState, !!hasNotSavedEdits, meta.title, meta.started);
  const fileName = replaceForWindowsPath(`${rawFileName}.zip`);
  const blob = new Blob([zipSync(zipConfig) as any]);
  if (startDownload)
    DG.Utils.download(fileName, blob);
  return [blob, zipConfig, fileName] as const;
}

function getExportName(
  state: PipelineState,
  hasNotSavedEdits: boolean,
  title?: string,
  started?: Dayjs.Dayjs,
  runError?: string,
) {
  const stateName = state.friendlyName ?? state.configId;
  const name = title ? `${stateName} - ${title}` : stateName;
  const fileName = (started && !hasNotSavedEdits) ? `${name} - ${started}` : (runError ? `${name} - failed` : `${name} - edited`);
  return fileName;
}


export function setDifference<T>(a: Set<T>, b: Set<T>) {
  return new Set(Array.from(a).filter((item) => !b.has(item)));
}
