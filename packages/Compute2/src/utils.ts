import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {
  isDynamicPipelineState,
  isFuncCallState,
  isStaticPipelineState,
  PipelineInstanceRuntimeData,
  PipelineState,
  PipelineStateDynamic,
  StepFunCallState,
  ViewAction,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {zipSync, Zippable} from 'fflate';
import {dfToViewerMapping, getStartedOrNull, replaceForWindowsPath, richFunctionViewReport, ValidationResult} from '@datagrok-libraries/compute-utils';
import {getCustomExports} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {DEFAULT_FLOAT_FORMAT} from '@datagrok-libraries/webcomponents-vue';
import {ConsistencyInfo, FuncCallStateInfo, MetaCallInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import type Dayjs from 'dayjs';
import {ExportCbInput, ViewersHook} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {BehaviorSubject} from 'rxjs';

export type NodeWithPath = {
  state: PipelineState,
  pathSegments: number[],
}

function hasSteps(state: PipelineState) {
  return isDynamicPipelineState(state) ||
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

// Returns a uuid guaranteed to exist in `tree`. Keeps the current selection if
// it is still alive; otherwise climbs the old positional path to the nearest
// surviving ancestor (handles a removed last child resolving out of bounds);
// otherwise falls back to root. Never returns a dead uuid.
export function resolveChosenUuid(
  currentUuid: string | undefined,
  tree: PipelineState,
  fallbackPath?: number[],
): string {
  if (currentUuid && findNodeWithPathByUuid(currentUuid, tree))
    return currentUuid;
  let path = fallbackPath ? [...fallbackPath] : [];
  while (path.length > 1) {
    const byPath = findTreeNodeByPath(path, tree);
    if (byPath)
      return byPath.state.uuid;
    path = path.slice(0, -1);
  }
  return tree.uuid;
}

export function findTreeNodeParrent(uuid: string, state: PipelineState): PipelineState | undefined {
  const notVisitedStates = [state];

  while (notVisitedStates.length > 0) {
    const currentState = notVisitedStates.pop()!;

    if (
      isDynamicPipelineState(currentState) ||
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
    if (state.uuid === uuid)
      return true;
    if (suitableForNavStep(state))
      prevUuid = state.uuid;
    return false;
  };
  return _findTreeNode([state], pred) ? findNodeWithPathByUuid(prevUuid, state) : undefined;
}

export function findNextStep(uuid: string, state: PipelineState): NodeWithPath | undefined {
  let found = false;
  const pred = (state: PipelineState) => {
    if (found && suitableForNavStep(state))
      return true;
    if (state.uuid === uuid)
      found = true;
    return false;
  };
  return _findTreeNode([state], pred);
}

export function findNextSubStep(state: PipelineState): NodeWithPath | undefined {
  return _findTreeNode([state], suitableForNavStep);
}

export type PipelineWithAdd = PipelineStateDynamic<StepFunCallState, PipelineInstanceRuntimeData>;

export const hasRunnableSteps = (data: PipelineState) =>
  isDynamicPipelineState(data) && !data.isReadonly && data.steps.length > 0;

export const hasAddControls = (data: PipelineState): data is PipelineWithAdd =>
  isDynamicPipelineState(data) && !data.isReadonly &&
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

export const hasAnyInconsistency = (consistencyStates?: Record<string, ConsistencyInfo>) =>
  Object.values(consistencyStates || {}).some((v) => v.inconsistent);

export const hasSubtreeAnyInconsistencies = (
  data: PipelineState,
  callStates: Record<string, FuncCallStateInfo | undefined>,
  consistencyStates: Record<string, Record<string, ConsistencyInfo> | undefined>,
) => {
  return _findTreeNode(
    [data],
    (state: PipelineState) => isFuncCallState(state) ?
      (!state.isReadonly && hasAnyInconsistency(consistencyStates[state.uuid]) && !callStates[state.uuid]?.pendingDependencies?.length) :
      false,
  );
};

export async function getViewers(call: DG.FuncCall, viewersHook?: ViewersHook, metaState?: Record<string, BehaviorSubject<any>>) {
  const mappings = await dfToViewerMapping(call);
  if (viewersHook) {
    for (const [ioName, viewers] of Object.entries(mappings ?? {})) {
      for (const viewer of viewers) {
        if (!viewer)
          continue;
        const meta = metaState?.[ioName]?.value;
        viewersHook(ioName, viewer.type, viewer, meta);
      }
    }
  }
  return mappings;
}

export async function reportTree(
  {
    startDownload,
    treeState,
    meta = {},
    callInfoStates,
    metaStates,
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
    metaStates?: Record<string, Record<string, BehaviorSubject<any>> | undefined>,
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
      const metaState = metaStates?.[state.uuid];
      const {isOutputOutdated, runError} = callInfo;
      const viewers = await getViewers(funcCall, state.viewersHook, metaState);

      const [blob, wb] = await richFunctionViewReport(
        'Excel',
        funcCall.func,
        funcCall,
        viewers,
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
          meta: metaState ?? {},
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

/** Resolves the custom export named `exportName` declared on the funcCall's function
 *  (via `meta.customExports`) and applies it, passing the funcCall through. */
export async function applyCustomExport(
  fc: DG.FuncCall,
  exportName: string,
  args: Record<string, unknown> = {},
): Promise<unknown> {
  const item = getCustomExports(fc.func).find((x) => x.name === exportName);
  if (!item)
    throw new Error(`No export named ${exportName} is defined for ${fc.func.nqName}`);
  return DG.Func.byName(item.function).apply({funcCall: fc, ...args});
}

export function applyDefaultGridFloatFormat(viewer: DG.Viewer | undefined, type: string) {
  if (!viewer || type !== DG.VIEWER.GRID) return;
  const grid = viewer as DG.Grid;
  for (let i = 0; i < grid.columns.length; i++) {
    const gc = grid.columns.byIndex(i);
    const col = gc?.column;
    if (!col || col.type !== DG.COLUMN_TYPE.FLOAT) continue;
    if (gc!.format) continue;
    if (col.tags?.['format']) continue;
    gc!.format = DEFAULT_FLOAT_FORMAT;
  }
}

// Guards result-dependent actions (export/save) shared by RFV and RFVApp. Returns true when the
// action may proceed; otherwise surfaces a shell message explaining why the results aren't ready.
export function canUseResults(
  state: {isRunning?: boolean, runError?: unknown, isOutputOutdated?: boolean} | undefined,
  action: string,
): boolean {
  if (state?.isRunning) {
    grok.shell.warning(`The model is still running — wait for it to finish before ${action}.`);
    return false;
  }
  if (state?.runError) {
    grok.shell.error(`The last run finished with an error — fix the inputs and rerun before ${action}.`);
    return false;
  }
  if (state?.isOutputOutdated) {
    grok.shell.warning(`Results are outdated — run the model to update them before ${action}.`);
    return false;
  }
  return true;
}
