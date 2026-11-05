import {
  isDynamicPipelineState,
  isFuncCallState,
  isStaticPipelineState,
  PipelineInstanceRuntimeData,
  PipelineState,
  PipelineStateDynamic,
  PipelineStateStatic,
  StepFunCallState,
  ViewAction,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import type {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import {ConsistencyInfo, FuncCallStateInfo} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {
  anyInput, inputHasError, inputHasWarning, inputsWhere, isChange, isFixableInconsistency, isInconsistent,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/run-summary';
import type {AugmentedStat, Status} from './types';

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

// A workflow with `compactView` that resolves to one script through a chain of one-child static pipelines
// (e.g. a root that only refs another provider) is rendered by TreeWizard in compact mode.
type SinglePipelineChain = PipelineStateStatic<StepFunCallState, PipelineInstanceRuntimeData>[];

export function resolveSingleStep(
  state: PipelineState,
): {step: StepFunCallState, chain: SinglePipelineChain} | undefined {
  const chain: SinglePipelineChain = [];
  let current = state;
  while (!isFuncCallState(current)) {
    if (!isStaticPipelineState(current) || current.isActionStep || current.steps.length !== 1)
      return undefined;
    chain.push(current);
    current = current.steps[0];
  }
  return chain[0]?.compactView ? {step: current, chain} : undefined;
}

export type PipelineWithAdd = PipelineStateDynamic<StepFunCallState, PipelineInstanceRuntimeData>;

export const hasAddControls = (data: PipelineState): data is PipelineWithAdd =>
  isDynamicPipelineState(data) && !data.isReadonly &&
    data.stepTypes.filter((item) => !item.disableUIAdding).length > 0;

const isItemActionAllowed = (
  stat: AugmentedStat, flag: 'disableUIDragging' | 'disableUIRemoving' | 'disableUIAdding',
) => !!stat.parent && !stat.parent.data.isReadonly && isDynamicPipelineState(stat.parent.data) &&
  !stat.parent.data.stepTypes.some((item) => item.configId === stat.data.configId && item[flag]);

export const isEachDraggable = (stat: AugmentedStat) => isItemActionAllowed(stat, 'disableUIDragging');
export const isDeletable = (stat: AugmentedStat) => isItemActionAllowed(stat, 'disableUIRemoving');
export const isDuplicable = (stat: AugmentedStat) =>
  !stat.data.isReadonly && isItemActionAllowed(stat, 'disableUIAdding');

export const couldBeSaved = (data: PipelineState) => !isFuncCallState(data) && !!data.nqName && !data.disableHistory;

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
  const globalActions = states.flatMap((state) =>
    state?.actions?.filter((action) => action.position === 'globalmenu') ?? []);
  return globalActions;
}

export const hasInconsistencies = (consistencyStates?: Record<string, ConsistencyInfo>) =>
  anyInput(consistencyStates, isFixableInconsistency);

export const hasAnyInconsistency = (consistencyStates?: Record<string, ConsistencyInfo>) =>
  anyInput(consistencyStates, isInconsistent);

type SubtreeStates = [
  data: PipelineState,
  callStates: Record<string, FuncCallStateInfo | undefined>,
  consistencyStates: Record<string, Record<string, ConsistencyInfo> | undefined>,
];

// the first writable step without pending dependencies that has an input matching cond
const findInconsistentStep = (
  [data, callStates, consistencyStates]: SubtreeStates, cond: (val: ConsistencyInfo) => boolean,
) =>
  _findTreeNode([data], (state: PipelineState) => isFuncCallState(state) && !state.isReadonly &&
    anyInput(consistencyStates[state.uuid], cond) && !callStates[state.uuid]?.pendingDependencies?.length);

export const hasSubtreeFixableInconsistencies = (...states: SubtreeStates) =>
  findInconsistentStep(states, isFixableInconsistency);

export const hasSubtreeAnyInconsistencies = (...states: SubtreeStates) => findInconsistentStep(states, isInconsistent);

export const statusToTooltip: Record<Status, string> = {
  [`next`]: `This step is available to run`,
  [`next warn`]: `This step is available to run, but has warnings`,
  [`next error`]: `This step needs user input`,
  ['pending']: 'This step has pending dependencies',
  ['pending executed']: 'This step has changed dependencies',
  ['running']: 'This step is running',
  ['succeeded']: 'This step is succeeded',
  ['succeeded info']: 'This step is succeeded with changes',
  ['succeeded warn']: 'This step is succeeded, but has warnings',
  ['succeeded inconsistent']: 'This step is succeeded, but has inconsistent inputs',
  ['failed']: 'Run failed',
};

/** The inputs behind a step status. */
export const statusInputs = (
  status: Status,
  validationsState?: Record<string, ValidationResult>,
  consistencyStates?: Record<string, ConsistencyInfo>,
): string[] => {
  switch (status) {
  case 'next error':
    return inputsWhere(validationsState, inputHasError);
  case 'next warn':
    return [...new Set([...inputsWhere(validationsState, inputHasWarning),
      ...inputsWhere(consistencyStates, isFixableInconsistency)])];
  case 'succeeded warn':
    return inputsWhere(validationsState, (val) => inputHasWarning(val) || inputHasError(val));
  case 'succeeded info':
    return inputsWhere(consistencyStates, isChange);
  case 'succeeded inconsistent':
    // everything the Update action resets, info inputs included
    return inputsWhere(consistencyStates, isInconsistent);
  default:
    return [];
  }
};

export const SELECTED_STEP_BACKGROUND = 'var(--grey-1)';
