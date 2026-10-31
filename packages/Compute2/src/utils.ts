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
  PipelineStateStatic,
  StepFunCallState,
  ViewAction,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {zipSync, Zippable} from 'fflate';
import {dfToViewerMapping, richFunctionViewReport} from '@datagrok-libraries/compute-utils/shared-utils/export';
import {getStartedOrNull, replaceForWindowsPath} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import type {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import type ExcelJS from 'exceljs';
import {getCustomExports} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {DEFAULT_FLOAT_FORMAT} from '@datagrok-libraries/webcomponents-vue';
import {ConsistencyInfo, FuncCallStateInfo, MetaCallInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import type Dayjs from 'dayjs';
import {ExportCbInput, ExportSummaryItem, ViewersHook} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {
  anyInput, getRunSummary, inputHasError, inputHasWarning, inputsWhere, isChange, isFixableInconsistency,
  isInconsistent, RunStates, statesToStatus, textTitle,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/run-summary';
import type {AugmentedStat, Status} from './components/TreeWizard/types';
import {BehaviorSubject} from 'rxjs';
import * as Vue from 'vue';

declare let ENABLE_VUE_DEV_TOOLS: any;

export function setVueAppOptions(app: Vue.App<any>) {
  app.config.compilerOptions.isCustomElement = (tag) => tag.startsWith('dg-') || tag === 'dock-spawn-ts';
  if (ENABLE_VUE_DEV_TOOLS)
    app.config.performance = true;
}

// fullScreen lifts the platform's cap on dialog content height, so a large dialog fills its size
export type DialogSize = {width: number, height: number, fullScreen?: boolean};

/** Mounts a Vue component into a modal platform dialog and resolves what the component reports. */
export function showVueDialog<T>(
  title: string,
  component: Vue.Component,
  props: Record<string, any>,
  size: DialogSize,
): Promise<T | null> {
  return new Promise((resolve) => {
    const host = ui.div([], {style: {width: '100%', height: '100%'}});
    // the component renders its own buttons
    const dialog = ui.dialog({title, showFooter: false});
    let settled = false;
    const app = Vue.createApp(component, {
      ...props,
      onDone: (result: T | null) => {
        settled = true;
        resolve(result);
        dialog.close();
      },
    });
    setVueAppOptions(app);
    dialog.add(host);
    dialog.onClose.subscribe(() => {
      app.unmount();
      if (!settled)
        resolve(null);
    });
    dialog.show({modal: true, resizable: true, ...size, center: true});
    app.mount(host);
  });
}

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
  const globalActions = states.flatMap(state => state?.actions?.filter(action => action.position === 'globalmenu') ?? []);
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

// export-time viewers are created ad hoc and never mounted; detach releases their dart side
export function disposeViewers(mapping: {[key: string]: (DG.Viewer | undefined)[]} | undefined) {
  for (const viewers of Object.values(mapping ?? {})) {
    for (const viewer of viewers)
      viewer?.detach();
  }
}

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

export const SUMMARY_FILE_NAME = '000_summary.xlsx';

const XLSX_BLOB_TYPE = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;charset=UTF-8';

const exportIndex = (idx: number) => String(idx + 1).padStart(3, '0');

const stepFileName = (state: StepFunCallState, idx: number, callInfo: FuncCallStateInfo, title?: string) => {
  const rawFileName = getExportName(
    state, callInfo.isOutputOutdated, title, getStartedOrNull(state.funcCall!), callInfo.runError);
  return `${exportIndex(idx)}_${replaceForWindowsPath(rawFileName)}.xlsx`;
};

const workflowDirName = (state: Exclude<PipelineState, StepFunCallState>, idx: number) =>
  `${exportIndex(idx)}_${replaceForWindowsPath(state.friendlyName ?? state.nqName ?? '')}`;

export function getExportSummary(treeState: PipelineState, states: RunStates): ExportSummaryItem[] {
  const files = new Map<string, Pick<ExportSummaryItem, 'path' | 'fileName'>>();
  // a step whose run failed to load has a call state but no call
  const callInfoStates = {...states.callInfoStates};
  const visit = (state: PipelineState, idx: number, path: string[]) => {
    if (isFuncCallState(state)) {
      const callInfo = callInfoStates[state.uuid];
      if (!state.funcCall)
        delete callInfoStates[state.uuid];
      const fileName = state.funcCall && callInfo ?
        stepFileName(state, idx, callInfo, textTitle(states.descriptions?.[state.uuid])) : undefined;
      files.set(state.uuid, {path, fileName});
      return;
    }
    const nPath = state === treeState ? [] : [...path, workflowDirName(state, idx)];
    files.set(state.uuid, {path: nPath});
    state.steps.forEach((step, stepIdx) => visit(step, stepIdx, nPath));
  };
  visit(treeState, 0, []);
  return getRunSummary(treeState, {...states, callInfoStates}).map((item) => ({...item, ...files.get(item.uuid)!}));
}

const SUMMARY_COLUMNS = [
  'Path', 'Type', 'Name', 'File', 'Status', 'Run error', 'Errors', 'Warnings', 'Notifications', 'Inconsistent inputs',
  'Steps', 'Not loaded', 'Failed', 'Outdated', 'With errors', 'With warnings', 'Inconsistent',
];

const summaryStatusText = (item: ExportSummaryItem) => {
  if (item.status)
    return statusToTooltip[item.status];
  return item.kind === 'step' ? 'Not loaded' : '';
};

export async function reportSummary(items: ExportSummaryItem[]) {
  await DG.Utils.loadJsCss(['/js/common/exceljs.min.js']);
  //@ts-ignore
  const wb = new window.ExcelJS.Workbook() as ExcelJS.Workbook;
  const sheet = wb.addWorksheet('Summary');
  sheet.addRow(SUMMARY_COLUMNS).font = {bold: true};
  for (const item of items) {
    const rollup = item.rollup;
    sheet.addRow([
      item.path.join('/'),
      item.kind === 'step' ? 'Step' : 'Workflow',
      item.title ? `${item.name} - ${item.title}` : item.name,
      item.fileName ?? '',
      summaryStatusText(item),
      item.runError ?? '',
      item.errors.join('\n'),
      item.warnings.join('\n'),
      item.notifications.join('\n'),
      item.inconsistentInputs.join(', '),
      ...(rollup ? [
        rollup.steps, rollup.notLoaded, rollup.failed, rollup.outdated,
        rollup.withErrors, rollup.withWarnings, rollup.inconsistent,
      ] : []),
    ]);
  }
  sheet.columns.forEach((column, idx) => {
    column.width = idx < 10 ? 30 : 12;
    column.alignment = {wrapText: true, vertical: 'top'};
  });
  const buffer = await wb.xlsx.writeBuffer();
  return [new Blob([buffer], {type: XLSX_BLOB_TYPE}), wb] as const;
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
    pipelineValidations,
    descriptions,
    hasNotSavedEdits,
    cb,
  }: {
    startDownload: boolean;
    treeState: PipelineState;
    meta?: MetaCallInfo;
    metaStates?: Record<string, Record<string, BehaviorSubject<any>> | undefined>,
    hasNotSavedEdits?: boolean;
    cb?: (input: ExportCbInput) => Promise<void>,
  } & RunStates) {
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
      ).finally(() => disposeViewers(viewers));

      const fileName = stepFileName(state, idx, callInfo, textTitle(description));
      const configKey = [...path, fileName].join('/')
      zipConfig[configKey] = [new Uint8Array(await blob.arrayBuffer()), { level: 0 }];
      if (cb) {
        await cb({
          fc: funcCall,
          wb,
          archive: zipConfig,
          path,
          fileName,
          status: statesToStatus(callInfo, validation, consistency),
          isOutputOutdated,
          runError,
          validation,
          consistency,
          meta: metaState ?? {},
          description
        });
      }
    } else {
      const nPath = state === treeState ? [] : [...path, workflowDirName(state, idx)];
      for (const [idx, stepState] of state.steps.entries()) {
        q.push({ state: stepState, idx, path: nPath })
      }
    }
  }

  const summary = getExportSummary(
    treeState, {callInfoStates, validationStates, consistencyStates, pipelineValidations, descriptions});
  const [summaryBlob] = await reportSummary(summary);
  zipConfig[SUMMARY_FILE_NAME] = [new Uint8Array(await summaryBlob.arrayBuffer()), {level: 0}];

  const rawFileName = getExportName(treeState, !!hasNotSavedEdits, meta.title, meta.started);
  const fileName = replaceForWindowsPath(`${rawFileName}.zip`);
  const blob = new Blob([zipSync(zipConfig) as any]);
  if (startDownload)
    DG.Utils.download(fileName, blob);
  return [blob, zipConfig, fileName, summary] as const;
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

export function pinView(view?: DG.ViewBase): void {
  if (view && !view.isPinned)
    view.pin();
}

// shared inline-style colors; tailwind arbitrary-value classes stay literal (the JIT needs them static)
export const STICKY_BAR_BACKGROUND = 'color-mix(in srgb, var(--white) 75%, transparent)';
export const SELECTED_STEP_BACKGROUND = 'var(--grey-1)';
