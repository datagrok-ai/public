import * as DG from 'datagrok-api/dg';
import {BehaviorSubject, defer} from 'rxjs';
import {concatMap, finalize} from 'rxjs/operators';
import {historyUtils} from '@datagrok-libraries/compute-utils/history-utils/src/history-utils';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {
  isFuncCallState, PipelineState,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {
  CustomExport, ExportCbInput, ExportUtils, ViewersHook,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {
  PipelineConfigurationProcessed,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {
  ConfigTraverseItem, isPipelineStaticConfig, isPipelineStepConfig,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-utils';
import {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import {
  ConsistencyInfo, MetaCallInfo,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {RunStates} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/run-summary';
import {
  loadSavedTree, metaCallInfo, readStates,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/headless';
import {
  RUN_ERROR_PATH, readOutputOutdated,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {
  callExcelReport, callFileName, exportedFile, getExportSummary, reportSummaryFile, reportTree,
} from './components/TreeWizard/tree-export';
import {findNodeWithPathByUuid, resolveSingleStep} from './components/TreeWizard/tree-utils';
import {runCompleted} from './components/RunStorage/run-storage';

/** An export as menus list it; `id` finds it again in a run of another version. */
export type ExportOption = {id: string, name: string};

// custom exports are named freely, so the default export takes a name they don't use
export const DEFAULT_EXPORT = '$default';
const DEFAULT_EXPORT_OPTION: ExportOption = {id: DEFAULT_EXPORT, name: 'Default Excel'};

/** Preferred exports first, in their declared order, then the default report and the other exports. */
const menuOrder = (defaults: ExportOption[], exports: (ExportOption & {preferred?: boolean})[]): ExportOption[] => [
  ...exports.filter((x) => x.preferred), ...defaults, ...exports.filter((x) => !x.preferred),
].map(({id, name}) => ({id, name}));

/** The exports of a function's runs, in the menu's order: the first is the default. */
export function callExports(func: DG.Func, hideDefault = false): ExportOption[] {
  const isReportEnabled = Utils.getFeature(Utils.getFeatures(func), 'export', true);
  return menuOrder(isReportEnabled && !hideDefault ? [DEFAULT_EXPORT_OPTION] : [],
    Utils.getCustomExports(func).filter((x) => x.function && x.name)
      .map((x) => ({id: x.name, name: x.name, preferred: x.preferred})));
}

export type CallExportStates = {
  viewersHook?: ViewersHook,
  metaState?: Record<string, BehaviorSubject<any>>,
  validation?: Record<string, ValidationResult>,
  consistency?: Record<string, ConsistencyInfo>,
  isOutputOutdated?: boolean,
  runError?: string,
  completed?: boolean,
};

/** Exports a call; custom exports come from `model`, which is a custom view's function for its runs. Every
 * export menu goes through here, so an export function always gets the same arguments. */
export async function exportCall(
  call: DG.FuncCall, exportId: string, states: CallExportStates = {}, model: DG.Func = call.func,
): Promise<DG.FileInfo> {
  if (exportId === DEFAULT_EXPORT) {
    const [file] = await callExcelReport(call, states);
    return file;
  }
  const item = Utils.getCustomExports(model).find((x) => x.name === exportId);
  if (!item)
    throw new Error(`No export named ${exportId} is defined for ${model.nqName}`);
  const result = await DG.Func.byName(item.function).apply({
    funcCall: call,
    startDownload: false,
    validationState: states.validation,
    consistencyState: states.consistency,
    isOutputOutdated: states.isOutputOutdated,
    runError: states.runError,
    viewersHook: states.viewersHook,
    metaState: states.metaState,
  });
  return exportedFile(result, exportId, callFileName(call, states.completed));
}

/** A custom view's run takes the viewers hook of its own function, which is not the view's function. */
const callHook = (call: DG.FuncCall, model?: DG.Func,
  hookOf: (func: DG.Func) => Promise<ViewersHook | undefined> = Utils.loadViewersHook) =>
  hookOf(model?.name === call.func.name ? model : call.func);

/** The Export menu of a custom view: the exports of its function, as bulk export runs them. */
export const viewExports = (model: DG.Func) => callExports(model).map(({id, name}) => ({
  name, export: async (call: DG.FuncCall) => exportCall(call, id, {viewersHook: await callHook(call, model)}, model),
}));

/** The exports of a workflow, in the menu's order: the first is the default. */
export const workflowExports = (node: {disableDefaultExport?: boolean, customExports?: CustomExport[]}) =>
  menuOrder(node.disableDefaultExport ? [] : [DEFAULT_EXPORT_OPTION],
    (node.customExports ?? []).map((x) => ({id: x.id, name: x.friendlyName ?? x.id, preferred: x.preferred})));

/** The exports of a workflow's runs, as the menus show them once a run is open: a compact view shows its step's. */
export function configExports(config: PipelineConfigurationProcessed): ExportOption[] {
  let node: ConfigTraverseItem = config;
  while (config.compactView && isPipelineStaticConfig(node) && !node.isActionStep && node.steps.length === 1)
    node = node.steps[0];
  return config.compactView && isPipelineStepConfig(node) ?
    callExports(DG.Func.byName(node.nqName), config.disableDefaultExport) : workflowExports(config);
}

export type TreeExportContext = {
  states: RunStates,
  metaStates?: Record<string, Record<string, BehaviorSubject<any>> | undefined>,
  meta?: MetaCallInfo,
  hasNotSavedEdits?: boolean,
  completed?: boolean,
};

const reportStates = ({states, metaStates, meta, hasNotSavedEdits, completed}: TreeExportContext) =>
  ({...states, metaStates, meta, hasNotSavedEdits, completed});

/** The states of a workflow step, found in `root` by its uuid, as its own export menu passes them. */
function stepStates(root: PipelineState, ctx: TreeExportContext, uuid: string): CallExportStates {
  const step = findNodeWithPathByUuid(uuid, root)?.state;
  const outdated = ctx.states.callInfoStates?.[uuid]?.isOutputOutdated;
  return {
    viewersHook: step && isFuncCallState(step) ? step.viewersHook : undefined,
    metaState: ctx.metaStates?.[uuid],
    validation: ctx.states.validationStates?.[uuid],
    consistency: ctx.states.consistencyStates?.[uuid],
    isOutputOutdated: outdated,
    runError: ctx.states.callInfoStates?.[uuid]?.runError,
    completed: outdated === undefined ? undefined : !outdated,
  };
}

/** The utils a workflow's custom export gets; a step is found in `root` by its uuid. */
export function treeExportUtils(root: PipelineState, ctx: TreeExportContext): ExportUtils {
  return {
    reportFuncCallExcel: (fc: DG.FuncCall, uuid: string) => callExcelReport(fc, stepStates(root, ctx, uuid)),
    reportStateExcel: (state: PipelineState, cb?: (input: ExportCbInput) => Promise<void>) =>
      reportTree({treeState: state, ...reportStates(ctx), cb}),
    getExportSummary: (state: PipelineState) => getExportSummary(state, ctx.states),
    reportSummaryExcel: (state: PipelineState) => reportSummaryFile(state, ctx.states, ctx),
    getFuncCallCustomExports: (fc: DG.FuncCall) => callExports(fc.func, true).map(({id}) => id),
    runFuncCallCustomExport: (fc: DG.FuncCall, uuid: string, exportName: string) =>
      exportCall(fc, exportName, stepStates(root, ctx, uuid)),
  };
}

export async function exportWorkflow(
  state: PipelineState, exportId: string, ctx: TreeExportContext,
): Promise<DG.FileInfo> {
  if (exportId === DEFAULT_EXPORT) {
    const [file] = await reportTree({treeState: state, ...reportStates(ctx)});
    return file;
  }
  const item = isFuncCallState(state) ? undefined : state.customExports?.find((x) => x.id === exportId);
  if (!item)
    throw new Error(`No export ${exportId} is defined for ${state.configId}`);
  return exportedFile(await item.handler(state, treeExportUtils(state, ctx), false), item.friendlyName ?? item.id);
}

async function exportTree(tree: StateTree, metaCall: DG.FuncCall, exportId: string) {
  const state = tree.toState();
  const states = readStates(tree);
  const completed = Object.values(states.callInfoStates ?? {}).every((info) => info?.isOutputOutdated === false);
  if (!completed)
    return null;
  const metaStates = Object.fromEntries(Object.entries(tree.getMeta()).map(([uuid, meta]) => [uuid, meta.value]));
  const ctx: TreeExportContext = {states, metaStates, meta: metaCallInfo(metaCall), completed};
  const single = resolveSingleStep(state);
  if (!single)
    return exportWorkflow(state, exportId, ctx);
  const {uuid, funcCall} = single.step;
  if (!funcCall)
    throw new Error(`The step of run ${metaCall.id} didn't load`);
  return exportCall(funcCall, exportId, stepStates(state, ctx, uuid));
}

/** Exports a saved workflow run as its menu does once the run is open; null when a step's outputs are outdated,
 * as the menu then offers nothing. */
export const exportSavedWorkflow = (runId: string, exportId: string) =>
  loadSavedTree(runId).pipe(concatMap(({tree, metaCall}) =>
    defer(() => exportTree(tree, metaCall, exportId)).pipe(finalize(() => tree.close())),
  )).toPromise() as Promise<DG.FileInfo | null>;

/** Exports a saved run of `model`; the report takes the viewers hook of the run's function, which is not the
 * model for a custom view. */
export async function exportSavedCall(runId: string, exportId: string, model?: DG.Func,
  hookOf: (func: DG.Func) => Promise<ViewersHook | undefined> = Utils.loadViewersHook) {
  const call = await historyUtils.loadRun(runId);
  return exportCall(call, exportId, {
    viewersHook: await callHook(call, model, hookOf),
    isOutputOutdated: readOutputOutdated(call),
    runError: call.options[RUN_ERROR_PATH] || undefined,
    completed: runCompleted(call),
  }, model);
}
