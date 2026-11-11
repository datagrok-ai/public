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
import {readOutputOutdated} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {
  applyCustomExport, callExcelReport, exportedFile, getExportSummary, reportSummary, reportTree,
} from './components/TreeWizard/tree-export';
import {findNodeWithPathByUuid, resolveSingleStep} from './components/TreeWizard/tree-utils';

/** An export as menus list it; `id` finds it again in a run of another version. */
export type ExportOption = {id: string, name: string};

// custom exports are named freely, so the default export takes a name they don't use
export const DEFAULT_EXPORT = '$default';
const DEFAULT_EXPORT_OPTION: ExportOption = {id: DEFAULT_EXPORT, name: 'Default Excel'};

const fileOf = async (fileName: string, blob: Blob) =>
  DG.FileInfo.fromBytes(fileName, new Uint8Array(await blob.arrayBuffer()));

/** The exports of a function's runs, in the menu's order: the first is the default. */
export function callExports(func: DG.Func, hideDefault = false): ExportOption[] {
  const isReportEnabled = Utils.getFeature(Utils.getFeatures(func), 'export', true);
  return [
    ...(isReportEnabled && !hideDefault ? [DEFAULT_EXPORT_OPTION] : []),
    ...Utils.getCustomExports(func).filter((x) => x.function && x.name).map((x) => ({id: x.name, name: x.name})),
  ];
}

export type CallExportStates = {
  viewersHook?: ViewersHook,
  metaState?: Record<string, BehaviorSubject<any>>,
  validation?: Record<string, ValidationResult>,
  consistency?: Record<string, ConsistencyInfo>,
  isOutputOutdated?: boolean,
  runError?: string,
};

export async function exportCall(
  call: DG.FuncCall, exportId: string, states: CallExportStates = {},
): Promise<DG.FileInfo> {
  if (exportId !== DEFAULT_EXPORT) {
    return applyCustomExport(call, exportId, {
      validationState: states.validation,
      consistencyState: states.consistency,
      isOutputOutdated: states.isOutputOutdated,
      runError: states.runError,
    });
  }
  const [blob] = await callExcelReport(call, states);
  return fileOf(`${call.func.nqName} - ${Utils.getStartedOrNull(call) ?? 'Not completed'}.xlsx`, blob);
}

export async function loadViewersHook(func: DG.Func): Promise<ViewersHook | undefined> {
  const maker = Utils.getViewersHook(func);
  return maker ? DG.Func.byName(maker).apply() : undefined;
}

/** The exports of a workflow, in the menu's order: the first is the default. */
export const workflowExports = (node: {disableDefaultExport?: boolean, customExports?: CustomExport[]}) => [
  ...(node.disableDefaultExport ? [] : [DEFAULT_EXPORT_OPTION]),
  ...(node.customExports ?? []).map((x): ExportOption => ({id: x.id, name: x.friendlyName ?? x.id})),
];

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
};

const reportStates = ({states, metaStates, meta, hasNotSavedEdits}: TreeExportContext) =>
  ({...states, metaStates, meta, hasNotSavedEdits});

/** The utils a workflow's custom export gets; a step is found in `root` by its uuid. */
export function treeExportUtils(root: PipelineState, ctx: TreeExportContext): ExportUtils {
  const {validationStates, consistencyStates, callInfoStates} = ctx.states;
  return {
    reportFuncCallExcel: (fc: DG.FuncCall, uuid: string) => {
      const step = findNodeWithPathByUuid(uuid, root)?.state;
      return callExcelReport(fc, {
        viewersHook: step && isFuncCallState(step) ? step.viewersHook : undefined,
        metaState: ctx.metaStates?.[uuid],
        validation: validationStates?.[uuid],
        consistency: consistencyStates?.[uuid],
      });
    },
    reportStateExcel: (state: PipelineState, cb?: (input: ExportCbInput) => Promise<void>) =>
      reportTree({treeState: state, ...reportStates(ctx), cb}),
    getExportSummary: (state: PipelineState) => getExportSummary(state, ctx.states),
    reportSummaryExcel: (state: PipelineState) => reportSummary(getExportSummary(state, ctx.states)),
    getFuncCallCustomExports: (fc: DG.FuncCall) => Utils.getCustomExports(fc.func).map((x) => x.name),
    runFuncCallCustomExport: (fc: DG.FuncCall, uuid: string, exportName: string) =>
      applyCustomExport(fc, exportName, {
        validationState: validationStates?.[uuid],
        consistencyState: consistencyStates?.[uuid],
        isOutputOutdated: callInfoStates?.[uuid]?.isOutputOutdated,
        runError: callInfoStates?.[uuid]?.runError,
      }),
  };
}

export async function exportWorkflow(
  state: PipelineState, exportId: string, ctx: TreeExportContext,
): Promise<DG.FileInfo> {
  if (exportId === DEFAULT_EXPORT) {
    const [blob, , fileName] = await reportTree({treeState: state, ...reportStates(ctx)});
    return fileOf(fileName, blob);
  }
  const item = isFuncCallState(state) ? undefined : state.customExports?.find((x) => x.id === exportId);
  if (!item)
    throw new Error(`No export ${exportId} is defined for ${state.configId}`);
  return exportedFile(await item.handler(state, treeExportUtils(state, ctx)), item.friendlyName ?? item.id);
}

async function exportTree(tree: StateTree, metaCall: DG.FuncCall, exportId: string) {
  const state = tree.toState();
  const states = readStates(tree);
  if (!Object.values(states.callInfoStates ?? {}).every((info) => info?.isOutputOutdated === false))
    return null;
  const metaStates = Object.fromEntries(Object.entries(tree.getMeta()).map(([uuid, meta]) => [uuid, meta.value]));
  const single = resolveSingleStep(state);
  if (!single)
    return exportWorkflow(state, exportId, {states, metaStates, meta: metaCallInfo(metaCall)});
  const {uuid, funcCall, viewersHook} = single.step;
  if (!funcCall)
    throw new Error(`The step of run ${metaCall.id} didn't load`);
  return exportCall(funcCall, exportId, {
    viewersHook, metaState: metaStates[uuid],
    validation: states.validationStates?.[uuid], consistency: states.consistencyStates?.[uuid],
    isOutputOutdated: states.callInfoStates?.[uuid]?.isOutputOutdated,
    runError: states.callInfoStates?.[uuid]?.runError,
  });
}

/** Exports a saved workflow run as its menu does once the run is open; null when a step's outputs are outdated,
 * as the menu then offers nothing. */
export const exportSavedWorkflow = (runId: string, exportId: string) =>
  loadSavedTree(runId).pipe(concatMap(({tree, metaCall}) =>
    defer(() => exportTree(tree, metaCall, exportId)).pipe(finalize(() => tree.close())),
  )).toPromise() as Promise<DG.FileInfo | null>;

export async function exportSavedCall(runId: string, exportId: string, viewersHook?: ViewersHook) {
  const call = await historyUtils.loadRun(runId);
  return exportCall(call, exportId, {viewersHook, isOutputOutdated: readOutputOutdated(call)});
}
