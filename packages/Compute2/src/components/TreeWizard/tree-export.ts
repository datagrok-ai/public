import * as DG from 'datagrok-api/dg';
import {
  isFuncCallState,
  PipelineState,
  StepFunCallState,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {zipSync, Zippable} from 'fflate';
import {richFunctionViewReport} from '@datagrok-libraries/compute-utils/shared-utils/export';
import {
  getCustomExports, getStartedOrNull, replaceForWindowsPath,
} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import type ExcelJS from 'exceljs';
import {FuncCallStateInfo, MetaCallInfo} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import type Dayjs from 'dayjs';
import {ExportCbInput, ExportSummaryItem} from
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {
  getRunSummary, RunStates, statesToStatus, textTitle,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/run-summary';
import {ViewersHook} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import {ConsistencyInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {BehaviorSubject} from 'rxjs';
import {disposeViewers, getViewers} from '../../utils';
import {statusToTooltip} from './tree-utils';

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

/** The default Excel report of a call, with its viewers built and hooked for it. */
export async function callExcelReport(call: DG.FuncCall, {viewersHook, metaState, validation, consistency}: {
  viewersHook?: ViewersHook,
  metaState?: Record<string, BehaviorSubject<any>>,
  validation?: Record<string, ValidationResult>,
  consistency?: Record<string, ConsistencyInfo>,
} = {}) {
  const viewers = await getViewers(call, viewersHook, metaState);
  return richFunctionViewReport('Excel', call.func, call, viewers, validation, consistency)
    .finally(() => disposeViewers(viewers));
}

export async function reportTree(
  {
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
    treeState: PipelineState;
    meta?: MetaCallInfo;
    metaStates?: Record<string, Record<string, BehaviorSubject<any>> | undefined>,
    hasNotSavedEdits?: boolean;
    cb?: (input: ExportCbInput) => Promise<void>,
  } & RunStates) {
  const zipConfig: Zippable = {};

  const q = [{state: treeState, idx: 0, path: [] as string[]}];

  while (q.length > 0) {
    const {state, idx, path} = q.shift()!;
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
      const [blob, wb] = await callExcelReport(funcCall,
        {viewersHook: state.viewersHook, metaState, validation, consistency});

      const fileName = stepFileName(state, idx, callInfo, textTitle(description));
      const configKey = [...path, fileName].join('/');
      const stepFile: Zippable[string] = [new Uint8Array(await blob.arrayBuffer()), {level: 0}];
      zipConfig[configKey] = stepFile;
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
          description,
        });
        // the callback may edit the workbook, or replace or remove the file itself
        if (zipConfig[configKey] === stepFile)
          zipConfig[configKey] = [new Uint8Array(await wb.xlsx.writeBuffer()), {level: 0}];
      }
    } else {
      const nPath = state === treeState ? [] : [...path, workflowDirName(state, idx)];
      for (const [idx, stepState] of state.steps.entries())
        q.push({state: stepState, idx, path: nPath});
    }
  }

  const summary = getExportSummary(
    treeState, {callInfoStates, validationStates, consistencyStates, pipelineValidations, descriptions});
  const [summaryBlob] = await reportSummary(summary);
  zipConfig[SUMMARY_FILE_NAME] = [new Uint8Array(await summaryBlob.arrayBuffer()), {level: 0}];

  const rawFileName = getExportName(treeState, !!hasNotSavedEdits, meta.title, meta.started);
  const fileName = replaceForWindowsPath(`${rawFileName}.zip`);
  const blob = new Blob([zipSync(zipConfig) as any]);
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


export function exportedFile(result: unknown, exportName: string): DG.FileInfo {
  if (!(result instanceof DG.FileInfo))
    throw new Error(`Export ${exportName} returned no file`);
  return result;
}

export const downloadFile = (file: DG.FileInfo) => DG.Utils.download(file.fileName, new Blob([file.data as any]));

/** Resolves the custom export named `exportName` declared on `exportsOf` (via `meta.customExports`), the
 *  funcCall's function unless a custom view declares them, and applies it, passing the funcCall through. */
export async function applyCustomExport(
  fc: DG.FuncCall,
  exportName: string,
  args: Record<string, unknown> = {},
  exportsOf: DG.Func = fc.func,
): Promise<DG.FileInfo> {
  const item = getCustomExports(exportsOf).find((x) => x.name === exportName);
  if (!item)
    throw new Error(`No export named ${exportName} is defined for ${exportsOf.nqName}`);
  return exportedFile(await DG.Func.byName(item.function).apply({funcCall: fc, ...args}), exportName);
}
