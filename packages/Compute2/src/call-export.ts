import * as DG from 'datagrok-api/dg';
import {richFunctionViewReport} from '@datagrok-libraries/compute-utils/shared-utils/export';
import {getStartedOrNull, replaceForWindowsPath} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {ViewersHook} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import {ConsistencyInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {disposeViewers, getViewers} from './utils';

export const blobFile = async (fileName: string, blob: Blob) =>
  DG.FileInfo.fromBytes(fileName, new Uint8Array(await blob.arrayBuffer()));

/** The default Excel report of a call, named after its run, with its viewers built and hooked for it. */
export async function callExcelReport(call: DG.FuncCall, {viewersHook, metaState, validation, consistency, completed}: {
  viewersHook?: ViewersHook,
  metaState?: Record<string, {value: any}>,
  validation?: Record<string, ValidationResult>,
  consistency?: Record<string, ConsistencyInfo>,
  completed?: boolean,
} = {}) {
  const viewers = await getViewers(call, viewersHook, metaState);
  const [blob, wb] = await richFunctionViewReport('Excel', call.func, call, viewers, validation, consistency)
    .finally(() => disposeViewers(viewers));
  return [await blobFile(`${callFileName(call, completed)}.xlsx`, blob), wb] as const;
}

/** Unknown completion, as for a custom view's run, keeps the run's date. */
export function callFileName(call: DG.FuncCall, completed?: boolean) {
  const started = completed === false ? undefined : getStartedOrNull(call);
  return replaceForWindowsPath(`${call.func.nqName} - ${started ?? 'Not completed'}`);
}

/** A legacy `meta.customExports` function returns a Blob, named `blobName` and `.xlsx` when it is Excel. */
export async function exportedFile(result: unknown, exportName: string, blobName?: string): Promise<DG.FileInfo> {
  if (blobName && result instanceof Blob) {
    const ext = result.type.includes('spreadsheetml.sheet') ? '.xlsx' : '';
    return blobFile(`${blobName}${ext}`, result);
  }
  if (!(result instanceof DG.FileInfo))
    throw new Error(`Export ${exportName} returned no file`);
  return result;
}

export const downloadFile = (file: DG.FileInfo) => DG.Utils.download(file.fileName, new Blob([file.data as any]));
