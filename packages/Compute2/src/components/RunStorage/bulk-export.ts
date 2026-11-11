import * as DG from 'datagrok-api/dg';
import {zipSync, Zippable} from 'fflate';
import {loadViewersHook, replaceForWindowsPath} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {ViewersHook} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {loadWorkflowConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/headless';
import {
  ExportOption, callExports, configExports, exportSavedCall, exportSavedWorkflow,
} from '../../run-exports';

/** A run to export, with the folders it goes in: the spaces down to it from the exported space. */
export type ExportRun = {run: DG.FuncCall, model: DG.Func, spaceId?: string, folders: string[]};

/** Workflow runs are loaded and exported as their workflow; the rest, step runs included, as a function call. */
export type ModelKind = 'call' | 'workflow';

export const modelKind = (model: DG.Func, stepHistory: boolean): ModelKind =>
  !stepHistory && model.options['editor'] === 'Compute2:TreeWizardEditor' ? 'workflow' : 'call';

/** The exports of a model's runs as the menus list them; a workflow's come from its current config. */
export async function modelExports(model: DG.Func, kind: ModelKind): Promise<ExportOption[]> {
  if (kind === 'call')
    return callExports(model);
  const config = await loadWorkflowConfig(model.options['provider'] ?? model.nqName).toPromise();
  return configExports(config!);
}

export type ExportJob = {item: ExportRun, kind: ModelKind, exportId: string};
export type ExportedFile = {folders: string[], file: DG.FileInfo};
export type FailedExport = {run: DG.FuncCall, error: string};

/** Exports runs one at a time; a run that fails, or turns out not completed once loaded, is reported and the rest
 * go on. Cancelling stops before the next run. A step's history passes the step's hook, which replaces the
 * function's own. */
export async function exportRuns(jobs: ExportJob[], {onProgress, isCancelled, viewersHook}: {
  onProgress: (done: number) => void, isCancelled: () => boolean, viewersHook?: ViewersHook,
}) {
  const files: ExportedFile[] = [];
  const failed: FailedExport[] = [];
  const notCompleted: DG.FuncCall[] = [];
  const hooks = new Map<string, Promise<ViewersHook | undefined>>();
  const hookOf = (func: DG.Func) => {
    if (viewersHook)
      return Promise.resolve(viewersHook);
    if (!hooks.has(func.nqName))
      hooks.set(func.nqName, loadViewersHook(func));
    return hooks.get(func.nqName)!;
  };
  for (const [idx, {item, kind, exportId}] of jobs.entries()) {
    if (isCancelled())
      break;
    try {
      const file = kind === 'workflow' ? await exportSavedWorkflow(item.run.id, exportId) :
        await exportSavedCall(item.run.id, exportId, item.model, hookOf);
      if (file)
        files.push({folders: item.folders, file});
      else
        notCompleted.push(item.run);
    } catch (e: any) {
      failed.push({run: item.run, error: String(e?.message ?? e)});
    }
    onProgress(idx + 1);
  }
  return {files, failed, notCompleted};
}

const splitName = (fileName: string) => {
  const dot = fileName.lastIndexOf('.');
  return dot > 0 ? [fileName.slice(0, dot), fileName.slice(dot)] : [fileName, ''];
};

/** Zips the files into their folders; a name taken in a folder, ignoring case, gets a number. */
export function zipExports(files: ExportedFile[]): Uint8Array {
  const zip: Zippable = {};
  const taken = new Set<string>();
  for (const {folders, file} of files) {
    // a folder named . or .. would point outside its place in the zip
    const dir = folders.map((name) => `${/^\.+$/.test(name.trim()) ? '_' : replaceForWindowsPath(name)}/`).join('');
    const [stem, ext] = splitName(replaceForWindowsPath(file.fileName));
    let path = `${dir}${stem}${ext}`;
    for (let n = 2; taken.has(path.toLowerCase()); n++)
      path = `${dir}${stem} (${n})${ext}`;
    taken.add(path.toLowerCase());
    zip[path] = [file.data, {level: 0}];
  }
  return zipSync(zip);
}
