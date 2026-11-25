import * as DG from 'datagrok-api/dg';
import {historyUtils} from '@datagrok-libraries/compute-utils/history-utils/src/history-utils';
import {openSaveRunDialog} from './SaveRunDialog';
import {openRunBrowserDialog} from './run-browser-dialog';
import {defaultRunName} from './run-storage';

/** The Open dialog of custom function views and history inputs: a saved run of the function, loaded, or null. */
export async function openRunDialog(func: DG.Func): Promise<DG.FuncCall | null> {
  const [run] = await openRunBrowserDialog({title: 'Open run', model: func, allowExport: false}) ?? [];
  return run ? historyUtils.loadRun(run.id) : null;
}

/** The Save dialog of custom function views: where and under which name to save the run, or null. */
export async function saveRunDialog(call: DG.FuncCall) {
  const {name, description, tags} = historyUtils.getRunMetadata(call);
  const result = await openSaveRunDialog({
    func: call.func, name: defaultRunName(call.func, call.author ? name : ''), description, tags,
  });
  return result ? {metadata: result.metadata, target: result.target} : null;
}
