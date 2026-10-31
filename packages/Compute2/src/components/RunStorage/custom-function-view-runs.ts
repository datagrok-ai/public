import * as DG from 'datagrok-api/dg';
import type {CustomFunctionView} from '@datagrok-libraries/compute-utils/function-views/src/custom-function-view';
import {historyUtils} from '@datagrok-libraries/compute-utils/history-utils/src/history-utils';
import {openSaveRunDialog} from './SaveRunDialog';
import {openRunBrowserDialog} from './run-browser-dialog';
import {defaultRunName} from './run-storage';

/** Gives a custom function view the Save and Open dialogs; returns the cleanup. Runs saved with a custom view
 * function open in that view again. */
export function attachRunStorage(view: CustomFunctionView, customView?: DG.Func) {
  view.runStorage = {
    customView: customView?.nqName,
    saveDialog: async (call: DG.FuncCall) => {
      const {name, description, tags} = historyUtils.getRunMetadata(call);
      const result = await openSaveRunDialog({
        func: call.func, name: defaultRunName(call.func, call.author ? name : ''), description, tags,
      });
      return result ? {metadata: result.metadata, target: result.target} : null;
    },
    openDialog: async (func: DG.Func) => {
      const [run] = await openRunBrowserDialog({title: 'Open run', model: func}) ?? [];
      return run ? historyUtils.loadRun(run.id) : null;
    },
  };
  return () => {
    view.runStorage = undefined;
  };
}
