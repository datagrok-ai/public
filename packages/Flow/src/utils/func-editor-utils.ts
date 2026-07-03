import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as rxjs from 'rxjs';
// this code handles the function editors for func nodes.
// prerequisit is that funccall needs to be configured with tables.

export const DEFAULT_EDITOR_SUPPORTED_TYPES = new Set([
  'int', 'num', 'double', 'qnum', 'datetime', 'dataframe', 'table', 'column', 'list', 'list<string>',
  'list<column>', 'string_list', 'column_list', 'string', 'bigint', 'dataframe_list', 'file', 'files', 'blob',
]);

export const EXPLICITLY_SUPPORTED_EDITABLE_FUNCTIONS = new Set([
  'core:AddNewColumn',
]);

export function shouldUseFunctionEditor(func: DG.Func) {
  if (EXPLICITLY_SUPPORTED_EDITABLE_FUNCTIONS.has(func.nqName))
    return true;
  if (!func.options.editor) // anything without editor can be edited on context panel
    return false;

  // anything with custom editor will be supported.
  return true;
}

export async function pollDialogCreation(timeout = 30_000): Promise<DG.Dialog | null> {
  return new Promise((res) => {
    let timeoutNum: any = null;
    const pollInterval = setInterval(() => {
      const cur = DG.Dialog.getOpenDialogs();
      const newD = cur[0];
      if (newD) {
        clearTimeout(timeoutNum);
        clearInterval(pollInterval);
        res(newD);
      }
    }, 100);
    timeoutNum = setTimeout(() => {
      clearInterval(pollInterval);
      res(null);
    }, timeout);
  });
}

export async function createFuncCallEditor(fc: DG.FuncCall): Promise<DG.FuncCall> {
  // here we expect that fc is generated from func.prepare and all needed parames are already passed to it,
  // especially table
  return new Promise<DG.FuncCall>( async (res) => {
    DG.Dialog.getOpenDialogs().forEach((d) => d.close());
    fc.edit();
    const d = await pollDialogCreation();
    if (!d)
      throw new Error('Could not find the dialog for function');
      // override the call methods, so that
    d.root.classList.add('d4-flow-function-funccall-editor'); // style for disabling table inputs
    let dialogSub: rxjs.Subscription | null = null;
    const sub = grok.events.onEvent('d4-before-run-action').subscribe((f: DG.FuncCall) => {
      if (f?.func === fc.func) {
        try {
          f.status = 'Canceled'; // this ,makes sure the funccall not be run,
        } catch (e) {
          console.error(e); // unsupported on current released version
        }
        // and all the parameters are saved to the funccall
        sub.unsubscribe();
        dialogSub?.unsubscribe();
        res(f);
      }
    });

    dialogSub = d.onClose.subscribe(() => {
      setTimeout(() => {
        dialogSub?.unsubscribe();
        sub.unsubscribe();
        res(fc);
      });
    });
  });
}
