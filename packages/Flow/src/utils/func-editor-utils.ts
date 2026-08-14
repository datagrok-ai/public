import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as rxjs from 'rxjs';
// Function-editor dialog plumbing for func nodes; the FuncCall must arrive pre-configured with its tables.

export const EXPLICITLY_SUPPORTED_EDITABLE_FUNCTIONS = new Set([
  'core:AddNewColumn',
]);

export function shouldUseFunctionEditor(func: DG.Func) {
  if (EXPLICITLY_SUPPORTED_EDITABLE_FUNCTIONS.has(func.nqName))
    return true;
  if (!func.options.editor) // anything without an editor can be edited on the context panel
    return false;

  return true;
}

/** Inputs that practically REQUIRE the function's own editor: the panel renders an inline
 *  pencil that opens it. Keyed `nqName:inputName`, case-insensitive. */
export const EDITOR_SHORTCUT_INPUTS = new Set([
  // NOT AddNewColumn's `expression`: its inline custom editor short-circuits the
  // DG-input path this pencil attaches to, so an entry for it would never render.
  'Chem:descriptorsDocker:selected',
  'Chem:deprotect:fragment',
].map((s) => s.toLowerCase()));

/** Implies {@link shouldUseFunctionEditor} — never true when the launcher would refuse to open the editor. */
export function hasEditorShortcut(func: DG.Func, inputName: string): boolean {
  try {
    if (!shouldUseFunctionEditor(func)) return false;
    return EDITOR_SHORTCUT_INPUTS.has(`${func.nqName}:${inputName}`.toLowerCase());
  } catch {
    return false; // Dart proxy access can throw — treat as no shortcut
  }
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

export async function createFuncCallEditor(
  fc: DG.FuncCall,
  opts?: {
    /** Return true to skip a `d4-before-run-action` event not from this dialog — the event fires
     *  for EVERY client funccall, and without the guard a concurrent Flow run of the same function
     *  gets canceled and resolves this round-trip with the wrong funccall. */
    ignoreEvent?: () => boolean;
  },
): Promise<DG.FuncCall> {
  return new Promise<DG.FuncCall>( async (res) => {
    DG.Dialog.getOpenDialogs().forEach((d) => d.close());
    fc.setAuxValue('forceEditParameters', true);
    fc.edit();
    const d = await pollDialogCreation();
    if (!d) {
      // Settle rather than throw: a throw inside an async executor is swallowed and the
      // caller's finally (which releases the autorun hold) would never run.
      console.warn('Flow: the function editor dialog never appeared — returning the call as-is');
      res(fc);
      return;
    }
    d.root.classList.add('d4-flow-function-funccall-editor'); // style for disabling table inputs
    let dialogSub: rxjs.Subscription | null = null;
    const sub = grok.events.onEvent('d4-before-run-action').subscribe((f: DG.FuncCall) => {
      if (opts?.ignoreEvent?.()) return; // someone else's funccall (e.g. a Flow run) — not this dialog
      if (f?.func === fc.func) {
        try {
          f.status = 'Canceled'; // prevents the funccall from actually running
        } catch (e) {
          console.error(e); // unsupported on current released version
        }
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
