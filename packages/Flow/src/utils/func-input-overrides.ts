/** Per-function input overrides — curated maps keyed by `Func.nqName` (same
 *  convention as `EXCLUDED_FUNC_NQNAMES`), consulted when a func node and its
 *  parameter form are built:
 *
 *  - {@link HIDDEN_FUNC_INPUTS} — inputs users should never see: no socket row
 *    on the node, no editor on the context panel. Purely visual — the slot,
 *    its seeded value, compilation, and creation-script import/emit stay
 *    untouched (data-sync scripts legitimately carry e.g.
 *    `subscribeOnChanges = true` and must round-trip). Only hide inputs with a
 *    declared default.
 *  - {@link CUSTOM_FUNC_INPUT_EDITORS} — a bespoke editor replacing the default
 *    one for a specific parameter (storage stays `inputValues[name]`, so
 *    compilation, required-input checks, and serialization are unaffected). */

import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {getParamDisplayName} from './dart-proxy-utils';

/** `{[func.nqName]: {[inputName]: true}}` — inputs hidden from the node and
 *  the context panel (visually only — see the module doc). Edit freely: flag
 *  an input to hide it everywhere. */
export const HIDDEN_FUNC_INPUTS: Record<string, Record<string, boolean>> = {
  // Reactivity/error plumbing for the formula engine — never user-relevant.
  'core:AddNewColumn': {subscribeOnChanges: true, errorBehavior: true},
  'core:AddNewColumnList': {subscribeOnChanges: true, errorBehavior: true},
};

/** The hidden-input names for a function, `∅` when none are registered
 *  (or the Dart proxy throws on `nqName`). */
export function hiddenInputsOf(func: DG.Func): ReadonlySet<string> {
  try {
    const map = HIDDEN_FUNC_INPUTS[func.nqName];
    return new Set(map ? Object.keys(map).filter((k) => map[k]) : []);
  } catch {
    return new Set();
  }
}

/** A custom editor instance for one parameter of one node. The panel assigns
 *  {@link onChanged} after construction; the editor must invoke it on every
 *  user edit with the value to store in `inputValues`. */
export interface CustomInputEditor {
  element: HTMLElement;
  getValue: () => unknown;
  /** Initialize the editor from the stored value (may be blank/null). */
  setValue: (v: unknown) => void;
  /** When present and false, the pending value is not stored. */
  isValid?: () => boolean;
  onChanged?: (v: unknown) => void;
}

export type CustomInputEditorFactory = (param: DG.Property) => CustomInputEditor;

/** `{[func.nqName]: {[inputName]: factory}}` — a factory per overridden input
 *  (an editor holds a live HTMLElement, so it must be built per node/render). */
export const CUSTOM_FUNC_INPUT_EDITORS: Record<string, Record<string, CustomInputEditorFactory>> = {
  'core:OpenFile': {fullPath: filePathEditor},
};

/** The registered custom-editor factory for a function input, or null. */
export function customEditorFor(func: DG.Func, inputName: string): CustomInputEditorFactory | null {
  try {
    return CUSTOM_FUNC_INPUT_EDITORS[func.nqName]?.[inputName] ?? null;
  } catch {
    return null;
  }
}

/** A file picker (`ui.input.file`) for string path parameters like OpenFile's
 *  `fullPath`. The stored value is the plain full-path string
 *  (`System:AppData/...`) — exactly what the parameter takes — read from the
 *  picked `FileInfo`. `setValue` rebuilds a `FileInfo` from the stored path
 *  (synchronous — no server round-trip) rather than routing through
 *  `stringValue`, whose async server resolution is what validates hand-typed
 *  paths. */
function filePathEditor(param: DG.Property): CustomInputEditor {
  const ed: CustomInputEditor = {} as CustomInputEditor;
  const input = ui.input.file(getParamDisplayName(param), {
    onValueChanged: (v) => ed.onChanged?.(v ? v.fullPath : ''),
    // temporary thing, remove local file opening. once we figure out how to handle local files, remove this
    onCreated: (a) => a.root.querySelector('.ui-input-options')?.querySelector('.fa-folder-open')?.remove?.()
  });
  ed.element = input.root;
  ed.getValue = (): unknown => {
    try {
      return input.value?.fullPath ?? '';
    } catch {
      return '';
    }
  };
  ed.setValue = (v): void => {
    try {
      if (v !== undefined && v !== null && String(v) !== '')
        input.value = DG.FileInfo.fromString(String(v), '');
    } catch {/* leave the editor blank */}
  };
  return ed;
}
