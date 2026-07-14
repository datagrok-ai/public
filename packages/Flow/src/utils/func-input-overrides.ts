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
 *    compilation, required-input checks, and serialization are unaffected).
 *  - {@link FUNC_WRAPPERS} — the node exposes a reshaped, Flow-friendly input
 *    list instead of the function's own awkward signature; at compile time the
 *    wrapper folds the resolved inputs back into the real arguments (e.g.
 *    AppendTables' unwirable `tables: dataframe_list` → two plain table
 *    sockets → `tables: [table1, table2]`). */

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

// ---------- function wrappers ----------

/** One input exposed by a {@link FuncWrapper} in place of the function's own. */
export interface WrappedFuncInput {
  name: string;
  /** DG property type of the exposed socket ('dataframe', 'string', …). */
  type: string;
  caption?: string;
  description?: string;
  /** Optional inputs don't gate the run (no "Needs input"). */
  optional?: boolean;
  defaultValue?: unknown;
}

/** Reshapes a function's signature into Flow-friendly node inputs. The node
 *  (sockets, pass-throughs, seeds, required checks, panel) is built from
 *  {@link inputs} exactly as if the function declared them; the compiler then
 *  runs the resolved input expressions through {@link mapInputs} to build the
 *  real `grok.functions.call` arguments. */
export interface FuncWrapper {
  inputs: WrappedFuncInput[];
  /** Resolved JS expressions of the exposed inputs (a key is absent when the
   *  input is unconnected and blank) → the function's real named arguments,
   *  also JS expressions. */
  mapInputs: (exposed: Record<string, string>) => Record<string, string>;
}

/** `{[func.nqName]: wrapper}` — functions kept in the catalog but exposed
 *  through reshaped inputs. Edit freely: add an entry to make an awkward
 *  signature wirable. */
export const FUNC_WRAPPERS: Record<string, FuncWrapper> = {
  // `tables: dataframe_list` is unwirable — nothing on a canvas outputs a
  // dataframe list. Expose two plain tables and fold them into the list.
  'core:AppendTables': {
    inputs: [
      {name: 'table1', type: 'dataframe', description: 'Defines the result columns'},
      {name: 'table2', type: 'dataframe', description: 'Appended below the first table'},
    ],
    mapInputs: (v) => v.table1 && v.table2 ? {tables: `[${v.table1}, ${v.table2}]`} : {} as Record<string, string>,
  },
};

/** The registered wrapper for a function, or null. */
export function funcWrapperOf(func: DG.Func): FuncWrapper | null {
  try {
    return FUNC_WRAPPERS[func.nqName] ?? null;
  } catch {
    return null;
  }
}

/** The exposed inputs as real `DG.Property` objects, so `FuncNode` builds
 *  sockets, seeds, and required checks through the same code path as real
 *  params. `nullable: false` unless declared optional — the base Property
 *  defaults to nullable (= optional) for non-strings, which would kill the
 *  "Needs input" gate on the exposed tables. */
export function wrapperProperties(wrapper: FuncWrapper): DG.Property[] {
  return wrapper.inputs.map((s) => DG.Property.fromOptions({
    name: s.name, type: s.type, caption: s.caption, description: s.description,
    nullable: s.optional === true,
    ...(s.defaultValue !== undefined ? {defaultValue: s.defaultValue as string} : {}),
  }));
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
