/** Per-function input overrides keyed by `Func.nqName` — hidden/panel-only inputs, hidden
 *  outputs, captions, custom editors, validators, and wrappers, consulted at node build time. */

import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {getParamDisplayName} from './dart-proxy-utils';
import {mpoColumnMappingEditor, mpoMappingRequirements} from '../panel/editors/mpo-mapping-editor';
import {columnFormulaEditor, expressionRequirements, rowConditionEditor} from '../panel/editors/expression-editor';
import {aggregationEditor} from '../panel/editors/aggregation-editor';
import {aggregationProblems, parseAggregations} from '../ops/data-ops';
// Type-only (erased at build time), so the scheme ↔ overrides pair avoids a runtime import cycle.
import type {FlowNode} from '../rete/scheme';

/** `{[func.nqName]: {[inputName]: true}}` — inputs hidden from the node AND the panel.
 *  Visual only (the data layer round-trips untouched); only hide inputs with a declared default. */
export const HIDDEN_FUNC_INPUTS: Record<string, Record<string, boolean>> = {
  // Reactivity/error plumbing for the formula engine — never user-relevant.
  'core:AddNewColumn': {subscribeOnChanges: true, errorBehavior: true},
  'core:AddNewColumnList': {subscribeOnChanges: true, errorBehavior: true},
  // `join` puts the converted column into the table and `overwrite` would replace the
  // source — both must stay at their defaults or the node's output carries nothing.
  'Chem:convertNotation': {join: true, overwrite: true},
  'Chem:recalculateCoords': {join: true},
};

/** Inputs with no socket row on the node body but a regular editor on the panel
 *  (same visual-only contract as {@link HIDDEN_FUNC_INPUTS}; a wired slot renders its row). */
export const PANEL_ONLY_FUNC_INPUTS: Record<string, Record<string, boolean>> = {
  'core:OpenFile': {fullPath: true, sheetName: true},
};

/** Declared outputs the node should not expose (visual only). */
export const HIDDEN_FUNC_OUTPUTS: Record<string, Record<string, boolean>> = {
  // Returns the NAMES of the appended columns; the result users want is the table pass-through.
  'Chem:runElementalAnalysis': {res: true},
};

/** Display captions for parameters whose declared names read as implementation;
 *  slot keys, stored values, and compilation stay on the real names. */
export const FUNC_INPUT_CAPTIONS: Record<string, Record<string, string>> = {
  'core:JoinTables': {
    table1: 'Left table', table2: 'Right table',
    keys1: 'Left keys', keys2: 'Right keys',
    values1: 'Columns from left', values2: 'Columns from right',
  },
};

export function inputCaptionOf(func: DG.Func, name: string): string | null {
  try {
    return FUNC_INPUT_CAPTIONS[func.nqName]?.[name] ?? null;
  } catch {
    return null;
  }
}

/** Inputs whose panel row only renders while the predicate holds; the slot and stored value stay untouched. */
export const CONDITIONAL_FUNC_INPUTS: Record<string, Record<string, (node: FlowNode) => boolean>> = {
  'core:OpenFile': {sheetName: (node) => /\.xlsx$/i.test(String(node.inputValues['fullPath'] ?? ''))},
};

/** Fails open (row shown) on any error. */
export function inputHiddenByCondition(func: DG.Func, name: string, node: FlowNode): boolean {
  try {
    const p = CONDITIONAL_FUNC_INPUTS[func.nqName]?.[name];
    return p ? !p(node) : false;
  } catch {
    return false;
  }
}

/** Inputs hidden EVERYWHERE — the panel filters by this; the node body additionally
 *  hides {@link nodeHiddenInputsOf}. */
export function hiddenInputsOf(func: DG.Func): ReadonlySet<string> {
  return hiddenNamesOf(HIDDEN_FUNC_INPUTS, func);
}

/** The input names the node body hides: fully hidden ones plus panel-only ones. */
export function nodeHiddenInputsOf(func: DG.Func): ReadonlySet<string> {
  const s = new Set(hiddenNamesOf(HIDDEN_FUNC_INPUTS, func));
  for (const n of hiddenNamesOf(PANEL_ONLY_FUNC_INPUTS, func)) s.add(n);
  return s;
}

export function hiddenOutputsOf(func: DG.Func): ReadonlySet<string> {
  return hiddenNamesOf(HIDDEN_FUNC_OUTPUTS, func);
}

function hiddenNamesOf(map: Record<string, Record<string, boolean>>, func: DG.Func): ReadonlySet<string> {
  try {
    const entry = map[func.nqName];
    return new Set(entry ? Object.keys(entry).filter((k) => entry[k]) : []);
  } catch {
    return new Set(); // Dart proxy access can throw
  }
}

/** A custom editor for one parameter of one node; the panel assigns {@link onChanged}
 *  after construction and the editor must invoke it on every user edit. */
export interface CustomInputEditor {
  element: HTMLElement;
  getValue: () => unknown;
  setValue: (v: unknown) => void;
  /** When present and false, the pending value is not stored. */
  isValid?: () => boolean;
  onChanged?: (v: unknown) => void;
  /** Release anything that outlives the DOM; called before the panel discards the editors. */
  detach?: () => void;
}

/** What a custom editor may know about the node it edits — passed in (rather than the
 *  editor reaching into the graph) so the factories stay unit-testable. */
export interface CustomEditorContext {
  /** Current stored value of a sibling input on the same node. */
  inputValue: (name: string) => unknown;
  /** Columns of the table wired to `tableParam`, or null. Never runs the flow. */
  columns: (tableParam: string) => DG.Column[] | null;
  /** The captured table itself, same never-runs-the-flow rule as {@link columns}. */
  table?: (tableParam: string) => DG.DataFrame | null;
  /** Distinguishes "connect a table" from "run the flow up to here". */
  isConnected?: (tableParam: string) => boolean;
  /** Materialize `tableParam` by running the slice up to its source — only ever from an
   *  explicit user action, never from a render. */
  produceTable?: (tableParam: string) => Promise<DG.DataFrame | null>;
  /** Re-run `cb` when a sibling input is edited — the panel does NOT re-render while focus is inside it. */
  watch: (inputName: string, cb: (value: unknown) => void) => void;
  /** The node being edited. Cache editor-computed verdicts keyed by node IDENTITY, never in
   *  `node.properties` — those serialize, and a background check landing would dirty the flow. */
  node?: FlowNode;
}

export type CustomInputEditorFactory = (param: DG.Property, ctx: CustomEditorContext) => CustomInputEditor;

/** A factory per overridden input (an editor holds a live HTMLElement, so it is built per node/render). */
export const CUSTOM_FUNC_INPUT_EDITORS: Record<string, Record<string, CustomInputEditorFactory>> = {
  'core:OpenFile': {fullPath: filePathEditor},
  'Chem:mpoScoreByProfile': {columnMapping: mpoColumnMappingEditor},
  // Takes precedence over the `EDITOR_SHORTCUT_INPUTS` pencil — the panel branches to a
  // custom editor before it ever builds a DG input to hang options on.
  'core:AddNewColumn': {expression: columnFormulaEditor},
  'Flow:filterRows': {condition: rowConditionEditor},
  'Flow:deleteRows': {condition: rowConditionEditor},
  'Flow:extractRows': {condition: rowConditionEditor},
  'Flow:selectRows': {condition: rowConditionEditor},
  'Flow:expressionToColumn': {expression: columnFormulaEditor},
  'Flow:aggregate': {aggregations: aggregationEditor},
};

/** Extra per-function readiness, returned as the labels the "Needs input" hint lists.
 *  Must be SYNCHRONOUS (runs on every render and in the run gate) and fail OPEN when it cannot decide. */
export type FuncNodeValidator = (node: FlowNode) => string[];

export const FUNC_NODE_VALIDATORS: Record<string, FuncNodeValidator> = {
  'Chem:mpoScoreByProfile': mpoMappingRequirements,
  // `[{type: 'avg', column: ''}]` is a non-blank string that still aggregates nothing,
  // so the generic blank check passes while the node isn't ready.
  'Flow:aggregate': aggregateRequirements,
  // The verdict is computed asynchronously by the formula editor — see
  // `expressionRequirements` for how it reaches this synchronous gate.
  'Flow:filterRows': expressionRequirements('condition'),
  'Flow:deleteRows': expressionRequirements('condition'),
  'Flow:extractRows': expressionRequirements('condition'),
  'Flow:selectRows': expressionRequirements('condition'),
  'Flow:expressionToColumn': expressionRequirements('expression'),
  'core:AddNewColumn': expressionRequirements('expression'),
};

/** Only the aggregation list is checked — empty group-by/pivot are legitimate. */
function aggregateRequirements(node: FlowNode): string[] {
  if (node.editorBridge?.isSocketConnected(node.id, 'input', 'aggregations'))
    return []; // fed by a wire — its value isn't in `inputValues` to inspect
  const problems = aggregationProblems(parseAggregations(node.inputValues['aggregations']));
  return problems.length === 0 ? [] : [`Aggregations — needs ${problems.join(', ')}`];
}

/** The registered validator, wrapped so it can never throw into a render. */
export function funcValidatorOf(func: DG.Func): FuncNodeValidator | undefined {
  let validator: FuncNodeValidator | undefined;
  try {
    validator = FUNC_NODE_VALIDATORS[func.nqName];
  } catch {
    return undefined; // nqName can throw on odd Dart proxies
  }
  if (!validator) return undefined;
  return (node) => {
    try {
      return validator!(node);
    } catch (e) {
      console.error(`Flow: readiness check failed for ${node.dgFuncName}`, e);
      return []; // fail open — never block a run over a broken check
    }
  };
}

export function customEditorFor(func: DG.Func, inputName: string): CustomInputEditorFactory | null {
  try {
    return CUSTOM_FUNC_INPUT_EDITORS[func.nqName]?.[inputName] ?? null;
  } catch {
    return null;
  }
}

/** One input exposed by a {@link FuncWrapper} in place of the function's own. */
export interface WrappedFuncInput {
  name: string;
  type: string;
  caption?: string;
  description?: string;
  /** Filters the column picker and routes a `string` slot to the semantic-type editor. */
  semType?: string;
  /** Optional inputs don't gate the run (no "Needs input"). */
  optional?: boolean;
  defaultValue?: unknown;
  choices?: string[];
}

/** Reshapes a function's signature into Flow-friendly node inputs; the compiler runs the
 *  resolved input expressions through {@link mapInputs} to build the real call arguments. */
export interface FuncWrapper {
  inputs: WrappedFuncInput[];
  /** Exposed-input JS expressions (a key is absent when unconnected and blank) →
   *  the function's real named arguments, also JS expressions. */
  mapInputs: (exposed: Record<string, string>) => Record<string, string>;
}

/** Functions kept in the catalog but exposed through reshaped inputs. */
export const FUNC_WRAPPERS: Record<string, FuncWrapper> = {
  // `tables: dataframe_list` is unwirable — expose two plain tables and fold them into the list.
  'core:AppendTables': {
    inputs: [
      {name: 'table1', type: 'dataframe', description: 'Defines the result columns'},
      {name: 'table2', type: 'dataframe', description: 'Appended below the first table'},
    ],
    mapInputs: (v) => v.table1 && v.table2 ? {tables: `[${v.table1}, ${v.table2}]`} : {} as Record<string, string>,
  },
};

export function funcWrapperOf(func: DG.Func): FuncWrapper | null {
  try {
    return FUNC_WRAPPERS[func.nqName] ?? null;
  } catch {
    return null;
  }
}

/** The parameter list a wrapped function PRESENTS. Every surface that renders or reasons
 *  about a node's parameters must go through this, or the node and the panel drift. */
export function effectiveFuncInputs(func: DG.Func): DG.Property[] {
  const wrapper = funcWrapperOf(func);
  return wrapper ? wrapperProperties(wrapper) : func.inputs;
}

/** The exposed inputs as real `DG.Property` objects. `nullable: false` unless declared optional —
 *  the base Property defaults to nullable, which would kill the "Needs input" gate. */
export function wrapperProperties(wrapper: FuncWrapper): DG.Property[] {
  return wrapper.inputs.map((s) => DG.Property.fromOptions({
    name: s.name, type: s.type, caption: s.caption, description: s.description,
    nullable: s.optional === true,
    ...(s.semType !== undefined ? {semType: s.semType} : {}),
    ...(s.choices !== undefined ? {choices: s.choices} : {}),
    ...(s.defaultValue !== undefined ? {defaultValue: s.defaultValue as string} : {}),
  }));
}

// NOTE: don't wrap column-only functions — a bare `column` has no table to resolve against;
// write a proper (table, column) twin instead (Chem:filterBySubstructure is the pattern).

/** File picker for string path parameters; stores the plain full-path string. `setValue` rebuilds
 *  a FileInfo synchronously rather than routing through `stringValue`, whose async server
 *  resolution is what validates hand-typed paths. */
function filePathEditor(param: DG.Property): CustomInputEditor {
  const ed: CustomInputEditor = {} as CustomInputEditor;
  let current = '';
  const input = ui.input.file(getParamDisplayName(param), {
    onValueChanged: (v) => {
      current = v ? v.fullPath : '';
      ed.onChanged?.(current);
    },
    // temporary thing, remove local file opening. once we figure out how to handle local files, remove this
    onCreated: (a) => a.root.querySelector('.ui-input-options')?.querySelector('.fa-folder-open')?.remove?.(),
  });
  // A hand-typed path must commit on change/Enter: DG's FileInput reports a value only once
  // the server resolves it, so a typed-but-unresolved path was silently dropped.
  const text = input.root.querySelector('input[type="text"], input:not([type])') as HTMLInputElement | null;
  text?.addEventListener('change', () => {
    const typed = text.value.trim();
    if (typed === current) return;
    current = typed;
    ed.onChanged?.(typed);
  });
  ed.element = input.root;
  ed.getValue = (): unknown => current;
  ed.setValue = (v): void => {
    current = v === undefined || v === null ? '' : String(v);
    try {
      if (current !== '')
        input.value = DG.FileInfo.fromString(current, '');
    } catch {/* leave the editor blank */}
  };
  return ed;
}
