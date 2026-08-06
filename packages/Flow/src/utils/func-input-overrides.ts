/** Per-function input overrides — curated maps keyed by `Func.nqName` (same
 *  convention as `INCLUDED_FUNC_NQNAMES`), consulted when a func node and its
 *  parameter form are built:
 *
 *  - {@link HIDDEN_FUNC_INPUTS} — inputs users should never see: no socket row
 *    on the node, no editor on the context panel. Purely visual — the slot,
 *    its seeded value, compilation, and creation-script import/emit stay
 *    untouched (data-sync scripts legitimately carry e.g.
 *    `subscribeOnChanges = true` and must round-trip). Only hide inputs with a
 *    declared default.
 *  - {@link PANEL_ONLY_FUNC_INPUTS} — inputs kept off the node body (no socket
 *    row, no pass-through row — unless wired) but still edited on the context
 *    panel. For configuration-like parameters (a file path, a sheet name) that
 *    would otherwise bloat the node card.
 *  - {@link HIDDEN_FUNC_OUTPUTS} — the mirror image: declared outputs that are
 *    bookkeeping rather than a result (e.g. the *names* of the columns a
 *    function appended in place). No socket, no strip binding — the real result
 *    is the mutated table, which the pass-through already carries.
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
import {mpoColumnMappingEditor, mpoMappingRequirements} from '../panel/editors/mpo-mapping-editor';
import {columnFormulaEditor, expressionRequirements, rowConditionEditor} from '../panel/editors/expression-editor';
import {aggregationEditor} from '../panel/editors/aggregation-editor';
import {aggregationProblems, parseAggregations} from '../ops/data-ops';
// Type-only (erased at build time), so the scheme ↔ overrides pair stays free
// of a runtime import cycle.
import type {FlowNode} from '../rete/scheme';

/** `{[func.nqName]: {[inputName]: true}}` — inputs hidden from the node and
 *  the context panel (visually only — see the module doc). Edit freely: flag
 *  an input to hide it everywhere. */
export const HIDDEN_FUNC_INPUTS: Record<string, Record<string, boolean>> = {
  // Reactivity/error plumbing for the formula engine — never user-relevant.
  'core:AddNewColumn': {subscribeOnChanges: true, errorBehavior: true},
  'core:AddNewColumnList': {subscribeOnChanges: true, errorBehavior: true},
  // `join` (default true) is what puts the converted column into the table, and
  // `overwrite` (default false) would replace the source column — on a canvas
  // both must stay at their defaults or the node's output carries nothing.
  'Chem:convertNotation': {join: true, overwrite: true},
  'Chem:recalculateCoords': {join: true},
};

/** `{[func.nqName]: {[inputName]: true}}` — inputs with no socket row on the
 *  node body but a regular editor on the context panel. Same visual-only
 *  contract as {@link HIDDEN_FUNC_INPUTS} (slot, seeded value, compilation,
 *  and script round-trips untouched; a wired slot renders its row), except the
 *  panel keeps the input editable. */
export const PANEL_ONLY_FUNC_INPUTS: Record<string, Record<string, boolean>> = {
  // A path and a sheet name are configuration, not data flow — editing them in
  // the panel keeps the node a compact "this file" card.
  'core:OpenFile': {fullPath: true, sheetName: true},
};

/** `{[func.nqName]: {[outputName]: true}}` — declared outputs the node should
 *  not expose. Same contract as {@link HIDDEN_FUNC_INPUTS}: visual only, so
 *  compilation and script round-trips are untouched. */
export const HIDDEN_FUNC_OUTPUTS: Record<string, Record<string, boolean>> = {
  // Returns the NAMES of the columns it appended to `table`, not data. The
  // result users want is the mutated table, carried by the `table` pass-through.
  'Chem:runElementalAnalysis': {res: true},
};

/** `{[func.nqName]: {[inputName]: caption}}` — display captions for
 *  parameters whose declared names read as implementation ("Keys1"). Purely
 *  cosmetic: slot keys, stored values, and compilation stay on the real
 *  parameter names. */
export const FUNC_INPUT_CAPTIONS: Record<string, Record<string, string>> = {
  'core:JoinTables': {
    table1: 'Left table', table2: 'Right table',
    keys1: 'Left keys', keys2: 'Right keys',
    values1: 'Columns from left', values2: 'Columns from right',
  },
};

/** The registered display caption for a function input, or null. */
export function inputCaptionOf(func: DG.Func, name: string): string | null {
  try {
    return FUNC_INPUT_CAPTIONS[func.nqName]?.[name] ?? null;
  } catch {
    return null;
  }
}

/** `{[func.nqName]: {[inputName]: predicate}}` — inputs whose panel row only
 *  renders while the predicate holds against the node's current values. The
 *  slot and stored value stay untouched (the compiler applies its own guard). */
export const CONDITIONAL_FUNC_INPUTS: Record<string, Record<string, (node: FlowNode) => boolean>> = {
  // A sheet name only means anything for an Excel file — for every other
  // extension the row is noise (and the compiler drops the value anyway).
  'core:OpenFile': {sheetName: (node) => /\.xlsx$/i.test(String(node.inputValues['fullPath'] ?? ''))},
};

/** True when a registered visibility predicate says the panel should NOT
 *  render this input right now. Fails open (row shown) on any error. */
export function inputHiddenByCondition(func: DG.Func, name: string, node: FlowNode): boolean {
  try {
    const p = CONDITIONAL_FUNC_INPUTS[func.nqName]?.[name];
    return p ? !p(node) : false;
  } catch {
    return false;
  }
}

/** The hidden-input names for a function, `∅` when none are registered
 *  (or the Dart proxy throws on `nqName`). Hidden EVERYWHERE — the panel
 *  filters by this set; the node body additionally hides
 *  {@link nodeHiddenInputsOf}. */
export function hiddenInputsOf(func: DG.Func): ReadonlySet<string> {
  return hiddenNamesOf(HIDDEN_FUNC_INPUTS, func);
}

/** The input names the NODE BODY hides: fully hidden ones plus panel-only
 *  ones. Feeds `FlowNode.hiddenInputs` (consumed by the node component). */
export function nodeHiddenInputsOf(func: DG.Func): ReadonlySet<string> {
  const s = new Set(hiddenNamesOf(HIDDEN_FUNC_INPUTS, func));
  for (const n of hiddenNamesOf(PANEL_ONLY_FUNC_INPUTS, func)) s.add(n);
  return s;
}

/** The hidden-output names for a function, `∅` when none are registered. */
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
  /** Release anything that outlives the DOM — a subscription, a hosted widget.
   *  Called by the panel before it throws the editors away (a new node shown, a
   *  re-render); an editor that only builds DOM needs no implementation. */
  detach?: () => void;
}

/** What a custom editor is allowed to know about the node it edits.
 *
 *  Some parameters can only be edited against the rest of the node: an MPO
 *  column mapping needs the chosen profile AND the columns of the upstream
 *  table. Passing this in (rather than letting an editor reach into the editor
 *  graph) keeps the factories pure enough to unit-test. */
export interface CustomEditorContext {
  /** Current stored value of a sibling input on the same node. */
  inputValue: (name: string) => unknown;
  /** Columns of the table wired to `tableParam`, or null when it isn't
   *  connected or hasn't been computed. Never runs the flow — an editor must
   *  degrade gracefully rather than trigger work on a panel render. */
  columns: (tableParam: string) => DG.Column[] | null;
  /** The captured table itself, under the same never-runs-the-flow rule as
   *  {@link columns} — for editors that need the data, not just the schema
   *  (the formula editor autocompletes and validates against a real table). */
  table?: (tableParam: string) => DG.DataFrame | null;
  /** Whether `tableParam` is wired at all. Distinguishes "connect a table"
   *  from "run the flow up to here", which are different instructions. */
  isConnected?: (tableParam: string) => boolean;
  /** Materialize `tableParam` by running the slice up to its source. **Only
   *  ever from an explicit user action** (a button) — never from a render, or
   *  opening the panel would start computing. */
  produceTable?: (tableParam: string) => Promise<DG.DataFrame | null>;
  /** Re-run `cb` whenever a sibling input is edited. Needed because the panel
   *  does NOT re-render while focus is inside it — which is exactly when a
   *  combo the editor depends on gets changed. */
  watch: (inputName: string, cb: (value: unknown) => void) => void;
  /** The node being edited — for editors that produce something a
   *  {@link FuncNodeValidator} has to read back later (a validity verdict only
   *  the editor can compute). Cache such answers keyed by node IDENTITY, never
   *  in `node.properties`: those serialize, and a background check landing
   *  would dirty a flow nobody edited. */
  node?: FlowNode;
}

export type CustomInputEditorFactory = (param: DG.Property, ctx: CustomEditorContext) => CustomInputEditor;

/** `{[func.nqName]: {[inputName]: factory}}` — a factory per overridden input
 *  (an editor holds a live HTMLElement, so it must be built per node/render). */
export const CUSTOM_FUNC_INPUT_EDITORS: Record<string, Record<string, CustomInputEditorFactory>> = {
  'core:OpenFile': {fullPath: filePathEditor},
  'Chem:mpoScoreByProfile': {columnMapping: mpoColumnMappingEditor},
  // The platform's own formula node gets the platform's own formula editor —
  // inline, instead of a bare text box with a pencil that opened the same thing
  // behind a dialog. Same editor Flow's `expressionToColumn` uses; it takes
  // precedence over the `EDITOR_SHORTCUT_INPUTS` pencil (the panel branches to
  // a custom editor before it ever builds a DG input to hang options on).
  'core:AddNewColumn': {expression: columnFormulaEditor},
  // Flow's own data operations: a row condition is a boolean formula and an
  // aggregation list is a list — neither is a thing to type into a text box.
  'Flow:filterRows': {condition: rowConditionEditor},
  'Flow:deleteRows': {condition: rowConditionEditor},
  'Flow:extractRows': {condition: rowConditionEditor},
  'Flow:selectRows': {condition: rowConditionEditor},
  'Flow:expressionToColumn': {expression: columnFormulaEditor},
  'Flow:aggregate': {aggregations: aggregationEditor},
};

// ---------- per-function readiness ----------

/** Extra, function-specific requirements a node must satisfy before it can run,
 *  returned as the labels the "Needs input" hint should list (empty = ready).
 *
 *  The generic checks only see whether a slot is wired or non-blank, which
 *  can't express "this JSON must cover every property of the chosen profile".
 *  Must be SYNCHRONOUS (it runs on every node render and in the run gate) and
 *  must fail OPEN when it cannot decide — blocking a run over a cold cache
 *  would be worse than letting the function report the problem itself. */
export type FuncNodeValidator = (node: FlowNode) => string[];

export const FUNC_NODE_VALIDATORS: Record<string, FuncNodeValidator> = {
  'Chem:mpoScoreByProfile': mpoMappingRequirements,
  // `aggregations` is a non-empty string the moment the editor writes one entry,
  // so the generic "is it blank" check passes while the list still reads
  // `[{type: 'avg', column: ''}]` — an aggregation over nothing.
  //
  // Nothing similar is needed for the blank-is-enough cases (`unpivot`'s
  // `mergeColumns`, `deleteColumns`' `columns`, `tagColumns`' `tag`): they are
  // declared `nullable: false`, so `FuncNode.requiredInputs` covers them —
  // and, unlike a validator, that check knows whether the slot is WIRED. A
  // validator reading `inputValues` alone would flag a slot fed by a Column
  // List Input node, which is why this one guards on the connection too.
  'Flow:aggregate': aggregateRequirements,
  // A condition that isn't a true/false expression, or names a column that
  // isn't there, is a run-time failure the formula editor already caught while
  // the user typed. Blocking here turns it into a hint instead of a red node.
  // The check itself is asynchronous (deciding a formula's type means
  // evaluating it) — see `expressionRequirements` for how the verdict reaches
  // this synchronous gate, and why it fails open when it hasn't.
  'Flow:filterRows': expressionRequirements('condition'),
  'Flow:deleteRows': expressionRequirements('condition'),
  'Flow:extractRows': expressionRequirements('condition'),
  'Flow:selectRows': expressionRequirements('condition'),
  'Flow:expressionToColumn': expressionRequirements('expression'),
  'core:AddNewColumn': expressionRequirements('expression'),
};

/** Aggregate is ready when every aggregation names something to aggregate.
 *  Group-by and pivot are legitimately empty (folding a whole table into one
 *  row is a real thing to want), so only the list is checked. */
function aggregateRequirements(node: FlowNode): string[] {
  if (node.editorBridge?.isSocketConnected(node.id, 'input', 'aggregations'))
    return []; // fed by a wire — its value isn't in `inputValues` to inspect
  const problems = aggregationProblems(parseAggregations(node.inputValues['aggregations']));
  return problems.length === 0 ? [] : [`Aggregations — needs ${problems.join(', ')}`];
}

/** The registered validator for a function, wrapped so it can never throw into
 *  a render, or undefined when the function has none. */
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
  /** Filters the column picker, and routes a `string` slot to the registered
   *  editor for that semantic type (e.g. Chem's molecule sketcher). */
  semType?: string;
  /** Optional inputs don't gate the run (no "Needs input"). */
  optional?: boolean;
  defaultValue?: unknown;
  choices?: string[];
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

/** The parameter list a wrapped function PRESENTS — the wrapper's exposed
 *  inputs when one is registered, else the function's own.
 *
 *  Every surface that renders or reasons about a node's parameters must go
 *  through this, or the node and the panel drift: the node builds its sockets
 *  from the wrapper while the panel would iterate the raw signature, showing
 *  different captions, no column picker (the raw signature may have no
 *  dataframe input to resolve against), and editors for parameters the node
 *  folds away. */
export function effectiveFuncInputs(func: DG.Func): DG.Property[] {
  const wrapper = funcWrapperOf(func);
  return wrapper ? wrapperProperties(wrapper) : func.inputs;
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
    ...(s.semType !== undefined ? {semType: s.semType} : {}),
    ...(s.choices !== undefined ? {choices: s.choices} : {}),
    ...(s.defaultValue !== undefined ? {defaultValue: s.defaultValue as string} : {}),
  }));
}

/** NOTE on wrapping column-only functions: don't. A function taking a bare
 *  `column` has no table for Flow to resolve the column against, and bolting a
 *  synthetic `table` input on with a wrapper only papers over it — the
 *  parameter the user fills isn't the parameter the function declares, and the
 *  semantic types and captions live in the wrong place. Write a proper twin
 *  taking `(table, column, …)` with its semTypes declared instead; Chem's
 *  `filterBySubstructure` / `similarityTo` / `diverseSubset` are the pattern. */

/** A file picker (`ui.input.file`) for string path parameters like OpenFile's
 *  `fullPath`. The stored value is the plain full-path string
 *  (`System:AppData/...`) — exactly what the parameter takes — read from the
 *  picked `FileInfo`. `setValue` rebuilds a `FileInfo` from the stored path
 *  (synchronous — no server round-trip) rather than routing through
 *  `stringValue`, whose async server resolution is what validates hand-typed
 *  paths. */
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
  // A hand-typed path must commit on change/Enter like every other string row.
  // DG's FileInput reports a value only once the server RESOLVES it, so a
  // typed-but-unresolved path was silently dropped: the panel reverted on the
  // next refresh and the following run reported success on the OLD file.
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
