/* The func-call editor: a form generated from a FuncCall's input params, writing through
   `setParamValue` and following `param.onChanged` back —
   the u2 counterpart of `DG.InputForm.forFuncCall`. Structural over {@link FuncCallLike}, so
   headless doubles and DG.FuncCall both fit; the platform itself is only reached by
   `editors: 'auto'`, through the global the bundler binds `datagrok-api/dg` to.
   Wave 2 adds the dynamic routes — async choices/suggestions through the `evalParam*` members
   and computed defaults written into the call (R6) — with `settled` as the readiness member.
   Wave 3 adds the table routes — dataframe/column/column_list fields with auto-fill, the
   default-column pick and dependent rewiring (param-tables.ts), all synchronous: the call holds
   OBJECTS, the inputs hold names. */
import {batch} from '../../core/signals.js';
import type {IProperty} from '../../core/property-like.js';
import {Input} from '../../core/input-base.js';
import type {InputOptions} from '../../core/input-base.js';
import {div} from '../../core/elements.js';
import {SuggestionList} from '../../core/suggestion-list.js';
import {Form} from '../../components/forms/form.js';
import {ChoiceInput, MultiChoiceInput} from '../../components/inputs/choice-input.js';
import {SuggestInput} from '../../components/inputs/suggest-input.js';
import {ObjectForm, inputForProperty, kindOf} from '../forms/object-form.js';
import type {FieldOverride, Kind} from '../forms/object-form.js';
import {Editors} from '../forms/editors.js';
import {ParamSource, ParamState} from './param-sources.js';
import type {ChoicesResult} from './param-sources.js';
import {tableInput} from '../inputs/pickers.js';
import {ColumnInput} from '../inputs/column-combo.js';
import {ColumnsInput} from '../inputs/columns.js';
import {TableBinding, columnPredicate, parentTableName, tableByName, preferredTable,
  asNamed, asTable, resolveTable, columnNamesOf, noMatchMessage,
  markAuto, unmarkAuto} from './param-tables.js';
import type {ColumnDependent} from './param-tables.js';
import type * as DG from 'datagrok-api/dg';

export interface FuncFormOptions {
  /** Two-way is the u2 default. `false` is honored for DG.InputForm contract
   * compatibility only: FuncCall edits then do not refresh the form. */
  twoWayBinding?: boolean;
  /** Suppresses computed-default evaluation and `propagateChoice` writes (R6). Dynamic choice
   * ITEMS still load and apply — items are not defaults (divergence #8). Literal defaults stay
   * display-only either way. */
  skipDefaultInit?: boolean;
  /** Suppresses the table auto-fill write (the current-or-first open table into a null
   * dataframe param, fpe:62-64) — the js-api option (forms.ts:226) reproduced. */
  skipTableAutoFill?: boolean;
  /** `false` hides the table-field roots (the `getEditor(condensed, showTableSelectors)`
   * contract, functions.ts:412); auto-fill and dependent wiring still run (fpe:73). */
  showTableSelectors?: boolean;
  condensed?: boolean;
  /** Fires on every non-echo input change — user or programmatic input write — after the
   * FuncCall has been updated. Mirrors DG.InputForm.onInputChanged. */
  onInputChanged?: (name: string, value: unknown) => void;
  /** Per-param editor overrides / options, keyed by param name (same shape as ObjectForm). */
  overrides?: Record<string, FieldOverride>;
  editors?: 'auto' | 'generated';
}

/** The structural slice funcForm reads off a func-call parameter — DG.FuncCallParam and the
 * headless doubles both satisfy it. */
export interface FuncCallParamLike {
  readonly name: string;
  /** Optional-default substitution included (`func_call_param.dart:78-88`). */
  readonly value: any;
  readonly property: IProperty & {options?: Record<string, any>, dart?: unknown};
  onChanged: {subscribe(fn: (p: FuncCallParamLike) => void): {unsubscribe(): void}};
}

// The eval members are OPTIONAL: DG.FuncCall has them; a call without them renders
// dynamic-routed params with a visible error state, never a TypeError.
export interface FuncCallLike {
  inputParams: {values(): Iterable<FuncCallParamLike>};
  setParamValue(name: string, value: any): void;
  evalParamChoices?(name: string): Promise<ChoicesResult>;
  evalParamSuggestions?(name: string, text: string): Promise<{items: string[],
    tooltips: Record<string, string>}>;
  evalParamDefault?(name: string): Promise<any>;
}

/** The per-field seam the dynamic sources plug into. Null kind marks a platform editor, which
 * reads and writes the param in its own native type. Orphaned: the current `source` carries no
 * param of this name — the field is disabled and its write path is off. */
export interface FuncField {
  param: FuncCallParamLike;
  input: Input<any>;
  kind: Kind | null;
  /** The dynamic route the field was built for; 'field' for every static editor, including a
   * dynamic-routed param an override or a registered editor claimed (no dynamic wiring then). */
  route: 'field' | 'choices' | 'suggestions' | 'multiChoices' | 'table' | 'column' | 'columns';
  orphaned: boolean;
  propagate: boolean;
  /** True while a choices result is being applied to the input — the batched flush of the
   * setItems write runs after {@link Input.system}'s counter has dropped, so the propagate
   * guard needs its own flag. */
  applying: boolean;
  lookup?: Record<string, Record<string, any>> | null;
  state?: ParamState;
  source?: ParamSource;
  /** W3 converters: the call holds OBJECTS (a raw name write comes back as a pending Resolve\*
   * FuncCall, P-W3-2), the input holds names. Installed only by the table/column/columns
   * routes. */
  fromParam?: (v: any) => any;
  toParam?: (v: any) => any;
  /** The table param a column/columns field resolves its objects against (FP-W3-2's
   * association). */
  parentName?: string;
  /** The `ColumnFilter.fromProp`-parity predicate, built once per column/columns field and
   * shared by its input and its {@link TableBinding} dependents. */
  filter?: (c: DG.Column) => boolean;
  /** The `auto` badge marking a value the binding guessed (table auto-fill, column auto-pick);
   * created once by {@link markAuto} and kept for re-marking after {@link unmarkAuto}. */
  autoBadge?: HTMLElement;
  /** Set by the first user edit of the field — the guess mark never returns after it. */
  userTouched?: boolean;
  /** Dismisses the field's transient notice and its interaction listeners. */
  noticeClear?: () => void;
}

/** The W1 scalar set; everything else lands in {@link FuncCallForm.unsupported}. */
const SCALARS = new Set(['string', 'int', 'double', 'float', 'num', 'bigint', 'qnum', 'bool',
  'datetime']);
/** The ib editor hints the property form honors (`input_base.dart:702-728`). The Dart FUNC form
 * errors on them (it evals them as function editors), so honoring them here is a documented W1
 * divergence, not a compatibility break. */
const HINTS = new Set(['textarea', 'password', 'switch', 'slider']);
const CONJUNCTIONS = new Set(['and', 'or', 'than', 'if', 'but', 'so', 'as', 'that']);

/** prop_gen's `camelCaseToWords` (`prop_gen_annotation.dart:89-112`): all-caps and already-spaced
 * names pass through, humps split, first word capitalized, conjunctions lowercased. */
function camelCaseToWords(name: string): string {
  if (name === name.toUpperCase() || name.includes(' '))
    return name;
  const words = name.match(/[A-Z]+(?![a-z])|[A-Z]?[^A-Z]+/g) ?? [name];
  return words
    .map((w, i) => {
      const word = CONJUNCTIONS.has(w.toLowerCase()) ? w.toLowerCase() : w;
      return i === 0 ? word.charAt(0).toUpperCase() + word.slice(1) : word;
    })
    .join(' ');
}

/** A {@link Form} over a FuncCall's input params, kept in sync with the call both ways: an edit
 * runs `setParamValue`, a `param.onChanged` refreshes the field (echo-suppressed by value in both
 * directions). It edits the call and nothing else — running it stays with the caller. */
export class FuncCallForm extends Form {
  private readonly _formOptions: FuncFormOptions;
  private readonly _fields: FuncField[] = [];
  private readonly _unsupported: string[] = [];
  private _call: FuncCallLike;
  private readonly _autoFilled = new Set<string>();
  private _paramSubs: (() => void)[] = [];
  private _tableBindings: TableBinding[] = [];
  private _refreshing = false;
  private _generation = 0;
  private _settled: Promise<void> = Promise.resolve();

  constructor(call: FuncCallLike, options: FuncFormOptions = {}) {
    super({condensed: options.condensed});
    this._formOptions = options;
    this._call = call;
    this.root.dataset.u2 = 'func-form';
    this._build();
    this._bind();
    // before _arm: the auto-pick writes are sync and must land before the default evals read
    this._bindTables();
    this._arm();
  }

  /** Resolves when the initial async population of the CURRENT source has landed (each dynamic
   * choices source ready-or-error + each computed default written-or-failed). Re-arms on
   * `source` rebind. Never rejects. The W6 flip shim awaits it to honor the Dart form's
   * awaits-initial-choice-loads contract (fpe:862-867). */
  get settled(): Promise<void> {
    return this._settled;
  }

  /** The FuncCall being edited. Setting it rebinds every field to the new call, refreshes
   * values from it, and re-subscribes onChanged — the DG.InputForm `source` contract. A field
   * whose param the new call does not carry is disabled and stops writing, until a later call
   * carries the name again. */
  get source(): FuncCallLike {
    return this._call;
  }

  set source(call: FuncCallLike) {
    this._call = call;
    this._generation++;
    for (const off of this._paramSubs) {
      off();
      this.scope.disown(off);
    }
    this._paramSubs = [];
    const params = new Map<string, FuncCallParamLike>();
    const tables: FuncCallParamLike[] = [];
    this._unsupported.length = 0;
    for (const param of call.inputParams.values()) {
      params.set(param.name, param);
      const route = FuncCallForm._route(param.property);
      if (route === 'table')
        tables.push(param);
      if (route !== 'skip' && !this._fields.some((f) => f.param.name === param.name))
        this._unsupported.push(param.name);
    }
    this._autoFill(tables);
    this._refreshing = true;
    try {
      batch(() => {
        for (const field of this._fields) {
          const param = params.get(field.param.name);
          const orphaned = param === undefined;
          if (orphaned !== field.orphaned) {
            field.orphaned = orphaned;
            this._applyEnabled(field);
          }
          if (param === undefined)
            continue;
          field.param = param;
          field.input.value.value = this._initialValue(field);
        }
      });
    } finally {
      this._refreshing = false;
    }
    this._bind();
    for (const field of this._fields) {
      field.state?.set(1, {kind: 'idle'});
      field.noticeClear?.();
      if (field.route === 'table')
        (this._autoFilled.has(field.param.name) ? markAuto : unmarkAuto)(field);
      if (field.route !== 'choices' && field.route !== 'multiChoices')
        continue;
      if (field.orphaned) {
        field.source?.dispose();
        field.source = undefined;
        field.state?.set(0, {kind: 'idle'});
      }
      else
        this._wireSource(field);
    }
    this._bindTables();
    this._arm();
  }

  /** The input editing a param, keyed by PARAM NAME (never by caption — a deliberate
   * divergence from the Dart form's caption-keyed lookups). */
  getInput(name: string): Input<any> | undefined {
    return this._fields.find((f) => f.param.name === name)?.input;
  }

  get isValid(): boolean {
    return this.validity.peek() === null;
  }

  validateInputs(): boolean {
    return this.isValid;
  }

  /** Param names of the current call the form cannot edit (non-scalar types, nested editors,
   * layout; after a `source` rebind, also any new-call param the form has no field for) —
   * the seam W6's per-capability fallback reads. */
  get unsupported(): ReadonlyArray<string> {
    return this._unsupported;
  }

  private _build(): void {
    const params = [...this._call.inputParams.values()];
    const routes = new Map<FuncCallParamLike, ReturnType<typeof FuncCallForm._route>>();
    const tables = new Map<string, FuncCallParamLike>();
    for (const param of params) {
      const route = FuncCallForm._route(param.property);
      routes.set(param, route);
      if (route === 'table')
        tables.set(param.name, param);
    }
    this._autoFill(tables.values());
    const categories: {name: string, entries: {param: FuncCallParamLike,
      route: FuncField['route'], parent?: FuncCallParamLike}[]}[] = [];
    for (const param of params) {
      let route = routes.get(param)!;
      if (route === 'skip')
        continue;
      let parent: FuncCallParamLike | undefined;
      if (route === 'column' || route === 'columns') {
        parent = tables.get(parentTableName(param.property) ?? '');
        // divergence #13: an unassociated column param is LISTED, where Dart silently renders
        // nothing for it (renderParam falls through every branch)
        if (parent === undefined)
          route = 'unsupported';
      }
      if (route === 'unsupported') {
        this._unsupported.push(param.name);
        continue;
      }
      const name = param.property.category ?? 'Misc';
      let category = categories.find((c) => c.name === name);
      if (category === undefined) {
        category = {name, entries: []};
        categories.push(category);
      }
      category.entries.push({param, route, parent});
    }
    const headers = categories.length > 1 ||
      (categories.length === 1 && categories[0].name !== 'Misc');
    for (const category of categories) {
      if (headers)
        this.addElement(div([category.name], 'u2-form-category'));
      for (const {param, route, parent} of category.entries)
        this._addField(param, route, parent);
    }
  }

  /** Dart's table auto-fill (fpe:62-64, minus the Dart-only `applicableTable`): a real write
   * into the call — the current table, else the first open one — before any field seeds its
   * initial value, so the field opens showing it. */
  private _autoFill(tables: Iterable<FuncCallParamLike>): void {
    this._autoFilled.clear();
    if (this._formOptions.skipTableAutoFill === true)
      return;
    for (const param of tables) {
      if (param.value != null)
        continue;
      const table = preferredTable();
      if (table != null) {
        this._call.setParamValue(param.name, table);
        this._autoFilled.add(param.name);
      }
    }
  }

  /** First match wins: SERVICE_PARAM skipped; `editorParam` mirrors unsupported; a present
   * `choices`/`suggestions` OPTION routes dynamic — the Dart discriminator (fpe:569-572, 645),
   * with choices winning over suggestions. Never `prop.choices`: annotation parsing leaves
   * garbage there for a dynamic source (scripting.dart:316), while the evaluator answers even
   * static list literals. A typed `prop.choices` with NO option keeps the W1 static path. */
  private static _route(prop: FuncCallParamLike['property']):
      'skip' | 'unsupported' | 'field' | 'choices' | 'suggestions' | 'multiChoices' |
      'table' | 'column' | 'columns' {
    const options = prop.options;
    const editor = options?.['editor'] ?? prop.editor ?? null;
    if (editor === 'none')
      return 'skip';
    if (options?.['editorParam'] != null)
      return 'unsupported';
    const type = prop.propertyType ?? prop.type ?? '';
    // guarded on a null editor: in Dart any non-null editor other than layout wins over the
    // type branches (fpe:452/497/519), so an editor-carrying dataframe/column param — and
    // `editor: columnsMap` (type map, P-W3-7) — stays unsupported
    if (editor === null) {
      if (type === 'dataframe')
        return 'table';
      if (type === 'column')
        return 'column';
      if (type === 'column_list')
        return 'columns';
    }
    if (type === 'string' && options?.['choices'] != null)
      return 'choices';
    if (type === 'string' && options?.['suggestions'] != null)
      return 'suggestions';
    if (type === 'list' && options?.['choices'] != null)
      return 'multiChoices';
    if (type === 'list')
      return prop.choices != null && prop.choices.length > 0 ? 'field' : 'unsupported';
    if (!SCALARS.has(type))
      return 'unsupported';
    return editor === null || HINTS.has(String(editor).toLowerCase()) ? 'field' : 'unsupported';
  }

  private _addField(param: FuncCallParamLike, routed: FuncField['route'],
    parent?: FuncCallParamLike): void {
    const prop = param.property;
    const {input: custom, ...rest} = this._formOptions.overrides?.[param.name] ?? {};
    const options: InputOptions<any> = {
      name: param.name,
      label: FuncCallForm._caption(prop, param.name),
      tooltipText: prop.description ?? undefined,
      nullable: prop.nullable,
      postfix: prop.units || prop.options?.['units'] || undefined,
      ...rest,
    };
    const registered = custom ? null : this.run(() => Editors.resolve(prop, options));
    // a custom editor owns its field entirely and gets no dynamic wiring (the W1 override contract)
    const route: FuncField['route'] = (custom ?? registered) != null ? 'field' : routed;
    const filter = route === 'column' || route === 'columns' ? columnPredicate(prop) : undefined;
    let input: Input<any>;
    let kind: Kind | null;
    if (route === 'choices') {
      input = this.run(() => new ChoiceInput({...options, items: []}));
      kind = 'choice';
    }
    else if (route === 'multiChoices') {
      input = this.run(() => new MultiChoiceInput({...options, items: []}));
      kind = 'list';
    }
    else if (route === 'suggestions') {
      input = this._suggestInput(param.name, options);
      kind = 'string';
    }
    else if (route === 'table') {
      const {label: caption, ...tableOptions} = options;
      input = this.run(() => tableInput(caption ?? param.name, tableOptions));
      kind = 'string';
    }
    else if (route === 'column') {
      input = this.run(() => new ColumnInput({...options,
        table: asTable(parent?.value), filter}));
      kind = 'string';
    }
    else if (route === 'columns') {
      input = this.run(() => new ColumnsInput({...options,
        table: asTable(parent?.value), filter}));
      kind = 'list';
    }
    else {
      const platform = (custom ?? registered) != null ? null :
        ObjectForm.platformInput(this, prop, this._call, this._formOptions.editors === 'auto');
      input = custom ?? registered ?? platform ??
        this.run(() => inputForProperty(prop, {...options, assumeWritable: true}));
      kind = platform !== null && input === platform ? null : kindOf(prop, true);
    }
    const field: FuncField = {param, input, kind, route, filter, orphaned: false,
      propagate: route === 'choices' && prop.options?.['propagateChoice'] === 'all',
      applying: false};
    if (route === 'table') {
      field.fromParam = (v) => asNamed(v)?.name ?? null;
      field.toParam = (name) => name == null ? null : tableByName(name);
    }
    else if (route === 'column' || route === 'columns') {
      field.parentName = parent!.name;
      const table = () => resolveTable(this._parentField(field)?.param.value);
      if (route === 'column') {
        field.fromParam = (v) => asNamed(v)?.name ?? null;
        field.toParam = (name) => name == null ? null : table()?.columns.byName(name) ?? null;
      }
      else {
        field.fromParam = (v) => columnNamesOf(v);
        field.toParam = (names: string[]) => {
          const t = table();
          return t == null ? [] : names.map((n) => t.columns.byName(n)).filter((c) => c != null);
        };
      }
    }
    if (kind !== null) {
      input.value.value = this._initialValue(field);
      // an empty column_list counts as empty too (isEmpty([]) is false; Dart's ColumnsInput
      // enforces a non-empty pick)
      if (prop.nullable === false)
        input.addValidator((value) => ObjectForm.isEmpty(value) ||
          (route === 'columns' && Array.isArray(value) && value.length === 0) ?
          this._requiredMessage(field) : null);
    }
    this.add(input);
    if (route === 'table' && this._formOptions.showTableSelectors === false)
      input.root.hidden = true;
    if (route === 'table' && this._autoFilled.has(param.name))
      markAuto(field);
    this._fields.push(field);
    if (route === 'choices' || route === 'multiChoices')
      this._wireSource(field);
    const descriptions = route === 'choices' ? prop.options?.['descriptions'] : null;
    // the SELECTED item's description as the input tooltip (choice_input.dart:116-117) — on the
    // editor area, so the label keeps the property description
    if (descriptions != null && typeof descriptions === 'object')
      this.effect(() => {
        const value = input.value.value;
        input.box.title = (value != null ? descriptions[value] : undefined) ??
          prop.description ?? '';
      });

    let initial = true;
    this.effect(() => {
      const value = input.value.value;
      if (initial) {
        initial = false;
        return;
      }
      if (this._refreshing || field.orphaned)
        return;
      // a badge over nothing claims a value that is gone (a pruned table, a no-match clear)
      if (value == null && field.autoBadge !== undefined)
        unmarkAuto(field);
      // a platform editor is bound by `forProperty` and writes the param itself
      if (field.kind !== null) {
        if (ObjectForm.same(value, this._read(field)))
          return;
        this._call.setParamValue(param.name, field.toParam !== undefined ?
          field.toParam(value) : FuncCallForm._paramValue(field.kind, value));
        if (field.propagate && !field.applying && !Input.isSystemWrite &&
            this._formOptions.skipDefaultInit !== true)
          this._propagate(field, value);
      }
      // under `twoWayBinding: false` no param.onChanged refresh runs to clear a failed default
      if (field.param.value != null)
        field.state?.set(1, {kind: 'idle'});
      this._formOptions.onInputChanged?.(param.name, value);
    });
  }

  private _suggestInput(name: string, options: InputOptions<any>): Input<string> {
    let tooltips: Record<string, string> = {};
    return this.run(() => new SuggestInput({
      ...options,
      minChars: 1,
      source: async (text: string, abort: AbortSignal) => {
        const call = this._call;
        if (typeof call.evalParamSuggestions !== 'function')
          throw new Error('call does not support suggestion evaluation');
        const r = await call.evalParamSuggestions(name, text);
        // only the response whose items become current keeps its tooltips — a stale landing
        // (superseded query) is dropped by the source and must not overwrite them
        if (!abort.aborted)
          tooltips = r.tooltips ?? {};
        return r.items;
      },
      // the second result column as the row tooltip (fpe:655)
      renderItem: (item: string) => {
        const el = SuggestionList.row('u2-typeahead-text', item);
        if (tooltips[item] != null)
          el.title = tooltips[item];
        return el;
      },
    }));
  }

  /** Builds (or rebuilds, on a `source` rebind) the dynamic-choices machinery of one field.
   * A call without `evalParamChoices` gets the visible error state instead — never a TypeError. */
  private _wireSource(field: FuncField): void {
    const state = this._stateOf(field);
    field.source?.dispose();
    field.source = undefined;
    if (typeof this._call.evalParamChoices !== 'function') {
      state.set(0, {kind: 'error',
        message: 'Couldn\'t load choices: call does not support choice evaluation',
        retry: () => this._wireSource(field)});
      return;
    }
    field.source = this.run(() => new ParamSource(this._call, field.param.name,
      (r) => this._applyItems(field, r), state));
  }

  private _applyItems(field: FuncField, r: ChoicesResult): void {
    field.lookup = r.lookup;
    // items always — including under `skipDefaultInit` (divergence #8): items are not defaults.
    // The setItems prune then flows into the call through the field effect (divergence #9),
    // whose batched flush runs before the microtask clears the flag.
    field.applying = true;
    (field.input as ChoiceInput | MultiChoiceInput).setItems(r.items);
    queueMicrotask(() => field.applying = false);
  }

  private _stateOf(field: FuncField): ParamState {
    let state = field.state;
    if (state === undefined) {
      const created = this.run(() => new ParamState());
      state = created;
      field.state = created;
      field.input.box.append(created.root);
      this.effect(() => {
        created.busy.value;
        this._applyEnabled(field);
      });
    }
    return state;
  }

  /** The one owner of `input.enabled`: orphaned × initial loading (a refresh never disables). */
  private _applyEnabled(field: FuncField): void {
    field.input.enabled = !field.orphaned && !(field.state?.busy.peek() ?? false);
    // the disable sweeps every button under root — the state element's Retry must stay clickable
    if (!field.input.enabled) {
      const retry = field.state?.root.querySelector<HTMLButtonElement>('.u2-param-source-retry');
      if (retry != null)
        retry.disabled = false;
    }
  }

  /** `propagateChoice: 'all'`: the picked item's lookup row lands in the sibling INPUT signals,
   * so sibling validation runs, their own effects write the call, and dependent sources
   * re-trigger — the Dart cascade (fpe:586-591). */
  private _propagate(field: FuncField, value: unknown): void {
    const row = value == null ? undefined : field.lookup?.[String(value)];
    if (row == null)
      return;
    for (const name of Object.keys(row)) {
      const sibling = this._fields.find((f) => f.param.name === name && !f.orphaned && f !== field);
      if (sibling === undefined)
        continue;
      sibling.input.value.value = sibling.kind === null ? row[name] :
        ObjectForm.coerce(sibling.kind, row[name]);
    }
  }

  /** {@link settled}'s producer, re-armed per bind: the computed-default evals plus each
   * source's initial, immediate run. */
  private _arm(): void {
    const starts: Promise<void>[] = [];
    for (const field of this._fields) {
      if (field.source !== undefined && !field.orphaned)
        starts.push(field.source.start());
    }
    this._settled = Promise.allSettled([this._initDefaults(), ...starts]).then(() => undefined);
  }

  /** Computed defaults (R6): evaluated once per bind, in parallel, each isolated onto its
   * field's state element; written through `setParamValue` so the param stream refreshes the
   * field. Never re-run on dependency changes — a computed default cannot see sibling params. */
  private _initDefaults(): Promise<void> {
    if (this._formOptions.skipDefaultInit === true)
      return Promise.resolve();
    const evals: Promise<void>[] = [];
    for (const param of this._call.inputParams.values()) {
      if (FuncCallForm._route(param.property) === 'skip')
        continue;
      if (param.property.options?.['default'] == null || param.value != null)
        continue;
      evals.push(this._evalDefault(param.name));
    }
    return Promise.allSettled(evals).then(() => undefined);
  }

  private async _evalDefault(name: string, retry = false): Promise<void> {
    const call = this._call;
    const generation = this._generation;
    const param = [...call.inputParams.values()].find((p) => p.name === name);
    const field = this._fields.find((f) => f.param.name === name && !f.orphaned);
    const state = field !== undefined ? this._stateOf(field) : undefined;
    try {
      if (typeof call.evalParamDefault !== 'function')
        throw new Error('call does not support default evaluation');
      state?.set(1, {kind: 'loading', refresh: retry});
      const value = await call.evalParamDefault(name);
      if (this._generation !== generation || this.scope.isDisposed)
        return;
      // the user may have filled the field while the eval was in flight — never overwrite
      if (param?.value == null)
        call.setParamValue(name, value);
      state?.set(1, {kind: 'idle'});
    } catch (e) {
      if (this._generation !== generation || this.scope.isDisposed)
        return;
      state?.set(1, {kind: 'error', message: 'Couldn\'t compute the default: ' +
        (e instanceof Error ? e.message : String(e)),
        retry: () => void this._evalDefault(name, true)});
    }
  }

  /** One subscription per param, owned by the form's scope and dropped on a `source` rebind. */
  private _bind(): void {
    if (this._formOptions.twoWayBinding === false)
      return;
    for (const field of this._fields) {
      if (field.orphaned)
        continue;
      const sub = field.param.onChanged.subscribe(() => this._refreshField(field));
      const off = () => sub.unsubscribe();
      this.own(off);
      this._paramSubs.push(off);
    }
  }

  private _refreshField(field: FuncField): void {
    // a param the user (or anyone) filled no longer misses its computed default
    if (field.param.value != null)
      field.state?.set(1, {kind: 'idle'});
    const value = this._read(field);
    if (ObjectForm.same(value, field.input.value.peek()))
      return;
    this._refreshing = true;
    try {
      field.input.value.value = value;
    } finally {
      this._refreshing = false;
    }
  }

  private _read(field: FuncField): any {
    if (field.fromParam !== undefined)
      return field.fromParam(field.param.value);
    const value = field.param.value;
    return field.kind === null ? value : ObjectForm.coerce(field.kind, value);
  }

  /** The param's value; a null one with a literal default shows the default — display-only, never
   * written into the FuncCall (the ApiTests `form without default initialization` contract).
   * W3 routes skip the default: a dataframe defaultValue is never serialized Dart-side
   * (func_param.dart:190). */
  private _initialValue(field: FuncField): any {
    const value = this._read(field);
    if (field.toParam !== undefined)
      return value;
    const defaultValue = field.param.property.defaultValue;
    return value == null && defaultValue != null && field.kind !== null ?
      ObjectForm.coerce(field.kind, defaultValue) : value;
  }

  /** The required-error text, cause-aware for the column routes: a null parent table suppresses
   * it — one root cause, one message, on the table field — and a table with no passing columns
   * names the real reason. {@link TableBinding} revalidates on every retarget, so the verdict
   * follows the table even under an unchanged (null) value. */
  private _requiredMessage(field: FuncField): string | null {
    if (field.route !== 'column' && field.route !== 'columns')
      return 'Value can\'t be empty';
    const table = resolveTable(this._parentField(field)?.param.value);
    if (table === null)
      return null;
    return table.columns.toList().some((c) => field.filter!(c)) ? 'Value can\'t be empty' :
      noMatchMessage(field.param.property, table.name);
  }

  /** By name alone, never by route: an override or a registered editor on the table param
   * demotes its route to 'field', and the dependents keep working through the param stream. */
  private _parentField(field: FuncField): FuncField | undefined {
    return this._fields.find((f) => f.param.name === field.parentName);
  }

  /** One {@link TableBinding} per bound table param — every 'table'-routed field, plus an
   * override-demoted one some dependent names — rebuilt whole on a `source` rebind. The binding
   * rides the PARAM stream, so a custom table editor still drives its dependents. A dependent
   * whose parent field is orphaned (or missing) orphans with it: its `toParam` would otherwise
   * resolve names against the OLD call's table. */
  private _bindTables(): void {
    for (const binding of this._tableBindings)
      binding.dispose();
    this._tableBindings = [];
    const dependents = new Map<FuncField, ColumnDependent[]>();
    for (const f of this._fields) {
      if (f.route !== 'column' && f.route !== 'columns')
        continue;
      const parent = this._parentField(f);
      const orphaned = f.orphaned || parent === undefined || parent.orphaned;
      if (orphaned !== f.orphaned) {
        f.orphaned = orphaned;
        this._applyEnabled(f);
      }
      if (orphaned)
        continue;
      const deps = dependents.get(parent!) ?? [];
      deps.push({kind: f.route, field: f, filter: f.filter!});
      dependents.set(parent!, deps);
    }
    for (const field of this._fields) {
      if (field.orphaned || (field.route !== 'table' && !dependents.has(field)))
        continue;
      const binding = this.run(() => new TableBinding(this._call, field.param,
        (name) => tableByName(name), dependents.get(field) ?? [],
        (dep, oldTableName) => this._clearedNotice(dep.field, oldTableName)));
      this._tableBindings.push(binding);
      binding.start();
    }
  }

  /** A table switch wiped a non-empty column selection (kept as the ruled `[]` write) — say so
   * on the field, until the user's next interaction with it. */
  private _clearedNotice(field: FuncField, oldTableName: string): void {
    const state = this._stateOf(field);
    state.notice(`Selection cleared — columns belonged to '${oldTableName}'`);
    if (field.noticeClear !== undefined)
      return;
    const root = field.input.root;
    const clear = () => {
      state.notice(null);
      root.removeEventListener('pointerdown', clear);
      root.removeEventListener('keydown', clear);
      field.noticeClear = undefined;
    };
    root.addEventListener('pointerdown', clear);
    root.addEventListener('keydown', clear);
    field.noticeClear = clear;
  }

  /** `camelCaseToWords` over the `friendlyName` when set, else over the name — the Dart rule
   * applies it to both (ib:679). */
  private static _caption(prop: IProperty, name: string): string {
    const friendly = prop.friendlyName;
    return camelCaseToWords(friendly != null && friendly !== '' ? friendly : name);
  }

  /** js-api's `toDart` marshals dayjs but passes a plain `Date` through to Dart, so a datetime
   * write goes out as dayjs — via the global the platform bundle installs (js-api dg.ts:56);
   * headless there is no platform and the Date itself is the doubles' native value. */
  private static _paramValue(kind: Kind, value: unknown): unknown {
    const dayjs = (globalThis as any).dayjs;
    return kind === 'datetime' && value instanceof Date && typeof dayjs === 'function' ?
      dayjs(value) : value;
  }
}

/** Generates a {@link FuncCallForm} over `call`'s input params — the u2 replacement surface for
 * `DG.InputForm.forFuncCall` (the W6 flip happens behind that name). */
export function funcForm(call: DG.FuncCall, options?: FuncFormOptions): FuncCallForm {
  return new FuncCallForm(call as unknown as FuncCallLike, options);
}
