/* The func-call editor: a form generated from a FuncCall's input params, writing through
   `setParamValue` and following `param.onChanged` back —
   the u2 counterpart of `DG.InputForm.forFuncCall`. Structural over {@link FuncCallLike}, so
   headless doubles and DG.FuncCall both fit; the platform itself is only reached by
   `editors: 'auto'`, through the global the bundler binds `datagrok-api/dg` to.
   Wave 2 adds the dynamic routes — async choices/suggestions through the `evalParam*` members
   and computed defaults written into the call (R6) — with `settled` as the readiness member. */
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
import type * as DG from 'datagrok-api/dg';

export interface FuncFormOptions {
  /** Two-way is the u2 default. `false` is honored for DG.InputForm contract
   * compatibility only: FuncCall edits then do not refresh the form. */
  twoWayBinding?: boolean;
  /** Suppresses computed-default evaluation and `propagateChoice` writes (R6). Dynamic choice
   * ITEMS still load and apply — items are not defaults (divergence #8). Literal defaults stay
   * display-only either way. */
  skipDefaultInit?: boolean;
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
interface FuncField {
  param: FuncCallParamLike;
  input: Input<any>;
  kind: Kind | null;
  /** The dynamic route the field was built for; 'field' for every static editor, including a
   * dynamic-routed param an override or a registered editor claimed (no dynamic wiring then). */
  route: 'field' | 'choices' | 'suggestions' | 'multiChoices';
  orphaned: boolean;
  propagate: boolean;
  /** True while a choices result is being applied to the input — the batched flush of the
   * setItems write runs after {@link Input.system}'s counter has dropped, so the propagate
   * guard needs its own flag. */
  applying: boolean;
  lookup?: Record<string, Record<string, any>> | null;
  state?: ParamState;
  source?: ParamSource;
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
  private _paramSubs: (() => void)[] = [];
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
    this._unsupported.length = 0;
    for (const param of call.inputParams.values()) {
      params.set(param.name, param);
      if (FuncCallForm._route(param.property) !== 'skip' &&
          !this._fields.some((f) => f.param.name === param.name))
        this._unsupported.push(param.name);
    }
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
    const categories: {name: string, params: {param: FuncCallParamLike, route: FuncField['route']}[]}[] = [];
    for (const param of this._call.inputParams.values()) {
      const route = FuncCallForm._route(param.property);
      if (route === 'skip')
        continue;
      if (route === 'unsupported') {
        this._unsupported.push(param.name);
        continue;
      }
      const name = param.property.category ?? 'Misc';
      let category = categories.find((c) => c.name === name);
      if (category === undefined) {
        category = {name, params: []};
        categories.push(category);
      }
      category.params.push({param, route});
    }
    const headers = categories.length > 1 ||
      (categories.length === 1 && categories[0].name !== 'Misc');
    for (const category of categories) {
      if (headers)
        this.addElement(div([category.name], 'u2-form-category'));
      for (const {param, route} of category.params)
        this._addField(param, route);
    }
  }

  /** First match wins: SERVICE_PARAM skipped; `editorParam` mirrors unsupported; a present
   * `choices`/`suggestions` OPTION routes dynamic — the Dart discriminator (fpe:569-572, 645),
   * with choices winning over suggestions. Never `prop.choices`: annotation parsing leaves
   * garbage there for a dynamic source (scripting.dart:316), while the evaluator answers even
   * static list literals. A typed `prop.choices` with NO option keeps the W1 static path. */
  private static _route(prop: FuncCallParamLike['property']):
      'skip' | 'unsupported' | 'field' | 'choices' | 'suggestions' | 'multiChoices' {
    const options = prop.options;
    const editor = options?.['editor'] ?? prop.editor ?? null;
    if (editor === 'none')
      return 'skip';
    if (options?.['editorParam'] != null)
      return 'unsupported';
    const type = prop.propertyType ?? prop.type ?? '';
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

  private _addField(param: FuncCallParamLike, routed: FuncField['route']): void {
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
    else {
      const platform = (custom ?? registered) != null ? null :
        ObjectForm.platformInput(this, prop, this._call, this._formOptions.editors === 'auto');
      input = custom ?? registered ?? platform ??
        this.run(() => inputForProperty(prop, {...options, assumeWritable: true}));
      kind = platform !== null && input === platform ? null : kindOf(prop, true);
    }
    const field: FuncField = {param, input, kind, route, orphaned: false,
      propagate: route === 'choices' && prop.options?.['propagateChoice'] === 'all',
      applying: false};
    if (kind !== null) {
      input.value.value = this._initialValue(field);
      if (prop.nullable === false)
        input.addValidator((value) => ObjectForm.isEmpty(value) ? 'Value can\'t be empty' : null);
    }
    this.add(input);
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
      // a platform editor is bound by `forProperty` and writes the param itself
      if (field.kind !== null) {
        if (ObjectForm.same(value, this._read(field)))
          return;
        this._call.setParamValue(param.name, FuncCallForm._paramValue(field.kind, value));
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
    const value = field.param.value;
    return field.kind === null ? value : ObjectForm.coerce(field.kind, value);
  }

  /** The param's value; a null one with a literal default shows the default — display-only, never
   * written into the FuncCall (the ApiTests `form without default initialization` contract). */
  private _initialValue(field: FuncField): any {
    const value = this._read(field);
    const defaultValue = field.param.property.defaultValue;
    return value == null && defaultValue != null && field.kind !== null ?
      ObjectForm.coerce(field.kind, defaultValue) : value;
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
