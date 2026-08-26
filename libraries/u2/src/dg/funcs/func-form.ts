/* The func-call editor, wave 1: a form generated from a FuncCall's input params, writing
   through `setParamValue` and following `param.onChanged` back —
   the u2 counterpart of `DG.InputForm.forFuncCall`. Structural over {@link FuncCallLike}, so
   headless doubles and DG.FuncCall both fit; the platform itself is only reached by
   `editors: 'auto'`, through the global the bundler binds `datagrok-api/dg` to. */
import {batch} from '../../core/signals.js';
import type {IProperty} from '../../core/property-like.js';
import {Input} from '../../core/input-base.js';
import type {InputOptions} from '../../core/input-base.js';
import {div} from '../../core/elements.js';
import {Form} from '../../components/forms/form.js';
import {ObjectForm, inputForProperty, kindOf} from '../forms/object-form.js';
import type {FieldOverride, Kind} from '../forms/object-form.js';
import {Editors} from '../forms/editors.js';
import type * as DG from 'datagrok-api/dg';

export interface FuncFormOptions {
  /** Two-way is the u2 default. `false` is honored for DG.InputForm contract
   * compatibility only: FuncCall edits then do not refresh the form. */
  twoWayBinding?: boolean;
  /** Reserved for computed defaults (W2). Accepted and stored now; W1 has no computed
   * defaults, and literal defaults are display-only, so it is a no-op in W1. */
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

export interface FuncCallLike {
  inputParams: {values(): Iterable<FuncCallParamLike>};
  setParamValue(name: string, value: any): void;
}

/** The per-field seam W2's dynamic sources plug into. Null kind marks a platform editor, which
 * reads and writes the param in its own native type. Orphaned: the current `source` carries no
 * param of this name — the field is disabled and its write path is off. */
interface FuncField {
  param: FuncCallParamLike;
  input: Input<any>;
  kind: Kind | null;
  orphaned: boolean;
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

  constructor(call: FuncCallLike, options: FuncFormOptions = {}) {
    super({condensed: options.condensed});
    this._formOptions = options;
    this._call = call;
    this.root.dataset.u2 = 'func-form';
    this._build();
    this._bind();
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
            field.input.enabled = !orphaned;
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
    const categories: {name: string, params: FuncCallParamLike[]}[] = [];
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
      category.params.push(param);
    }
    const headers = categories.length > 1 ||
      (categories.length === 1 && categories[0].name !== 'Misc');
    for (const category of categories) {
      if (headers)
        this.addElement(div([category.name], 'u2-form-category'));
      for (const param of category.params)
        this._addField(param);
    }
  }

  /** First match wins: SERVICE_PARAM skipped; nested/layout editors, `editorParam` mirrors
   * and non-scalar types unsupported; the four ib hints route like the property form's. */
  private static _route(prop: FuncCallParamLike['property']): 'skip' | 'unsupported' | 'field' {
    const options = prop.options;
    const editor = options?.['editor'] ?? prop.editor ?? null;
    if (editor === 'none')
      return 'skip';
    if (options?.['editorParam'] != null)
      return 'unsupported';
    if (!SCALARS.has(prop.propertyType ?? prop.type ?? ''))
      return 'unsupported';
    return editor === null || HINTS.has(String(editor).toLowerCase()) ? 'field' : 'unsupported';
  }

  private _addField(param: FuncCallParamLike): void {
    const prop = param.property;
    const {input: custom, ...rest} = this._formOptions.overrides?.[param.name] ?? {};
    const kind = kindOf(prop, true);
    const options: InputOptions<any> = {
      name: param.name,
      label: FuncCallForm._caption(prop, param.name),
      tooltipText: prop.description ?? undefined,
      nullable: prop.nullable,
      postfix: prop.units || prop.options?.['units'] || undefined,
      ...rest,
    };
    const registered = custom ? null : this.run(() => Editors.resolve(prop, options));
    const platform = (custom ?? registered) != null ? null :
      ObjectForm.platformInput(this, prop, this._call, this._formOptions.editors === 'auto');
    const input = custom ?? registered ?? platform ??
      this.run(() => inputForProperty(prop, {...options, assumeWritable: true}));
    const native = platform !== null && input === platform;
    const field: FuncField = {param, input, kind: native ? null : kind, orphaned: false};
    if (!native) {
      input.value.value = this._initialValue(field);
      if (prop.nullable === false)
        input.addValidator((value) => ObjectForm.isEmpty(value) ? 'Value can\'t be empty' : null);
    }
    this.add(input);
    this._fields.push(field);

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
      }
      this._formOptions.onInputChanged?.(param.name, value);
    });
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
