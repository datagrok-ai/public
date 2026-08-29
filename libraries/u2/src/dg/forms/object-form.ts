/* Schema-driven forms — GOAL.md's priority #1: a whole form derived from Property metadata, with
   per-field overrides. Everything the generator reads off a property is declared by IProperty,
   which DG.Property satisfies structurally; the platform itself is only reached by `editors:
   'auto'`, through the global the bundler binds `datagrok-api/dg` to. */
import {batch} from '../../core/signals.js';
import type {IProperty} from '../../core/property-like.js';
import {Input, InputOptions} from '../../core/input-base.js';
import {text} from '../../core/text.js';
import {Form} from '../../components/forms/form.js';
import type {FormLayout} from '../../components/forms/form.js';
import {TextInput} from '../../components/inputs/text-input.js';
import {NumberInput, NumberInputOptions} from '../../components/inputs/number-input.js';
import {BigIntInput} from '../../components/inputs/bigint-input.js';
import {QNumInput} from '../../components/inputs/qnum-input.js';
import {BoolInput} from '../../components/inputs/bool-input.js';
import {ChoiceInput, MultiChoiceInput} from '../../components/inputs/choice-input.js';
import {DateTimeInput} from '../../components/inputs/date-input.js';
import {TextArea} from '../../components/inputs/text-input.js';
import {RadioInput} from '../../components/inputs/radio-input.js';
import {ColorInput} from '../../components/inputs/color-input.js';
import {FontInput} from '../../components/inputs/font-input.js';
import {IconInput} from '../../components/inputs/icon-input.js';
import {ImageInput} from '../../components/inputs/image-input.js';
import {TagsInput} from '../../components/inputs/tags-input.js';
import {SliderInput} from '../../components/inputs/slider-input.js';
import {ListInput} from '../../components/inputs/list-input.js';
import {MapInput} from '../../components/inputs/map-input.js';
import {fromDartInput} from '../inputs/from-dart-input.js';
import {Editors} from './editors.js';

type DgApi = typeof import('datagrok-api/dg');

export interface FieldOverride extends Partial<InputOptions<any>> {
  /** Replaces the generated input entirely — a semType-specific editor, say. The form wires it to
   * the property and lays it out, but never owns it: dispose it with whatever created it. */
  input?: Input<any>;
}

export interface ObjectFormOptions {
  /** Property names to show, in this order; everything, in metadata order, by default. */
  include?: string[];
  exclude?: string[];
  /** Per-field overrides, keyed by property name. */
  overrides?: Record<string, FieldOverride>;
  /** `auto` asks the platform for the editor of every real `DG.Property` (molecule, column, date,
   * color…), generating one only where it has none; `generated` (the default) never asks. */
  editors?: 'auto' | 'generated';
  layout?: FormLayout;
  /** The DG.InputForm-contract flag; in the platform condensed IS the tall view, so it maps to
   * `layout: 'tall'` where {@link layout} is not given. */
  condensed?: boolean;
  onChanged?: (name: string, value: unknown) => void;
}

/** What {@link objectForm} enumerates properties off. JS-implemented `DG.Widget` subclasses are
 * ones; `DG.Entity` is NOT, and Dart-backed objects need their dart handle as the target — see
 * src/dg/README.md. */
export interface PropertySource {
  getProperties(): IProperty[];
}

/** The value type a field is read and written in — what the editor is built from is decided
 * separately, by {@link routeFor}. */
export type Kind = 'string' | 'int' | 'float' | 'bigint' | 'qnum' | 'bool' | 'choice' | 'datetime' |
  'list' | 'map' | 'file' | 'readonly';

interface Field {
  prop: IProperty;
  /** Null for a platform editor, which reads and writes the property in its own native type. */
  kind: Kind | null;
  input: Input<any>;
}

export type InputFactory = (prop: IProperty, options: InputOptions<any>) => Input<any>;

/** Editors the router routes to but cannot import: this module has to load without the platform
 * (its own tests, `editors` and `from-dart-input` run headless), while `FileInput` reaches for
 * `grok.dapi` at import time. `src/dg/index.ts` registers them; a kind nobody registered falls
 * back to a read-only text field. */
export class PlatformInputs {
  private static readonly _factories: Record<string, InputFactory> = {};

  /** Returns the unregister function, as {@link Editors.register} does. */
  static register(kind: string, factory: InputFactory): () => void {
    PlatformInputs._factories[kind] = factory;
    return () => {
      if (PlatformInputs._factories[kind] === factory)
        delete PlatformInputs._factories[kind];
    };
  }

  static create(kind: string, prop: IProperty, options: InputOptions<any>): Input<any> | null {
    const factory = PlatformInputs._factories[kind];
    return factory ? factory(prop, options) : null;
  }
}

/** `InputType` → editor, the u2 half of the platform's `inputFactories` map
 * (`input_base.dart:580-607`); an input type u2 has no editor for falls through to the rest of
 * the routing, where Dart would reach for a JS-registered input instead. */
const BY_INPUT_TYPE: Record<string, InputFactory> = {
  Text: (prop, options) => new TextInput(options),
  Search: (prop, options) => new TextInput({...options, search: true}),
  TextArea: (prop, options) => new TextArea(options),
  Int: (prop, options) => new NumberInput({...options, ...numberOptions(prop, 'int', options)}),
  Float: (prop, options) => new NumberInput({...options, ...numberOptions(prop, 'float', options)}),
  BigInt: (prop, options) => new BigIntInput(options),
  QNum: (prop, options) => new QNumInput({...options, format: formatter(prop)}),
  Slider: (prop, options) => sliderFor(prop, options),
  Bool: (prop, options) => new BoolInput(options),
  Switch: (prop, options) => new BoolInput({...options, switch: true}),
  Date: (prop, options) => new DateTimeInput(options),
  Color: (prop, options) => new ColorInput(options),
  Font: (prop, options) => new FontInput(options),
  Image: (prop, options) => new ImageInput(options),
  Icon: (prop, options) => new IconInput(options),
  Radio: (prop, options) => new RadioInput({...options, items: prop.choices ?? []}),
  Choice: (prop, options) => new ChoiceInput({...options, items: prop.choices ?? []}),
  MultiChoice: (prop, options) => new MultiChoiceInput({...options, items: prop.choices ?? []}),
  Tags: (prop, options) => new TagsInput<string>({...options, items: prop.choices ?? undefined,
    allowNew: true}),
  List: (prop, options) => new ListInput(options),
  Map: (prop, options) => new MapInput(options),
  File: (prop, options) => fileInputFor(prop, options),
};

/** `editor` → editor (`input_base.dart:702-728`), matched case-insensitively: Dart compares against
 * lower-case 'textarea', while `Property.propertyOptions` writes the `InputType` spelling
 * ('TextArea') into its own `description` option (`property.ts:326`). */
const BY_EDITOR: Record<string, InputFactory> = {
  textarea: (prop, options) => new TextArea(options),
  password: (prop, options) => new TextInput({...options, password: true}),
  switch: (prop, options) => new BoolInput({...options, switch: true}),
  slider: (prop, options) => sliderFor(prop, options),
};

/** The property types each `editor` hint is honored for, mirroring where `_forProperty` reaches its
 * branch: switch under `pt == Types.BOOL` (`input_base.dart:702`), textarea and password under
 * `pt == Types.STRING` (`:725-728`), and the slider (`:705`) for everything the two bool branches
 * above it did not already take. A hint the type does not accept is ignored, and the type's own
 * editor is built — `{editor: 'textarea', type: 'int'}` is an int box on both sides. */
const EDITOR_TYPES: Record<string, (type: string | null | undefined) => boolean> = {
  textarea: (type) => type === 'string',
  password: (type) => type === 'string',
  switch: (type) => type === 'bool',
  slider: (type) => type !== 'bool',
};

export function kindOf(prop: IProperty, assumeWritable = false): Kind {
  if (!prop.set && !assumeWritable)
    return 'readonly';
  const type = prop.propertyType ?? prop.type;
  if (type === 'list')
    return 'list';
  if (type === 'map')
    return 'map';
  if (type === 'file' || type === 'blob')
    return 'file';
  if (prop.choices != null && prop.choices.length > 0)
    return 'choice';
  switch (type) {
    case 'string':
      return 'string';
    case 'datetime':
      return 'datetime';
    case 'int':
      return 'int';
    case 'bigint':
      return 'bigint';
    case 'qnum':
      return 'qnum';
    case 'double':
    case 'float':
    case 'num':
      return 'float';
    case 'bool':
      return 'bool';
    default:
      return 'readonly';
  }
}

/** The editor the property's own hints ask for, before its type gets a say. */
function routeFor(prop: IProperty, options: InputOptions<any>): Input<any> | null {
  const byInputType = prop.inputType ? BY_INPUT_TYPE[prop.inputType] : undefined;
  if (byInputType)
    return byInputType(prop, options);
  const editor = prop.editor?.toLowerCase();
  const byEditor = editor ? BY_EDITOR[editor] : undefined;
  return byEditor && EDITOR_TYPES[editor!](prop.propertyType ?? prop.type) ?
    byEditor(prop, options) : null;
}

function sliderFor(prop: IProperty, options: InputOptions<any>): Input<any> {
  return new SliderInput({...options, min: finite(prop.min) ?? 0, max: finite(prop.max) ?? 100,
    step: finite(prop.step)});
}

function fileInputFor(prop: IProperty, options: InputOptions<any>): Input<any> {
  return PlatformInputs.create('file', prop, options) ?? readonlyText(options);
}

function readonlyText(options: InputOptions<any>): Input<any> {
  const input = new TextInput(options);
  input.enabled = false;
  return input;
}

function finite(value: number | null | undefined): number | undefined {
  return typeof value === 'number' && isFinite(value) ? value : undefined;
}

/** `DG.format` where the platform is loaded (the same global `editors: 'auto'` reaches for);
 * without it the input keeps its own precision rules. */
function formatter(prop: IProperty): ((value: number) => string) | undefined {
  const format = prop.format;
  const dgFormat = (globalThis as any).DG?.format as ((x: number, f: string) => string) | undefined;
  if (format == null || format === '' || typeof dgFormat !== 'function')
    return undefined;
  return (value) => dgFormat(value, format);
}

/** What `NumberInput.bindProperty` applies Dart-side (`number_input.dart:119-129`): bounds, step
 * and format from the property, a clicker on bounded ints, a slider on floats or on explicit
 * `showSlider`, units as the postfix. */
function numberOptions(prop: IProperty, kind: 'int' | 'float',
  options: InputOptions<any>): NumberInputOptions {
  const min = finite(prop.min);
  const max = finite(prop.max);
  return {
    mode: kind === 'int' ? 'int' : 'float',
    min, max,
    step: finite(prop.step),
    slider: prop.showSlider === true || kind === 'float',
    clicker: prop.showPlusMinus ?? (kind === 'int' && min !== undefined && max !== undefined),
    format: formatter(prop),
    postfix: options.postfix ?? prop.units ?? undefined,
  };
}

export interface PropertyInputOptions extends InputOptions<any> {
  /** Builds the editor the property's type asks for even when the property has no setter. For a
   * value editor (`dartInputFor`) the value belongs to the platform proxy, not to the property,
   * and a `FuncParam` never carries a `set` — without this every func-param dialog would get a
   * disabled text box. */
  assumeWritable?: boolean;
}

/** The editor {@link propertyForm} generates for a property, with the platform's own metadata
 * mapping applied. Exported for value-editor builders (`dartInputFor`), which learn the property
 * only after the platform binds it; a null property — nothing bound — gets a plain text editor. */
export function inputForProperty(prop: IProperty | null,
  options: PropertyInputOptions = {}): Input<any> {
  const {assumeWritable, ...rest} = options;
  if (prop == null)
    return new TextInput(rest);
  const kind = kindOf(prop, assumeWritable);
  const merged: InputOptions<any> = {nullable: prop.nullable,
    tooltipText: prop.description ?? undefined, ...rest};
  if (kind === 'readonly')
    return readonlyText(merged);
  const routed = routeFor(prop, merged);
  if (routed != null)
    return routed;
  switch (kind) {
    case 'bool':
      return new BoolInput(merged);
    case 'choice':
      return new ChoiceInput({...merged, items: prop.choices ?? []});
    case 'int':
    case 'float':
      return new NumberInput({...merged, ...numberOptions(prop, kind, merged)});
    case 'bigint':
      return new BigIntInput(merged);
    case 'qnum':
      return new QNumInput({...merged, format: formatter(prop)});
    case 'datetime':
      return new DateTimeInput(merged);
    case 'list':
      return prop.choices != null && prop.choices.length > 0 ?
        new MultiChoiceInput({...merged, items: prop.choices}) : new ListInput(merged);
    case 'map':
      return new MapInput(merged);
    case 'file':
      return fileInputFor(prop, merged);
    default:
      return new TextInput(merged);
  }
}

/** A {@link Form} generated from property metadata and kept in sync with {@link target} both ways:
 * an edit runs `prop.set`, {@link refresh} re-reads `prop.get`. It mutates the object and nothing
 * else — persisting it (`grok.dapi.*.save`, a domain insert) stays with the caller. */
export class ObjectForm extends Form {
  readonly target: object;

  private readonly _fields: Field[] = [];
  private readonly _onChanged: ((name: string, value: unknown) => void) | undefined;
  private readonly _auto: boolean;
  private _refreshing = false;

  constructor(props: IProperty[], target: object, options: ObjectFormOptions = {}) {
    super({layout: options.layout ?? (options.condensed ? 'tall' : undefined)});
    this.target = target;
    this._onChanged = options.onChanged;
    this._auto = options.editors === 'auto';
    this.root.dataset.u2 = 'object-form';
    for (const prop of ObjectForm._select(props, options))
      this._addField(prop, options.overrides?.[prop.name!] ?? {});
  }

  /** The properties the form renders, in layout order. */
  get properties(): ReadonlyArray<IProperty> {
    return this._fields.map((f) => f.prop);
  }

  /** The input editing a property — generated, override-supplied, or the platform's own editor
   * wrapped by `fromDartInput`. Keyed by property name regardless of the input's caption. */
  input(name: string): Input<any> | undefined {
    return this._fields.find((f) => f.prop.name === name)?.input;
  }

  /** Re-reads every field off {@link target} — after a save, an external edit, a server refresh.
   * Echo-suppressed: nothing is written back and `onChanged` does not fire. */
  refresh(): void {
    // a plain call, so this flag still holds when the batch flushes the effects it triggered —
    // unlike a write made inside an effect, which is flushed only after that effect returns
    this._refreshing = true;
    try {
      batch(() => {
        for (const field of this._fields)
          field.input.value.value = this._read(field.prop, field.kind);
      });
    } finally {
      this._refreshing = false;
    }
  }

  private _addField(prop: IProperty, override: FieldOverride): void {
    const kind = kindOf(prop);
    const {input: custom, ...rest} = override;
    const options: InputOptions<any> = {
      name: prop.name,
      label: prop.caption ?? prop.friendlyName ?? prop.name,
      tooltipText: prop.description ?? undefined,
      nullable: prop.nullable,
      ...rest,
    };
    const registered = custom ? null : this.run(() => Editors.resolve(prop, options));
    const platform = (custom ?? registered) != null ? null :
      ObjectForm.platformInput(this, prop, this.target, this._auto);
    const input = custom ?? registered ?? platform ??
      this.run(() => inputForProperty(prop, options));
    const native = input === platform;
    if (!native) {
      input.value.value = this._read(prop, kind);
      if (prop.nullable === false)
        input.addValidator((value) => ObjectForm.isEmpty(value) ? 'Value can\'t be empty' : null);
    }
    this.add(input);
    this._fields.push({prop, kind: native ? null : kind, input});

    // a platform editor is bound to the target by `forProperty` and writes the property itself
    if (native) {
      let initial = true;
      this.effect(() => {
        const value = input.value.value;
        if (initial)
          initial = false;
        else if (!this._refreshing)
          this._onChanged?.(prop.name!, value);
      });
      return;
    }

    const set = prop.set;
    if (kind === 'readonly' || !set)
      return;
    // echo suppression by comparison, not by a flag: a write made inside an effect is flushed
    // after that effect returns, so a transient flag would already be down (as in PropertyGrid)
    this.effect(() => {
      const value = input.value.value;
      if (ObjectForm.same(value, this._read(prop, kind)))
        return;
      set(this.target, value);
      this._onChanged?.(prop.name!, value);
    });
  }

  /** The platform's own editor for a real `DG.Property` under `auto`, or null wherever it has
   * none — a property u2 declared itself, an older core, an editor that refuses the property.
   * The wrapper is owned by `form`'s scope, as its generated inputs are. */
  static platformInput(form: Form, prop: IProperty, target: object,
    auto: boolean): Input<any> | null {
    const dg = (globalThis as any).DG as DgApi | undefined;
    if (!auto || dg == null || (prop as {dart?: unknown}).dart == null)
      return null;
    try {
      const input = dg.InputBase.forProperty(prop as any, target);
      return input == null ? null : form.run(() => fromDartInput(input, prop.name));
    } catch {
      return null;
    }
  }

  private _read(prop: IProperty, kind: Kind | null): any {
    const value = prop.get ? prop.get(this.target) : undefined;
    return kind === null ? value : ObjectForm.coerce(kind, value);
  }

  private static _select(props: IProperty[], options: ObjectFormOptions): IProperty[] {
    let selected = props;
    if (options.include) {
      const byName = new Map(props.map((p) => [p.name, p]));
      selected = options.include.map((name) => byName.get(name)).filter((p): p is IProperty => p !== undefined);
    }
    const exclude = options.exclude;
    return exclude ? selected.filter((p) => !exclude.includes(p.name!)) : selected;
  }

  static coerce(kind: Kind, value: unknown): any {
    switch (kind) {
      case 'bool':
        return value === true;
      case 'int':
      case 'float':
      case 'qnum': {
        if (typeof value === 'number')
          return value;
        const parsed = value === null || value === undefined || value === '' ? NaN : Number(value);
        return isFinite(parsed) ? parsed : null;
      }
      // js-api marshals a Dart `BigInt` to a JS one and back (wrappers_impl.ts:88,133); a property
      // that hands over the digits as text is read here just as well
      case 'bigint':
        return typeof value === 'bigint' ? value :
          (value === null || value === undefined ? null : BigIntInput.parse(String(value)) ?? null);
      case 'choice':
        return value === null || value === undefined ? null : String(value);
      case 'datetime': {
        if (value === null || value === undefined || value === '')
          return null;
        if (value instanceof Date)
          return value;
        // dayjs values (Dart-bound datetime props) expose toDate()
        const dayjsLike = value as {toDate?: () => Date};
        const date = typeof dayjsLike.toDate === 'function' ? dayjsLike.toDate() : new Date(value as string | number);
        return isNaN(date.getTime()) ? null : date;
      }
      case 'list':
        return Array.isArray(value) ? value :
          (value === null || value === undefined || value === '' ? [] : ListInput.parse(String(value)));
      case 'map':
        return typeof value === 'object' && value !== null && !Array.isArray(value) ? value : {};
      // the platform's own value (a FileInfo), handed to the platform's own editor untouched
      case 'file':
        return value ?? null;
      default:
        return text(value);
    }
  }

  /** Identity for echo suppression. `Object.is` decides the scalars, but {@link coerce} builds a
   * fresh `Date`, list or record on every read of an empty or foreign-shaped value, so identity
   * alone would report an edit — and fire `onChanged` — the moment the field is constructed. */
  static same(a: unknown, b: unknown): boolean {
    if (a instanceof Date && b instanceof Date)
      return a.getTime() === b.getTime();
    if (Array.isArray(a) && Array.isArray(b))
      return a.length === b.length && a.every((item, i) => Object.is(item, b[i]));
    if (ObjectForm._isRecord(a) && ObjectForm._isRecord(b)) {
      const keys = Object.keys(a);
      return keys.length === Object.keys(b).length && keys.every((k) => Object.is(a[k], b[k]));
    }
    return Object.is(a, b);
  }

  private static _isRecord(x: unknown): x is Record<string, unknown> {
    return typeof x === 'object' && x !== null && !Array.isArray(x) && !(x instanceof Date);
  }

  static isEmpty(value: unknown): boolean {
    return value === null || value === undefined || value === '';
  }
}

/** Generates a form over `props`, editing `source` in place. The primary schema-driven surface:
 * `props` come from wherever the caller has them — `DG.Property.fromOptions`, a viewer's
 * `getProperties()`, `grok.dapi.domains` row properties, `grok.dapi.entities.getProperties()`. */
export function propertyForm(props: IProperty[], source: object, options?: ObjectFormOptions): ObjectForm {
  return new ObjectForm(props, source, options);
}

/** {@link propertyForm} over an object that enumerates its own properties AND is itself the
 * get/set target — a JS-implemented `DG.Widget` subclass. A `DG.Viewer` is not one: its properties
 * are defined over its LOOK (`grok_Viewer_Get_Look(viewer.dart)`), so both the viewer and its dart
 * handle throw as the target — use `viewerSettings(viewer)` (`dg/viewers/viewers.ts`), which
 * hands `propertyForm` the look. `DG.Entity` cannot enumerate at all. */
export function objectForm(source: PropertySource, options?: ObjectFormOptions): ObjectForm {
  return propertyForm(source.getProperties(), source, options);
}
