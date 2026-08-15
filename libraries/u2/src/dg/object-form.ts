/* Schema-driven forms — GOAL.md's priority #1: a whole form derived from Property metadata, with
   per-field overrides. Everything the generator reads off a property is declared by PropertyLike,
   which DG.Property satisfies structurally; the platform itself is only reached by `editors:
   'auto'`, through the global the bundler binds `datagrok-api/dg` to. */
import {batch} from '../core/signals.js';
import {Input, InputOptions} from '../core/input-base.js';
import {Form} from '../components/form.js';
import {TextInput} from '../components/text-input.js';
import {NumberInput} from '../components/number-input.js';
import {BoolInput} from '../components/bool-input.js';
import {ChoiceInput} from '../components/choice-input.js';
import {fromDartInput} from './from-dart-input.js';
import {Editors} from './editors.js';

type DgApi = typeof import('datagrok-api/dg');

/** The property metadata a form is generated from. `DG.Property` IS a `PropertyLike`; so is any
 * descriptor that names the field and can read and write it (`type` is read when `propertyType`
 * is absent, as on a `DG.IProperty` literal). A property with no `set` renders read-only. */
export interface PropertyLike {
  name: string;
  caption?: string | null;
  friendlyName?: string | null;
  propertyType?: string | null;
  type?: string | null;
  semType?: string | null;
  description?: string | null;
  choices?: string[] | null;
  nullable?: boolean;
  min?: number | null;
  max?: number | null;
  get?: (source: any) => any;
  set?: (source: any, value: any) => void;
}

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
  condensed?: boolean;
  onChanged?: (name: string, value: unknown) => void;
}

/** What {@link objectForm} enumerates properties off. JS-implemented `DG.Widget` subclasses are
 * ones; `DG.Entity` is NOT, and Dart-backed objects need their dart handle as the target — see
 * src/dg/README.md. */
export interface PropertySource {
  getProperties(): PropertyLike[];
}

type Kind = 'string' | 'int' | 'float' | 'bool' | 'choice' | 'readonly';

interface Field {
  prop: PropertyLike;
  /** Null for a platform editor, which reads and writes the property in its own native type. */
  kind: Kind | null;
  input: Input<any>;
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

  constructor(props: PropertyLike[], target: object, options: ObjectFormOptions = {}) {
    super({condensed: options.condensed});
    this.target = target;
    this._onChanged = options.onChanged;
    this._auto = options.editors === 'auto';
    this.root.dataset.u2 = 'object-form';
    for (const prop of ObjectForm._select(props, options))
      this._addField(prop, options.overrides?.[prop.name] ?? {});
  }

  /** The properties the form renders, in layout order. */
  get properties(): ReadonlyArray<PropertyLike> {
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

  private _addField(prop: PropertyLike, override: FieldOverride): void {
    const kind = ObjectForm._kindOf(prop);
    const {input: custom, ...rest} = override;
    const options: InputOptions<any> = {
      name: prop.name,
      label: prop.caption ?? prop.friendlyName ?? prop.name,
      tooltipText: prop.description ?? undefined,
      nullable: prop.nullable,
      ...rest,
    };
    const registered = custom ? null : this.run(() => Editors.resolve(prop, options));
    const platform = (custom ?? registered) != null ? null : this._platformInput(prop);
    const input = custom ?? registered ?? platform ??
      this.run(() => ObjectForm._createInput(prop, kind, options));
    const native = input === platform;
    if (!native) {
      input.value.value = this._read(prop, kind);
      if (prop.nullable === false)
        input.addValidator((value) => ObjectForm._isEmpty(value) ? 'Value can\'t be empty' : null);
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
          this._onChanged?.(prop.name, value);
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
      if (Object.is(value, this._read(prop, kind)))
        return;
      set(this.target, value);
      this._onChanged?.(prop.name, value);
    });
  }

  /** The platform's own editor for a real `DG.Property`, or null wherever it has none — a property
   * u2 declared itself, an older core, an editor that refuses the property. */
  private _platformInput(prop: PropertyLike): Input<any> | null {
    const dg = (globalThis as any).DG as DgApi | undefined;
    if (!this._auto || dg == null || (prop as {dart?: unknown}).dart == null)
      return null;
    try {
      const input = dg.InputBase.forProperty(prop as any, this.target);
      return input == null ? null : this.run(() => fromDartInput(input, prop.name));
    } catch {
      return null;
    }
  }

  private _read(prop: PropertyLike, kind: Kind | null): any {
    const value = prop.get ? prop.get(this.target) : undefined;
    return kind === null ? value : ObjectForm._coerce(kind, value);
  }

  private static _select(props: PropertyLike[], options: ObjectFormOptions): PropertyLike[] {
    let selected = props;
    if (options.include) {
      const byName = new Map(props.map((p) => [p.name, p]));
      selected = options.include.map((name) => byName.get(name)).filter((p): p is PropertyLike => p !== undefined);
    }
    const exclude = options.exclude;
    return exclude ? selected.filter((p) => !exclude.includes(p.name)) : selected;
  }

  private static _kindOf(prop: PropertyLike): Kind {
    if (!prop.set)
      return 'readonly';
    if (prop.choices != null && prop.choices.length > 0)
      return 'choice';
    switch (prop.propertyType ?? prop.type) {
      case 'string':
      case 'datetime': // the date picker is the one v1 control still deferred (PLAN.md P4)
        return 'string';
      case 'int':
      case 'bigint':
        return 'int';
      case 'double':
      case 'float':
      case 'num':
      case 'qnum':
        return 'float';
      case 'bool':
        return 'bool';
      default:
        return 'readonly';
    }
  }

  private static _createInput(prop: PropertyLike, kind: Kind, options: InputOptions<any>): Input<any> {
    switch (kind) {
      case 'bool':
        return new BoolInput(options);
      case 'choice':
        return new ChoiceInput({...options, items: prop.choices ?? []});
      case 'int':
      case 'float':
        return new NumberInput({...options, mode: kind === 'int' ? 'int' : 'float',
          min: ObjectForm._number(prop.min), max: ObjectForm._number(prop.max)});
      case 'readonly': {
        const input = new TextInput(options);
        input.enabled = false;
        return input;
      }
      default:
        return new TextInput(options);
    }
  }

  private static _coerce(kind: Kind, value: unknown): any {
    switch (kind) {
      case 'bool':
        return value === true;
      case 'int':
      case 'float': {
        if (typeof value === 'number')
          return value;
        const parsed = value === null || value === undefined || value === '' ? NaN : Number(value);
        return isFinite(parsed) ? parsed : null;
      }
      case 'choice':
        return value === null || value === undefined ? null : String(value);
      default:
        return value === null || value === undefined ? '' : String(value);
    }
  }

  private static _number(value: number | null | undefined): number | undefined {
    return typeof value === 'number' && isFinite(value) ? value : undefined;
  }

  private static _isEmpty(value: unknown): boolean {
    return value === null || value === undefined || value === '';
  }
}

/** Generates a form over `props`, editing `source` in place. The primary schema-driven surface:
 * `props` come from wherever the caller has them — `DG.Property.fromOptions`, a viewer's
 * `getProperties()`, `grok.dapi.domains` row properties, `grok.dapi.entities.getProperties()`. */
export function propertyForm(props: PropertyLike[], source: object, options?: ObjectFormOptions): ObjectForm {
  return new ObjectForm(props, source, options);
}

/** {@link propertyForm} over an object that enumerates its own properties AND is itself the
 * get/set target — a JS-implemented `DG.Widget` subclass. A Dart-backed one (`DG.Viewer`,
 * `DG.ViewBase`) reads and writes through its dart handle:
 * `propertyForm(viewer.getProperties(), viewer.dart)`. `DG.Entity` cannot enumerate at all. */
export function objectForm(source: PropertySource, options?: ObjectFormOptions): ObjectForm {
  return propertyForm(source.getProperties(), source, options);
}
