/* The editor of property metadata itself — what a property IS, not what it holds. Its field
   vocabulary is the platform's own `Property.propertyOptions` map (js-api `property.ts:322-348`),
   mirrored here (the map only exists on a loaded platform, and this control builds headless);
   `applicableTo` decides which fields a given property type shows, exactly as the map declares it.

   The value is a plain `IProperty` record, never a live `DG.Property`: `name` and `type` of a
   bound property are not safely mutable, so the caller rebuilds instead
   (`DG.Property.fromOptions`) — see plan.md §5. */
import type {IProperty} from 'datagrok-api/dg';
import {Control} from '../core/component.js';
import {Scope} from '../core/scope.js';
import {signal, Signal} from '../core/signals.js';
import {ObjectForm, propertyForm} from './object-form.js';
import type {PropertyLike} from '../core/property-like.js';

/** The types the editor offers, and the ones the convergence matrix covers (plan.md §"full input
 * matrix"). */
const TYPES = ['string', 'int', 'double', 'num', 'bigint', 'qnum', 'bool', 'datetime', 'list',
  'map', 'file'];

/** What `applicableTo: TYPE.NUMERICAL` covers here: the types whose editor actually reads
 * min/max/step/showSlider/showPlusMinus. `bigint` and `qnum` have their own editors, which take
 * no bounds. */
const NUMERICAL = ['int', 'double', 'num', 'float'];

interface FieldMeta extends PropertyLike {
  /** Property types this field applies to; shown for every type when absent. */
  applicableTo?: string;
}

/** Always shown: what the property is, and the hints that pick its editor. */
const IDENTITY: FieldMeta[] = [
  {name: 'name', type: 'string', nullable: false},
  {name: 'type', type: 'string', nullable: false, choices: TYPES,
    description: 'Property data type, such as "int" or "string".'},
  {name: 'friendlyName', type: 'string', caption: 'Friendly name',
    description: 'Custom field friendly name shown in [PropertyGrid]'},
  {name: 'description', type: 'string', editor: 'TextArea', description: 'Property description'},
  {name: 'nullable', type: 'bool',
    description: 'Whether an empty value is allowed. This is used by validators.'},
  {name: 'semType', type: 'string', caption: 'Semantic type', description: 'Semantic type'},
  {name: 'inputType', type: 'string', caption: 'Input type', description: 'Property input type'},
  {name: 'editor', type: 'string', description: 'Custom editor (such as slider or text area)'},
  {name: 'units', type: 'string', description: 'Units of measurement. See also: [postfix]'},
  {name: 'format', type: 'string', description: 'Value format, such as "0.00"'},
  {name: 'category', type: 'string', description: 'Corresponding category on the context panel'},
];

/** Shown for the property types they apply to. */
const TYPED: FieldMeta[] = [
  {name: 'min', type: 'double', applicableTo: 'numerical',
    description: 'Minimum value. Applicable to numerical properties only'},
  {name: 'max', type: 'double', applicableTo: 'numerical',
    description: 'Maximum value. Applicable to numerical properties only'},
  {name: 'step', type: 'double', applicableTo: 'numerical',
    description: 'Step to be used in a slider. Only applies to numerical properties.'},
  {name: 'showSlider', type: 'bool', applicableTo: 'numerical', caption: 'Show slider',
    description: 'Whether a slider appears next to the number input. Applies to numerical columns only.'},
  {name: 'showPlusMinus', type: 'bool', applicableTo: 'numerical', caption: 'Show plus/minus',
    description: 'Whether a plus/minus clicker appears next to the number input. Applies to numerical columns only.'},
  // 'string_list' platform-side; declared as a list so it edits as the CSV `ListInput`
  {name: 'choices', type: 'list', applicableTo: 'string',
    description: 'List of choices. Applicable to string properties only'},
];

/** The identity fields whose edit changes which typed fields apply. */
const REBUILDS = ['type', 'editor', 'inputType'];

export interface PropertyEditorOptions {
  /** The WHOLE record on every edit — the caller rebuilds whatever hangs off the property
   * (`DG.Property.fromOptions`) rather than patching it. */
  onChanged?: (options: IProperty) => void;
  condensed?: boolean;
  /** Extra validation per field name, for the rules only the caller knows — a name already taken by
   * another property, say. Re-applied every time the fields are rebuilt, so a re-target keeps it. */
  validators?: Record<string, (value: any) => string | null>;
}

/** Edits the metadata of one property: an identity section that never changes, and a
 * type-dependent one rebuilt whenever the type (or an editor hint) does. Re-target it with
 * {@link setTarget} instead of building a second one — the context panel shows one editor at a
 * time, and a re-target leaks nothing. */
export class PropertyEditor extends Control {
  /** The edited record, replaced by a fresh copy on every edit. */
  readonly options: Signal<IProperty>;

  private readonly _typed = document.createElement('div');
  private readonly _onChanged: ((options: IProperty) => void) | undefined;
  private readonly _condensed: boolean | undefined;
  private readonly _validators: Record<string, (value: any) => string | null>;
  private _record: IProperty = {};
  private _scope: Scope | undefined;
  private _typedScope: Scope | undefined;
  private _rebuilding = false;

  constructor(target: IProperty, options: PropertyEditorOptions = {}) {
    super();
    this._onChanged = options.onChanged;
    this._condensed = options.condensed;
    this._validators = options.validators ?? {};
    this._record = {...target};
    this.options = signal<IProperty>({...this._record});
    this.root.classList.add('u2-property-editor');
    this.root.dataset.u2 = 'property-editor';
    this._typed.className = 'u2-property-editor-typed';
    this.own(() => {
      this._typedScope?.dispose();
      this._scope?.dispose();
    });
    this._render();
  }

  /** Edits another record from now on. The caller's object is never touched: the editor works on
   * a copy and reports it whole.
   *
   * The record and the {@link options} signal change at once; the fields are rebuilt on the next
   * microtask. A caller that re-targets from its own `onChanged` handler runs inside the effect the
   * current generation owns, and rebuilding there would dispose that effect while it is running.
   * Several calls in one turn therefore cost one rebuild, of the last target. */
  setTarget(target: IProperty): void {
    this._record = {...target};
    this.options.value = {...this._record};
    if (this._rebuilding)
      return;
    this._rebuilding = true;
    queueMicrotask(() => {
      this._rebuilding = false;
      if (!this.scope.isDisposed)
        this._render();
    });
  }

  private _render(): void {
    this._typedScope?.dispose();
    this._scope?.dispose();
    const scope = new Scope();
    this._scope = scope;
    this.root.textContent = '';
    Scope.runWith(scope, () => {
      this.root.append(this._form(IDENTITY).root, this._typed);
      this._renderTyped();
    });
  }

  private _renderTyped(): void {
    this._typedScope?.dispose();
    const scope = new Scope();
    this._typedScope = scope;
    this._typed.textContent = '';
    const type = this._record.type ?? '';
    const fields = TYPED.filter((field) => PropertyEditor._applies(field, type));
    Scope.runWith(scope, () => this._typed.append(this._form(fields).root));
  }

  private _form(fields: FieldMeta[]): ObjectForm {
    // the name is what the caller keys the property by, and renaming is not something a reader can
    // do halfway: reported per keystroke, 'dose' → 'name' renames through 'n', 'na' and 'nam' too
    const form = propertyForm(fields.map((field) => PropertyEditor._propertyFor(field)), this._record,
      {condensed: this._condensed, onChanged: (name) => this._edited(name),
        overrides: {name: {commitOn: 'change'}}});
    for (const field of fields) {
      const validator = this._validators[field.name];
      if (validator)
        form.input(field.name)!.addValidator(validator);
    }
    return form;
  }

  private _edited(name: string): void {
    const record = {...this._record};
    this.options.value = record;
    this._onChanged?.(record);
    if (REBUILDS.includes(name))
      this._renderTyped();
  }

  private static _applies(field: FieldMeta, type: string): boolean {
    if (field.applicableTo === undefined)
      return true;
    return field.applicableTo === 'numerical' ? NUMERICAL.includes(type) : field.applicableTo === type;
  }

  private static _propertyFor(field: FieldMeta): PropertyLike {
    const name = field.name;
    return {caption: PropertyEditor._caption(name), ...field,
      get: (x: any) => x[name], set: (x: any, value: any) => x[name] = value};
  }

  /** Sentence case off the field name — `showSlider` reads as 'Show slider'. A field whose label
   * needs more than that (`semType`) carries its own caption. */
  private static _caption(name: string): string {
    const words = name.replace(/([A-Z])/g, ' $1').toLowerCase();
    return words.charAt(0).toUpperCase() + words.slice(1);
  }
}

/** Editor over a property-options record; every edit reports the whole record. */
export function propertyEditor(target: IProperty, options?: PropertyEditorOptions): PropertyEditor {
  return new PropertyEditor(target, options);
}
