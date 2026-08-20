/* u2 inputs as platform value editors (PLAN.md D8, "InputBase-conformant"). */
import * as DG from 'datagrok-api/dg';
import {Input} from '../core/input-base.js';
import {Scope} from '../core/scope.js';
import {BigIntInput} from '../components/inputs/bigint-input.js';
import {PlatformInput} from './from-dart-input.js';

const DATA_TYPES: {[id: string]: string} = {
  'text-input': DG.TYPE.STRING,
  'text-area': DG.TYPE.STRING,
  'choice-input': DG.TYPE.STRING,
  'number-input': DG.TYPE.FLOAT,
  'slider-input': DG.TYPE.FLOAT,
  'qnum-input': DG.TYPE.QNUM,
  'bigint-input': DG.TYPE.BIG_INT,
  'bool-input': DG.TYPE.BOOL,
  'color-input': DG.TYPE.STRING,
  'radio-input': DG.TYPE.STRING,
  'multi-choice-input': DG.TYPE.OBJECT,
  'list-input': DG.TYPE.OBJECT,
  'map-input': DG.TYPE.OBJECT,
  'file-input': DG.TYPE.FILE,
  // the platform's own type for a list of files (ddt `utils.dart:162`); DG.TYPE has no constant
  'files-input': 'files',
  'rsa-input': DG.TYPE.STRING,
  // the pickers are name-valued (pickers.ts), so a column comes over as its name and the
  // collections as plain arrays and maps of names — not as `column`/`column_list`
  'column-picker': DG.TYPE.STRING,
  'columns-input': DG.TYPE.OBJECT,
  'columns-map-input': DG.TYPE.OBJECT,
  'aggregated-columns-input': DG.TYPE.OBJECT,
};

export interface DartInputOptions {
  /** Overrides the type inferred from the u2 input; a `DG.TYPE` constant. */
  dataType?: string;
}

/** Builds the editor for the property the platform bound; null wherever it bound none. */
export type PropertyInputBuilder<T> = (prop: DG.Property | null) => Input<T>;

/** A u2 `Input` seen by the platform as a `DG.JsInputBase`: dialogs, forms, property editors and
 * `role: 'valueEditor'` registrations take it wherever they take a `DG.InputBase`.
 *
 * `getInput()` returns the u2 **editor box** — the editor together with its options rail, so units
 * and side icons cross the bridge with it — and not the input's root: the Dart proxy
 * (`JsInputProxy`) builds its own `ui-input-root` with a caption label and appends `getInput()` to
 * it, so handing over the root would nest two input roots and show two captions. The consequences,
 * all
 * deliberate: the caption comes from the platform (seeded here from the u2 label, then owned by
 * the dialog/property that binds the input), validation messages are rendered by the platform
 * (u2's own error line stays in the unmounted u2 root), and `enabled`/`readOnly` must be set on
 * this object rather than on the u2 input, whose setter looks for controls under its own root.
 *
 * Built from a builder ({@link dartInputFor}) the u2 input appears on the first `getInput()`
 * instead of in the constructor; everything below tolerates that gap. */
export class DartInput<T> extends DG.JsInputBase<T> {
  private readonly _options: DartInputOptions;
  private readonly _build: PropertyInputBuilder<T> | null;
  private readonly _scope = new Scope();
  private _input: Input<T> | null = null;
  private _editor: HTMLElement | null = null;
  private _pending: T[] | null = null;
  private _programmatic = false;
  private _validated: unknown = null;
  private _queued = false;
  private _userEdit = false;

  constructor(input: Input<T> | PropertyInputBuilder<T>, options: DartInputOptions = {}) {
    super();
    this._options = options;
    this._build = typeof input === 'function' ? input : null;
    this._syncValidator();
    if (this._build === null) {
      const eager = input as Input<T>;
      this.caption = eager.root.querySelector('.u2-input-label')?.textContent ?? '';
      this._attach(eager);
    }
  }

  get inputType(): string { return `u2-${this._id}`; }

  get dataType(): string {
    return this._options.dataType ?? DATA_TYPES[this._id] ?? this._inferDataType();
  }

  /** The wrapped u2 input — its `value` signal stays the source of truth. Null until a deferred
   * editor is built. */
  get u2(): Input<T> | null { return this._input; }

  getInput(): HTMLElement {
    this._syncValidator();
    // the platform builds from a Timer callback, where no scope is ambient — without this the
    // deferred input's effects would belong to nobody and outlive the editor
    if (this._editor === null) {
      const prop = this._boundProperty();
      this._attach(Scope.runWith(this._scope, () => this._build!(prop)));
      DartInput._dropPostfix(this._editor!, prop);
    }
    return this._editor!;
  }

  getValue(): T {
    if (this._input !== null)
      return this._input.value.peek();
    return (this._pending === null ? null : this._pending[0]) as T;
  }

  /** Programmatic set: fires `onChanged` only, and nothing at all when the value is unchanged.
   * Before a deferred editor exists the value waits for it, as the Dart proxy's `_tempValue` does
   * (`input_base_js_proxy.dart:45-51`). */
  setValue(value: T): void {
    if (this._input === null) {
      this._pending = [value];
      return;
    }
    this._programmatic = true;
    try {
      this._input.value.value = value;
    } finally {
      this._programmatic = false;
    }
  }

  getStringValue(): string {
    const value = this.getValue();
    if (value === null || value === undefined)
      return '';
    return typeof value === 'object' ? JSON.stringify(value) : String(value);
  }

  setStringValue(value: string): void {
    this.setValue(this._parse(value));
  }

  /** Severs the bridge. An input passed in ({@link asDartInput}) stays alive — it was never owned;
   * one built here ({@link dartInputFor}) is disposed with the bridge, since nothing else holds it.
   * Disposing the u2 input severs the bridge too, so this is only needed the other way round. */
  detach(): void {
    this._scope.dispose();
  }

  private get _id(): string {
    return this._input?.root.dataset.u2 ?? 'input';
  }

  private _attach(input: Input<T>): void {
    this._input = input;
    this._editor = input.box;
    this.root.dataset.u2 = this._id;
    if (this._pending !== null) {
      const [value] = this._pending;
      this._pending = null;
      this.setValue(value);
    }

    let initial = true;
    this._scope.effect(() => {
      input.value.value; // subscribes the effect to user and programmatic edits alike
      if (initial) {
        initial = false;
        return;
      }
      this._notify(!this._programmatic && !Input.isSystemWrite);
    });
    input.own(() => this._scope.dispose());
  }

  /** The platform hears about the edit from a microtask, never from inside the value effect. A
   * subscriber that writes any signal back — the other half of a form syncing itself, a picker
   * refreshed from an event — would land in the middle of the running effect round, and the signal
   * runtime aborts the round as a cycle, killing the effect for good. Deferring keeps the contract
   * intact (`fireInput` before `fireChanged`; a programmatic write fires `fireChanged` alone) and
   * reports one pair per edit, however many signal writes that edit took. */
  private _notify(userEdit: boolean): void {
    this._userEdit = this._userEdit || userEdit;
    if (this._queued)
      return;
    this._queued = true;
    queueMicrotask(() => {
      this._queued = false;
      const user = this._userEdit;
      this._userEdit = false;
      if (this._scope.isDisposed)
        return;
      if (user)
        this.fireInput();
      this.fireChanged();
    });
  }

  /** Validators live on the dart handle, and the platform can re-home the JS input on a second
   * one: the proxy `JsInputBase` builds is dropped when a `valueEditor` func resolves, and
   * `jsi.dart` is repointed at the proxy the form actually binds
   * (`input_base_js_proxy.dart:18-25,27-38`). The validator added in the constructor would then
   * gate nothing — no invalid state, no disabled OK — so every `getInput()` puts it on whatever
   * handle the platform holds now. */
  private _syncValidator(): void {
    const dart = this.dart;
    if (dart === this._validated)
      return;
    this._validated = dart;
    this.addValidator(() => this._input?.validity.peek() ?? null);
  }

  /** The property the platform bound to this editor before asking for its element — the proxy
   * binds first and runs `getInput()` from a `Timer` (`input_base_js_proxy.dart:18-25`), so the
   * builder sees everything `forProperty` applied. Null outside that path. */
  private _boundProperty(): DG.Property | null {
    try {
      return this.property ?? null;
    } catch {
      return null;
    }
  }

  /** The Dart host renders the bound property's units itself, into the root it wraps this editor
   * in (`input_base.dart:219-222`), so the rail's own postfix would say the same thing twice. Only
   * this path — a standalone or `propertyForm` input has no host and keeps it. */
  private static _dropPostfix(editor: HTMLElement, prop: DG.Property | null): void {
    if (!prop || !prop.units)
      return;
    const postfix = editor.querySelector('.u2-input-postfix');
    if (postfix === null)
      return;
    const rail = postfix.parentElement!;
    postfix.remove();
    if (rail.children.length === 0)
      rail.remove();
  }

  private _inferDataType(): string {
    if (this._input === null)
      return DG.TYPE.OBJECT;
    const value = this._input.value.peek();
    if (typeof value === 'string')
      return DG.TYPE.STRING;
    if (typeof value === 'number')
      return DG.TYPE.FLOAT;
    if (typeof value === 'boolean')
      return DG.TYPE.BOOL;
    return DG.TYPE.OBJECT;
  }

  private _parse(text: string): T {
    const type = this.dataType;
    if (type === DG.TYPE.STRING)
      return text as unknown as T;
    const trimmed = text.trim();
    if (type === DG.TYPE.BOOL)
      return (trimmed.toLowerCase() === 'true') as unknown as T;
    if (trimmed === '')
      return null as unknown as T;
    if (type === DG.TYPE.BIG_INT)
      return (BigIntInput.parse(trimmed) ?? null) as unknown as T;
    if (type === DG.TYPE.FLOAT || type === DG.TYPE.INT || type === DG.TYPE.QNUM) {
      const parsed = Number(trimmed);
      return (isNaN(parsed) ? null : parsed) as unknown as T;
    }
    if (type === DG.TYPE.OBJECT)
      return JSON.parse(trimmed) as T;
    return text as unknown as T;
  }
}

/** Wraps a u2 input as a platform value editor. See {@link DartInput} for what the platform owns
 * (caption, validation display, enabled state) and what stays with the u2 input (the value). */
export function asDartInput<T>(input: Input<T>, options: DartInputOptions = {}): DartInput<T> {
  // wrapping a wrapper back would build a second root and a second caption; the platform input the
  // caller already has is the better object in every respect (INTEROP.md: wrap exactly once)
  if (input instanceof PlatformInput)
    throw new Error('This input already wraps a DG.InputBase; use `input.dgInput` directly.');
  return new DartInput<T>(input, options);
}

/** A value editor that builds its u2 input on the first `getInput()`, when the Dart proxy is
 * already bound to the property — so `build` reads min/max/step/format/units/showSlider/
 * showPlusMinus off it and the editor carries the same metadata the platform's own would
 * (`inputForProperty` is that mapping). This is what `role: 'valueEditor'` registrations return;
 * {@link asDartInput} stays the wrapper for an input that already exists. */
export function dartInputFor<T>(build: PropertyInputBuilder<T>,
  options: DartInputOptions = {}): DartInput<T> {
  return new DartInput<T>(build, options);
}
