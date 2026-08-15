/* u2 inputs as platform value editors (PLAN.md D8, "InputBase-conformant"). */
import * as DG from 'datagrok-api/dg';
import {Input} from '../core/input-base.js';
import {Scope} from '../core/scope.js';
import {PlatformInput} from './from-dart-input.js';

const DATA_TYPES: {[id: string]: string} = {
  'text-input': DG.TYPE.STRING,
  'text-area': DG.TYPE.STRING,
  'choice-input': DG.TYPE.STRING,
  'number-input': DG.TYPE.FLOAT,
  'bool-input': DG.TYPE.BOOL,
  'multi-choice-input': DG.TYPE.OBJECT,
};

export interface DartInputOptions {
  /** Overrides the type inferred from the u2 input; a `DG.TYPE` constant. */
  dataType?: string;
}

/** A u2 `Input` seen by the platform as a `DG.JsInputBase`: dialogs, forms, property editors and
 * `role: 'valueEditor'` registrations take it wherever they take a `DG.InputBase`.
 *
 * `getInput()` returns the u2 **editor**, not the input's root: the Dart proxy (`JsInputProxy`)
 * builds its own `ui-input-root` with a caption label and appends `getInput()` to it, so handing
 * over the root would nest two input roots and show two captions. The consequences, all
 * deliberate: the caption comes from the platform (seeded here from the u2 label, then owned by
 * the dialog/property that binds the input), validation messages are rendered by the platform
 * (u2's own error line stays in the unmounted u2 root), and `enabled`/`readOnly` must be set on
 * this object rather than on the u2 input, whose setter looks for controls under its own root. */
export class DartInput<T> extends DG.JsInputBase<T> {
  private readonly _input: Input<T>;
  private readonly _editor: HTMLElement;
  private readonly _inputType: string;
  private readonly _dataType: string;
  private readonly _scope = new Scope();
  private _programmatic = false;

  constructor(input: Input<T>, options: DartInputOptions = {}) {
    super();
    const id = input.root.dataset.u2 ?? 'input';
    this._input = input;
    this._editor = input.root.querySelector<HTMLElement>('.u2-input-editor') ?? input.root;
    this._inputType = `u2-${id}`;
    this._dataType = options.dataType ?? DATA_TYPES[id] ?? this._inferDataType();
    this.caption = input.root.querySelector('.u2-input-label')?.textContent ?? '';
    this.root.dataset.u2 = id;
    this.addValidator(() => this._input.validity.peek());

    let initial = true;
    this._scope.effect(() => {
      this._input.value.value; // subscribes the effect to user and programmatic edits alike
      if (initial) {
        initial = false;
        return;
      }
      if (!this._programmatic)
        this.fireInput();
      this.fireChanged();
    });
    input.own(() => this._scope.dispose());
  }

  get inputType(): string { return this._inputType; }

  get dataType(): string { return this._dataType; }

  /** The wrapped u2 input — its `value` signal stays the source of truth. */
  get u2(): Input<T> { return this._input; }

  getInput(): HTMLElement { return this._editor; }

  getValue(): T { return this._input.value.peek(); }

  /** Programmatic set: fires `onChanged` only, and nothing at all when the value is unchanged. */
  setValue(value: T): void {
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

  /** Severs the bridge; the wrapped u2 input is never owned by it and stays alive. Disposing the
   * u2 input does the same, so this is only needed when the input outlives the bridge. */
  detach(): void {
    this._scope.dispose();
  }

  private _inferDataType(): string {
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
    if (this._dataType === DG.TYPE.STRING)
      return text as unknown as T;
    const trimmed = text.trim();
    if (this._dataType === DG.TYPE.BOOL)
      return (trimmed.toLowerCase() === 'true') as unknown as T;
    if (trimmed === '')
      return null as unknown as T;
    if (this._dataType === DG.TYPE.FLOAT || this._dataType === DG.TYPE.INT) {
      const parsed = Number(trimmed);
      return (isNaN(parsed) ? null : parsed) as unknown as T;
    }
    if (this._dataType === DG.TYPE.OBJECT)
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
