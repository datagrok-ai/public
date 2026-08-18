import {signal, Signal, ReadonlySignal} from './signals.js';
import {Component} from './component.js';
import {label, span} from './elements.js';

export interface InputOptions<T> {
  label?: string;
  /** Stable key for forms and serialization; falls back to {@link label}. */
  name?: string;
  value?: T;
  /** Two-way external binding; without it the input owns its signal. */
  bind?: Signal<T>;
  /** Compact editor-only variant: no label, one inline-control-height row. */
  inline?: boolean;
  nullable?: boolean;
  onChanged?: (value: T) => void;
  tooltipText?: string;
}

export type Validator<T> = (value: T) => string | null;

/** Label + editor + validity, shared by every u2 input. Subclasses provide the editor element
 * and the concrete `data-u2` id; everything below is the same for all of them. */
export abstract class Input<T, O extends InputOptions<T> = InputOptions<T>> extends Component {
  readonly value: Signal<T>;
  readonly validity: ReadonlySignal<string | null>;
  readonly name: string | undefined;

  protected readonly options: O;
  private readonly _validators: Validator<T>[] = [];
  private readonly _validity = signal<string | null>(null);
  private readonly _error = span('', 'u2-input-error');
  private readonly _editor: HTMLElement;
  private _enabled = true;

  constructor(options: O, defaultValue: T) {
    super();
    this.options = options;
    this.value = options.bind ?? signal(options.value ?? defaultValue);
    this.validity = this._validity;
    this.name = options.name ?? options.label;

    this.root.classList.add('u2-input-root');
    if (options.inline)
      this.root.classList.add('u2-input-inline');
    if (options.tooltipText)
      this.root.title = options.tooltipText;
    if (!options.inline && options.label !== undefined)
      this.root.append(label(options.label, 'u2-input-label'));

    this._editor = this.run(() => this.createEditor());
    this._editor.classList.add('u2-input-editor');
    this.root.append(this._editor, this._error);

    this.effect(() => this._validate(this.value.value));
    this.effect(() => this._showValidity(this._validity.value));
    const onChanged = options.onChanged;
    if (onChanged) {
      let initial = true;
      this.effect(() => {
        const value = this.value.value;
        if (initial)
          initial = false;
        else
          onChanged(value);
      });
    }
  }

  get nullable(): boolean {
    return this.options.nullable ?? true;
  }

  get enabled(): boolean {
    return this._enabled;
  }

  set enabled(x: boolean) {
    this._enabled = x;
    this.root.classList.toggle('u2-input-disabled', !x);
    const controls = this.root.querySelectorAll('input, textarea, button, select');
    for (let i = 0; i < controls.length; i++)
      (controls[i] as HTMLInputElement).disabled = !x;
  }

  addValidator(v: Validator<T>): void {
    this._validators.push(v);
    this._validate(this.value.peek());
  }

  protected abstract createEditor(): HTMLElement;

  private _validate(value: T): void {
    let message: string | null = null;
    for (const validator of this._validators) {
      message = validator(value);
      if (message)
        break;
    }
    this._validity.value = message;
  }

  private _showValidity(message: string | null): void {
    this._editor.classList.toggle('u2-invalid', message !== null);
    this._error.textContent = message ?? '';
    if (message === null)
      this._editor.removeAttribute('title');
    else
      this._editor.title = message;
  }
}
