import {signal, Signal, ReadonlySignal} from './signals.js';
import {bindTypeOf} from './widget-like.js';
import type {BindProp, BindSource} from './widget-like.js';
import {Control} from './component.js';
import {div, label, span} from './elements.js';

/** An input chrome option: a literal applied once, or a signal followed live (UB-10). */
export type LiveOption<T> = T | ReadonlySignal<T>;

/** The literal a chrome option carries — undefined where it is a signal or absent. */
export function labelText(x?: LiveOption<string>): string | undefined {
  return typeof x === 'string' ? x : undefined;
}

export interface InputOptions<T> {
  label?: LiveOption<string>;
  /** Stable key for forms and serialization; falls back to {@link label} where it is a string. */
  name?: string;
  value?: T;
  /** Two-way external binding; without it the input owns its signal. */
  bind?: Signal<T>;
  /** Compact editor-only variant: no label, one inline-control-height row. */
  inline?: boolean;
  nullable?: boolean;
  /** `false` grays the input out and blocks edits; live through the {@link Input.enabled} accessor. */
  enabled?: LiveOption<boolean>;
  /** Units or any short suffix shown after the editor (Dart `InputBase.addPostfix`). */
  postfix?: LiveOption<string>;
  onChanged?: (value: T) => void;
  tooltipText?: LiveOption<string>;
  /** When an edit reaches {@link Input.value}: on every keystroke (`input`, the default) or on the
   * commit — blur or Enter — as Dart's `fireChanged` does (`change`). Only the free-text editors
   * have a gap between the two; every other editor commits as it edits. Worth asking for where a
   * consumer acts on the value irreversibly, as a rename does. */
  commitOn?: 'input' | 'change';
}

export type Validator<T> = (value: T) => string | null;

/** Label + editor + options rail + validity, shared by every u2 input. Subclasses provide the
 * editor element and the concrete `data-u2` id; everything below is the same for all of them.
 *
 * The rail is the platform's `InputBase.options` (`input_base.dart:439-442`): units and every
 * side icon sit right of the editor on one 28px row whose underline continues the editor's, with
 * the postfix pinned leftmost (`ui.css:814-822` gives it `order: -1`). Editor and rail live in one
 * `u2-input-box` so the pair travels as a single element through the platform bridge, while the
 * `.u2-input-editor` selector keeps pointing at the editor alone. Hover affordances (the number
 * chrome, the text clear icon) hang off `.u2-input-root:hover` — hovering anywhere on the input,
 * label included, as `d4.css:130-134` does — with `.u2-input-box:hover` covering the bridged case,
 * where only the box is mounted and the label belongs to the platform. */
export abstract class Input<T, O extends InputOptions<T> = InputOptions<T>> extends Control {
  readonly value: Signal<T>;
  readonly validity: ReadonlySignal<string | null>;

  protected readonly options: O;
  private static _systemWrites = 0;
  private readonly _validators: Validator<T>[] = [];
  private readonly _validity = signal<string | null>(null);
  private readonly _error = span('', 'u2-input-error');
  private readonly _editor: HTMLElement;
  private readonly _inputBox: HTMLElement;
  private readonly _optionsRail: HTMLElement;
  private _labelEl: HTMLElement | undefined;
  private _postfixEl: HTMLElement | undefined;
  private _liveBinds: Map<string, Signal<unknown>> | undefined;
  private _enabled = true;

  constructor(options: O, defaultValue: T) {
    super();
    this.options = options;
    this.value = options.bind ?? signal(options.value ?? defaultValue);
    this.validity = this._validity;
    this.name = options.name ?? labelText(options.label);

    this.root.classList.add('u2-input-root');
    if (options.inline)
      this.root.classList.add('u2-input-inline');

    this._inputBox = div([], 'u2-input-box');
    this._optionsRail = div([], 'u2-input-options');
    this._editor = this.run(() => this.createEditor());
    this._editor.classList.add('u2-input-editor');
    // prepended: a subclass that filled the rail from createEditor already attached it
    this._inputBox.prepend(this._editor);
    this.root.append(this._inputBox, this._error);
    this.liveOption('label', options.label, (x) => this.label = x);
    this.liveOption('tooltipText', options.tooltipText, (x) => this.tooltipText = x);
    this.liveOption('postfix', options.postfix, (x) => {
      if (!this._postfixEl) {
        if (!x)
          return; // an empty postfix never materializes the span or its rail
        this.addOptions(this._postfixEl = span('', 'u2-input-postfix'));
      }
      this._postfixEl.textContent = x;
    });
    // `enabled: true` runs no sweep — internals a subclass disabled in createEditor must survive
    this.liveOption('enabled', options.enabled, (x) => {
      if (x === false)
        this.enabled = false;
      else if (!this._enabled)
        this.enabled = true;
    });

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

  /** Applies a value write the input makes on its own behalf — pruning a selection whose item
   * vanished, reconciling a list the platform rebuilt — so consumers that tell user edits from the
   * rest see it as one of the rest (the outbound bridge fires `onChanged` alone for it, as Dart
   * does in `multi_choice_input.dart:44-56`). */
  static system<T>(write: () => T): T {
    Input._systemWrites++;
    try {
      return write();
    } finally {
      Input._systemWrites--;
    }
  }

  /** True while a {@link system} write is running. */
  static get isSystemWrite(): boolean {
    return Input._systemWrites > 0;
  }

  get nullable(): boolean {
    return this.options.nullable ?? true;
  }

  get label(): string {
    return this._labelEl?.textContent ?? labelText(this.options.label) ?? '';
  }

  set label(x: string) {
    if (this._labelEl)
      this._labelEl.textContent = x;
    else if (!this.options.inline) {
      this._labelEl = label(x, 'u2-input-label');
      this.root.prepend(this._labelEl);
    }
  }

  get tooltipText(): string {
    return this.root.title;
  }

  set tooltipText(x: string) {
    this.root.title = x;
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

  /** Re-runs the validators against the current value — for a validator whose verdict depends
   * on state beyond the value itself, after that state moved under an unchanged value. */
  revalidate(): void {
    this._validate(this.value.peek());
  }

  /** The editor together with its options rail — what the platform bridge hands to a Dart form,
   * and what a host that wants the control without the label should mount. */
  get box(): HTMLElement {
    return this._inputBox;
  }

  /** Applies a chrome option now and — when it arrived as a signal — on every change. One-way,
   * owned by the input's scope; the signal answers {@link bindStep} so a deferred link wires it.
   * A literal installs NO effect: the plain setter stays the sole owner afterwards. */
  protected liveOption<T2>(name: string, value: LiveOption<T2> | undefined,
    apply: (x: T2) => void): void {
    if (value === undefined)
      return;
    if (!(value instanceof Signal)) {
      apply(value as T2);
      return;
    }
    (this._liveBinds ??= new Map()).set(name, value as unknown as Signal<unknown>);
    this.effect(() => {
      const x = (value as ReadonlySignal<T2>).value;
      if (x !== undefined)
        apply(x); // a forward-ref proxy starts undefined; skip until linked
    });
  }

  bindStep(name: string): Signal<unknown> | BindSource | null {
    return this._liveBinds?.get(name) ?? super.bindStep(name);
  }

  bindProps(): BindProp[] {
    const props = super.bindProps();
    for (const [name, sig] of this._liveBinds ?? []) {
      if (!props.some((p) => p.name === name))
        props.push({name, type: bindTypeOf(sig.peek()), writable: false}); // one-way (UB-10)
    }
    return props;
  }

  protected abstract createEditor(): HTMLElement;

  /** Puts an icon, a swatch or any other side control on the options rail, as Dart's
   * `InputBase.addOptions` does. The rail materializes with its first child, so an input with
   * nothing to put there keeps the platform's bare `[label][editor]` row. */
  protected addOptions(el: HTMLElement): void {
    this._optionsRail.append(el);
    if (this._optionsRail.parentElement === null)
      this._inputBox.append(this._optionsRail);
  }

  /** Re-applies the disabled state to controls a subclass rebuilt after construction — a fresh
   * `<input>` is enabled until told otherwise. */
  protected refreshEnabled(): void {
    if (!this._enabled)
      this.enabled = false;
  }

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
