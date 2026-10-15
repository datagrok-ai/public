import {signal, batch, Signal, ReadonlySignal, untracked} from './signals.js';
import {Scope} from './scope.js';
import {bindTypeOf, ITEMS_CAP} from './widget-like.js';
import type {BindProp, BindSource, ObservableLike, IItemState, IItemsState} from './widget-like.js';
import {Control} from './component.js';
import {div, label, span} from './elements.js';
import {EMPTY, isEmpty} from './text.js';
import type {IWidgetStatus} from './widget-like.js';

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
  /** Runs after each change of the value, owned by whoever built the input (the ambient scope at
   * construction; the input where there was none or once it is gone): what it opens survives the
   * next change. */
  onChanged?: (value: T, input: Input<T>) => void;
  tooltipText?: LiveOption<string>;
  /** When an edit reaches {@link Input.value}: on every keystroke (`input`, the default) or on the
   * commit — blur or Enter — as Dart's `fireChanged` does (`change`). The free-text editors (text,
   * number, bigint, qnum) have a gap between the two, and the slider commits on release; every
   * other editor commits as it edits. Worth asking for where a consumer acts on the value
   * irreversibly, as a rename does. */
  commitOn?: 'input' | 'change';
  /** `false` hides the input (the root's `hidden`); live. A hidden input never fails its form. */
  visible?: LiveOption<boolean>;
  /** Shows the value and lets it be selected and copied, but blocks edits; live. */
  readOnly?: LiveOption<boolean>;
  /** An empty value ({@link Input.isEmptyValue}) is invalid; live. `nullable: false` implies it. The
   * verdict shows once the input was edited or left, or its form validated ({@link Input.touch}). */
  required?: LiveOption<boolean>;
  /** Enter in the editor: handled here, never seen by the form or dialog around the input. */
  onEnter?: (e: KeyboardEvent, input: Input<T>) => void;
  /** Escape in the editor: handled here, never seen by the dialog around the input. */
  onEscape?: (e: KeyboardEvent, input: Input<T>) => void;
  /** The editor fills the row, as the text editors do in a wide form. */
  fill?: boolean;
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
  private readonly _enabled = signal(true);
  private readonly _shown = signal(true);
  private readonly _required = signal(false);
  /** The verdict is the required check's: it waits for {@link touch}. */
  private readonly _missing = signal(false);
  private readonly _touched = signal(false);
  private _readOnly = false;
  private _guarded = false;
  private _touchWired = false;
  private _onChanged: ((value: T, input: Input<T>) => void) | undefined;
  private _watching = false;
  /** The value {@link follow} last wrote, until the change watch has skipped it. */
  private _followed: {value: T} | undefined;

  constructor(options: O, defaultValue: T) {
    super();
    this.options = options;
    this.value = options.bind ?? signal(options.value ?? defaultValue);
    this.validity = this._validity;
    // an explicit name is a deliberate identity and stamps data-u2-name; the label fallback is a
    // form key only — captions get reworded, so they never become automation ids
    if (options.name !== undefined)
      this.name = options.name;
    else
      this.deriveName(labelText(options.label));

    this.root.classList.add('u2-input-root');
    if (options.inline)
      this.root.classList.add('u2-input-inline');
    if (options.fill)
      this.root.classList.add('u2-input-fill');

    this._inputBox = div([], 'u2-input-box');
    this._optionsRail = div([], 'u2-input-options');
    this._optionsRail.dataset.u2Part = 'options';
    this._error.dataset.u2Part = 'error';
    this._editor = this.runInScope(() => this.createEditor());
    this._editor.classList.add('u2-input-editor');
    this._editor.dataset.u2Part = 'editor';
    // prepended: a subclass that filled the rail from createEditor already attached it
    this._inputBox.prepend(this._editor);
    this.root.append(this._inputBox, this._error);
    const target = this.focusTarget;
    if (!target.id)
      target.id = `${this.rootId}-editor`;
    this._error.id = `${this.rootId}-error`;
    if (options.inline && options.name !== undefined)
      target.setAttribute('aria-label', options.name);
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
      else if (!this._enabled.peek())
        this.enabled = true;
    });
    this.liveOption('visible', options.visible, (x) => this.visible = x);
    this.liveOption('readOnly', options.readOnly, (x) => this.readOnly = x);
    this.liveOption('required', options.required ?? (options.nullable === false ? true : undefined),
      (x) => this.required = x);

    this.effect(() => this._validate(this.value.value));
    this.effect(() => this._showValidity(this._validity.value !== null && this._missing.value &&
      !this._touched.value ? null : this._validity.value));
    if (options.onChanged) {
      this._onChanged = Input.bindOwned(this.scope, options.onChanged) as (value: T, input: Input<T>) => void;
      this._watchChanges();
    }
    const {onEnter, onEscape} = options;
    if (onEnter || onEscape) {
      const keys: Record<string, ((e: KeyboardEvent, input: Input<T>) => void) | undefined> = {
        Enter: onEnter && Input.bindOwned(this.scope, onEnter),
        Escape: onEscape && Input.bindOwned(this.scope, onEscape)};
      this.scope.listen(target, 'keydown', (e) => {
        const handler = keys[(e as KeyboardEvent).key];
        if (handler === undefined || (e as KeyboardEvent).isComposing)
          return;
        e.preventDefault();
        e.stopPropagation();
        handler(e as KeyboardEvent, this);
      });
    }
  }

  /** True for `null`, `undefined`, `''`, an empty array and a plain object without keys — what
   * `required` refuses. `false`, `0` and a `Date` are values. */
  static isEmpty(value: unknown): boolean {
    if (isEmpty(value) || Array.isArray(value) && value.length === 0)
      return true;
    if (typeof value !== 'object')
      return false;
    const proto = Object.getPrototypeOf(value);
    return (proto === Object.prototype || proto === null) && Object.keys(value as object).length === 0;
  }

  /** `fn` run as {@link Scope.bind} runs it, in the scope building the caller — or in `own` where none
   * is, or once that one is gone: a dialog adopts its content and keeps it on screen past its builder,
   * and the callbacks must not go quiet. */
  static bindOwned<A extends unknown[], R>(own: Scope, fn: (...args: A) => R): (...args: A) => R | undefined {
    const owner = Scope.ambient;
    const fallback = own.bind(fn);
    if (owner === undefined || owner === own)
      return fallback;
    const bound = owner.bind(fn);
    return (...args: A) => owner.isDisposed ? fallback(...args) : bound(...args);
  }

  /** The items block of a control listing its choices: all of them, or the first {@link ITEMS_CAP}. */
  static itemsState(list: IItemState[], all?: boolean): IItemsState {
    const truncated = !all && list.length > ITEMS_CAP;
    return {count: list.length, list: truncated ? list.slice(0, ITEMS_CAP) : list, truncated};
  }

  /** One watch of the value for `onChanged` and the `change` event: the initial value and a
   * {@link follow} write are not changes. Built on first need. */
  private _watchChanges(): void {
    if (this._watching)
      return;
    this._watching = true;
    let initial = true;
    this.effect(() => {
      const value = this.value.value;
      const followed = this._followed;
      this._followed = undefined;
      if (initial)
        initial = false;
      else if (followed === undefined || !Object.is(followed.value, value)) {
        untracked(() => {
          this._onChanged?.(value, this);
          this.fireEvent('change', value);
        });
      }
    });
  }

  /** `change` carries the new value (`input.on('change', (value, input) => …)`), alongside the
   * DOM `change` events the root sees; every other id as {@link Control.onEvent}. */
  onEvent(eventId: string | null = null): ObservableLike<unknown> {
    if (eventId === null || eventId === 'change')
      this._watchChanges();
    return super.onEvent(eventId);
  }

  /** Keeps the value on `read` (the model field the input edits): re-read whenever a signal `read`
   * reads changes, written only when it differs, never echoed back through `onChanged`. A panel
   * built once over a live model thus shows edits made elsewhere without a rebuild, and the field
   * being typed in keeps its focus and caret. One-way; owned by the input. */
  follow(read: () => T): this {
    return this.effect(() => {
      const value = read();
      if (Object.is(value, this.value.peek()))
        return;
      this._followed = {value};
      Input.system(() => this.value.value = value);
    });
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

  /** Creates the label part on first set, tied to the {@link focusTarget}: `for` when it is labelable,
   * else its `aria-labelledby`. An inline input shows no label; the text becomes its `aria-label`. */
  set label(x: string) {
    if (this._labelEl)
      this._labelEl.textContent = x;
    else if (this.options.inline)
      this.focusTarget.setAttribute('aria-label', x);
    else {
      this._labelEl = label(x, 'u2-input-label');
      this._labelEl.dataset.u2Part = 'label';
      this._labelEl.id = `${this.rootId}-label`;
      const target = this.focusTarget;
      if (target.matches('input, textarea, select, button'))
        this._labelEl.setAttribute('for', target.id);
      else
        target.setAttribute('aria-labelledby', this._labelEl.id);
      this.root.prepend(this._labelEl);
    }
  }

  /** The element keyboard focus lands on — the native editor, or the element in the editor carrying
   * the role or the tab stop; it carries the id, `aria-invalid` and `aria-describedby`. */
  protected get focusTarget(): HTMLElement {
    const editor = this._editor;
    if (editor.matches('input, textarea, select') || editor.hasAttribute('role'))
      return editor;
    return editor.querySelector<HTMLElement>('input, textarea, select, [role="combobox"], [role="listbox"], ' +
      '[role="slider"], [role="spinbutton"], [role="switch"], [contenteditable="true"], ' +
      '[tabindex]:not([tabindex="-1"])') ?? editor;
  }

  /** The base state plus the value, raw (JSON-projected) and as displayed, the validity — `error` is
   * the message shown, so a required verdict waiting for {@link touch} reads `valid: false` with
   * `error: null` — and the read-only and required flags. */
  getState(): ReturnType<Control['getState']> {
    const error = this._validity.peek();
    const placeholder = this._placeholder();
    return {...super.getState(), value: Control.jsonValue(this.value.peek()) ?? null, text: this.displayText(),
      valid: error === null, error: this._error.textContent || null, readOnly: this._readOnly,
      required: this._required.peek(), ...placeholder ? {placeholder} : {}};
  }

  get visible(): boolean {
    return super.visible;
  }

  set visible(x: boolean) {
    super.visible = x;
    this._shown.value = x;
  }

  /** Whether a form counts this input's verdict: visible — no `hidden` on the root or above it —
   * and enabled. Reads signals, so a computed over it follows {@link visible} and {@link enabled};
   * an ancestor hidden by other means is seen on the next read. */
  get validates(): boolean {
    return this._shown.value && this._enabled.value && this.root.closest('[hidden]') === null;
  }

  get readOnly(): boolean {
    return this._readOnly;
  }

  /** Native `readOnly` on the text editors (the value stays selectable), `aria-readonly` on the
   * {@link focusTarget}, and every other editing gesture — a click, a key, the wheel — stopped at
   * the root. */
  set readOnly(x: boolean) {
    this._readOnly = x;
    this.root.classList.toggle('u2-input-readonly', x);
    if (x)
      this.focusTarget.setAttribute('aria-readonly', 'true');
    else
      this.focusTarget.removeAttribute('aria-readonly');
    for (const el of this.root.querySelectorAll<HTMLInputElement>('input, textarea')) {
      if (Input._typed(el))
        el.readOnly = x;
    }
    if (x && !this._guarded) {
      this._guarded = true;
      for (const type of ['keydown', 'mousedown', 'click', 'input', 'change'])
        this.scope.listen(this.root, type, (e) => this._guard(e), {capture: true});
      this.scope.listen(this.root, 'wheel', (e) => this._guard(e), {capture: true, passive: true});
    }
  }

  get required(): boolean {
    return this._required.peek();
  }

  /** `aria-required` on the {@link focusTarget}, and an empty value fails the input once it is
   * {@link touch}ed. */
  set required(x: boolean) {
    this._required.value = x;
    if (x)
      this.focusTarget.setAttribute('aria-required', 'true');
    else
      this.focusTarget.removeAttribute('aria-required');
    if (x && !this._touchWired) {
      this._touchWired = true;
      for (const type of ['focusout', 'input', 'change'])
        this.scope.listen(this.root, type, () => this.touch());
    }
  }

  /** Shows the verdicts that wait for the user — a required field left empty. The input touches
   * itself when edited or left; {@link Form.validate} touches every input it counts. */
  touch(): void {
    this._touched.value = true;
  }

  /** Focuses the editor — the {@link focusTarget}, or the first focusable control inside it. */
  focus(): void {
    const target = this.focusTarget;
    const focusable = target.matches(Input._FOCUSABLE) ? target :
      target.querySelector<HTMLElement>(Input._FOCUSABLE) ?? target;
    focusable.focus();
  }

  private static readonly _FOCUSABLE = 'input, textarea, select, button, [tabindex]';

  /** Input types without free text: native `readOnly` does not apply to them. */
  private static readonly _UNTYPED = ['checkbox', 'radio', 'range', 'button', 'submit', 'reset', 'color', 'file',
    'image'];

  /** A free-text field — native `readOnly` keeps its value from edits and leaves the caret free. */
  private static _typed(el: EventTarget | null): boolean {
    const tag = (el as Element | null)?.tagName;
    return tag === 'TEXTAREA' || tag === 'INPUT' && !Input._UNTYPED.includes((el as HTMLInputElement).type);
  }

  /** Read-only: a free-text field keeps the caret, selection, Tab, Escape and Enter; every other
   * gesture — a click on a box, a stepping key or the wheel on a number, a toggle key — ends here. */
  private _guard(e: Event): void {
    if (!this._readOnly)
      return;
    const typed = Input._typed(e.target);
    const key = (e as KeyboardEvent).key;
    const free = e.type === 'keydown' ?
      key === 'Tab' || key === 'Escape' || key === 'Enter' ||
        typed && (!(key === 'ArrowUp' || key === 'ArrowDown') || (e.target as Element).tagName === 'TEXTAREA') :
      typed && e.type !== 'wheel';
    if (free)
      return;
    if (e.type !== 'wheel')
      e.preventDefault();
    e.stopImmediatePropagation();
  }

  /** The value as the user reads it: the chosen option's caption, the text in the editor, a checkbox's
   * 'true'/'false'. */
  protected displayText(): string {
    const target = this.focusTarget;
    if (target.tagName === 'SELECT') {
      const select = target as HTMLSelectElement;
      const option = [...select.querySelectorAll('option')].find((o) => o.value === select.value);
      return option === undefined ? '' : Control.text(option);
    }
    if (target.tagName === 'INPUT' && ['checkbox', 'radio'].includes((target as HTMLInputElement).type))
      return String((target as HTMLInputElement).checked);
    if (target.tagName === 'INPUT' || target.tagName === 'TEXTAREA')
      return (target as HTMLInputElement).value;
    return target.getAttribute('aria-checked') ?? Control.text(this._editor);
  }

  /** Label part, then the focus target's `aria-label`/`aria-labelledby`, then the placeholder. */
  protected accessibleLabel(): string | undefined {
    const own = this._labelEl === undefined ? '' : Control.text(this._labelEl);
    return own || Control.labelOf(this.focusTarget) || super.accessibleLabel() || this._placeholder();
  }

  private _placeholder(): string | undefined {
    const target = this.focusTarget as HTMLInputElement;
    return target.placeholder || target.getAttribute('placeholder') || undefined;
  }

  /** The base status plus the structural parts every input shares — the same elements the
   * `data-u2-part` attributes mark, so pixel-driving and state-reading automation agree. */
  getWidgetStatus(): IWidgetStatus {
    const status = super.getWidgetStatus();
    status.parts = {editor: this._editor, error: this._error};
    if (this._labelEl)
      status.parts.label = this._labelEl;
    if (this._optionsRail.parentElement !== null)
      status.parts.options = this._optionsRail;
    return status;
  }

  get tooltipText(): string {
    return this.root.title;
  }

  set tooltipText(x: string) {
    this.root.title = x;
  }

  get enabled(): boolean {
    return this._enabled.peek();
  }

  set enabled(x: boolean) {
    this._enabled.value = x;
    this.root.classList.toggle('u2-input-disabled', !x);
    if (x)
      this.root.removeAttribute('aria-disabled');
    else
      this.root.setAttribute('aria-disabled', 'true');
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

  /** What `required` refuses: {@link isEmpty}, unless the value has its own nothing — a filter
   * without conditions. */
  protected isEmptyValue(value: T): boolean {
    return Input.isEmpty(value);
  }

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
    if (!this._enabled.peek())
      this.enabled = false;
    if (this._readOnly)
      this.readOnly = true;
  }

  /** The added validators first; the required check speaks only when they have nothing to say. */
  private _validate(value: T): void {
    let message: string | null = null;
    for (const validator of this._validators) {
      message = validator(value);
      if (message)
        break;
    }
    const missing = !message && this._required.value && this.isEmptyValue(value);
    batch(() => {
      this._missing.value = missing;
      this._validity.value = missing ? EMPTY : message || null;
    });
  }

  private _showValidity(message: string | null): void {
    this._editor.classList.toggle('u2-invalid', message !== null);
    this._error.textContent = message ?? '';
    const target = this.focusTarget;
    if (message === null) {
      this._editor.removeAttribute('title');
      target.removeAttribute('aria-invalid');
      target.removeAttribute('aria-describedby');
    } else {
      this._editor.title = message;
      target.setAttribute('aria-invalid', 'true');
      target.setAttribute('aria-describedby', this._error.id);
    }
  }
}
