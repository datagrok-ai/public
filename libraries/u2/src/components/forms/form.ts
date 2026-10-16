import {batch, computed, signal, untracked, ReadonlySignal, Signal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {Input, LiveOption} from '../../core/input-base.js';
import type {IItemsState} from '../../core/widget-like.js';

export type FormLayout = 'auto' | 'normal' | 'wide' | 'tall';

export interface FormOptions {
  /** `'auto'` (the default) picks by fit, following resizes: `'wide'` (caption left, editors
   * filling the row) while the label column plus `--dg-form-min-editor-width` fit, `'tall'`
   * (caption above) when they don't, with a 10px hysteresis; dialog-hosted and nested forms
   * never switch and keep the plain caption-left look, as the platform's do. Either way the
   * editors keep the widths `ui.form` gives them: text fills the row, a number is 150px, a choice
   * as wide as its longest option plus 30px. `'normal'`: caption left, those platform widths,
   * never switches. `'wide'`: caption left, every editor fills the row (`ui-form-wide`).
   * `'tall'`: caption above on its own small grey line, editor full-width below; bool stays
   * inline. */
  layout?: FormLayout;
  /** Caption text alignment in caption-left layouts; `'right'` (the platform default) or
   * `'left'`. Live: a signal is followed as a pure CSS class toggle. */
  captionAlign?: LiveOption<'right' | 'left'>;
  /** A fixed label column, in px: nothing is measured. For forms that must line up with something
   * else, or rows with one-letter captions (`R1`) that want a narrow column. */
  labelWidth?: number;
  /** One label column for several forms — the tabs of an editor, a header form over a body form: each
   * measures its own captions and all of them take the widest, so switching tabs moves no input. */
  labelColumn?: LabelColumn;
  /** Values read as text: an editor's underline shows while it is hovered, focused, empty or wrong,
   * not at rest — for an inspector or a docked editor, where a column of rules is noise. Any other
   * host gets the same look from the `u2-quiet` class. */
  quiet?: boolean;
}

export interface FormColumnsOptions {
  /** The narrowest a column may get, caption and editor together, in px; 200 by default. */
  minWidth?: number;
  /** The most columns a row holds; 2 by default. */
  columns?: number;
  /** Places the columns elsewhere, as {@link Form.add}'s `host` does. */
  host?: HTMLElement;
}

/** A label column shared by several forms ({@link FormOptions.labelColumn}): as wide as the widest
 * caption of any of them, a hidden one included (measured off its text). */
export class LabelColumn {
  /** The column's width in px; 0 until a member measured. */
  readonly width: ReadonlySignal<number>;
  private readonly _width = signal(0);
  private readonly _members = new Map<Form, number>();

  constructor() {
    this.width = this._width;
  }

  /** @internal A member's own widest caption. */
  report(form: Form, width: number): void {
    this._members.set(form, width);
    this._update();
  }

  /** @internal A disposed member leaves. */
  leave(form: Form): void {
    this._members.delete(form);
    this._update();
  }

  private _update(): void {
    let width = 0;
    for (const w of this._members.values())
      width = Math.max(width, w);
    this._width.value = width;
  }
}

/** Platform chrome shares the aligned label column: a bridged input brings its own `ui-` label. */
const LABELS = '.u2-input-label, .ui-input-label';

const _warnedUnkeyed = new Set<string>();
let _measure: CanvasRenderingContext2D | null | undefined;

/** Vertical input layout with an aligned label column and aggregated validity. The form does not
 * own the inputs it lays out: each one is already adopted by the scope that created it. */
export class Form extends Control {
  /** The first verdict of a counted input (visible and enabled) — required ones from the start. */
  readonly validity: ReadonlySignal<string | null>;
  /** The first message a counted input shows ({@link Input.error}): none on a fresh form. */
  readonly error: ReadonlySignal<string | null>;
  /** {@link validity} with the input's caption: "Name: Value can't be empty" — what a disabled OK or
   * Save says. */
  readonly problem: ReadonlySignal<string | null>;

  private readonly _rows = document.createElement('div');
  private readonly _inputs: Input<any>[] = [];
  private readonly _first = signal<ReadonlySignal<string | null>>(computed(() => null));
  private readonly _firstShown = signal<ReadonlySignal<string | null>>(computed(() => null));
  private readonly _firstNamed = signal<ReadonlySignal<string | null>>(computed(() => null));
  private readonly _layout: FormLayout;
  private readonly _labelWidth: number | undefined;
  private readonly _column: LabelColumn | undefined;
  private _buttons: HTMLElement | undefined;
  private _frame = 0;
  private _queued = false;
  private _tallBeforeMeasure = false;
  private static _groups = 0;
  private _minEditor: number | undefined;
  private readonly _labelMemo = new WeakMap<HTMLElement, number>();
  /** Re-fits the choice widths when a select's options change. */
  private _choices: MutationObserver | undefined;

  constructor(options: FormOptions = {}) {
    super();
    this.validity = computed(() => this._first.value.value);
    this.error = computed(() => this._firstShown.value.value);
    this.problem = computed(() => this._firstNamed.value.value);

    this.root.classList.add('u2-form');
    this.root.dataset.u2 = 'form';
    this._layout = options.layout ?? 'auto';
    this._labelWidth = options.labelWidth;
    if (this._labelWidth !== undefined)
      this.root.style.setProperty('--u2-form-label-width', `${this._labelWidth}px`);
    const column = this._column = options.labelColumn;
    if (column !== undefined) {
      this.own(() => column.leave(this));
      // a wider caption elsewhere moves this form's editors too, and may tip it into the tall layout
      this.effect(() => {
        const width = column.width.value;
        if (width > 0)
          untracked(() => this._applyLabelWidth(width));
      });
    }
    // the auto form's roomy state is wide, but its editors keep the platform widths ui.form gives
    if (this._layout === 'auto')
      this.root.classList.add('u2-form-auto');
    if (this._layout === 'wide')
      this.root.classList.add('u2-form-wide');
    if (this._layout === 'tall')
      this.root.classList.add('u2-form-tall');
    if (options.quiet)
      this.root.classList.add('u2-quiet');
    const align = options.captionAlign;
    if (align instanceof Signal)
      this.effect(() => this.root.classList.toggle('u2-form-captions-left', align.value === 'left'));
    else if (align === 'left')
      this.root.classList.add('u2-form-captions-left');
    this._rows.className = 'u2-form-rows';
    this.root.append(this._rows);

    const onKeyDown = (e: KeyboardEvent) => this._onKeyDown(e);
    this.root.addEventListener('keydown', onKeyDown);
    this.own(() => this.root.removeEventListener('keydown', onKeyDown));
    this.own(() => {
      if (this._frame)
        cancelAnimationFrame(this._frame);
      this._choices?.disconnect();
    });
    // any layout: a form built hidden (a collapsed section, a closed pencil) measures its labels at
    // 0 px and must align again once it shows; the auto layout also re-picks its state on a resize
    if (this._layout !== 'tall' && typeof ResizeObserver !== 'undefined') {
      const observer = new ResizeObserver(() => this._scheduleFrame());
      observer.observe(this.root);
      this.own(() => observer.disconnect());
    }
  }

  get inputs(): ReadonlyArray<Input<any>> {
    return this._inputs;
  }

  /** `host` lays the row out elsewhere — a Section body, say — while the input still joins the
   * form's label column, validity and Enter navigation. */
  add(input: Input<any>, host?: HTMLElement): Form {
    // dev-mode nudge (once per input kind): a keyless input has no getValues entry and no
    // automation id — nothing can address it
    if (input.name === undefined) {
      const kind = input.root.dataset.u2 ?? 'input';
      if (!_warnedUnkeyed.has(kind)) {
        _warnedUnkeyed.add(kind);
        console.warn(`u2: a ${kind} with neither name nor label was added to a form — ` +
          'it has no form key and no automation id; give it a name');
      }
    }
    this._inputs.push(input);
    // an input built in a row that goes away (an each() key, a rebuilt section) leaves with it
    input.own(() => this.remove(input));
    (host ?? this._rows).append(input.root);
    const select = input.root.dataset.u2 === 'choice-input' ? input.root.querySelector('select') : null;
    if (select !== null && typeof MutationObserver !== 'undefined')
      (this._choices ??= new MutationObserver(() => this._scheduleFrame())).observe(select, {childList: true});
    this._rebuildValidity();
    this._scheduleAlign();
    return this;
  }

  /** Takes `input` out of the form — its row, the label column, the validity and the Enter order.
   * The input stays its builder's to dispose; a disposed input leaves on its own. */
  remove(input: Input<any>): Form {
    const at = this._inputs.indexOf(input);
    if (at < 0)
      return this;
    this._inputs.splice(at, 1);
    input.root.remove();
    if (this.scope.isDisposed)
      return this;
    this._rebuildValidity();
    this._scheduleAlign();
    return this;
  }

  addAll(inputs: Input<any>[]): Form {
    for (const input of inputs)
      this.add(input);
    return this;
  }

  /** One caption in the label column over several inputs side by side — Colors: Line, Background,
   * Text. The inputs join the form's validity, values and Enter order; their own captions sit beside
   * their editors (after them with `captionsAfter`, as a checklist's do) and size no column. `host`
   * places the row elsewhere, as {@link add} does. */
  addGroup(caption: string, inputs: Input<any>[], options?: {captionsAfter?: boolean, host?: HTMLElement}): Form {
    const label = document.createElement('span');
    label.className = 'u2-input-label u2-form-group-label';
    label.id = `${this.rootId}-group-${++Form._groups}`;
    label.textContent = caption;
    const items = document.createElement('div');
    items.className = 'u2-form-group-items';
    const row = document.createElement('div');
    row.className = 'u2-form-group';
    row.setAttribute('role', 'group');
    row.setAttribute('aria-labelledby', label.id);
    row.append(label, items);
    (options?.host ?? this._rows).append(row);
    for (const input of inputs) {
      input.root.classList.add('u2-input-in-group');
      input.root.classList.toggle('u2-input-caption-after', options?.captionsAfter === true);
      this.add(input, items);
    }
    return this;
  }

  /** The inputs in equal columns, filled row by row — Symbol | Name, Polymer | Type. A column is never
   * narrower than `minWidth`: a row with no room for two holds one, so a narrow dock stacks them rather
   * than squeezing them. Every column keeps the form's label column, so captions and editors line up
   * down the rows and across them. The inputs join the form as {@link add}'s do. */
  addColumns(inputs: Input<any>[], options?: FormColumnsOptions): Form {
    const grid = document.createElement('div');
    grid.className = 'u2-form-columns';
    grid.dataset.u2Part = 'columns';
    if (options?.minWidth !== undefined)
      grid.style.setProperty('--u2-form-column-min', `${options.minWidth}px`);
    if (options?.columns !== undefined)
      grid.style.setProperty('--u2-form-columns', String(options.columns));
    (options?.host ?? this._rows).append(grid);
    for (const input of inputs)
      this.add(input, grid);
    return this;
  }

  /** A non-input row — a button, a heading, a plain element — kept in document order with the
   * input rows instead of trailing the form. */
  addElement(el: HTMLElement): Form {
    this._rows.append(el);
    return this;
  }

  /** The form's button row, under its rows: at the end by default, at the start (`'start'`), or
   * (`'editors'`) starting where the editors do, as `ui.buttonsInput` places the buttons of a platform
   * form. `host` puts the row elsewhere — a footer pinned under a form that scrolls — and its buttons
   * keep the platform's form-button look there. */
  addButtons(build: (row: HTMLElement) => void, options?: {align?: 'end' | 'editors' | 'start',
    host?: HTMLElement}): Form {
    if (!this._buttons) {
      this._buttons = document.createElement('div');
      this._buttons.className = 'u2-form-buttons';
      (options?.host ?? this.root).append(this._buttons);
    }
    if (options?.align !== undefined) {
      this._buttons.classList.toggle('u2-form-buttons-editors', options.align === 'editors');
      this._buttons.classList.toggle('u2-form-buttons-start', options.align === 'start');
    }
    build(this._buttons);
    return this;
  }

  /** Focus to the first input — where a form starts, and where it returns after Save; with
   * `where`, the first input it accepts, else the first. */
  focusFirst(where?: (input: Input<any>) => boolean): void {
    const input = (where === undefined ? undefined : this._inputs.find(where)) ?? this._inputs[0];
    input?.focus();
  }

  /** Shows every verdict the form counts — inputs that are visible and enabled — and focuses the
   * first failing one. */
  validate(): boolean {
    let failing: Input<any> | undefined;
    for (const input of this._inputs) {
      if (!input.validates)
        continue;
      input.touch();
      if (failing === undefined && input.validity.peek() !== null)
        failing = input;
    }
    failing?.focus();
    return failing === undefined;
  }

  /** The base state plus the form's verdict (the first message of a counted input) and one item per
   * input: key = its form key, label = its caption, `invalid`, `disabled`. */
  getState(): ReturnType<Control['getState']> {
    const error = this.validity.peek();
    return {...super.getState(), valid: error === null, error};
  }

  getItems(opts?: {all?: boolean}): IItemsState {
    const list = this._inputs.map((input, index) => ({key: input.name ?? `#${index}`, label: input.label, index,
      ...input.validity.peek() !== null ? {invalid: true} : {}, ...input.enabled ? {} : {disabled: true},
      rendered: input.root.isConnected}));
    return Input.itemsState(list, opts?.all);
  }

  getValues(): Record<string, unknown> {
    const values: Record<string, unknown> = {};
    for (const input of this._inputs) {
      if (input.name)
        values[input.name] = input.value.peek();
    }
    return values;
  }

  setValues(values: Record<string, unknown>): void {
    batch(() => {
      for (const input of this._inputs) {
        if (input.name && input.name in values)
          input.value.value = values[input.name];
      }
    });
  }

  private _rebuildValidity(): void {
    const inputs = this._inputs.slice();
    const first = (verdict: (input: Input<any>) => string | null) => computed(() => {
      for (const input of inputs) {
        const message = input.validates ? verdict(input) : null;
        if (message)
          return message;
      }
      return null;
    });
    batch(() => {
      this._first.value = first((input) => input.validity.value);
      this._firstShown.value = first((input) => input.error.value);
      this._firstNamed.value = first((input) => {
        const message = input.validity.value;
        return message !== null && input.label ? `${input.label}: ${message}` : message;
      });
    });
  }

  private _onKeyDown(e: KeyboardEvent): void {
    if (e.key !== 'Enter' || e.defaultPrevented || e.target instanceof HTMLTextAreaElement)
      return;
    const target = e.target as Node;
    const index = this._inputs.findIndex((input) => input.root.contains(target));
    if (index < 0 || index === this._inputs.length - 1)
      return;
    e.preventDefault();
    this._inputs[index + 1].focus();
  }

  /** Rows added while building align in a microtask: once the building code has run, before the
   * browser paints — a form built and attached in one go never shows a misaligned frame. */
  private _scheduleAlign(): void {
    if (this._queued)
      return;
    this._queued = true;
    queueMicrotask(() => {
      this._queued = false;
      this._alignLabels();
    });
  }

  /** A resize (or a form shown after being built hidden) aligns in the next frame: re-laying out
   * inside the observer's own callback would resize what it observes. */
  private _scheduleFrame(): void {
    if (this._frame)
      return;
    this._frame = requestAnimationFrame(() => {
      this._frame = 0;
      this._alignLabels();
    });
  }

  /** Lays the label column out now: measures the captions (unless the width is fixed), shares them
   * with the {@link LabelColumn}, fits the choices and picks the auto layout. Runs on its own after
   * every change; call it where a test or a host must see the final layout synchronously. */
  layoutNow(): void {
    this._alignLabels();
  }

  private _alignLabels(): void {
    if (this._layout === 'tall' || this.scope.isDisposed)
      return;
    if (this._labelWidth !== undefined) {
      this._applyLabelWidth(this._labelWidth);
      return;
    }
    const width = this._measureLabels();
    this._column?.report(this, width);
    this._applyLabelWidth(this._column?.width.peek() ?? width);
  }

  // natural label widths only read back with the alignment (and the tall skin) dropped —
  // min-width clamps them, and a tall label spans the full row
  private _measureLabels(): number {
    this._tallBeforeMeasure ||= this.root.classList.contains('u2-form-tall');
    this.root.classList.remove('u2-form-tall');
    this.root.style.setProperty('--u2-form-label-width', '0px');
    // a form nested in one gone tall would otherwise read its captions at the tall skin's full row
    this.root.classList.add('u2-form-measuring');
    const captions: HTMLElement[] = [...this.root.querySelectorAll<HTMLElement>('.u2-form-group-label')];
    for (const input of this._inputs) {
      // a checklist row's caption follows its box, and a grouped input's sits beside it: no column
      if (!input.root.classList.contains('u2-input-caption-after') && !input.root.classList.contains('u2-input-in-group'))
        captions.push(...input.root.querySelectorAll<HTMLElement>(LABELS));
    }
    let width = 0;
    for (const caption of captions) {
      // a hidden label (collapsed section, an inactive tab) measures 0 — its remembered width
      // stands in, else its text measured in its font, so a hidden form still sizes a shared
      // column; the fraction rounds up, as the platform's column does (js-api ui.ts getLabelWidth)
      let w = Math.max(caption.offsetWidth, Math.ceil(caption.getBoundingClientRect().width));
      if (w)
        this._labelMemo.set(caption, w);
      else
        w = this._labelMemo.get(caption) ?? Form._textWidth(caption);
      width = Math.max(width, w);
    }
    this.root.classList.remove('u2-form-measuring');
    return width;
  }

  private _applyLabelWidth(width: number): void {
    if (this.scope.isDisposed)
      return;
    const wasTall = this.root.classList.contains('u2-form-tall') || this._tallBeforeMeasure;
    this._tallBeforeMeasure = false;
    this.root.classList.remove('u2-form-tall');
    if (width)
      this.root.style.setProperty('--u2-form-label-width', `${width}px`);
    else
      this.root.style.removeProperty('--u2-form-label-width');
    Form.fitChoices(this.root, this._inputs.map((input) => input.root));
    if (this._layout !== 'auto')
      return;
    // the platform exempts dialog forms from the auto switch (js-api ui.ts:1437-1441), and a
    // form nested inside another form leaves the switching to the outer one (ui.ts:1419-1424);
    // exempt forms keep the plain caption-left look — the platform dialog form is the aligned
    // skin, never wide — so both auto states come off
    if (this.root.closest('.u2-dialog') !== null ||
        this.root.parentElement?.closest('.u2-form') != null) {
      this.root.classList.remove('u2-form-wide');
      return;
    }
    // clientWidth MUST be read after the tall class came off: the flip can move scrollbars and
    // paddings, and deciding from the tall-state width would compare against the wrong geometry
    const next = Form._pickLayout(this.root.clientWidth, width, this._minEditorWidth(),
      wasTall ? 'tall' : this.root.classList.contains('u2-form-wide') ? 'wide' : 'normal');
    this.root.classList.toggle('u2-form-tall', next === 'tall');
    this.root.classList.toggle('u2-form-wide', next === 'wide');
  }

  /** A caption's width from its text, in its own font — for one not rendered (an inactive tab). */
  private static _textWidth(label: HTMLElement): number {
    const text = label.textContent ?? '';
    if (text === '' || !label.isConnected)
      return 0;
    _measure ??= document.createElement('canvas').getContext?.('2d') ?? null;
    if (_measure === null)
      return 0;
    const style = getComputedStyle(label);
    _measure.font = style.font || `${style.fontSize} ${style.fontFamily}`;
    let width = _measure.measureText(text).width + (parseFloat(style.paddingLeft) || 0) +
      (parseFloat(style.paddingRight) || 0);
    if (label.parentElement?.classList.contains('u2-input-required'))
      width += _measure.measureText('*').width + (parseFloat(getComputedStyle(label, '::after').marginLeft) || 0);
    return Math.ceil(width);
  }

  /** The platform's fit rule (js-api ui.ts:1443-1448): tall when the label column leaves less
   * than the minimum editor width, WIDE when it leaves 10px more than that — the hysteresis —
   * and the current state everywhere in between or when the host is hidden (width 0). The
   * roomy state is wide by user ruling (2026-08-28); `'normal'` only ever ARRIVES as the
   * never-measured or exempt state — the fit test never returns it. */
  static _pickLayout(width: number, labelWidth: number, minEditor: number,
    current: 'normal' | 'wide' | 'tall'): 'normal' | 'wide' | 'tall' {
    if (width <= 0)
      return current;
    if (width - labelWidth < minEditor)
      return 'tall';
    if (width - labelWidth > minEditor + 10)
      return 'wide';
    return current;
  }

  // the fallback is never cached: a form measured before its stylesheet applied (detached,
  // display:none boot) would otherwise pin 150 and ignore the token forever
  private _minEditorWidth(): number {
    if (this._minEditor === undefined) {
      const v = parseFloat(getComputedStyle(this.root).getPropertyValue('--dg-form-min-editor-width'));
      if (!Number.isFinite(v))
        return 150;
      this._minEditor = v;
    }
    return this._minEditor;
  }

  /** The platform's choice width (js-api ui.ts getInputsMinWidths) for the choice inputs among
   * `roots`: the longest option's text, at least 40px, plus 30px for the chevron, at most 80% of
   * the window. Measured off a probe in `host`'s own font — a form, or a dialog body that lays
   * inputs out as one; a detached host keeps the stylesheet's width. */
  static fitChoices(host: HTMLElement, roots: Iterable<HTMLElement>): void {
    if (!host.isConnected)
      return;
    const probe = document.createElement('span');
    probe.style.cssText = 'position: absolute; visibility: hidden; white-space: pre; padding: 0; margin: 0;';
    host.append(probe);
    for (const root of roots) {
      if (root.dataset.u2 !== 'choice-input')
        continue;
      let width = 40;
      for (const option of root.querySelectorAll('option')) {
        probe.textContent = option.textContent;
        width = Math.max(width, Math.ceil(probe.getBoundingClientRect().width));
      }
      root.style.setProperty('--u2-choice-width', `${Math.min(width + 30, innerWidth * 0.8)}px`);
    }
    probe.remove();
  }
}
