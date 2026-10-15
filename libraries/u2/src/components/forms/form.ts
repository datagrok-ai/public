import {batch, computed, signal, ReadonlySignal, Signal} from '../../core/signals.js';
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
}

/** Platform chrome shares the aligned label column: a bridged input brings its own `ui-` label. */
const LABELS = '.u2-input-label, .ui-input-label';

const _warnedUnkeyed = new Set<string>();

/** Vertical input layout with an aligned label column and aggregated validity. The form does not
 * own the inputs it lays out: each one is already adopted by the scope that created it. */
export class Form extends Control {
  readonly validity: ReadonlySignal<string | null>;

  private readonly _rows = document.createElement('div');
  private readonly _inputs: Input<any>[] = [];
  private readonly _first = signal<ReadonlySignal<string | null>>(computed(() => null));
  private readonly _layout: FormLayout;
  private _buttons: HTMLElement | undefined;
  private _frame = 0;
  private _minEditor: number | undefined;
  private readonly _labelMemo = new WeakMap<HTMLElement, number>();
  /** Re-fits the choice widths when a select's options change. */
  private _choices: MutationObserver | undefined;

  constructor(options: FormOptions = {}) {
    super();
    this.validity = computed(() => this._first.value.value);

    this.root.classList.add('u2-form');
    this.root.dataset.u2 = 'form';
    this._layout = options.layout ?? 'auto';
    // the auto form's roomy state is wide, but its editors keep the platform widths ui.form gives
    if (this._layout === 'auto')
      this.root.classList.add('u2-form-auto');
    if (this._layout === 'wide')
      this.root.classList.add('u2-form-wide');
    if (this._layout === 'tall')
      this.root.classList.add('u2-form-tall');
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
      const observer = new ResizeObserver(() => this._scheduleAlign());
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
      (this._choices ??= new MutationObserver(() => this._scheduleAlign())).observe(select, {childList: true});
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

  /** A non-input row — a button, a heading, a plain element — kept in document order with the
   * input rows instead of trailing the form. */
  addElement(el: HTMLElement): Form {
    this._rows.append(el);
    return this;
  }

  /** The form's button row, under its rows: at the end by default, or (`align: 'editors'`) starting
   * where the editors do, as `ui.buttonsInput` places the buttons of a platform form. */
  addButtons(build: (row: HTMLElement) => void, options?: {align?: 'end' | 'editors'}): Form {
    if (!this._buttons) {
      this._buttons = document.createElement('div');
      this._buttons.className = 'u2-form-buttons';
      this.root.append(this._buttons);
    }
    if (options?.align !== undefined)
      this._buttons.classList.toggle('u2-form-buttons-editors', options.align === 'editors');
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
    this._first.value = computed(() => {
      for (const input of inputs) {
        const message = input.validates ? input.validity.value : null;
        if (message)
          return message;
      }
      return null;
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

  private _scheduleAlign(): void {
    if (this._frame)
      return;
    this._frame = requestAnimationFrame(() => {
      this._frame = 0;
      this._alignLabels();
    });
  }

  // natural label widths only read back with the alignment (and the tall skin) dropped —
  // min-width clamps them, and a tall label spans the full row
  private _alignLabels(): void {
    if (this._layout === 'tall')
      return;
    const wasTall = this.root.classList.contains('u2-form-tall');
    this.root.classList.remove('u2-form-tall');
    this.root.style.setProperty('--u2-form-label-width', '0px');
    let width = 0;
    for (const input of this._inputs) {
      // a checklist row's caption follows its box and reserves no column
      if (input.root.classList.contains('u2-input-caption-after'))
        continue;
      const labels = input.root.querySelectorAll<HTMLElement>(LABELS);
      for (let i = 0; i < labels.length; i++) {
        // a hidden label (collapsed section, expression-hidden row) measures 0 — its remembered
        // width stands in, so folding a section never reflows the shared column; the fraction
        // rounds up, as the platform's column does (js-api ui.ts getLabelWidth)
        let w = Math.max(labels[i].offsetWidth, Math.ceil(labels[i].getBoundingClientRect().width));
        if (w)
          this._labelMemo.set(labels[i], w);
        else
          w = this._labelMemo.get(labels[i]) ?? 0;
        width = Math.max(width, w);
      }
    }
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
