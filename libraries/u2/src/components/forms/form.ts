import {batch, computed, signal, ReadonlySignal, Signal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {Input, LiveOption} from '../../core/input-base.js';

export type FormLayout = 'auto' | 'normal' | 'wide' | 'tall';

export interface FormOptions {
  /** `'auto'` (the default) picks by fit, following resizes: `'wide'` (caption left, editors
   * filling the row) while the label column plus `--dg-form-min-editor-width` fit, `'tall'`
   * (caption above) when they don't, with a 10px hysteresis; dialog-hosted and nested forms
   * never switch and keep the plain caption-left look, as the platform's do. `'normal'`:
   * caption left, platform-pinned editor widths, never switches. `'wide'`: caption left,
   * editors fill the row. `'tall'`: caption above on its own small grey line, editor
   * full-width below; bool stays inline. */
  layout?: FormLayout;
  /** Caption text alignment in caption-left layouts; `'right'` (the platform default) or
   * `'left'`. Live: a signal is followed as a pure CSS class toggle. */
  captionAlign?: LiveOption<'right' | 'left'>;
}

const FOCUSABLE = 'input, textarea, select, button, [tabindex]';
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

  constructor(options: FormOptions = {}) {
    super();
    this.validity = computed(() => this._first.value.value);

    this.root.classList.add('u2-form');
    this.root.dataset.u2 = 'form';
    this._layout = options.layout ?? 'auto';
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
    });
    if (this._layout === 'auto' && typeof ResizeObserver !== 'undefined') {
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
    (host ?? this._rows).append(input.root);
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

  addButtons(build: (row: HTMLElement) => void): Form {
    if (!this._buttons) {
      this._buttons = document.createElement('div');
      this._buttons.className = 'u2-form-buttons';
      this.root.append(this._buttons);
    }
    build(this._buttons);
    return this;
  }

  validate(): boolean {
    for (const input of this._inputs) {
      if (input.validity.peek() !== null) {
        Form._focus(input);
        return false;
      }
    }
    return true;
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
        const message = input.validity.value;
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
    Form._focus(this._inputs[index + 1]);
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
      const labels = input.root.querySelectorAll<HTMLElement>(LABELS);
      for (let i = 0; i < labels.length; i++) {
        // a hidden label (collapsed section, expression-hidden row) measures 0 — its remembered
        // width stands in, so folding a section never reflows the shared column
        let w = labels[i].offsetWidth;
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

  private static _focus(input: Input<any>): void {
    const editor = input.root.querySelector('.u2-input-editor') as HTMLElement | null;
    if (!editor)
      return;
    const focusable = editor.matches(FOCUSABLE) ? editor : editor.querySelector<HTMLElement>(FOCUSABLE);
    (focusable ?? editor).focus();
  }
}
