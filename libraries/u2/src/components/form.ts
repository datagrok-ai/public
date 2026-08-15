import {batch, computed, signal, ReadonlySignal} from '../core/signals.js';
import {Component} from '../core/component.js';
import {Input} from '../core/input-base.js';

export interface FormOptions {
  condensed?: boolean;
  wide?: boolean;
}

const FOCUSABLE = 'input, textarea, select, button, [tabindex]';
/** Platform chrome shares the aligned label column: a bridged input brings its own `ui-` label. */
const LABELS = '.u2-input-label, .ui-input-label';

/** Vertical input layout with an aligned label column and aggregated validity. The form does not
 * own the inputs it lays out: each one is already adopted by the scope that created it. */
export class Form extends Component {
  readonly validity: ReadonlySignal<string | null>;

  private readonly _rows = document.createElement('div');
  private readonly _inputs: Input<any>[] = [];
  private readonly _first = signal<ReadonlySignal<string | null>>(computed(() => null));
  private _buttons: HTMLElement | undefined;
  private _frame = 0;

  constructor(options: FormOptions = {}) {
    super();
    this.validity = computed(() => this._first.value.value);

    this.root.classList.add('u2-form');
    this.root.dataset.u2 = 'form';
    if (options.condensed)
      this.root.classList.add('u2-form-condensed');
    if (options.wide)
      this.root.classList.add('u2-form-wide');
    this._rows.className = 'u2-form-rows';
    this.root.append(this._rows);

    const onKeyDown = (e: KeyboardEvent) => this._onKeyDown(e);
    this.root.addEventListener('keydown', onKeyDown);
    this.own(() => this.root.removeEventListener('keydown', onKeyDown));
    this.own(() => {
      if (this._frame)
        cancelAnimationFrame(this._frame);
    });
  }

  get inputs(): ReadonlyArray<Input<any>> {
    return this._inputs;
  }

  add(input: Input<any>): Form {
    this._inputs.push(input);
    this._rows.append(input.root);
    this._rebuildValidity();
    this._scheduleAlign();
    return this;
  }

  addAll(inputs: Input<any>[]): Form {
    for (const input of inputs)
      this.add(input);
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

  // natural label widths only read back with the alignment dropped — min-width clamps them
  private _alignLabels(): void {
    this.root.style.setProperty('--u2-form-label-width', '0px');
    let width = 0;
    for (const input of this._inputs) {
      const labels = input.root.querySelectorAll<HTMLElement>(LABELS);
      for (let i = 0; i < labels.length; i++)
        width = Math.max(width, labels[i].offsetWidth);
    }
    if (width)
      this.root.style.setProperty('--u2-form-label-width', `${width}px`);
    else
      this.root.style.removeProperty('--u2-form-label-width');
  }

  private static _focus(input: Input<any>): void {
    const editor = input.root.querySelector('.u2-input-editor') as HTMLElement | null;
    if (!editor)
      return;
    const focusable = editor.matches(FOCUSABLE) ? editor : editor.querySelector<HTMLElement>(FOCUSABLE);
    (focusable ?? editor).focus();
  }
}
