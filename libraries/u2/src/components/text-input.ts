import {Input, InputOptions} from '../core/input-base.js';
import {bindValue} from '../core/bind.js';
import {icon} from './icon.js';

export interface TextInputOptions extends InputOptions<string> {
  /** Magnifier affordance plus a clear button, revealed on hover as the platform's is
   * (`text_input.dart:131-141` wires it through `htmlMouseOverVisibility`). */
  search?: boolean;
  placeholder?: string;
  /** Masked field with an eye toggle. The platform's `PASSWORD_PLACEHOLDER` focus/blur dance
   * (`text_input.dart:54-73`) is a credentials convention and stays out of the core input. */
  password?: boolean;
  /** Width follows the text between {@link minWidth} and {@link maxWidth}, as Dart's autoResize. */
  autoResize?: boolean;
  minWidth?: number;
  maxWidth?: number;
}

/** Single-line text field. The clear icon of the search variant sits absolutely at the right edge
 * of the field and only shows while the pointer is over the input, as the platform's does
 * (`text_input.dart:131-141`, `ui.css:540-546`); u2 additionally keeps it out of an empty field,
 * where there is nothing to clear. */
export class TextInput extends Input<string, TextInputOptions> {
  private _input!: HTMLInputElement;

  constructor(options: TextInputOptions = {}) {
    super(options, '');
    this.root.dataset.u2 = 'text-input';
  }

  protected createEditor(): HTMLElement {
    const input = document.createElement('input');
    this._input = input;
    input.type = this.options.password ? 'password' : 'text';
    input.autocomplete = 'off';
    input.placeholder = this.options.placeholder ?? '';
    bindValue(this.scope, input, this.value, this.options.commitOn);
    if (this.options.autoResize)
      this._autoResize();
    if (this.options.password)
      return this._passwordEditor();
    if (!this.options.search)
      return input;

    const editor = document.createElement('div');
    editor.className = 'u2-input-search';
    const clear = document.createElement('button');
    clear.type = 'button';
    clear.className = 'u2-input-clear';
    clear.textContent = '✕';
    clear.setAttribute('aria-label', 'Clear');
    const onClear = () => {
      this.value.value = '';
      input.focus();
    };
    clear.addEventListener('click', onClear);
    this.own(() => clear.removeEventListener('click', onClear));
    this.effect(() => clear.hidden = this.value.value === '');
    editor.append(input, clear);
    return editor;
  }

  private _passwordEditor(): HTMLElement {
    const editor = document.createElement('div');
    editor.className = 'u2-input-password';
    const eye = document.createElement('button');
    eye.type = 'button';
    eye.className = 'u2-input-eye';
    eye.setAttribute('aria-label', 'Show password');
    const glyph = icon('eye');
    eye.append(glyph);
    const onClick = () => {
      const masked = this._input.type === 'password';
      this._input.type = masked ? 'text' : 'password';
      glyph.classList.toggle('fa-eye', !masked);
      glyph.classList.toggle('fa-eye-slash', masked);
      eye.setAttribute('aria-label', masked ? 'Hide password' : 'Show password');
    };
    eye.addEventListener('click', onClick);
    this.own(() => eye.removeEventListener('click', onClick));
    editor.append(this._input, eye);
    return editor;
  }

  /** Off-screen span measuring the text the field carries — the platform's own trick
   * (`text_input.dart:218-239`); an empty field is measured by its placeholder. */
  private _autoResize(): void {
    const measure = document.createElement('span');
    measure.className = 'u2-input-measure';
    measure.setAttribute('aria-hidden', 'true');
    this.root.append(measure);
    this.effect(() => {
      const text = this.value.value;
      measure.textContent = text === '' ? this._input.placeholder : text;
      const min = this.options.minWidth ?? 100;
      const max = this.options.maxWidth ?? 300;
      this._input.style.width = `${Math.min(Math.max(measure.offsetWidth + 10, min), max)}px`;
    });
  }
}

export interface TextAreaOptions extends InputOptions<string> {
  placeholder?: string;
  autoGrow?: boolean;
}

export class TextArea extends Input<string, TextAreaOptions> {
  private _area!: HTMLTextAreaElement;

  constructor(options: TextAreaOptions = {}) {
    super(options, '');
    this.root.dataset.u2 = 'text-area';
  }

  protected createEditor(): HTMLElement {
    const area = document.createElement('textarea');
    this._area = area;
    area.placeholder = this.options.placeholder ?? '';
    bindValue(this.scope, area, this.value, this.options.commitOn);
    if (this.options.autoGrow)
      this.effect(() => this._grow(this.value.value));
    return area;
  }

  // scrollHeight only reads back while laid out; a detached textarea keeps its CSS height.
  private _grow(text: string): void {
    const area = this._area;
    if (!area.isConnected)
      return;
    area.style.height = 'auto';
    area.style.height = text ? `${area.scrollHeight}px` : '';
  }
}
