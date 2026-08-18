import {Input, InputOptions} from '../core/input-base.js';
import {bindValue} from '../core/bind.js';

export interface TextInputOptions extends InputOptions<string> {
  /** Magnifier affordance plus a clear button that shows once there is something to clear. */
  search?: boolean;
  placeholder?: string;
}

export class TextInput extends Input<string, TextInputOptions> {
  private _input!: HTMLInputElement;

  constructor(options: TextInputOptions = {}) {
    super(options, '');
    this.root.dataset.u2 = 'text-input';
  }

  protected createEditor(): HTMLElement {
    const input = document.createElement('input');
    this._input = input;
    input.type = 'text';
    input.autocomplete = 'off';
    input.placeholder = this.options.placeholder ?? '';
    bindValue(this.scope, input, this.value);
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
    bindValue(this.scope, area, this.value);
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
