import {Input, InputOptions} from '../../core/input-base.js';
import {button} from '../../core/elements.js';
import {icon} from '../display/icon.js';

/** Splits on commas outside double quotes — the same lookahead Dart's `ListInput` uses. */
const SPLIT = /,(?=(?:[^"]*"[^"]*")*(?![^"]*"))/;

export interface ListInputOptions extends InputOptions<string[]> {
  placeholder?: string;
}

/** Comma-separated list editor, the port of Dart's `ListInput`: the same quote-aware split, and
 * the same expand toggle — on the options rail, where the platform puts it
 * (`list_input.dart:40-53`) — that swaps the one-line field for a textarea. Both fields exist all
 * along and carry the same text, so expanding never loses what is being typed. */
export class ListInput extends Input<string[], ListInputOptions> {
  private _field!: HTMLInputElement;
  private _area!: HTMLTextAreaElement;
  private _toggle!: HTMLElement;
  private _glyph!: HTMLElement;
  private _expanded!: boolean;

  constructor(options: ListInputOptions = {}) {
    super(options, []);
    this.root.dataset.u2 = 'list-input';
  }

  /** `a,"b,c",d` → `['a', 'b,c', 'd']`; blank text is the empty list. */
  static parse(text: string): string[] {
    const line = text.replace(/[\r\n]/g, '');
    if (line === '')
      return [];
    return line.split(SPLIT).map((part) => {
      const item = part.trim();
      return item.startsWith('"') ? item.substring(1, item.length - 1) : item;
    });
  }

  /** Quotes the items that carry a comma, so {@link parse} gets them back whole. */
  static join(items: string[]): string {
    return items.map((item) => `${item}`.includes(',') ? `"${item}"` : `${item}`).join(',');
  }

  get expanded(): boolean {
    return this._expanded;
  }

  set expanded(x: boolean) {
    this._expanded = x;
    this._field.hidden = x;
    this._area.hidden = !x;
    this._glyph.classList.toggle('fa-expand', !x);
    this._glyph.classList.toggle('fa-compress', x);
    this._toggle.setAttribute('aria-label', x ? 'Collapse' : 'Expand');
  }

  protected createEditor(): HTMLElement {
    this._expanded = false;
    const placeholder = this.options.placeholder ?? 'Comma-separated values';
    const field = document.createElement('input');
    this._field = field;
    field.type = 'text';
    field.autocomplete = 'off';
    field.className = 'u2-list-text';
    field.placeholder = placeholder;

    const area = document.createElement('textarea');
    this._area = area;
    area.className = 'u2-list-text u2-list-area';
    area.placeholder = placeholder;
    area.hidden = true;

    let echo = false;
    this.effect(() => {
      const value = this.value.value;
      if (echo)
        return;
      const text = ListInput.join(value);
      field.value = text;
      area.value = text;
    });
    const onInput = (source: HTMLInputElement | HTMLTextAreaElement) => {
      const text = source.value;
      echo = true;
      try {
        this.value.value = ListInput.parse(text);
      } finally {
        echo = false;
      }
      (source === field ? area : field).value = text;
    };
    this._listen(field, 'input', () => onInput(field));
    this._listen(area, 'input', () => onInput(area));

    this._glyph = icon('expand');
    this._toggle = button(this._glyph, () => this.expanded = !this._expanded);
    this._toggle.classList.add('u2-icon-btn', 'u2-list-toggle');
    this._toggle.setAttribute('aria-label', 'Expand');

    const editor = document.createElement('div');
    editor.className = 'u2-list';
    editor.append(field, area);
    this.addOptions(this._toggle);
    return editor;
  }

  private _listen(el: HTMLElement, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }
}
