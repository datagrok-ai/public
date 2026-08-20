/* Font editor over the platform's strict font string (`font_input.dart`): size box with a preset
   popup, B / I toggles and a family select, all composing back into
   `"<weight> <style> <size>px \"<family>\""` — generic families stay lowercase and unquoted, as
   the platform writes them. Text that does not parse leaves the controls at their defaults, which
   is what the platform does too; the value is normalized on the next edit. */
import {Input, InputOptions} from '../core/input-base.js';
import {Overlay, OVERLAY_CLOSE_EVENT} from '../core/overlay.js';
import {icon} from './icon.js';

/** Written unquoted and lowercase by {@link FontInput.compose} (`font_input.dart:41-42`). */
const GENERIC = ['monospace', 'serif', 'sans-serif', 'cursive', 'fantasy', 'system-ui'];
/** The subset the strict format accepts without quotes (`font_input.dart:40`). */
const PARSED_GENERIC = ['monospace', 'serif', 'cursive'];
const SIZES = [8, 9, 10, 11, 12, 14, 18, 24, 36, 48, 72, 96];
const FAMILIES = ['Monospace', 'Arial', 'Open Sans', 'Roboto'];
const DEFAULT_SIZE = 10;
const DEFAULT_FAMILY = 'Roboto';
const INT = /^-?\d+$/;

export interface Font {
  bold: boolean;
  italic: boolean;
  size: number;
  family: string;
}

export interface FontInputOptions extends InputOptions<string> {
  families?: string[];
  /** Sizes offered by the popup under the size box. */
  sizes?: number[];
}

export class FontInput extends Input<string, FontInputOptions> {
  private static _seq = 0;

  private _id!: string;
  private _size!: HTMLInputElement;
  private _family!: HTMLSelectElement;
  private _bold!: HTMLElement;
  private _italic!: HTMLElement;
  private _popup!: HTMLElement;
  private _font!: Font;
  private _closeOverlay: (() => void) | undefined;

  constructor(options: FontInputOptions = {}) {
    super(options, FontInput.compose({bold: false, italic: false, size: DEFAULT_SIZE,
      family: DEFAULT_FAMILY}));
    this.root.dataset.u2 = 'font-input';
  }

  /** The strict format, or null when the string is not one (`font_input.dart:45-88`). */
  static parse(text: string): Font | null {
    if (!text)
      return null;
    const trimmed = text.trim();
    const first = trimmed.indexOf('"');
    const last = trimmed.lastIndexOf('"');
    let parts: string[];
    let family: string;
    if (first >= 0 && last > first) {
      family = trimmed.substring(first + 1, last);
      parts = trimmed.substring(0, first).trim().split(' ').filter((p) => p);
      if (parts.length !== 3)
        return null;
    } else {
      parts = trimmed.split(' ').filter((p) => p);
      if (parts.length !== 4)
        return null;
      family = parts[3].toLowerCase();
      if (!PARSED_GENERIC.includes(family))
        return null;
    }
    const weight = parts[0].toLowerCase();
    const style = parts[1].toLowerCase();
    const size = parts[2].toLowerCase();
    if ((weight !== 'bold' && weight !== 'normal') || (style !== 'italic' && style !== 'normal'))
      return null;
    if (!size.endsWith('px'))
      return null;
    const px = size.substring(0, size.length - 2);
    if (!INT.test(px))
      return null;
    return {bold: weight === 'bold', italic: style === 'italic', size: Number(px), family};
  }

  static compose(font: Font): string {
    const weight = font.bold ? 'bold' : 'normal';
    const style = font.italic ? 'italic' : 'normal';
    const family = font.family === '' ? DEFAULT_FAMILY : font.family;
    const size = `${font.size}px`;
    return GENERIC.includes(family.toLowerCase()) ?
      `${weight} ${style} ${size} ${family.toLowerCase()}` :
      `${weight} ${style} ${size} "${family}"`;
  }

  protected createEditor(): HTMLElement {
    this._id = `u2-font-${++FontInput._seq}`;
    this._font = {bold: false, italic: false, size: DEFAULT_SIZE, family: DEFAULT_FAMILY};

    const size = document.createElement('input');
    this._size = size;
    size.type = 'text';
    size.autocomplete = 'off';
    size.className = 'u2-font-size';
    size.setAttribute('aria-label', 'Font size');
    size.setAttribute('aria-controls', `${this._id}-sizes`);

    const minus = this._step(-1, 'Decrease font size');
    const plus = this._step(1, 'Increase font size');
    this._bold = this._toggle('B', 'Bold', 'u2-font-bold');
    this._italic = this._toggle('I', 'Italic', 'u2-font-italic');

    const family = document.createElement('select');
    this._family = family;
    family.className = 'u2-font-family';
    family.setAttribute('aria-label', 'Font family');
    for (const name of this.options.families ?? FAMILIES)
      family.append(new Option(name, name));

    this._popup = document.createElement('div');
    this._popup.className = 'u2-font-sizes';
    this._popup.id = `${this._id}-sizes`;
    this._popup.setAttribute('role', 'listbox');
    for (const value of this.options.sizes ?? SIZES) {
      const row = document.createElement('div');
      row.className = 'u2-font-size-option';
      row.dataset.size = String(value);
      row.setAttribute('role', 'option');
      row.textContent = String(value);
      this._popup.append(row);
    }

    let echo = false;
    this.effect(() => {
      const value = this.value.value;
      if (!echo)
        this._apply(value);
    });
    const commit = () => {
      echo = true;
      try {
        this.value.value = FontInput.compose(this._font);
      } finally {
        echo = false;
      }
    };
    this._listen(size, 'input', () => {
      const typed = size.value.trim();
      this._font.size = INT.test(typed) ? Number(typed) : DEFAULT_SIZE;
      commit();
    });
    // mousedown, not click: it beats the blur that closes the popup, so a click on an open box
    // toggles it shut instead of closing and reopening (`font_input.dart:230-234`)
    this._listen(size, 'mousedown', () => this._toggleSizes());
    this._listen(size, 'blur', () => this._closeSizes());
    this._listen(family, 'change', () => {
      this._font.family = family.value;
      commit();
    });
    const step = (delta: number) => {
      const next = this._font.size + delta;
      if (next < 0)
        return;
      this._font.size = next;
      this._size.value = String(next);
      commit();
    };
    this._listen(minus, 'click', () => step(-1));
    this._listen(plus, 'click', () => step(1));
    this._listen(this._bold, 'click', () => {
      this._font.bold = !this._font.bold;
      this._paint();
      commit();
    });
    this._listen(this._italic, 'click', () => {
      this._font.italic = !this._font.italic;
      this._paint();
      commit();
    });
    // pointerdown, not click: the size box blurs on click, and its blur closes the popup first
    this._listen(this._popup, 'pointerdown', (e) => {
      const row = (e.target as HTMLElement).closest('.u2-font-size-option') as HTMLElement | null;
      if (!row)
        return;
      e.preventDefault();
      this._font.size = Number(row.dataset.size);
      this._size.value = row.dataset.size as string;
      this._closeSizes();
      commit();
    });
    this._listen(this._popup, OVERLAY_CLOSE_EVENT, () => this._closeOverlay = undefined);

    const editor = document.createElement('div');
    editor.className = 'u2-font';
    editor.append(size, minus, plus, this._bold, this._italic, family);
    return editor;
  }

  /** The − / + pair the platform puts beside the size box (`font_input.dart:189-203`): ±1px, and
   * a step that would go negative is refused. */
  private _step(delta: number, label: string): HTMLElement {
    const el = document.createElement('button');
    el.type = 'button';
    el.className = 'u2-font-step';
    el.append(icon(delta < 0 ? 'minus' : 'plus'));
    el.setAttribute('aria-label', label);
    el.title = label;
    return el;
  }

  private _toggle(text: string, label: string, cls: string): HTMLElement {
    const el = document.createElement('button');
    el.type = 'button';
    el.className = `u2-font-toggle ${cls}`;
    el.textContent = text;
    el.setAttribute('aria-label', label);
    el.setAttribute('aria-pressed', 'false');
    return el;
  }

  private _apply(value: string): void {
    const parsed = FontInput.parse(value);
    this._font = {
      bold: parsed?.bold ?? false,
      italic: parsed?.italic ?? false,
      size: parsed?.size ?? DEFAULT_SIZE,
      family: parsed && parsed.family !== '' ? parsed.family : DEFAULT_FAMILY,
    };
    this._size.value = parsed === null ? '' : String(parsed.size);
    this._selectFamily(this._font.family);
    this._paint();
  }

  /** A family the select does not offer joins it as a hidden option, so an unknown family
   * survives the round-trip instead of silently becoming the first item (`font_input.dart:151`). */
  private _selectFamily(family: string): void {
    const options = Array.from(this._family.querySelectorAll('option')) as HTMLOptionElement[];
    if (!options.some((option) => option.value === family)) {
      const option = new Option(family, family);
      option.hidden = true;
      this._family.append(option);
    }
    this._family.value = family;
  }

  private _paint(): void {
    this._paintToggle(this._bold, this._font.bold);
    this._paintToggle(this._italic, this._font.italic);
  }

  private _paintToggle(el: HTMLElement, on: boolean): void {
    el.classList.toggle('u2-font-toggle-active', on);
    el.setAttribute('aria-pressed', String(on));
  }

  private _toggleSizes(): void {
    if (this._closeOverlay)
      this._closeSizes();
    else
      this._closeOverlay = Overlay.show(this._size, this._popup, this.scope);
  }

  private _closeSizes(): void {
    const close = this._closeOverlay;
    this._closeOverlay = undefined;
    if (close)
      close();
  }

  private _listen(el: EventTarget, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }
}
