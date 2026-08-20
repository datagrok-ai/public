import {Input, InputOptions} from '../core/input-base.js';
import {signal, Signal} from '../core/signals.js';
import {icon} from './icon.js';

export type MapInputOptions = InputOptions<Record<string, string>>;

/** Key/value rows with per-row add and remove, the port of Dart's `MapInput`. Two deliberate
 * divergences: duplicate keys report through the standard validity signal instead of Dart's
 * silent `validate()` override, and a row whose key is empty never reaches the value. */
export class MapInput extends Input<Record<string, string>, MapInputOptions> {
  private _rows!: [string, string][];
  private _body!: HTMLElement;
  private _dup!: Signal<boolean>;
  private _echo!: boolean;

  constructor(options: MapInputOptions = {}) {
    super(options, {});
    this.root.dataset.u2 = 'map-input';
  }

  protected createEditor(): HTMLElement {
    this._dup = signal(false);
    this._echo = false;
    const body = document.createElement('div');
    this._body = body;
    body.className = 'u2-map';
    this.effect(() => {
      const value = this.value.value;
      if (!this._echo)
        this._load(value);
    });
    this._listen(body, 'input', (e) => this._onInput(e.target as HTMLElement));
    this._listen(body, 'click', (e) => this._onClick(e.target as HTMLElement));
    // registered here, before the base builds its validation effect, so the effect tracks `_dup`
    this.addValidator(() => this._dup.value ? 'Duplicate keys are not allowed' : null);
    return body;
  }

  private _load(value: Record<string, string>): void {
    this._rows = Object.keys(value).map((key): [string, string] => [key, value[key]]);
    if (this._rows.length === 0)
      this._rows.push(['', '']);
    this._render();
  }

  private _render(): void {
    this._body.textContent = '';
    for (let i = 0; i < this._rows.length; i++)
      this._body.append(this._row(i));
    this._mark();
    this.refreshEnabled();
  }

  private _row(index: number): HTMLElement {
    const row = document.createElement('div');
    row.className = 'u2-map-row';
    const actions = document.createElement('div');
    actions.className = 'u2-map-actions';
    actions.append(this._button(index, 'add', 'plus', 'Add row'),
      this._button(index, 'remove', 'minus', 'Remove row'));
    row.append(this._cell(index, 0, 'Key'), this._cell(index, 1, 'Value'), actions);
    return row;
  }

  private _cell(index: number, part: 0 | 1, placeholder: string): HTMLInputElement {
    const cell = document.createElement('input');
    cell.type = 'text';
    cell.autocomplete = 'off';
    cell.className = part === 0 ? 'u2-map-key' : 'u2-map-value';
    cell.placeholder = placeholder;
    cell.value = this._rows[index][part];
    cell.dataset.row = String(index);
    cell.dataset.part = String(part);
    return cell;
  }

  private _button(index: number, action: string, glyph: string, label: string): HTMLElement {
    const el = document.createElement('button');
    el.type = 'button';
    el.className = 'u2-map-btn';
    el.append(icon(glyph));
    el.dataset.row = String(index);
    el.dataset.act = action;
    el.setAttribute('aria-label', label);
    el.title = label;
    return el;
  }

  private _onInput(target: HTMLElement): void {
    const part = target.dataset.part;
    if (part === undefined)
      return;
    this._rows[Number(target.dataset.row)][Number(part)] = (target as HTMLInputElement).value;
    this._commit();
  }

  private _onClick(target: HTMLElement): void {
    const el = target.closest('.u2-map-btn') as HTMLElement | null;
    if (!el)
      return;
    const index = Number(el.dataset.row);
    if (el.dataset.act === 'add')
      this._rows.splice(index + 1, 0, ['', '']);
    else
      this._rows.splice(index, 1);
    if (this._rows.length === 0)
      this._rows.push(['', '']);
    this._render();
    this._commit();
  }

  /** Rewrites the value from the rows, and only when it actually changed — adding or clearing an
   * empty row must not read as an edit to whoever listens on the signal. */
  private _commit(): void {
    const value: Record<string, string> = {};
    for (const [key, item] of this._rows) {
      if (key !== '')
        value[key] = item;
    }
    this._mark();
    if (MapInput._same(this.value.peek(), value))
      return;
    this._echo = true;
    try {
      this.value.value = value;
    } finally {
      this._echo = false;
    }
  }

  private _mark(): void {
    const seen = new Set<string>();
    const dups = new Set<string>();
    for (const [key] of this._rows) {
      if (key === '')
        continue;
      if (seen.has(key))
        dups.add(key);
      else
        seen.add(key);
    }
    const dup = dups.size > 0;
    // _render reaches here from the value effect, which must not be restarted by a write of its own
    if (this._dup.peek() !== dup)
      this._dup.value = dup;
    const cells = this._body.querySelectorAll('.u2-map-key');
    for (let i = 0; i < cells.length; i++)
      cells[i].classList.toggle('u2-invalid', dups.has(this._rows[i][0]));
  }

  private static _same(a: Record<string, string>, b: Record<string, string>): boolean {
    const keys = Object.keys(a);
    if (keys.length !== Object.keys(b).length)
      return false;
    for (const key of keys) {
      if (a[key] !== b[key])
        return false;
    }
    return true;
  }

  private _listen(el: HTMLElement, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }
}
