import {Input, InputOptions} from '../../core/input-base.js';
import {Scope} from '../../core/scope.js';
import {Overlay, OVERLAY_CLOSE_EVENT} from '../../core/overlay.js';
import {div, divV, span} from '../../core/elements.js';
import {icon} from '../display/icon.js';
import {ICON_NAMES, BRAND_ICON_NAMES} from '../display/icon-names.js';
import {VirtualGrid} from '../collections/grid.js';
import {TextInput} from './text-input.js';

export interface IconInputOptions extends InputOptions<string> {
  /** The icons on offer; every Font Awesome name the platform ships by default. */
  names?: readonly string[];
}

const ALL = [...ICON_NAMES, ...BRAND_ICON_NAMES];
const GRID_KEYS = new Set(['ArrowLeft', 'ArrowRight', 'ArrowUp', 'ArrowDown', 'PageUp', 'PageDown', 'Enter']);

/** A Font Awesome icon name, picked from a searchable grid in an anchored popup. A pick commits
 * and closes; Esc, an outside click or the anchor leaving the DOM close without a write. */
export class IconInput extends Input<string, IconInputOptions> {
  private _control!: HTMLElement;
  private _popup?: Scope;

  constructor(options: IconInputOptions = {}) {
    super(options, '');
    this.root.dataset.u2 = 'icon-input';
  }

  protected createEditor(): HTMLElement {
    const preview = span('', 'u2-icon-input-preview');
    const name = span('', 'u2-icon-input-name');
    const control = this._control = div([preview, name], 'u2-icon-input');
    control.tabIndex = 0;
    for (const [k, v] of [['role', 'button'], ['aria-haspopup', 'listbox'], ['aria-expanded', 'false']])
      control.setAttribute(k, v);
    this.effect(() => {
      const value = this.value.value;
      preview.replaceChildren(...(value ? [icon(value)] : []));
      name.textContent = value || '—';
    });
    control.addEventListener('click', () => this._toggle());
    control.addEventListener('keydown', (e) => {
      const clear = (e.key === 'Backspace' || e.key === 'Delete') && this.nullable;
      if (clear)
        this.value.value = '';
      else if (e.key === 'Enter' || e.key === ' ')
        this._toggle();
      else
        return;
      e.preventDefault();
    });
    this.own(() => this._close());
    return control;
  }

  private _toggle(): void {
    if (!this.enabled) // the editor is a div, so the base's disabled sweep never reaches it
      return;
    if (this._popup)
      return this._close();
    const scope = this._popup = new Scope();
    Scope.runWith(scope, () => {
      const names = this.options.names ?? ALL;
      const search = new TextInput({placeholder: 'Search…', search: true, inline: true});
      const grid = new VirtualGrid<string>({cellWidth: 28, cellHeight: 28,
        render: (n, _i, cell) => { cell.title = n; return icon(n); },
        onActivate: (n) => { this.value.value = n; this._close(); }});
      scope.effect(() => {
        const q = search.value.value.trim().toLowerCase().replace(/\s+/g, '-');
        grid.setItems(q ? names.filter((n) => n.includes(q)) : names);
      });
      const popup = divV([search, grid], 'u2-icon-input-popup');
      popup.addEventListener('keydown', (e) => {
        if (GRID_KEYS.has(e.key) && !grid.root.contains(e.target as Node))
          grid.onKeyDown(e);
      });
      popup.addEventListener(OVERLAY_CLOSE_EVENT, () => this._close());
      Overlay.show(this._control, popup, scope);
      this._control.setAttribute('aria-expanded', 'true');
      search.root.querySelector<HTMLElement>('input')?.focus();
    });
  }

  private _close(): void {
    this._popup?.dispose();
    this._popup = undefined;
    this._control.setAttribute('aria-expanded', 'false');
  }
}
