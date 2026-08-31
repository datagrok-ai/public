/* A function picked by its namespace-qualified name. The input stays one row — the platform icon
   and the friendly name of whatever the value resolves to — and the whole FunctionsBrowser opens
   in an anchored popup for picking. The value is always the nqName string, never the `DG.Func`
   (the pickers.ts convention), which is what the `FunctionName` semantic type stores. */
import * as DG from 'datagrok-api/dg';
import {Input, InputOptions} from '../../core/input-base.js';
import {Scope} from '../../core/scope.js';
import {Overlay, OVERLAY_CLOSE_EVENT} from '../../core/overlay.js';
import {div, span} from '../../core/elements.js';
import {functionsBrowser, funcIcon} from '../entities/functions-browser.js';
import {Editors} from '../forms/editors.js';

export interface FunctionInputOptions extends InputOptions<string> {
  filter?: (f: DG.Func) => boolean;
  /** Only functions with one scalar (or dynamic) output. */
  scalarOnly?: boolean;
  /** Packages whose functions are left out of the popup. */
  ignorePackages?: string[];
  placeholder?: string;
}

/** A namespace-qualified function name (`Chem:SmilesToMw`), picked from the FunctionsBrowser in an
 * anchored popup. A row click, Enter or double-click commits and closes; Esc, an outside click or
 * the anchor leaving the DOM close without a write. The editor renders the platform icon and
 * friendly name where the value resolves in the client registry, the raw string where it doesn't. */
export class FunctionInput extends Input<string, FunctionInputOptions> {
  private _control!: HTMLElement;
  private _popup?: Scope;

  constructor(options: FunctionInputOptions = {}) {
    super(options, '');
    this.root.dataset.u2 = 'function-input';
  }

  protected createEditor(): HTMLElement {
    const icon = span('', 'u2-function-input-icon');
    const name = span('', 'u2-function-input-name');
    const control = this._control = div([icon, name], 'u2-function-input');
    control.tabIndex = 0;
    for (const [k, v] of [['role', 'button'], ['aria-haspopup', 'listbox'], ['aria-expanded', 'false']])
      control.setAttribute(k, v);
    this.effect(() => {
      const value = this.value.value;
      const func = value === '' ? undefined :
        DG.Func.find({}).find((f) => f.nqName === value);
      const glyph = func ? funcIcon(func) : null;
      icon.replaceChildren(...(glyph ? [glyph] : []));
      icon.hidden = glyph === null;
      name.textContent = value === '' ? (this.options.placeholder ?? 'Pick a function…') :
        (func?.friendlyName || value);
      name.classList.toggle('u2-function-input-empty', value === '');
      name.title = value; // not the editor's own title — the base owns that for validity messages
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
      const commit = (picked: string): void => {
        this.value.value = picked;
        this._close();
      };
      const browser = functionsBrowser({
        showTags: false, showRunButton: false, showSignature: true,
        setCurrentObject: false,
        filter: this.options.filter,
        scalarOnly: this.options.scalarOnly,
        ignorePackages: this.options.ignorePackages,
        onActivate: (item) => commit(item.name),
      });
      const popup = div([browser.root], 'u2-function-input-popup');
      // one click picks: the row carries the qualified name it renders, so the commit never
      // depends on the browser's selection effects having flushed
      popup.addEventListener('click', (e) => {
        const picked = (e.target as Element | null)?.closest('[data-u2-func]')
          ?.getAttribute('data-u2-func');
        if (picked)
          commit(picked);
      });
      popup.addEventListener(OVERLAY_CLOSE_EVENT, () => this._close());
      Overlay.show(this._control, popup, scope);
      this._control.setAttribute('aria-expanded', 'true');
      browser.root.querySelector<HTMLElement>('input')?.focus();
    });
  }

  private _close(): void {
    this._popup?.dispose();
    this._popup = undefined;
    this._control.setAttribute('aria-expanded', 'false');
  }
}

export function functionInput(label: string, options: FunctionInputOptions = {}): FunctionInput {
  return new FunctionInput({...options, label});
}

Editors.register({
  match: (prop) => prop.semType === 'FunctionName',
  create: (prop, options) => new FunctionInput(options),
});
