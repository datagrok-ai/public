/* A saved run of a function, picked from its history in an anchored popup. The value is always
   the FuncCall id string, never the object (the pickers.ts convention); `functionName` dictates
   whose runs the popup lists. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {signal, Signal, ReadonlySignal} from '../../core/signals.js';
import {Input, InputOptions, LiveOption} from '../../core/input-base.js';
import {Scope} from '../../core/scope.js';
import {div, span} from '../../core/elements.js';
import {FuncCallHistoryBrowser} from '../entities/func-call-history-browser.js';
import {Editors} from '../forms/editors.js';

export interface FuncCallInputOptions extends InputOptions<string> {
  /** Whose runs the popup lists. Live: a change while open re-filters in place. */
  functionName?: LiveOption<string>;
  placeholder?: string;
}

/** The saved-run label: function name plus the local start time. */
function callLabel(call: DG.FuncCall): string {
  const name = call.func?.friendlyName || call.func?.name || 'Function call';
  const started = call.started?.isValid() ? call.started.format('MMM D, YYYY HH:mm') : '';
  return started === '' ? name : `${name} — ${started}`;
}

/** A FuncCall id, picked from the FuncCallHistoryBrowser in an anchored popup. A row click,
 * Enter or double-click commits and closes; Esc, an outside click or the anchor leaving the DOM
 * close without a write. The editor shows the picked run's label where the call is known (from
 * the pick, or a background `find` by id), the raw id until then. */
export class FuncCallInput extends Input<string, FuncCallInputOptions> {
  readonly functionName: Signal<string>;

  private _control!: HTMLElement;
  private _popup?: Scope;
  // assigned in createEditor — the base constructor runs it before subclass field initializers
  private _picked!: Signal<DG.FuncCall | undefined>;
  private _resolveGen!: number;

  constructor(options: FuncCallInputOptions = {}) {
    super(options, '');
    this.root.dataset.u2 = 'func-call-input';
    const value = options.functionName;
    this.functionName = value instanceof Signal ?
      this._followName(value) : signal((value as string | undefined) ?? '');
  }

  protected createEditor(): HTMLElement {
    this._picked = signal<DG.FuncCall | undefined>(undefined);
    this._resolveGen = 0;
    const name = span('', 'u2-func-call-input-name');
    const control = this._control = div([name], 'u2-func-call-input');
    control.tabIndex = 0;
    for (const [k, v] of [['role', 'button'], ['aria-haspopup', 'listbox'], ['aria-expanded', 'false']])
      control.setAttribute(k, v);
    this.effect(() => {
      const value = this.value.value;
      const picked = this._picked.value;
      if (value !== '' && picked?.id !== value)
        this._resolve(value);
      name.textContent = value === '' ? (this.options.placeholder ?? 'Pick a run…') :
        picked?.id === value ? callLabel(picked) : value;
      name.classList.toggle('u2-func-call-input-empty', value === '');
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
    const scope = this._popup = FuncCallHistoryBrowser.popup(this._control,
      'u2-func-call-input-popup',
      {functionName: this.functionName}, // the same signal — the popup follows retargets live
      (call) => {
        this._picked.value = call;
        this.value.value = call.id;
      });
    scope.own(() => {
      this._popup = undefined;
      this._control.setAttribute('aria-expanded', 'false');
    });
    this._control.setAttribute('aria-expanded', 'true');
  }

  private _close(): void {
    this._popup?.dispose();
  }

  private _resolve(id: string): void {
    const gen = ++this._resolveGen;
    void grok.dapi.functions.calls.allPackageVersions()
      .include('func.package,func.params,inputs').find(id)
      .then((call) => {
        if (call != null && gen === this._resolveGen && this.value.peek() === id)
          this._picked.value = call;
      })
      .catch(() => {}); // an unknown id keeps showing raw — the value stands
  }

  private _followName(value: ReadonlySignal<string>): Signal<string> {
    const own = signal((value.peek() as string | undefined) ?? '');
    this.effect(() => {
      const x = value.value;
      if (x !== undefined)
        own.value = x; // a forward-ref proxy starts undefined; skip until linked
    });
    return own;
  }
}

export function funcCallInput(label: string, options: FuncCallInputOptions = {}): FuncCallInput {
  return new FuncCallInput({...options, label});
}

Editors.register({
  match: (prop) => prop.semType === 'FuncCall',
  create: (prop, options) => new FuncCallInput({...options,
    functionName: ((prop as {options?: Record<string, any>}).options?.['funcName'] as string | undefined) ?? ''}),
});
