/* Access grid: who may do what, as a small table of principals × capabilities — one checkbox per
   cell, a remove button per row, a picker at the bottom that adds a principal. Rows the caller
   marks inherited sit above the editable ones, read-only, with a note of where they come from.
   The value is the editable rows alone; a capability the caller locks (Edit and Delete without a
   writable binding) stays visible but cannot be granted. */
import {Input, InputOptions, LiveOption} from '../../core/input-base.js';
import {span} from '../../core/elements.js';
import type {ChoiceItem} from '../inputs/choice-input.js';

export interface AccessRow {
  principal: string;
  can: Record<string, boolean>;
}

/** A row shown for information: granted elsewhere (`from`), never edited here. */
export interface InheritedAccessRow extends AccessRow {
  from: string;
}

export type AccessCapabilityItem = string | {name: string, label: string};

export interface AccessGridOptions extends InputOptions<AccessRow[]> {
  /** The columns, in order: a name, or a name with the header it gets. */
  capabilities: AccessCapabilityItem[];
  /** What the add picker offers; a principal already in a row is left out of it. */
  principals: ChoiceItem[];
  /** The first column's header ('Group' by default). */
  principalLabel?: string;
  /** Rows granted elsewhere, shown muted above the editable ones. */
  inherited?: LiveOption<InheritedAccessRow[]>;
  /** Capabilities that cannot be granted now; their boxes stay, disabled. */
  locked?: LiveOption<string[]>;
  /** The picker's placeholder ('+ add group or user…' by default). */
  addText?: string;
}

function itemValue(item: ChoiceItem): string {
  return typeof item === 'string' ? item : item.value;
}

function itemLabel(item: ChoiceItem): string {
  return typeof item === 'string' ? item : item.label;
}

export class AccessGrid extends Input<AccessRow[], AccessGridOptions> {
  private _body!: HTMLElement;
  private _picker!: HTMLSelectElement;
  private _capabilities!: {name: string, label: string}[];
  private _principals!: ChoiceItem[];
  private _inherited!: InheritedAccessRow[];
  private _locked!: Set<string>;

  constructor(options: AccessGridOptions) {
    super(options, []);
    this.root.dataset.u2 = 'access-grid';
    this.liveOption('inherited', options.inherited, (rows) => {
      this._inherited = rows;
      this._render();
    });
    this.liveOption('locked', options.locked, (names) => {
      this._locked = new Set(names);
      this._render();
    });
  }

  protected createEditor(): HTMLElement {
    this._capabilities = this.options.capabilities.map((c) =>
      typeof c === 'string' ? {name: c, label: c.charAt(0).toUpperCase() + c.slice(1)} : c);
    this._principals = this.options.principals;
    const inherited = this.options.inherited;
    const locked = this.options.locked;
    this._inherited = Array.isArray(inherited) ? inherited : [];
    this._locked = new Set(Array.isArray(locked) ? locked : []);
    const table = document.createElement('table');
    table.className = 'u2-access-grid';
    const head = document.createElement('tr');
    head.append(AccessGrid._cell('th', this.options.principalLabel ?? 'Group'));
    for (const c of this._capabilities)
      head.append(AccessGrid._cell('th', c.label, 'u2-access-grid-cap'));
    head.append(AccessGrid._cell('th', ''));
    const thead = document.createElement('thead');
    thead.append(head);
    this._body = document.createElement('tbody');
    const cell = AccessGrid._cell('td', '');
    cell.colSpan = this._capabilities.length + 2;
    const foot = document.createElement('tr');
    foot.append(cell);
    const tfoot = document.createElement('tfoot');
    tfoot.append(foot);
    table.append(thead, this._body, tfoot);
    this._picker = document.createElement('select');
    this._picker.className = 'u2-access-grid-add';
    this._picker.setAttribute('aria-label', this.options.addText ?? '+ add group or user…');
    cell.append(this._picker);
    this._listen(this._picker, 'change', () => this._add(this._picker.value));
    this._listen(this._body, 'change', (e) => this._onChange(e.target as HTMLInputElement));
    this._listen(this._body, 'click', (e) => this._onClick(e.target as HTMLElement));
    this.effect(() => {
      this.value.value;
      this._render();
    });
    return table;
  }

  private _render(): void {
    this._body.textContent = '';
    for (const row of this._inherited)
      this._body.append(this._row(row, row.from));
    for (const row of this.value.peek())
      this._body.append(this._row(row, null));
    this._fillPicker();
    this.refreshEnabled();
  }

  private _row(row: AccessRow, from: string | null): HTMLElement {
    const tr = document.createElement('tr');
    tr.className = from === null ? 'u2-access-grid-row' : 'u2-access-grid-row u2-access-grid-inherited';
    tr.dataset.principal = row.principal;
    const who = AccessGrid._cell('td', '');
    tr.append(who);
    who.append(span(row.principal, 'u2-access-grid-principal'));
    if (from !== null && from !== '')
      who.append(span(from, 'u2-access-grid-from'));
    for (const c of this._capabilities) {
      const box = document.createElement('input');
      box.type = 'checkbox';
      box.className = 'u2-input-checkbox u2-access-grid-check';
      box.checked = row.can[c.name] === true;
      box.disabled = from !== null || this._locked.has(c.name);
      box.dataset.capability = c.name;
      box.setAttribute('aria-label', `${c.label} for ${row.principal}`);
      const cell = AccessGrid._cell('td', '', 'u2-access-grid-cap');
      cell.append(box);
      tr.append(cell);
    }
    const last = AccessGrid._cell('td', '');
    tr.append(last);
    if (from === null) {
      const remove = document.createElement('button');
      remove.type = 'button';
      remove.className = 'u2-access-grid-remove';
      remove.textContent = '✕';
      remove.setAttribute('aria-label', `Remove ${row.principal}`);
      last.append(remove);
    }
    return tr;
  }

  private _fillPicker(): void {
    const taken = new Set([...this.value.peek(), ...this._inherited].map((r) => r.principal));
    const picker = this._picker;
    picker.textContent = '';
    picker.append(new Option(this.options.addText ?? '+ add group or user…', ''));
    for (const item of this._principals) {
      if (!taken.has(itemValue(item)))
        picker.append(new Option(itemLabel(item), itemValue(item)));
    }
    picker.value = '';
  }

  private _add(principal: string): void {
    if (principal === '' || this.value.peek().some((r) => r.principal === principal))
      return;
    const can: Record<string, boolean> = {};
    for (const c of this._capabilities)
      can[c.name] = c === this._capabilities[0];
    this.value.value = [...this.value.peek(), {principal, can}];
  }

  private _onChange(target: HTMLInputElement): void {
    const capability = target.dataset.capability;
    const principal = (target.closest('tr') as HTMLElement | null)?.dataset.principal;
    if (capability === undefined || principal === undefined)
      return;
    this.value.value = this.value.peek().map((r) =>
      r.principal === principal ? {...r, can: {...r.can, [capability]: target.checked}} : r);
  }

  private _onClick(target: HTMLElement): void {
    if (!target.closest('.u2-access-grid-remove'))
      return;
    const principal = (target.closest('tr') as HTMLElement | null)?.dataset.principal;
    this.value.value = this.value.peek().filter((r) => r.principal !== principal);
  }

  private _listen(el: HTMLElement, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }

  private static _cell(tag: 'th' | 'td', text: string, cls?: string): HTMLTableCellElement {
    const el = document.createElement(tag);
    el.textContent = text;
    if (cls !== undefined)
      el.className = cls;
    return el;
  }
}
