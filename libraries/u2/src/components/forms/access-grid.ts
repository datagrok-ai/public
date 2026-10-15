/* Access grid: who may do what, as a small table of principals × capabilities — one checkbox per
   cell, a remove button per row, a picker at the bottom that adds a principal. Rows the caller
   marks inherited sit above the editable ones, read-only, with a note of where they come from.
   The value is the editable rows alone; a capability the caller locks (Edit and Delete without a
   writable binding) stays visible but cannot be granted. Rows are keyed by principal and patched
   in place, so the focused box keeps its focus across any change. */
import {batch, signal, Signal, untracked} from '../../core/signals.js';
import {Input, InputOptions, LiveOption} from '../../core/input-base.js';
import {span} from '../../core/elements.js';
import {each} from '../../core/region.js';
import type {ChoiceItem} from '../inputs/choice-input.js';
import type {IItemsState} from '../../core/widget-like.js';

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
  /** What the add picker offers and how a principal is labelled; a principal already in a row is left
   * out of the picker. Live. */
  principals: LiveOption<ChoiceItem[]>;
  /** The first column's header ('Group' by default). */
  principalLabel?: string;
  /** Rows granted elsewhere, shown muted above the editable ones. */
  inherited?: LiveOption<InheritedAccessRow[]>;
  /** Capabilities that cannot be granted now; their boxes stay, disabled. */
  locked?: LiveOption<string[]>;
  /** The add select's placeholder ('+ add group or user…' by default); ignored with `picker`. */
  addText?: string;
  /** The add control in place of the select — a look-up over a directory, say; every pick is
   * handed to `add`, and a principal picked that way joins the labels. */
  picker?: (add: (item: {value: string, label: string}) => void) => HTMLElement;
}

/** What a row shows: the region's item, rebuilt on every change and patched onto the kept row. */
interface GridRow {
  key: string;
  row: AccessRow;
  label: string;
  /** Where an inherited row comes from; null on an editable one. */
  from: string | null;
  enabled: boolean;
  locked: string[];
}

function itemValue(item: ChoiceItem): string {
  return typeof item === 'string' ? item : item.value;
}

function itemLabel(item: ChoiceItem): string {
  return typeof item === 'string' ? item : item.label;
}

/** Keys that move a closed native select's choice — Windows and Linux fire `change` for each at once. */
const NAV_KEYS = ['ArrowUp', 'ArrowDown', 'Home', 'End', 'PageUp', 'PageDown'];

export class AccessGrid extends Input<AccessRow[], AccessGridOptions> {
  private declare _body: HTMLElement;
  private declare _capabilities: {name: string, label: string}[];
  private declare _principals: Signal<ChoiceItem[]>;
  private declare _inherited: Signal<InheritedAccessRow[]>;
  private declare _locked: Signal<string[]>;
  /** Bumped after every `enabled` sweep, which re-enables every box: the rows re-apply their own. */
  private declare _sweeps: Signal<number>;
  /** A key is moving the add select's choice: the `change` it fires in its task waits for Enter. */
  private declare _keyNav: boolean;

  constructor(options: AccessGridOptions) {
    super(options, []);
    this.root.dataset.u2 = 'access-grid';
    this.liveOption('principals', options.principals, (items) => this.setPrincipals(items ?? []));
    this.liveOption('inherited', options.inherited, (rows) => this._inherited.value = rows);
    this.liveOption('locked', options.locked, (names) => this._locked.value = names);
  }

  get principals(): ChoiceItem[] {
    return this._principals.peek();
  }

  set principals(x: ChoiceItem[]) {
    this.setPrincipals(x);
  }

  /** The first column's header. */
  get principalLabel(): string {
    return this.options.principalLabel ?? 'Group';
  }

  set principalLabel(x: string) {
    this.options.principalLabel = x ?? undefined;
    this.root.querySelector('.u2-access-grid th')!.textContent = this.principalLabel;
  }

  /** The add picker's placeholder and accessible name. */
  get addText(): string {
    return this.options.addText ?? '+ add group or user…';
  }

  set addText(x: string) {
    this.options.addText = x ?? undefined;
    const select = this.root.querySelector<HTMLSelectElement>('.u2-access-grid-add');
    if (select === null)
      return;
    select.setAttribute('aria-label', this.addText);
    if (select.options.length > 0)
      select.options[0].text = this.addText;
  }

  /** Replaces what the add picker offers and the labels rows are shown by; the rows stay. */
  setPrincipals(items: ChoiceItem[]): void {
    this._principals.value = items;
  }

  get enabled(): boolean {
    return super.enabled;
  }

  set enabled(x: boolean) {
    super.enabled = x;
    this._sweeps.value = this._sweeps.peek() + 1;
  }

  /** One item per row, keyed by principal (an inherited row by `inherited:<principal>`), labelled as
   * shown; inherited rows are disabled. */
  getItems(opts?: {all?: boolean}): IItemsState {
    return Input.itemsState(untracked(() => this._rows()).map((r, index) =>
      ({key: r.key, label: r.label, index, rendered: true, ...r.from === null ? {} : {disabled: true}})), opts?.all);
  }

  /** The editable rows' principals, as shown. */
  protected displayText(): string {
    return untracked(() => this._rows()).filter((r) => r.from === null).map((r) => r.label).join(', ');
  }

  protected createEditor(): HTMLElement {
    const {principals, inherited, locked} = this.options;
    this._capabilities = this.options.capabilities.map((c) =>
      typeof c === 'string' ? {name: c, label: c.charAt(0).toUpperCase() + c.slice(1)} : c);
    this._principals = signal(Array.isArray(principals) ? principals : []);
    this._inherited = signal(Array.isArray(inherited) ? inherited : []);
    this._locked = signal(Array.isArray(locked) ? locked : []);
    this._sweeps = signal(0);
    this._keyNav = false;
    const table = document.createElement('table');
    table.className = 'u2-access-grid';
    table.setAttribute('role', 'grid');
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
    const picker = this.options.picker;
    if (picker !== undefined)
      cell.append(picker((item) => this._add(item.value, item.label)));
    else {
      const select = this._select();
      cell.append(select);
      this.effect(() => this._fillPicker(select));
    }
    this.scope.listen(this._body, 'change', (e) => this._onChange(e.target as HTMLInputElement));
    this.scope.listen(this._body, 'click', (e) => this._onClick(e.target as HTMLElement));
    each(this._body, () => this._rows(), (r) => r.key, (r) => this._row(r), {update: (tr, r) => this._patch(tr, r)});
    return table;
  }

  private _rows(): GridRow[] {
    this._sweeps.value;
    const enabled = this.enabled;
    const locked = this._locked.value;
    const principals = this._principals.value;
    const label = (principal: string) => {
      const item = principals.find((p) => itemValue(p) === principal);
      return item === undefined ? principal : itemLabel(item);
    };
    return [
      ...this._inherited.value.map((row) => ({key: `inherited:${row.principal}`, row, label: label(row.principal),
        from: row.from, enabled, locked})),
      ...(this.value.value ?? []).map((row) => ({key: row.principal, row, label: label(row.principal), from: null,
        enabled, locked})),
    ];
  }

  private _row(r: GridRow): HTMLElement {
    const tr = document.createElement('tr');
    tr.className = r.from === null ? 'u2-access-grid-row' : 'u2-access-grid-row u2-access-grid-inherited';
    tr.setAttribute('role', 'row');
    tr.dataset.u2Key = r.row.principal;
    const who = AccessGrid._cell('td', '');
    const name = span(r.label, 'u2-access-grid-principal');
    name.dataset.u2Part = 'item-label';
    who.append(name);
    if (r.from !== null && r.from !== '')
      who.append(span(r.from, 'u2-access-grid-from'));
    tr.append(who);
    for (const c of this._capabilities) {
      const box = document.createElement('input');
      box.type = 'checkbox';
      box.className = 'u2-input-checkbox u2-access-grid-check';
      box.dataset.u2Key = c.name;
      const cell = AccessGrid._cell('td', '', 'u2-access-grid-cap');
      cell.append(box);
      tr.append(cell);
    }
    const last = AccessGrid._cell('td', '');
    tr.append(last);
    if (r.from === null) {
      const remove = document.createElement('button');
      remove.type = 'button';
      remove.className = 'u2-access-grid-remove';
      remove.dataset.u2Key = 'remove';
      remove.textContent = '✕';
      last.append(remove);
    }
    this._patch(tr, r);
    return tr;
  }

  /** Applies a row's state to its kept element; a changed `from` note rebuilds it. */
  private _patch(tr: HTMLElement, r: GridRow): boolean {
    if ((tr.querySelector('.u2-access-grid-from')?.textContent ?? '') !== (r.from ?? ''))
      return false;
    tr.querySelector('.u2-access-grid-principal')!.textContent = r.label;
    for (const box of tr.querySelectorAll<HTMLInputElement>('.u2-access-grid-check')) {
      const capability = this._capabilities.find((c) => c.name === box.dataset.u2Key)!;
      const locked = r.locked.includes(capability.name);
      // a locked capability is not granted whatever the row holds: the box shows what applies
      box.checked = r.row.can[capability.name] === true && !locked;
      box.disabled = !r.enabled || r.from !== null || locked;
      box.setAttribute('aria-label', `${capability.label} for ${r.label}`);
    }
    const remove = tr.querySelector<HTMLButtonElement>('.u2-access-grid-remove');
    if (remove !== null) {
      remove.disabled = !r.enabled;
      remove.setAttribute('aria-label', `Remove ${r.label}`);
    }
    return true;
  }

  /** The native add select: a pick from the list adds at once; a key that moves the closed select's
   * choice (an arrow, Home/End, Page Up/Down, a letter) only moves it, and Enter adds. Leaving it puts
   * the placeholder back. A key's `change` comes in the key's own task; a pick from the list opened by
   * a key (macOS) comes later and adds. */
  private _select(): HTMLSelectElement {
    const select = document.createElement('select');
    select.className = 'u2-access-grid-add';
    select.setAttribute('aria-label', this.options.addText ?? '+ add group or user…');
    const commit = () => {
      this._keyNav = false;
      this._add(select.value);
    };
    this.scope.listen(select, 'keydown', (e) => {
      const key = (e as KeyboardEvent).key;
      if (key === 'Enter' && select.value !== '') {
        e.preventDefault();
        e.stopPropagation();
        commit();
      } else if ((NAV_KEYS.includes(key) || key.length === 1) && !(e as KeyboardEvent).altKey) {
        this._keyNav = true;
        this.scope.setTimeout(() => this._keyNav = false, 0);
      }
    });
    this.scope.listen(select, 'change', () => {
      if (!this._keyNav)
        commit();
    });
    this.scope.listen(select, 'blur', () => {
      this._keyNav = false;
      select.value = '';
    });
    return select;
  }

  private _fillPicker(picker: HTMLSelectElement): void {
    const taken = new Set([...this.value.value ?? [], ...this._inherited.value].map((r) => r.principal));
    const principals = this._principals.value;
    picker.textContent = '';
    picker.append(new Option(this.options.addText ?? '+ add group or user…', ''));
    for (const item of principals) {
      if (!taken.has(itemValue(item)))
        picker.append(new Option(itemLabel(item), itemValue(item)));
    }
    picker.value = '';
  }

  /** Adds an editable row granting the first capability not locked. */
  private _add(principal: string, label?: string): void {
    const rows = this.value.peek() ?? [];
    if (principal === '' || !this.enabled || rows.some((r) => r.principal === principal))
      return;
    const locked = this._locked.peek();
    const first = this._capabilities.find((c) => !locked.includes(c.name));
    const can: Record<string, boolean> = {};
    for (const c of this._capabilities)
      can[c.name] = c === first;
    batch(() => {
      if (label !== undefined && !this._principals.peek().some((p) => itemValue(p) === principal))
        this._principals.value = [...this._principals.peek(), {value: principal, label}];
      this.value.value = [...rows, {principal, can}];
    });
  }

  /** A box that may not change — a disabled grid, an inherited row, a locked capability — springs back. */
  private _onChange(target: HTMLInputElement): void {
    const capability = target.dataset.u2Key;
    const tr = target.closest('tr') as HTMLElement | null;
    const principal = tr?.dataset.u2Key;
    if (capability === undefined || principal === undefined)
      return;
    const fixed = tr!.classList.contains('u2-access-grid-inherited') || this._locked.peek().includes(capability);
    if (!this.enabled || fixed) {
      target.checked = !target.checked;
      return;
    }
    this.value.value = (this.value.peek() ?? []).map((r) =>
      r.principal === principal ? {...r, can: {...r.can, [capability]: target.checked}} : r);
  }

  /** Focus on a removed row goes to the next row's remove, else the previous one's, else the add control. */
  private _onClick(target: HTMLElement): void {
    const remove = target.closest<HTMLElement>('.u2-access-grid-remove');
    if (remove === null || !this.enabled)
      return;
    const tr = remove.closest('tr') as HTMLElement;
    if (tr.contains(document.activeElement)) {
      const removes = [...this._body.querySelectorAll<HTMLElement>('.u2-access-grid-remove')];
      const at = removes.indexOf(remove);
      const add = this._body.parentElement!.querySelector('tfoot')!.querySelector<HTMLElement>(
        'input, select, button, [tabindex]');
      (removes[at + 1] ?? removes[at - 1] ?? add)?.focus();
    }
    this.value.value = (this.value.peek() ?? []).filter((r) => r.principal !== tr.dataset.u2Key);
  }

  private static _cell(tag: 'th' | 'td', text: string, cls?: string): HTMLTableCellElement {
    const el = document.createElement(tag);
    el.textContent = text;
    if (cls !== undefined)
      el.className = cls;
    return el;
  }
}
