import {batch, computed, signal, untracked, Signal, ReadonlySignal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {Scope} from '../../core/scope.js';
import {Input} from '../../core/input-base.js';
import {each} from '../../core/region.js';
import {camelCaseToWords, text as asText} from '../../core/text.js';
import type {IItemsState} from '../../core/widget-like.js';
import {TextInput} from '../inputs/text-input.js';
import {NumberInput} from '../inputs/number-input.js';
import {BoolInput} from '../inputs/bool-input.js';
import {ChoiceInput} from '../inputs/choice-input.js';
import {ListInput} from '../inputs/list-input.js';

export interface PropDescriptor {
  name: string;
  /** Platform TYPE strings, as the registry and `IProperty` use them; `choice` is the one
   * addition — a string constrained to {@link choices}. */
  type: 'string' | 'int' | 'double' | 'bool' | 'choice' | 'string_list';
  /** Rows sharing a category sit under one collapsible header; uncategorized rows come first. */
  category?: string;
  choices?: string[];
  min?: number;
  max?: number;
  description?: string;
  readonly?: boolean;
}

export interface PropChange {
  name: string;
  value: unknown;
}

export interface PropertyGridOptions {
  /** An edit made in the grid (not a write to {@link PropertyGrid.values}): the property and its new
   * value. Also the `change` event (`grid.on('change', (change, grid) => …)`). */
  onChanged?: (change: PropChange, grid: PropertyGrid) => void;
}

let gridCount = 0;

/** Compact two-column editor over the u2 inputs: a name column and an inline editor per row.
 * The value record is replaced (never mutated) on every edit, so consumers can bind to it. Rows are
 * keyed by property name: each follows its own value only, so a keystroke re-runs one row. */
export class PropertyGrid extends Control {
  readonly values: Signal<Record<string, unknown>> = signal<Record<string, unknown>>({});

  private readonly _idPrefix = `u2-propgrid-${++gridCount}`;
  private readonly _onChange: PropertyGridOptions['onChanged'];
  private readonly _collapsed = new Map<string, Signal<boolean>>();
  private readonly _groups = signal<readonly [string, PropDescriptor[]][]>([]);
  /** The descriptor each kept row was built for, so a changed one is told from a moved one. */
  private readonly _shown = new WeakMap<HTMLElement, PropDescriptor>();
  private _groupId = 0;

  constructor(options: PropertyGridOptions = {}) {
    super();
    this._onChange = options.onChanged &&
      Input.bindOwned(options.onChanged, {owner: this.scope}) as PropertyGridOptions['onChanged'];
    this.root.classList.add('u2-propgrid');
    this.root.dataset.u2 = 'property-grid';
    this.root.setAttribute('role', 'grid');
    // a kept category follows its rows itself, so a new descriptor list never rebuilds it
    this.build(() => each(this.root, this._groups, ([category]) => `category:${category}`,
      ([category]) => this._group(category), {update: () => true}));
  }

  /** The row descriptors shown; a new list keeps the current values and every unchanged row. */
  get properties(): PropDescriptor[] {
    return this._groups.peek().flatMap(([, props]) => props);
  }

  set properties(x: PropDescriptor[]) {
    this.setProperties(x ?? [], this.values.peek());
  }

  /** Shows `props` over `values`. A row whose property keeps its name and descriptor keeps its
   * editor (and focus); a changed descriptor rebuilds that row, a gone one disposes it. Category
   * collapse state is kept by name. */
  setProperties(props: PropDescriptor[], values: Record<string, unknown>): void {
    batch(() => {
      this.values.value = {...values};
      this._groups.value = PropertyGrid._group(props);
    });
  }

  /** One item per row: key = the property name, label = the caption the row shows. */
  getItems(opts?: {all?: boolean}): IItemsState {
    const rendered = new Set<string>();
    for (const row of this.root.querySelectorAll<HTMLElement>('.u2-propgrid-row'))
      rendered.add(row.dataset.u2Item!);
    const list = this._groups.peek().flatMap(([, props]) => props).map((prop, index) => ({key: prop.name,
      label: camelCaseToWords(prop.name), index, ...prop.readonly ? {disabled: true} : {}, rendered: rendered.has(prop.name)}));
    return Input.itemsState(list, opts?.all);
  }

  /** A category: its header (when it has a name) and its rows, keyed by property name. */
  private _group(category: string): HTMLElement {
    const section = document.createElement('div');
    section.className = 'u2-propgrid-section';
    section.setAttribute('role', 'rowgroup');
    const group = document.createElement('div');
    group.className = 'u2-propgrid-group';
    if (category) {
      group.id = `${this._idPrefix}-${this._groupId++}`;
      section.append(this._createHeader(category, group, Scope.ambient!));
    }
    section.append(group);
    const props = computed(() => this._groups.value.find(([c]) => c === category)?.[1] ?? []);
    each(group, props, (prop) => prop.name, (prop) => this._createRow(prop, Scope.ambient!), {update: (el, prop) => {
      const kept = PropertyGrid._sameRow(this._shown.get(el), prop);
      if (kept)
        this._shown.set(el, prop);
      return kept;
    }});
    return section;
  }

  private _createHeader(category: string, group: HTMLElement, scope: Scope): HTMLElement {
    let collapsed = this._collapsed.get(category);
    if (!collapsed) {
      collapsed = signal(false);
      this._collapsed.set(category, collapsed);
    }
    const state = collapsed;

    const header = document.createElement('div');
    header.className = 'u2-propgrid-category';
    header.dataset.u2Part = 'category';
    header.tabIndex = 0;
    header.title = category;
    header.setAttribute('role', 'button');
    header.setAttribute('aria-controls', group.id);
    const chevron = document.createElement('span');
    chevron.className = 'u2-propgrid-chevron';
    chevron.textContent = '▸';
    const title = document.createElement('span');
    title.className = 'u2-propgrid-category-title';
    title.textContent = category;
    header.append(chevron, title);

    scope.effect(() => {
      const on = state.value;
      header.setAttribute('aria-expanded', String(!on));
      group.style.display = on ? 'none' : '';
    });
    const toggle = () => state.value = !state.peek();
    scope.listen(header, 'click', toggle);
    scope.listen(header, 'keydown', (e) => {
      const key = (e as KeyboardEvent).key;
      if (key !== 'Enter' && key !== ' ')
        return;
      e.preventDefault();
      toggle();
    });
    return header;
  }

  /** A row in its own scope: it follows its own value through a computed, so an edit elsewhere in
   * the record never wakes it. */
  private _createRow(prop: PropDescriptor, scope: Scope): HTMLElement {
    const row = document.createElement('div');
    this._shown.set(row, prop);
    row.className = 'u2-propgrid-row';
    row.setAttribute('role', 'row');
    // the caption is reworded for reading; automation keeps addressing the row by its property name
    row.dataset.u2Name = prop.name;
    const name = document.createElement('div');
    name.className = 'u2-propgrid-name';
    name.setAttribute('role', 'rowheader');
    name.dataset.u2Part = 'item-label';
    // the platform's caption for a property without a friendly name (`showLegend` → Show Legend)
    name.textContent = camelCaseToWords(prop.name);
    if (prop.description)
      name.title = prop.description;
    const cell = document.createElement('div');
    cell.className = 'u2-propgrid-editor';
    cell.setAttribute('role', 'gridcell');
    row.append(name, cell);
    const value = computed(() => this.values.value[prop.name]);

    if (prop.readonly) {
      row.classList.add('u2-propgrid-row-readonly');
      row.setAttribute('aria-readonly', 'true');
      const text = document.createElement('span');
      text.className = 'u2-propgrid-value';
      cell.append(text);
      scope.effect(() => text.textContent = asText(value.value));
      return row;
    }

    const input = this._createInput(prop);
    cell.append(input.root);
    this._bridge(prop, input, value, scope);
    return row;
  }

  private _createInput(prop: PropDescriptor): Input<any> {
    switch (prop.type) {
      case 'bool':
        return new BoolInput({inline: true});
      case 'choice':
        return new ChoiceInput({inline: true, items: prop.choices ?? []});
      case 'int':
      case 'double':
        return new NumberInput({inline: true, mode: prop.type === 'int' ? 'int' : 'float',
          min: prop.min, max: prop.max});
      case 'string_list':
        return new ListInput({inline: true, placeholder: prop.description});
      default:
        return new TextInput({inline: true});
    }
  }

  // Echo suppression compares values instead of raising a flag: a write made inside an effect is
  // flushed after that effect returns, so a transient flag would already be down by then.
  private _bridge(prop: PropDescriptor, input: Input<any>, value: ReadonlySignal<unknown>, scope: Scope): void {
    scope.effect(() => input.value = PropertyGrid._coerce(prop, value.value));
    scope.effect(() => {
      const edited = input.value;
      const record = this.values.peek();
      if (PropertyGrid.same(prop, PropertyGrid._coerce(prop, record[prop.name]), edited))
        return;
      this.values.value = {...record, [prop.name]: edited};
      const change = {name: prop.name, value: edited};
      untracked(() => {
        this._onChange?.(change, this);
        this.fireEvent('change', change);
      });
    });
  }

  /** Whether a kept row can show `next` as it is: the same editor over the same constraints. */
  private static _sameRow(prev: PropDescriptor | undefined, next: PropDescriptor): boolean {
    return prev === next || prev !== undefined && prev.type === next.type && !!prev.readonly === !!next.readonly &&
      prev.min === next.min && prev.max === next.max && prev.description === next.description &&
      (prev.choices ?? []).join('\n') === (next.choices ?? []).join('\n');
  }

  private static _group(props: PropDescriptor[]): [string, PropDescriptor[]][] {
    const groups = new Map<string, PropDescriptor[]>();
    groups.set('', []);
    for (const prop of props) {
      const key = prop.category ?? '';
      let list = groups.get(key);
      if (!list) {
        list = [];
        groups.set(key, list);
      }
      list.push(prop);
    }
    return [...groups].filter(([, list]) => list.length > 0);
  }

  private static _coerce(prop: PropDescriptor, value: unknown): any {
    switch (prop.type) {
      case 'bool':
        return value === true;
      case 'int':
      case 'double': {
        if (typeof value === 'number')
          return value;
        const parsed = value === null || value === undefined || value === '' ? NaN : Number(value);
        return isFinite(parsed) ? parsed : null;
      }
      case 'choice':
        return value === null || value === undefined ? null : String(value);
      case 'string_list':
        return Array.isArray(value) && value.every((item) => typeof item === 'string') ? value : [];
      default:
        return asText(value);
    }
  }

  /** A list editor answers a fresh array for the same list, so the echo check above compares
   * element-wise there — identity would report every rebuild as an edit and loop. Public because
   * anything comparing a row's old and new value has to compare it the way the grid does. */
  static same(prop: PropDescriptor, a: unknown, b: unknown): boolean {
    if (prop.type !== 'string_list' || !Array.isArray(a) || !Array.isArray(b))
      return a === b;
    return a.length === b.length && a.every((item, i) => item === b[i]);
  }
}
