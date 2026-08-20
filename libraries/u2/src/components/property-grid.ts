import {signal, Signal, ReadonlySignal} from '../core/signals.js';
import {Scope} from '../core/scope.js';
import {Control} from '../core/component.js';
import {Input} from '../core/input-base.js';
import {TextInput} from './text-input.js';
import {NumberInput} from './number-input.js';
import {BoolInput} from './bool-input.js';
import {ChoiceInput} from './choice-input.js';

export interface PropDescriptor {
  name: string;
  /** Platform TYPE strings, as the registry and `PropertyLike` use them; `choice` is the one
   * addition — a string constrained to {@link choices}. */
  type: 'string' | 'int' | 'double' | 'bool' | 'choice';
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

let gridCount = 0;

/** Compact two-column editor over the u2 inputs: a name column and an inline editor per row.
 * The value record is replaced (never mutated) on every edit, so consumers can bind to it. */
export class PropertyGrid extends Control {
  readonly values: Signal<Record<string, unknown>> = signal<Record<string, unknown>>({});
  readonly onChanged: ReadonlySignal<PropChange | null>;

  private readonly _idPrefix = `u2-propgrid-${++gridCount}`;
  private readonly _lastChange = signal<PropChange | null>(null);
  private readonly _collapsed = new Map<string, Signal<boolean>>();
  private _rows = new Scope();
  private _groupId = 0;

  constructor() {
    super();
    this.onChanged = this._lastChange;
    this.root.classList.add('u2-propgrid');
    this.root.dataset.u2 = 'property-grid';
    this.root.setAttribute('role', 'group');
    this.own(() => this._rows.dispose());
  }

  /** Rebuilds every row under a fresh scope: the previous editors and their bridges are disposed,
   * so nothing stale writes into the new value record. Category collapse state is kept by name. */
  setProperties(props: PropDescriptor[], values: Record<string, unknown>): void {
    this._rows.dispose();
    this._rows = new Scope();
    this._groupId = 0;
    this.root.textContent = '';
    this._lastChange.value = null;
    this.values.value = {...values};
    Scope.runWith(this._rows, () => {
      for (const [category, group] of PropertyGrid._group(props))
        this._addGroup(category, group);
    });
  }

  private _addGroup(category: string, props: PropDescriptor[]): void {
    const group = document.createElement('div');
    group.className = 'u2-propgrid-group';
    if (category) {
      group.id = `${this._idPrefix}-${this._groupId++}`;
      this.root.append(this._createHeader(category, group));
    }
    for (const prop of props)
      group.append(this._createRow(prop));
    this.root.append(group);
  }

  private _createHeader(category: string, group: HTMLElement): HTMLElement {
    let collapsed = this._collapsed.get(category);
    if (!collapsed) {
      collapsed = signal(false);
      this._collapsed.set(category, collapsed);
    }
    const state = collapsed;

    const header = document.createElement('div');
    header.className = 'u2-propgrid-category';
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

    this._rows.effect(() => {
      const on = state.value;
      header.setAttribute('aria-expanded', String(!on));
      group.style.display = on ? 'none' : '';
    });
    const toggle = () => state.value = !state.peek();
    this._listen(header, 'click', toggle);
    this._listen(header, 'keydown', (e) => {
      const key = (e as KeyboardEvent).key;
      if (key !== 'Enter' && key !== ' ')
        return;
      e.preventDefault();
      toggle();
    });
    return header;
  }

  private _createRow(prop: PropDescriptor): HTMLElement {
    const row = document.createElement('div');
    row.className = 'u2-propgrid-row';
    const name = document.createElement('div');
    name.className = 'u2-propgrid-name';
    name.textContent = prop.name;
    name.title = prop.description ?? prop.name;
    const cell = document.createElement('div');
    cell.className = 'u2-propgrid-editor';
    row.append(name, cell);

    if (prop.readonly) {
      row.classList.add('u2-propgrid-row-readonly');
      const text = document.createElement('span');
      text.className = 'u2-propgrid-value';
      cell.append(text);
      this._rows.effect(() => text.textContent = PropertyGrid._format(this.values.value[prop.name]));
      return row;
    }

    const input = this._createInput(prop);
    cell.append(input.root);
    this._bridge(prop, input);
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
      default:
        return new TextInput({inline: true});
    }
  }

  // Echo suppression compares values instead of raising a flag: a write made inside an effect is
  // flushed after that effect returns, so a transient flag would already be down by then.
  private _bridge(prop: PropDescriptor, input: Input<any>): void {
    this._rows.effect(() => input.value.value = PropertyGrid._coerce(prop, this.values.value[prop.name]));
    this._rows.effect(() => {
      const value = input.value.value;
      const record = this.values.peek();
      if (PropertyGrid._coerce(prop, record[prop.name]) === value)
        return;
      this.values.value = {...record, [prop.name]: value};
      this._lastChange.value = {name: prop.name, value};
    });
  }

  private _listen(el: HTMLElement, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this._rows.own(() => el.removeEventListener(type, handler));
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
      default:
        return value === null || value === undefined ? '' : String(value);
    }
  }

  private static _format(value: unknown): string {
    return value === null || value === undefined ? '' : String(value);
  }
}
