/* The components tray (Q1): the strip along the bottom of the design surface, one chip per
   non-visual component. A component has no element on the canvas, so its chip is where it is
   selected, acted on — and addressed by automation, which is why `data-u2-name` lands here. The `+`
   chip is the guided way in (Q7): the palette drag stays the low-level one, and both end in exactly
   one `add-component`. */
import * as grok from 'datagrok-api/grok';
import {Component, Control} from '../../core/component.js';
import {Scope} from '../../core/scope.js';
import {div} from '../../core/elements.js';
import {iconButton} from '../../components/actions/buttons.js';
import {ChoiceInput} from '../../components/inputs/choice-input.js';
import {Dialog} from '../../components/containers/dialog.js';
import {Menu} from '../../components/navigation/menu.js';
import {NumberInput} from '../../components/inputs/number-input.js';
import {TextInput} from '../../components/inputs/text-input.js';
import {nameForTag, seedNode} from '../../spec/editor.js';
import type {ComponentMeta} from '../../spec/registry.js';
import type {SpecInstance, SpecNode} from '../../spec/spec.js';
import {backends} from '../../sources/backends.js';
import {COLLECTIONS} from '../../sources/registrations.js';
import {funcPicker, paramValues} from './func-picker.js';
import type {FuncPick} from './func-picker.js';
import {SpecNodeRef, nodeLabel} from './node-ref.js';
import {sourceStatus, statusText} from './source-status.js';

export interface TrayOptions {
  onSelect: (node: SpecNode) => void;
  onContext: (node: SpecNode, x: number, y: number) => void;
  /** The `+` menu's insert: the view names the node — uniqueness spans both trees — and applies
   * one `add-component` patch. */
  onAdd: (base: string, seed: (name: string) => SpecNode) => void;
}

const HINT = 'Data sources';
const FUNC = 'u2-func-source';
const TABLE = 'u2-table-source';
const ENTITY = 'u2-entity-source';
const STATE = 'u2-state';
/** DATA_BINDING §8: a nudge, not a gate — a table opened by hand is not a durable dependency. */
const BY_NAME = 'references the open table by name; prefer a query or project table for a durable form';

/** A seeded source node: the tag's defaults, the wizard's own props on top, and its binds. */
export function sourceNode(meta: ComponentMeta | undefined, tag: string, name: string,
  props: Record<string, unknown>, binds?: Record<string, string>): SpecNode {
  const node = seedNode(meta, tag, name);
  const merged = {...node.props, ...props};
  if (Object.keys(merged).length > 0)
    node.props = merged;
  if (binds !== undefined && Object.keys(binds).length > 0)
    node.bind = binds;
  return node;
}

/** The func-source a picked function becomes: its literal arguments as `params`, its bound ones as
 * the dotted `params.<name>` keys. No `designData` — the source computes it from the function's
 * kind, and seeding it here would freeze a choice made before a function was picked. */
export function funcSourceNode(meta: ComponentMeta | undefined, name: string, pick: FuncPick): SpecNode {
  const props: Record<string, unknown> = {func: pick.name};
  const params = paramValues(pick.args);
  if (Object.keys(params).length > 0)
    props.params = params;
  const binds: Record<string, string> = {};
  for (const [param, path] of Object.entries(pick.binds))
    binds[`params.${param}`] = path;
  return sourceNode(meta, FUNC, name, props, binds);
}

export class Tray extends Control {
  private readonly _chips = new Map<HTMLElement, SpecNode>();
  private readonly _plus: HTMLButtonElement;
  private _instance: SpecInstance | undefined;
  /** The chips' live effects: one generation per {@link update}, since every chip is a new element
   * over a component the last patch may have rebuilt. */
  private _pulse: Scope | undefined;

  constructor(private readonly _options: TrayOptions) {
    super();
    this.own(() => this._pulse?.dispose());
    this.root.classList.add('u2-designer-tray');
    this.root.dataset.u2 = 'tray';
    this._plus = iconButton('plus', () => this._menu(), {tooltip: 'Add a data source'});
    this._plus.classList.add('u2-designer-tray-add');
    const chip = (e: Event): SpecNode | undefined =>
      this._chips.get((e.target as Element).closest('.u2-designer-chip') as HTMLElement);
    this._listen('click', (e) => {
      const node = chip(e);
      if (node)
        _options.onSelect(node);
    });
    this._listen('contextmenu', (e) => {
      const me = e as MouseEvent;
      me.preventDefault();
      me.stopPropagation();
      const node = chip(e);
      if (node)
        _options.onContext(node, me.clientX, me.clientY);
    });
  }

  /** Rebuilt from the spec on every selection and every structural patch — a component is rebuilt
   * as a new object, so nothing here may outlive the state it was drawn from. */
  update(instance: SpecInstance | undefined, selected: SpecNode | null): void {
    this._instance = instance;
    this._chips.clear();
    this._pulse?.dispose();
    this._pulse = undefined;
    this.root.textContent = '';
    this.root.append(this._plus);
    const components = instance?.spec.components ?? [];
    if (components.length === 0) {
      this.root.append(div([HINT], 'u2-designer-tray-hint'));
      return;
    }
    const pulse = new Scope();
    this._pulse = pulse;
    for (const node of components) {
      const chip = div([], 'u2-designer-chip');
      const ref = new SpecNodeRef(instance!, node);
      chip.textContent = nodeLabel(node);
      chip.title = ref.error() ?? node.tag;
      if (node.name !== undefined)
        chip.dataset.u2Name = node.name;
      chip.classList.toggle('u2-designer-chip-broken', ref.error() !== null);
      chip.classList.toggle('u2-designer-chip-selected', node === selected);
      this._chips.set(chip, node);
      this.root.append(chip);
      Tray._pulseOf(pulse, chip, ref);
    }
  }

  /** A chip follows the source it stands for: loading and failure are the two things a form full
   * of empty fields cannot tell apart on its own. */
  private static _pulseOf(scope: Scope, chip: HTMLElement, ref: SpecNodeRef): void {
    const built = ref.built();
    if (!Component.is(built))
      return;
    const label = ref.error() ?? ref.node.tag;
    scope.effect(() => {
      const at = sourceStatus(built);
      chip.classList.toggle('u2-designer-chip-loading', at?.state === 'loading');
      chip.classList.toggle('u2-designer-chip-error', at?.error != null);
      chip.title = at === null ? label : `${label} — ${statusText(at)}`;
    });
  }

  /** On while a non-visual palette item is being dragged: wherever it is dropped, it lands here. */
  setDropTarget(on: boolean): void {
    this.root.classList.toggle('u2-designer-tray-drop', on);
  }

  /** The four ways in (Q7): the tables the user already has open, a function or query, a server
   * collection, and a value the form owns. */
  private _menu(): void {
    const instance = this._instance;
    if (instance === undefined)
      return;
    const menu = new Menu();
    const tables = backends.workspace?.tableNames() ?? [];
    menu.group('Open tables');
    if (tables.length === 0)
      menu.item('No open tables', () => {}, {enabled: false});
    for (const table of tables)
      menu.item(table, () => this._addTable(instance, table));
    menu.endGroup();
    menu.item('Function or query…', () => this._addFunc(instance));
    menu.item('Entity collection…', () => this._addEntities(instance));
    menu.item('State variable', () => this._insert(instance, STATE, nameForTag(STATE), {}));
    menu.show({anchor: this._plus});
  }

  private _addTable(instance: SpecInstance, table: string): void {
    this._insert(instance, TABLE, specName(table, nameForTag(TABLE)), {table});
    grok.shell.info(`${table} ${BY_NAME}.`);
  }

  private _addFunc(instance: SpecInstance): void {
    this.run(() => funcPicker({title: 'Data source: function or query', instance,
      onPick: (pick) => this._options.onAdd(nameForTag(FUNC), (name) =>
        funcSourceNode(instance.registry.get(FUNC), name, pick))}));
  }

  private _addEntities(instance: SpecInstance): void {
    this.run(() => {
      const collection = new ChoiceInput({label: 'Collection', items: COLLECTIONS, value: COLLECTIONS[0]});
      const filter = new TextInput({label: 'Filter', placeholder: 'A smart-search filter'});
      const pageSize = new NumberInput({label: 'Page size', mode: 'int', value: 20});
      const dialog = new Dialog('Data source: entity collection');
      dialog.add(collection).add(filter).add(pageSize)
        .onOK(() => {
          const props: Record<string, unknown> = {entity: collection.value.peek() ?? COLLECTIONS[0],
            pageSize: pageSize.value.peek() ?? 20};
          const text = filter.value.peek().trim();
          if (text !== '')
            props.filter = text;
          dialog.dispose();
          this._insert(instance, ENTITY, nameForTag(ENTITY), props);
        })
        .onCancel(() => dialog.dispose())
        .show({modal: true, width: 360});
    });
  }

  private _insert(instance: SpecInstance, tag: string, base: string,
    props: Record<string, unknown>): void {
    this._options.onAdd(base, (name) => sourceNode(instance.registry.get(tag), tag, name, props));
  }

  private _listen(type: string, handler: (e: Event) => void): void {
    this.root.addEventListener(type, handler);
    this.own(() => this.root.removeEventListener(type, handler));
  }
}

/** A table name as a spec name — what a `$.` path can address; the tag's own base where nothing
 * usable is left of it. */
export function specName(text: string, fallback: string): string {
  const parts = text.split(/[^A-Za-z0-9_$]+/).filter((part) => part !== '');
  const name = parts.map((part, i) =>
    i === 0 ? part : part.charAt(0).toUpperCase() + part.slice(1)).join('');
  return name === '' || /^\d/.test(name) ? fallback : name;
}
