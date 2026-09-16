/* Reference pickers with an `Input<string | null>` face — the value is the id a ref column
   holds, the box shows the name. `PickInput` is the face over any `TypeAhead`; `domains.pick` is
   the one over a domain table, querying it by its name column through the domain seam, so it
   works over the memory backend as it does over the server. */
import * as grok from 'datagrok-api/grok';
import {Input, InputOptions, labelText} from '../../core/input-base.js';
import {div, span} from '../../core/elements.js';
import type {ObjectRenderer} from '../../core/object-renderer.js';
import {text} from '../../core/text.js';
import {TypeAhead} from '../../components/inputs/typeahead.js';
import {Filters} from '../../core/filter/index.js';
import type {DomainCondition, DomainConditionNode, DomainConditionTree, FilterSchema} from '../../core/filter/index.js';
import {backends} from '../../sources/backends.js';
import type {DomainQueryLike, DomainTableInfoLike, DomainTableLike} from '../../sources/domain-backend.js';
import {DomainSource} from '../../sources/domain-source.js';
import {Rows} from '../../sources/rows-like.js';
import {DgDomainBackend} from './backend.js';

export interface PickItem {
  id: string;
  name: string;
}

export interface PickInputOptions<T> extends InputOptions<string | null> {
  /** Builds the type-ahead the input is a face of. */
  typeAhead: () => TypeAhead<T>;
  idOf: (item: T) => string;
  /** The item behind a value written from outside — the id a form loaded; null when unknown. */
  resolve: (id: string) => Promise<T | null>;
}

export class PickInput<T> extends Input<string | null, PickInputOptions<T>> {
  private _typeAhead!: TypeAhead<T>;
  private _hint!: HTMLElement;

  constructor(options: PickInputOptions<T>) {
    super(options, null);
    this.root.dataset.u2 = 'pick-input';
    this.root.append(this._hint);
  }

  get typeAhead(): TypeAhead<T> {
    return this._typeAhead;
  }

  // runs from the base constructor, so the coupling state lives here rather than in fields the
  // subclass initializers would reset afterwards
  protected createEditor(): HTMLElement {
    const {typeAhead, idOf, resolve} = this.options;
    const picker = typeAhead();
    this._typeAhead = picker;
    // the last pick is held through typing, so blur can put its name back
    let held: T | null = null;
    let resolving: string | null = null;
    let gen = 0;
    const lookUp = async (id: string) => {
      const mine = ++gen;
      resolving = id;
      try {
        const item = await resolve(id);
        if (mine === gen && this.value.peek() === id && item !== null)
          picker.selected.value = item;
      } catch {
        // an unknown id keeps the value and shows nothing
      } finally {
        if (mine === gen)
          resolving = null;
      }
    };
    // a form seeds a string field with '' for an empty cell: empty either way, and never rewritten
    // to null on the input's own account — that write would dirty the row on render
    const empty = (id: string | null) => id === null || id === '';
    this.effect(() => {
      const id = this.value.value;
      const current = picker.selected.peek();
      if ((current === null ? null : idOf(current)) === id)
        return;
      if (empty(id))
        picker.selected.value = null;
      else
        void lookUp(id!);
    });
    this.effect(() => {
      const item = picker.selected.value;
      // a value still awaiting its item is not a clear
      if (item === null && resolving !== null)
        return;
      const id = item === null ? null : idOf(item);
      const current = this.value.peek();
      // a selection dropped by typing keeps the value; one dropped by clearing the text does not
      if (id !== null)
        held = item;
      else if (picker.text.peek() !== '')
        return;
      else
        held = null;
      if (id === null ? !empty(current) : id !== current)
        this.value.value = id;
    });
    const input = picker.root.querySelector('input')!;
    const hint = span('', 'u2-pick-hint');
    hint.dataset.u2Part = 'hint';
    hint.hidden = true;
    this._hint = hint;
    // text that is not a pick never survives a blur: the held pick's name comes back (and a line
    // says so until the next focus), else the box is cleared — text that looks chosen is worse
    // than none (an Enter still awaiting its candidates is left to land)
    const onBlur = () => {
      if (picker.isPickPending || picker.selected.peek() !== null || picker.text.peek() === '')
        return;
      if (held !== null) {
        const typed = picker.text.peek();
        picker.selected.value = held;
        hint.textContent = `No match for "${typed}" — kept ${picker.text.peek()}`;
        hint.hidden = false;
      } else
        picker.resync();
    };
    const onFocus = () => hint.hidden = true;
    input.addEventListener('blur', onBlur);
    input.addEventListener('focus', onFocus);
    this.own(() => {
      input.removeEventListener('blur', onBlur);
      input.removeEventListener('focus', onFocus);
    });
    const clear = document.createElement('button');
    clear.type = 'button';
    clear.className = 'u2-input-clear';
    clear.textContent = '✕';
    clear.setAttribute('aria-label', 'Clear');
    const onClear = () => {
      held = null;
      picker.selected.value = null;
      picker.resync();
      this.value.value = null;
      input.focus();
    };
    clear.addEventListener('click', onClear);
    this.own(() => clear.removeEventListener('click', onClear));
    this.effect(() => clear.hidden = empty(this.value.value) && picker.text.value === '');
    return div([picker.root, clear], 'u2-pick-box');
  }
}

export interface DomainPickOptions extends InputOptions<string | null> {
  /** Narrows the candidates: a smart-filter string or a condition tree, AND-ed with the search.
   * A `$name` in the string is a sibling's value, bound through {@link params} at every search. */
  filter?: string | DomainCondition | DomainConditionTree;
  /** The values `$name`s bind to — the row being edited; a `$name` left empty means no
   * candidates and a "Pick a <sibling> first" hint. */
  params?: () => Record<string, unknown>;
  /** The schema the `$name`s belong to — the hint's captions. */
  siblings?: FilterSchema;
  placeholder?: string;
  /** Row presentation; the name as text by default. */
  renderer?: ObjectRenderer<PickItem>;
  /** How many candidates one search shows (default 10). */
  limit?: number;
  debounceMs?: number;
}

export class DomainPick extends PickInput<PickItem> {
  /** The target's plural name once a search has seen its table — the empty row's wording. */
  private static readonly _plural = new Map<string, string>();

  constructor(readonly table: string, options: DomainPickOptions = {}) {
    const {filter, params, siblings, placeholder, renderer, limit, debounceMs, ...rest} = options;
    const listItem = renderer?.listItem?.bind(renderer);
    let unbound: string[] = [];
    const caption = (name: string) => siblings?.properties.find((p) => p.name === name)?.friendlyName ?? name;
    super({
      ...rest,
      typeAhead: () => new TypeAhead<PickItem>({
        source: (query) => {
          const bound = DomainPick.bind(filter, params?.());
          unbound = bound.unbound;
          return unbound.length > 0 ? Promise.resolve([]) :
            DomainPick.search(table, query, {filter: bound.filter, limit});
        },
        itemText: (item) => item.name,
        render: listItem,
        placeholder: placeholder ?? `${labelText(options.label) ?? table.split('.').pop()}…`,
        debounceMs,
        openOnFocus: true,
        // an empty table is not a dead end: say what there is none of yet
        emptyText: (query) => unbound.length > 0 ? `Pick a ${caption(unbound[0]).toLowerCase()} first` :
          query.trim() === '' ?
            `No ${(DomainPick._plural.get(table) ?? table.split('.').pop()!).toLowerCase()} yet` : 'No matches',
      }),
      idOf: (item) => item.id,
      resolve: (id) => DomainPick.resolve(table, id),
    });
    this.root.dataset.u2 = 'domain-pick';
  }

  /** The candidates for `query`: the rows whose name column contains it, in the table's own order
   * (as seeded — statuses, priorities; a declared order is a phase-4 `views` item). */
  static async search(table: string, query: string,
    options: {filter?: DomainPickOptions['filter'], limit?: number} = {}): Promise<PickItem[]> {
    const t = await DomainPick.table(table);
    DomainPick._plural.set(table, t.info.pluralName);
    const nameColumn = DomainPick.nameColumn(t.info);
    const rows = await t.query({
      filter: DomainPick.filter(nameColumn, query, options.filter),
      limit: options.limit ?? 10,
    });
    return rows.map((row) => DomainPick.item(row, nameColumn));
  }

  /** The item behind an id; a draft id is the draft's caption from the live source holding it,
   * no query issued; an id the table does not answer is shown as itself. Over the platform the
   * registry's resolver answers it — batched over 30 ms and cached, so a formful of refs is one
   * call; every other backend reads the row. */
  static async resolve(table: string, id: string): Promise<PickItem> {
    const draft = Rows.isDraft(id) ? DomainSource.draftOf(id) : undefined;
    if (draft !== undefined) {
      const info = draft.source.schema.info;
      const name = text(draft.row[DomainPick.nameColumn(info)]);
      return {id, name: name === '' ? `New ${info.singularName.toLowerCase() || 'row'}` : name};
    }
    if (backends.domain instanceof DgDomainBackend) {
      const name = (await grok.dapi.domains.registry.resolveNames(table, [id]))[id];
      return {id, name: name || id};
    }
    const t = await DomainPick.table(table);
    const [row] = await t.query({filter: {property: 'id', operator: '=', value: id}, limit: 1});
    return row === undefined ? {id, name: id} : DomainPick.item(row, DomainPick.nameColumn(t.info));
  }

  /** The declared name column, else the first business-key column, else the id. */
  static nameColumn(info: DomainTableInfoLike): string {
    return info.nameColumn ?? info.businessKey[0] ?? 'id';
  }

  /** A string filter with its `$name`s bound to `params` (an empty sibling counts as absent), as a
   * tree; the names still unbound, so the picker can say which sibling to pick first. */
  static bind(extra: DomainPickOptions['filter'], params?: Record<string, unknown>):
    {filter?: DomainCondition | DomainConditionTree, unbound: string[]} {
    if (extra === undefined || (typeof extra === 'string' && extra.trim() === ''))
      return {unbound: []};
    if (typeof extra !== 'string')
      return {filter: extra, unbound: []};
    const values = Object.fromEntries(Object.entries(params ?? {})
      .filter(([, v]) => v !== null && v !== undefined && v !== ''));
    const root = Filters.bind(Filters.parse(extra).root, values);
    const unbound: string[] = [];
    Filters.walk(root, (n) => {
      if (Filters.isGroup(n))
        return;
      for (const v of Array.isArray(n.value) ? n.value : [n.value]) {
        if (Filters.isParam(v) && !unbound.includes(v.param))
          unbound.push(v.param);
      }
    });
    return {filter: unbound.length > 0 ? undefined : Filters.toDomainTree(root), unbound};
  }

  static filter(nameColumn: string, query: string,
    extra?: DomainPickOptions['filter']): DomainQueryLike['filter'] {
    const nodes: DomainConditionNode[] = [];
    const q = query.trim();
    if (q !== '')
      nodes.push({property: nameColumn, operator: 'like', value: `%${Filters.escapeLike(q)}%`});
    const more = typeof extra === 'string' ?
      (extra.trim() === '' ? undefined : Filters.toDomainTree(Filters.parse(extra).root)) : extra;
    if (more !== undefined) {
      if (nodes.length > 0)
        nodes.push('and');
      nodes.push(more);
    }
    return nodes.length === 0 ? undefined : nodes;
  }

  static item(row: Record<string, unknown>, nameColumn: string): PickItem {
    const name = text(row[nameColumn]);
    return {id: String(row.id), name: name === '' ? String(row.id) : name};
  }

  private static table(address: string): Promise<DomainTableLike> {
    const backend = backends.domain;
    if (backend === undefined)
      return Promise.reject(new Error('no platform backend for domain tables'));
    return backend.table(address);
  }
}
