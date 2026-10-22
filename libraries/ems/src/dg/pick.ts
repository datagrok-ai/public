/* Reference pickers with an `Input<string | null>` face — the value is the id a ref column
   holds, the box shows the name. `PickInput` is the face over any `TypeAhead`; `domains.pick` is
   the one over a domain table, querying it by its name column through the domain seam, so it
   works over the memory backend as it does over the server. */
import * as DG from 'datagrok-api/dg';
import {emsBackends} from '../core/domain-backend.js';
import {Input, InputOptions, labelText, div, span, text, TypeAhead, Filters, Rows,
  iconButton} from '@datagrok-libraries/u2';
import type {ObjectRenderer, DomainCondition, DomainConditionNode, DomainConditionTree, FilterSchema,
  IProperty} from '@datagrok-libraries/u2';
import type {DomainQueryLike, DomainTableLike} from '../core/domain-backend.js';
import {DomainSource} from '../core/domain-source.js';
import {DomainApp} from './app.js';
import {DomainTable} from './facade.js';

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

/** Names {@link DomainPick.lookUp} asks the target for in one query. */
const LOOKUP_CHUNK = 100;

export class DomainPick extends PickInput<PickItem> {
  /** What a reference no row answers for reads as — the uuid it holds would read as its value. */
  static readonly UNKNOWN = '(unknown)';
  /** How an app opens a row of one of its tables that no open view shows — its own view at the row;
   * false leaves the row to the platform's entity view ({@link open}). */
  static opener: ((table: string, id: string) => boolean) | undefined;
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
          query.trim() !== '' ? 'No matches' :
            `No ${DomainPick._plural.get(table) ?? table.split('.').pop()!.replace(/_/g, ' ').toLowerCase()} yet`,
      }),
      idOf: (item) => item.id,
      resolve: (id) => DomainPick.resolve(table, id),
    });
    this.root.dataset.u2 = 'domain-pick';
    const open = this.runInScope(() =>
      iconButton('external-link-alt', () => DomainPick.open(table, this.value.peek()!), {tooltip: 'Open'}));
    open.dataset.u2Part = 'open';
    this.addOptions(open);
    this.effect(() => open.hidden = !DomainPick.opens(table, this.value.value));
  }

  /** Opens the row a reference holds: the app open over its table brings it forward, else the
   * {@link opener} the app set, else the platform's entity view. After the click that asked for it:
   * the view a click lands in makes itself current as the click bubbles up to it. */
  static open(table: string, id: string): void {
    setTimeout(() => {
      if (DomainApp.activate(DomainApp.baseOf(table), id) || DomainPick.opener?.(table, id) === true)
        return;
      const handler = new DG.DomainObjectHandler(table);
      handler.openRow(handler.rowFrom({id}));
    });
  }

  /** Whether a reference holds a row there is to open — a saved one — and something here opens it: an app open
   * over its table, or the {@link opener} an app set. */
  static opens(table: string, id: string | null): boolean {
    return id !== null && id !== '' && !Rows.isDraft(id) &&
      (DomainApp.baseOf(table) !== undefined || DomainPick.opener !== undefined);
  }

  /** The candidates for `query`: the rows whose name column contains it, in the table's own order
   * (as seeded — statuses, priorities; a declared order is a phase-4 `views` item). A name column
   * that is a formatted number is matched by the digits a typed id carries; text with none is
   * searched for over the table's other searchable columns, and with none of those finds nothing. */
  static async search(table: string, query: string,
    options: {filter?: DomainPickOptions['filter'], limit?: number} = {}): Promise<PickItem[]> {
    const t = await DomainPick.table(table);
    DomainPick._plural.set(table, DomainSource.noun(t.info));
    const nameColumn = DomainSource.captionColumn(t.info);
    const numeric = DomainSource.formattedNumber(t.properties, nameColumn) !== undefined;
    const search = numeric && DomainPick.digits(query) === '' ? query.trim() : '';
    if (search !== '' && t.info.searchableColumns.every((c) => c === nameColumn))
      return [];
    const rows = await t.query({
      filter: DomainPick.filter(nameColumn, query, options.filter, numeric),
      ...(search === '' ? {} : {search}),
      limit: options.limit ?? 10,
    });
    return rows.map((row) => DomainPick.item(row, nameColumn, t.properties));
  }

  /** The item behind an id; a draft id is the draft's caption from the live source holding it,
   * no query issued; an id the table does not answer is {@link UNKNOWN}. A backend that declares
   * `resolveNames` answers it there — over the platform that is the registry's resolver, batched
   * and cached, so a formful of refs is one call; every other backend reads the row. */
  static async resolve(table: string, id: string): Promise<PickItem> {
    const draft = Rows.isDraft(id) ? DomainSource.draftOf(id) : undefined;
    if (draft !== undefined) {
      const info = draft.source.schema.info;
      const name = text(draft.row[DomainSource.captionColumn(info)]);
      return {id, name: name === '' ? `New ${DomainSource.noun(info, false)}` : name};
    }
    const backend = emsBackends.domain;
    if (backend?.resolveNames !== undefined) {
      const name = (await backend.resolveNames(table, [id]))[id];
      if (!name)
        return {id, name: DomainPick.UNKNOWN};
      const format = await DomainPick.nameFormat(table);
      return {id, name: format === null ? name : format(name)};
    }
    const t = await DomainPick.table(table);
    const nameColumn = DomainSource.captionColumn(t.info);
    const [row] = await t.query({filter: {property: 'id', operator: '=', value: id}, limit: 1});
    return row === undefined ? {id, name: DomainPick.UNKNOWN} : DomainPick.item(row, nameColumn, t.properties);
  }

  /** How a name the server answered bare reads in `table`: the registry's resolver and a query's
   * `~caption_<column>` give a numbered name as its digits, and the display format is the table's
   * name column's. Null where that column is not a formatted number, or the table cannot be read;
   * the backend reads a table's schema once per address. */
  static nameFormat(table: string): Promise<((name: string) => string) | null> {
    return DomainPick.table(table).then((t) => {
      const nameColumn = DomainSource.captionColumn(t.info);
      if (DomainSource.formattedNumber(t.properties, nameColumn) === undefined)
        return null;
      return (name: string) =>
        /^\d+$/.test(name) ? DomainTable.captionText(t.properties, nameColumn, Number(name)) : name;
    }, () => null);
  }

  /** The ids of the rows `names` name in `table`, matched exactly on its name column — a number
   * by its value, a formatted one by the digits it carries; a name no row holds is absent. */
  static async lookUp(table: string, names: readonly string[]): Promise<Map<string, string>> {
    const t = await DomainPick.table(table);
    const nameColumn = DomainSource.captionColumn(t.info);
    const prop = t.properties.find((p) => p.name === nameColumn);
    const formatted = DomainSource.formattedNumber(t.properties, nameColumn) !== undefined;
    const number = prop !== undefined && DomainSource.isNumber(prop);
    const byValue = new Map<unknown, string[]>();
    for (const name of names) {
      const typed = formatted ? DomainPick.digits(name) : name.trim();
      const value = number ? Number(typed) : typed;
      if (typed === '' || Number.isNaN(value))
        continue;
      byValue.set(value, [...(byValue.get(value) ?? []), name]);
    }
    const values = [...byValue.keys()];
    const queries: Promise<Record<string, unknown>[]>[] = [];
    for (let at = 0; at < values.length; at += LOOKUP_CHUNK) {
      const chunk = values.slice(at, at + LOOKUP_CHUNK);
      queries.push(t.query({filter: {property: nameColumn, operator: '=', value: chunk},
        columns: [nameColumn], limit: chunk.length}));
    }
    const ids = new Map<string, string>();
    for (const rows of await Promise.all(queries)) {
      for (const row of rows) {
        for (const name of byValue.get(row[nameColumn]) ?? [])
          ids.set(name, String(row.id));
      }
    }
    return ids;
  }

  static isRowRef(prop: IProperty): boolean {
    return DomainTable.isReference(prop) && prop.semType !== 'User' && prop.semType !== 'Group';
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

  /** The digits a typed id carries — `DG-000123` and `123` both mean 123; '' where it has none. */
  static digits(query: string): string {
    return /\d+/.exec(query)?.[0] ?? '';
  }

  static filter(nameColumn: string, query: string,
    extra?: DomainPickOptions['filter'], numeric = false): DomainQueryLike['filter'] {
    const nodes: DomainConditionNode[] = [];
    const q = numeric ? DomainPick.digits(query) : query.trim();
    if (q !== '') {
      nodes.push(numeric ? {property: nameColumn, operator: '=', value: Number(q)} :
        {property: nameColumn, operator: 'like', value: `%${Filters.escapeLike(q)}%`});
    }
    const more = typeof extra === 'string' ?
      (extra.trim() === '' ? undefined : Filters.toDomainTree(Filters.parse(extra).root)) : extra;
    if (more !== undefined) {
      if (nodes.length > 0)
        nodes.push('and');
      nodes.push(more);
    }
    return nodes.length === 0 ? undefined : nodes;
  }

  static item(row: Record<string, unknown>, nameColumn: string, properties: IProperty[]): PickItem {
    const name = DomainTable.captionText(properties, nameColumn, row[nameColumn]);
    return {id: String(row.id), name: name === '' ? String(row.id) : name};
  }

  private static table(address: string): Promise<DomainTableLike> {
    const backend = emsBackends.domain;
    if (backend === undefined)
      return Promise.reject(new Error('no platform backend for domain tables'));
    return backend.table(address);
  }
}
