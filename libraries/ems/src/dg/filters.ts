/* `domains.filters` — the query box (or the builder) over a source's query: the table's own
   filter schema, the platform's facet values where there is a platform, two-way with
   `source.query`. A change the USER makes goes through the gate (STATE-CONTRACT H6) — it answers
   at once while there is nothing to lose and refuses while a batch is being written back; cancel
   puts the previous filter back. A query set in code is the source's to skip. */
import {domainFilterSchema} from './filter-schema.js';
import {emsBackends} from '../core/domain-backend.js';
import {Control, computed, signal, ReadonlySignal, Filters, FilterQueryInput,
  FilterBuilder} from '@datagrok-libraries/u2';
import type {FilterGroup, FilterProblem, FilterProperty, FilterSchema, FilterValueItem} from '@datagrok-libraries/u2';
import type {DomainSource} from '../core/domain-source.js';
import {confirmDiscard} from '../core/session.js';
import {DomainPick} from './pick.js';
import {DgDomainBackend, SYSTEM_COLUMNS} from './backend.js';
import {DomainTable} from './facade.js';

export type DomainFiltersMode = 'query' | 'builder';

/** The subtree operator, offered only where the hierarchy it walks exists. */
const UNDER = 'under';
/** What a `'basic'` backend (`support.filters`) runs on no column at all; the datetime `!=` and
 * the bool null tests it also lacks go by kind (ARCHITECTURE §5.2 of the external bindings). */
/** What a `'basic'` backend is not offered: the operators it refuses by name, and the two whose
 * emitted payload is a nested group it cannot run — `fuzzy` (two clauses under an `or`) and
 * `between` (a pair under an `and`, on every kind). */
const NOT_BASIC = new Set([UNDER, 'matches', '!matches', '!like', 'fuzzy', 'between']);
const NULL_TESTS = new Set(['is null', 'is not null']);
/** How many rows of the target table an `under` value slot offers at once. */
const PICK_LIMIT = 50;
const EMPTY: Promise<FilterValueItem[]> = Promise.resolve([]);

export interface DomainFiltersOptions {
  /** The one-line query box with completion (default), or the row-per-condition builder. */
  mode?: DomainFiltersMode;
  placeholder?: string;
}

export class DomainFilters extends Control {
  readonly mode: DomainFiltersMode;
  /** The input, once the table's schema is known. */
  readonly input: ReadonlySignal<FilterQueryInput | FilterBuilder | null>;
  /** What is wrong with the filter in the box while it is one the rows do NOT answer — a refused
   * commit standing over the applied query, or a builder row whose value is wrong; null while the two
   * agree, the query in force having been refused included (the list's error says what happened to
   * that one). */
  readonly problem: ReadonlySignal<string | null>;

  /** Per source: how a query a control wrote bound is spelled in the box — `DomainApp.presets`
   * registers `assignee = "<id>"` → `assignee = $me`, so the switch does not paste a uuid there. */
  private static readonly _shownAs = new WeakMap<DomainSource, Map<string, string>>();

  private readonly _input = signal<FilterQueryInput | FilterBuilder | null>(null);
  /** The tree the source's query stands for — what a cancelled change goes back to. */
  private _accepted: FilterGroup = Filters.group('and');
  /** The query text the grammar or the schema refuses, kept in the box beside the tree it did not
   * make; null while the query and the tree agree. */
  private _raw: string | null = null;
  private _rawProblems: FilterProblem[] = [];
  private _syncing = false;
  private _complete: (tree: FilterGroup) => FilterGroup = (tree) => tree;

  constructor(readonly source: DomainSource, private readonly _options: DomainFiltersOptions = {}) {
    super();
    this.mode = _options.mode ?? 'query';
    this.input = this._input;
    this.problem = computed(() => {
      const input = this._input.value;
      if (input instanceof FilterBuilder)
        return DomainFilters.wrong(input);
      if (!(input instanceof FilterQueryInput))
        return null;
      const message = input.problems.value[0]?.message ?? null;
      return message === null || input.query.value === input.raw.value ? null : message;
    });
    this.root.classList.add('u2-domain-filters');
    this.root.dataset.u2 = 'domain-filters';
    let live = true;
    this.own(() => live = false);
    let building = false;
    this.effect(() => {
      source.state.value;
      if (building || this._input.peek() !== null || source.schema.properties.length === 0)
        return;
      building = true;
      void this._schema().then((schema) => {
        if (live)
          this._build(schema);
      }).catch((e) => source.fail(e)).finally(() => building = false);
    });
  }

  /** The table's schema with the platform's values behind it; over another backend (the gallery,
   * the tests) the properties alone. Either way the offer is narrowed the same way. */
  private async _schema(): Promise<FilterSchema> {
    const own = {...this.source.schema, properties: DomainFilters.ordered(this.source.schema.properties)};
    const operators = await this._operators();
    if (!(emsBackends.domain instanceof DgDomainBackend))
      return {...own, operators, values: (prop, query, signal, ctx) => this._under(prop, query, ctx) ?? EMPTY};
    const {values, resolveRef, renderer} = await domainFilterSchema(this.source.table);
    return {...own, operators, resolveRef, renderer, values: (prop, query, signal, ctx) =>
      this._under(prop, query, ctx) ?? values!(prop, query, signal, ctx).then(DomainFilters.held)};
  }

  /** The table's own columns in the schema's order, the system ones after them: a new condition
   * starts on a column people fill in, not on the row id. */
  static ordered<T extends FilterProperty>(properties: T[]): T[] {
    const system = (p: T) => SYSTEM_COLUMNS.some(([name]) => name === p.name);
    return [...properties.filter((p) => !system(p)), ...properties.filter(system)];
  }

  /** The values to pick from, without the empty cell — the null tests are how to ask for that; a
   * reference's arrives spelled as the id "null". */
  static held(items: FilterValueItem[]): FilterValueItem[] {
    return items.filter(({value}) => value !== null && value !== undefined && value !== '' &&
      !(Filters.isRef(value) && value.id === 'null'));
  }

  /** The candidates for an `under` term: the ROWS of the hierarchy it walks — every location,
   * not the handful the column happens to hold — read through the same search `u2-domain-pick`
   * uses, so a caption is shown and the id is what lands in the query. Null for anything else,
   * which leaves the schema's own values in force. */
  private _under(prop: FilterProperty, query: string,
    ctx?: {operator?: string}): Promise<FilterValueItem[]> | null {
    if (ctx?.operator !== UNDER)
      return null;
    const address = prop.ref ?? (prop.name === 'id' ? this.source.table : undefined);
    if (address === undefined)
      return null;
    return DomainPick.search(address, query, {limit: PICK_LIMIT}).then((items) => items.map((item) => ({
      value: {type: address, id: item.id, name: item.name}, label: item.name})));
  }

  /** `under` is registered for every ref and string column — the operator registry cannot know
   * which target is a hierarchy, and "is under" on a column the server would refuse is a dead end.
   * Offered here on a ref column whose target IS one, and on `id` when this table is; the value
   * slot's candidates are the target's rows, which the schema's `values` already answers. A
   * `'basic'` backend is offered nothing it cannot run ({@link NOT_BASIC}); whatever is typed
   * regardless is still refused by name in the status line. */
  private async _operators(): Promise<FilterSchema['operators']> {
    const backend = emsBackends.domain;
    const own = this.source.schema;
    const hierarchies = new Set<string>();
    let basic = false;
    if (backend !== undefined) {
      const refs = new Set(own.properties.map((p) => p.ref).filter((ref): ref is string => ref !== undefined));
      await Promise.all([...refs].map(async (address) => {
        // a target the caller cannot read is simply not offered a subtree filter
        const target = await backend.table(address).catch(() => undefined);
        if (target?.info.hierarchy === true)
          hierarchies.add(address);
      }));
      basic = (await backend.table(this.source.table)).support.filters === 'basic';
    }
    const applies = (prop: FilterProperty) => prop.ref !== undefined ? hierarchies.has(prop.ref) :
      prop.name === 'id' && own.info.hierarchy === true;
    const runs = (prop: FilterProperty, id: string) => {
      const kind = Filters.kindOf(prop);
      return !NOT_BASIC.has(id) && !(kind === Filters.KIND.DATE_TIME && id === '!=') &&
        !(kind === Filters.KIND.BOOL && NULL_TESTS.has(id));
    };
    // a host operator is offered where its table answers it — the registry behind it is the page's, and
    // every table keeps the core operators there; a column its table answers one for offers that set
    // alone, beside the null tests
    const answered = (id: string) => DomainTable.of(this.source)?.operators.all.some((o) => o.id === id) === true;
    return (prop, offered) => {
      const hosted = offered.some((o) => o.host === true && answered(o.id));
      const out = (applies(prop) ? offered : offered.filter((o) => o.id !== UNDER))
        .filter((o) => o.host === true ? answered(o.id) : !hosted || NULL_TESTS.has(o.id));
      return basic ? out.filter((o) => runs(prop, o.id)) : out;
    };
  }

  private _build(schema: FilterSchema): void {
    const source = this.source;
    const input = this.runInScope(() => this.mode === 'builder' ?
      new FilterBuilder({schema, target: 'domain', inline: true}) :
      new FilterQueryInput({schema, target: 'domain', inline: true, name: 'filters',
        placeholder: this._options.placeholder ?? 'Filter…',
        display: (tree) => DomainFilters._shownAs.get(source)?.get(Filters.format(tree)) ?? null}));
    const box = input instanceof FilterQueryInput ? input : null;
    // the builder holds rows still being written: what is applied is the rows that are complete
    if (box === null)
      this._complete = (tree) => Filters.pruneInvalid(tree, schema, 'domain');
    this.effect(() => {
      const q = source.query.value;
      const parsed = typeof q === 'string' && q.trim() !== '' ? Filters.parse(q, schema, 'domain') : null;
      this._accepted = parsed === null ? DomainFilters.treeOf(q, schema) : parsed.root;
      // a query the grammar or the schema refuses would leave the box empty and the user with no
      // way back but the address bar: the text stays, with its problems, for them to fix
      this._raw = parsed !== null && parsed.problems.length > 0 ? q as string : null;
      // the same wording a typed refusal gets: a `?q=` a link carried fails the same way
      this._rawProblems = parsed === null ? [] :
        FilterQueryInput.worded(q as string, schema, parsed.problems);
      this._write(input, this._accepted);
      box?.showText(this._raw, this._rawProblems);
    });
    this.effect(() => {
      const tree = this._complete(input.value.value);
      const raw = box?.raw.value ?? null;
      // a row whose value is wrong — not one still being written — holds the filter applied before it,
      // as a refused text holds it in the box: the rows never answer some other question meanwhile
      if (this._syncing || (raw === this._raw && Filters.equals(tree, this._accepted)) ||
          (input instanceof FilterBuilder && DomainFilters.wrong(input) !== null))
        return;
      void confirmDiscard(source.session, {action: 'change the filter'}).then((ok) => {
        if (ok)
          this._apply(tree);
        else {
          this._write(input, this._accepted);
          box?.showText(this._raw, this._rawProblems);
        }
      });
    });
    this.root.replaceChildren(input.root);
    this._input.value = input;
  }

  /** The first problem of a builder row that is wrong rather than unfinished — a mistyped or
   * out-of-range number, not a value still empty. */
  static wrong(builder: FilterBuilder): string | null {
    return builder.problems.value.find((p) => p.code !== 'missing-value')?.message ?? null;
  }

  /** The tree goes to the source in the form its query has: text stays text. */
  private _apply(tree: FilterGroup): void {
    const q = this.source.query.peek();
    this._accepted = tree;
    this._raw = null;
    this._rawProblems = [];
    this.source.query.value = typeof q === 'string' ? Filters.format(tree) : tree;
  }

  private _write(input: FilterQueryInput | FilterBuilder, tree: FilterGroup): void {
    if (Filters.equals(tree, this._complete(input.value.peek())))
      return;
    this._syncing = true;
    try {
      input.value.value = tree;
    } finally {
      this._syncing = false;
    }
  }

  /** Spells `bound` as `unbound` in the box over `source` — the form a control wrote the query
   * from, kept for the user to read and edit. */
  static showAs(source: DomainSource, bound: string, unbound: string): void {
    let map = DomainFilters._shownAs.get(source);
    if (map === undefined)
      DomainFilters._shownAs.set(source, map = new Map<string, string>());
    map.set(bound, unbound);
  }

  static treeOf(query: string | FilterGroup, schema?: FilterSchema): FilterGroup {
    if (typeof query !== 'string')
      return query;
    return query.trim() === '' ? Filters.group('and') : Filters.parse(query, schema).root;
  }
}
