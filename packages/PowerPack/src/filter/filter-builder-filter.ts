/* The `role: filter` host for u2's FilterBuilder: one filter over the whole frame (the column it
   was picked for, if any, stays its `columnName`), evaluated to a BitSet on every tree change and
   AND-ed into `df.filter` collaboratively. Only the complete rows are evaluated — a row being
   edited never un-filters the frame. The Dart proxy stamps `type`/`column`/`active` on the saved
   state; `model` (the typed tree as JSON) and `query` (canonical string) are ours. */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';
import {signal, Scope, FilterBuilder, FilterQueryInput, Filters, OperatorRegistry} from '@datagrok-libraries/u2';
import type {FilterGroup as FilterTree, FilterOperatorSet} from '@datagrok-libraries/u2';
import {FilterSchemas, toBitSet} from '@datagrok-libraries/u2/src/dg/index.js';
import type {DataFrameFilterSchema} from '@datagrok-libraries/u2/src/dg/index.js';

const NOT_FILTERING = 'Not filtering: no complete condition';

/** The operator sets of the `meta.role: filterOperators` package functions (Chem's Molecule
 * operators): discovered, validated and registered once per page. */
class OperatorSets {
  private static _discovery: Promise<void> | null = null;
  static done = false;

  static discover(): Promise<void> {
    return OperatorSets._discovery ??= OperatorSets._register()
      .catch((e: any) => console.warn(`Filter builder: operator discovery failed: ${e?.message ?? e}`))
      .finally(() => OperatorSets.done = true);
  }

  private static async _register(): Promise<void> {
    for (const f of DG.Func.find({meta: {role: 'filterOperators'}})) {
      try {
        const set: FilterOperatorSet = await f.apply();
        const problems = OperatorRegistry.checkSet(set);
        if (problems.length > 0)
          throw new Error(problems.join('; '));
        Filters.operators.registerSet(set);
      } catch (e: any) {
        console.warn(`Filter builder: operators from ${f.nqName} skipped: ${e?.message ?? e}`);
      }
    }
  }
}

export class FilterBuilderFilter extends DG.Filter {
  builder: FilterBuilder | null = null;
  queryInput: FilterQueryInput | null = null;
  bitset: DG.BitSet | null = null;
  /** The line under the rows — the complete rows' query, or why nothing filters yet; the panel's
   * own header summary says the same, so it is off unless the state asks. */
  showStatus = false;
  private _tree = signal<FilterTree>(Filters.group('and'));
  private _scope: Scope | null = null;
  private _schema: DataFrameFilterSchema | null = null;
  private _showQuery = false;
  private _subs: Subscription[] = [];
  private _run = 0;
  private _abort: AbortController | null = null;
  private _loader = ui.loader();
  private _status = ui.divText('', 'power-pack-filter-builder-status');

  constructor() {
    super();
    this.root.classList.add('power-pack-filter-builder');
    this._loader.style.display = 'none';
  }

  get caption(): string { return 'Filter builder'; }

  /** The complete rows; explicit about rows that do not filter yet. */
  get filterSummary(): string {
    const valid = this._valid();
    return Filters.count(valid) > 0 ? Filters.format(valid) :
      Filters.count(this._tree.peek()) > 0 ? NOT_FILTERING : '';
  }

  get isFiltering(): boolean { return super.isFiltering && this._hasQuery(); }
  get isReadyToApplyFilter(): boolean { return this.bitset != null; }

  attach(dataFrame: DG.DataFrame): void {
    super.attach(dataFrame);
    this._schema = FilterSchemas.forDataFrame(dataFrame);
    this._mount();
    // the columns change under the filter: a semantic type detected after the snapshot (Chem
    // answers asynchronously), a tag set by hand, a column added or renamed
    for (const event of [dataFrame.onSemanticTypeDetected, dataFrame.onMetadataChanged, dataFrame.onColumnsChanged])
      this._subs.push(event.subscribe(() => this._refreshSchema()));
    // rows built before the provider sets landed derive their operators again
    if (!OperatorSets.done) {
      OperatorSets.discover().then(() => {
        if (this._scope)
          this._mount();
      });
    }
  }

  applyFilter(): void {
    const df = this.dataFrame;
    if (!df || !this.bitset)
      return;
    if (this.bitset.length !== df.filter.length) {
      this._recompute(this._tree.peek());
      return;
    }
    df.filter.and(this.bitset, false);
    df.rows.addFilterState(this.saveState());
  }

  saveState(): any {
    const tree = this._tree.peek();
    return {...super.saveState(), model: Filters.toJson(tree), query: Filters.format(tree), showStatus: this.showStatus};
  }

  /** The model is lossless (a semType operator such as Chem's `Contains` has no grammar spelling);
   * the query is the form states written by hand come in. An operator not registered yet stays on
   * its row and is validated again when the schema refreshes. */
  applyState(state: any): void {
    super.applyState(state);
    this.showStatus = state.showStatus ?? this.showStatus;
    this._status.hidden = !this.showStatus;
    if (state.model)
      this._tree.value = Filters.fromJson(state.model);
    else if (state.query)
      this._tree.value = Filters.parse(state.query, this._schema ?? undefined, 'dataframe').root;
    else if (state.columnName && Filters.count(this._tree.peek()) === 0)
      this._seed(state.columnName);
    this.columnName = state.columnName ?? '';
  }

  detach(): void {
    this._run++;
    this._abort?.abort();
    for (const sub of this._subs.splice(0))
      sub.unsubscribe();
    super.detach();
    this.bitset = null;
    this._scope?.dispose();
    this._scope = null;
    this.builder = null;
    this.queryInput = null;
  }

  /** The editors over the schema as it is now; the tree — the rows — is theirs to show. */
  private _mount(): void {
    this._scope?.dispose();
    const scope = this._scope = new Scope();
    const common = {schema: this._schema!, bind: this._tree, inline: true, target: 'dataframe' as const};
    const builder = this.builder = Scope.runWith(scope, () => new FilterBuilder(common));
    const queryInput = this.queryInput = Scope.runWith(scope, () => new FilterQueryInput(common));
    const show = () => {
      builder.root.hidden = this._showQuery;
      queryInput.root.hidden = !this._showQuery;
    };
    show();
    const toggle = ui.iconFA('code', () => {
      this._showQuery = !this._showQuery;
      show();
    }, 'Switch between the builder and the query text');
    toggle.dataset.u2Part = 'query-toggle';
    this._status.hidden = !this.showStatus;
    ui.empty(this.root);
    this.root.append(builder.root, queryInput.root,
      ui.divH([toggle, this._loader, this._status], 'power-pack-filter-builder-bar'));
    scope.effect(() => this._recompute(this._tree.value));
  }

  /** Rebuilds the editors when a column's shape changed: the rows keep their property names, the
   * operators and value editors derive again — a molecule column gains its sketcher. */
  private _refreshSchema(): void {
    const schema = this._schema;
    if (!schema || !this.dataFrame)
      return;
    const shape = (props: {name: string, type?: string, semType?: string}[]) =>
      props.map((p) => `${p.name}\t${p.type}\t${p.semType ?? ''}`).join('\n');
    const before = shape(schema.properties);
    schema.refresh();
    if (shape(schema.properties) !== before)
      this._mount();
  }

  /** The column-picker entry path: the picked column seeds the first condition, the filter
   * itself stays columnless. */
  private _seed(columnName: string): void {
    const prop = this._schema?.properties.find((p) => p.name === columnName);
    const op = prop ? Filters.operators.for(prop)[0] : undefined;
    if (prop && op)
      this._tree.value = Filters.group('and', [Filters.cond(prop.name, op.id)]);
  }

  /** The tree without its incomplete rows — what gets evaluated. */
  private _valid(): FilterTree {
    const tree = this._tree.peek();
    return this._schema ? Filters.pruneInvalid(tree, this._schema, 'dataframe') : Filters.group('and');
  }

  private _hasQuery(): boolean {
    return this.builder != null && Filters.count(this._valid()) > 0;
  }

  private async _recompute(tree: FilterTree): Promise<void> {
    const run = ++this._run;
    this._abort?.abort();
    const had = this.bitset != null;
    this.bitset = null;
    const df = this.dataFrame;
    if (!df)
      return;
    const valid = this._valid();
    this._status.textContent = this.filterSummary;
    if (!this._hasQuery()) {
      if (had)
        df.rows.requestFilter();
      return;
    }
    const abort = this._abort = new AbortController();
    this._loader.style.display = '';
    try {
      const bitset = await toBitSet(df, valid, {signal: abort.signal});
      if (run !== this._run)
        return;
      this.bitset = bitset;
      df.rows.requestFilter();
    } catch (e: any) {
      if (run === this._run)
        grok.shell.warning(`Filter builder: ${e?.message ?? e}`);
    } finally {
      if (run === this._run)
        this._loader.style.display = 'none';
    }
  }
}
