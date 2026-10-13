import {Subscription} from 'rxjs';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Aggregation, AGGREGATIONS, EnumeratorConfig, PropagatedColumns} from './config';
import {OutputRow} from './enumerate';
import {getStringColumn, NO_FILTER_TAG, nonBlankRows} from './shared';

/** A qnum packs its qualifier into the value and a bigint reads back as a BigInt, so only these add up. */
function isAggregatable(col: DG.Column): boolean {
  return col.type === DG.COLUMN_TYPE.INT || col.type === DG.COLUMN_TYPE.FLOAT;
}

/** Empty when any part is missing: a partial total would pass a range filter as if it were whole. */
function aggregate(values: (number | null)[], agg: Aggregation): number | null {
  if (values.length === 0 || values.some((v) => v == null)) return null;
  const nums = values as number[];
  if (agg === 'multiply') return nums.reduce((a, b) => a * b, 1);
  if (agg === 'min') return Math.min(...nums);
  if (agg === 'max') return Math.max(...nums);
  const sum = nums.reduce((a, b) => a + b, 0);
  return agg === 'sum' ? sum : sum / nums.length;
}

/** Per-step copies of any other type become text: a bool column cannot hold an empty cell, so it would
 * read "false" where a route has no such step, and a bigint column does not take plain numbers. */
const KEPT_TYPES: string[] =
  [DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.FLOAT, DG.COLUMN_TYPE.STRING, DG.COLUMN_TYPE.DATE_TIME];

interface PropagationTable {
  /** Holds every column in `picked`. */
  df: DG.DataFrame;
  picked: PropagatedColumns;
  /** The source row of the i-th template, building block or reagent the enumerator was given. */
  rows: number[];
}

export interface PropagationSnapshot {
  templates: PropagationTable;
  buildingBlocks: PropagationTable;
  reagents: PropagationTable;
  missing: string[];
}

/** Copied when a run starts: the enumerator works on its own copy of the inputs, so reading the tables
 * once it finishes would misalign every value after a row edited meanwhile. */
export function snapshotPropagation(
  config: EnumeratorConfig, tDf: DG.DataFrame, bDf: DG.DataFrame, rDf: DG.DataFrame | null,
): PropagationSnapshot {
  const en = config.enumeration;
  const missing: string[] = [];
  const take = (df: DG.DataFrame, keyCol: string, picked: PropagatedColumns, file: string): PropagationTable => {
    const present: string[] = [];
    for (const n of Object.keys(picked)) {
      if (df.col(n)) present.push(n);
      else missing.push(`Propagated column "${n}" is not in the ${file} file, so it was left out.`);
    }
    return present.length === 0 ? {df, picked: {}, rows: []} : {
      df: df.clone(null, present.map((n) => df.col(n)!.name)),
      picked: Object.fromEntries(present.map((n) => [n, picked[n]])),
      // A reagents file without its SMILES column yields no reagents, so no row is read.
      rows: df.col(keyCol) ? nonBlankRows(getStringColumn(df, keyCol)) : [],
    };
  };
  return {
    templates: take(tDf, en.smarts_col, en.template_propagated_columns, 'reaction templates'),
    buildingBlocks: take(bDf, en.bb_smiles_column, en.bb_propagated_columns, 'building blocks'),
    reagents: rDf ? take(rDf, en.reagent_smiles_column, en.reagent_propagated_columns, 'reagents') :
      {df: DG.DataFrame.create(), picked: {}, rows: []},
    missing,
  };
}

interface PropagationSource extends PropagationTable {
  prefix: string;
  unit: string;
  nth: (k: number) => string;
  /** The rows of `df` an output row drew on, in route order; null for one with no row. */
  sourceRows: (r: OutputRow) => (number | null)[];
}

function columnsFor(rows: OutputRow[], s: PropagationSource): DG.Column[] {
  if (Object.keys(s.picked).length === 0) return [];
  const perRow = rows.map(s.sourceRows);
  const width = perRow.reduce((m, r) => Math.max(m, r.length), 0);
  const out: DG.Column[] = [];
  for (const [name, aggregations] of Object.entries(s.picked)) {
    const src = s.df.col(name)!;
    const kept = KEPT_TYPES.includes(src.type);
    // get() hands back the type's null sentinel (-2147483648 for int), which would add up as a number.
    const all = Array.from({length: src.length}, (_, i) => src.isNone(i) ? null : kept ? src.get(i) : src.getString(i));
    const values = perRow.map((r) => r.map((row) => row == null ? null : all[row]));
    for (const aggregation of isAggregatable(src) ? new Set(aggregations) : []) {
      const c = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, `${s.prefix}_${name}_${aggregation}`,
        values.map((v) => aggregate(v as (number | null)[], aggregation)));
      c.setTag(DG.TAGS.DESCRIPTION,
        `"${name}" of every ${s.unit} of the route, combined by ${aggregation}; empty when any of them has no value.`);
      out.push(c);
    }
    for (let k = 0; k < width; k++) {
      const c = DG.Column.fromList(kept ? src.type : DG.COLUMN_TYPE.STRING, `${s.prefix}_${k + 1}_${name}`,
        values.map((v) => v[k] ?? null));
      if (kept && src.semType) c.semType = src.semType;
      if (kept && src.meta.format) c.meta.format = src.meta.format;
      c.setTag(DG.TAGS.DESCRIPTION, `"${name}" of ${s.nth(k + 1)}.`);
      c.setTag(NO_FILTER_TAG, 'true');
      out.push(c);
    }
  }
  return out;
}

export function propagatedColumns(rows: OutputRow[], snapshot: PropagationSnapshot): DG.Column[] {
  const {templates, buildingBlocks, reagents} = snapshot;
  return [
    ...columnsFor(rows, {...templates, prefix: 'reaction', unit: 'step',
      nth: (k) => `the reaction template step ${k} ran`,
      sourceRows: (r) => r.steps.map((s) => templates.rows[s.templateIndex])}),
    ...columnsFor(rows, {...buildingBlocks, prefix: 'bb', unit: 'building block',
      nth: (k) => `building block ${k} of the route, counting step by step in reactant order`,
      sourceRows: (r) => r.steps.flatMap((s) =>
        s.buildingBlocks.map((p) => p == null ? null : buildingBlocks.rows[p]))}),
    ...columnsFor(rows, {...reagents, prefix: 'reagent', unit: 'reagent',
      nth: (k) => `reagent ${k} of the route, counting step by step in reactant order`,
      sourceRows: (r) => r.steps.flatMap((s) => s.reagents.map((p) => p == null ? null : reagents.rows[p]))}),
  ];
}

export interface PropagatedColumnsPickerOpts {
  tableInput: DG.InputBase<DG.DataFrame | null>;
  /** Shown between the table input and the picker. */
  columnInputs: DG.InputBase<unknown>[];
  /** Shown after the aggregation rows, in the same form. */
  trailingInputs?: DG.InputBase<unknown>[];
  initial: PropagatedColumns;
  tooltip: string;
  onChanged: () => void;
  subs: Subscription[];
}

/** A column picker plus a row of aggregation checkboxes per picked numeric column. `root` holds the whole
 * form, the table and column inputs included: the platform sizes the label column per form. */
export class PropagatedColumnsPicker {
  readonly columnsInput: DG.InputBase<DG.Column[]>;
  readonly root: HTMLElement = ui.div([]);
  /** Source of truth, keyed by the shown table's own column names. Keeps picks the table lacks or cannot
   * aggregate, for a table or a renamed or added column that has them. */
  private picks: PropagatedColumns;
  private tableSubs: Subscription[] = [];
  private aggregationInputs = new Map<string, DG.InputBase<Aggregation[] | null>>();
  private aggregationSubs: Subscription[] = [];

  constructor(private readonly opts: PropagatedColumnsPickerOpts) {
    this.picks = {...opts.initial};
    // Ticked by show(): the `checked` creation option ticks nothing.
    this.columnsInput = ui.input.columns('Propagate columns (optional)', {nullable: true});
    this.columnsInput.setTooltip(opts.tooltip);
    this.bindTable();
    opts.subs.push(
      // onInput, not onChanged: show() sets these inputs in code, and only a user's change updates picks.
      this.columnsInput.onInput.subscribe(() => {
        this.takeCheckedColumns();
        opts.onChanged();
      }),
      // A subset clone or a re-loaded file gets the same picks by name.
      opts.tableInput.onChanged.subscribe(() => this.bindTable()),
      new Subscription(() => [...this.tableSubs, ...this.aggregationSubs].forEach((s) => s.unsubscribe())),
    );
  }

  get value(): PropagatedColumns {
    return {...this.picks};
  }

  set value(picked: PropagatedColumns) {
    this.picks = {...picked};
    this.show();
  }

  private get table(): DG.DataFrame | null {
    return this.opts.tableInput.value;
  }

  private bindTable(): void {
    this.tableSubs.forEach((s) => s.unsubscribe());
    const t = this.table;
    this.tableSubs = t ? [
      // A column renamed or added to match a remembered pick gets ticked.
      t.onColumnNameChanged.subscribe((e: DG.EventData) => {
        const {oldName, newName} = e.args as {oldName: string; newName: string};
        this.picks = Object.fromEntries(
          Object.entries(this.picks).map(([n, a]) => [n === oldName ? newName : n, a]));
        this.show();
      }),
      t.onColumnsChanged.subscribe(() => this.show()),
    ] : [];
    if (t) ui.input.setColumnsInputTable(this.columnsInput, t);
    // An optional file can be cleared, leaving no columns to offer.
    this.columnsInput.enabled = t != null;
    this.show();
  }

  private show(): void {
    const t = this.table;
    // A table finds its columns ignoring case, so a saved `mw` resolves to `MW` and takes its name.
    if (t) this.picks = Object.fromEntries(Object.entries(this.picks).map(([n, a]) => [t.col(n)?.name ?? n, a]));
    this.columnsInput.value = t ?
      Object.keys(this.picks).map((n) => t.col(n)).filter((c): c is DG.Column => c != null) : [];
    this.renderAggregations();
  }

  private takeCheckedColumns(): void {
    const t = this.table;
    const checked = (this.columnsInput.value ?? []).map((c) => c.name);
    const next: PropagatedColumns = {};
    // Unchecking drops a pick with its aggregations; a pick the shown table lacks cannot be unchecked.
    for (const [name, aggs] of Object.entries(this.picks))
      if (checked.includes(name) || !t?.col(name)) next[name] = aggs;
    for (const name of checked) next[name] ??= [];
    this.picks = next;
    this.renderAggregations();
  }

  /** Rebuilds the form only when the aggregation rows change: a rebuild moves every input into a new
   * form element, which costs the keyboard focus. */
  private renderAggregations(): void {
    const names = (this.columnsInput.value ?? []).filter(isAggregatable).map((c) => c.name);
    const shown = [...this.aggregationInputs.keys()];
    if (names.length === shown.length && names.every((n, i) => n === shown[i]) && this.root.firstChild) {
      for (const n of names) this.aggregationInputs.get(n)!.value = this.picks[n] ?? [];
      return;
    }
    this.aggregationSubs.forEach((s) => s.unsubscribe());
    this.aggregationSubs = [];
    this.aggregationInputs.clear();
    for (const name of names) {
      const input = ui.input.multiChoice<Aggregation>(`${name} aggregation`, {items: [...AGGREGATIONS],
        value: this.picks[name] ?? []});
      input.setTooltip(`Each ticked box adds one column combining the route's "${name}" values, to filter on: ` +
        'sum, multiply (yields as fractions: 0.8 × 0.5 = 0.4), avg (the average), min (e.g. the worst step\'s ' +
        'yield) or max (e.g. the most expensive building block). A column is empty when any of the values is ' +
        'missing. Tick none for the individual columns only.');
      input.root.classList.add('chem-enum-aggregation');
      this.aggregationSubs.push(input.onInput.subscribe(() => {
        this.picks[name] = input.value ?? [];
        this.opts.onChanged();
      }));
      this.aggregationInputs.set(name, input);
    }
    const focused = document.activeElement;
    this.root.replaceChildren(ui.form([this.opts.tableInput, ...this.opts.columnInputs, this.columnsInput,
      ...this.aggregationInputs.values(), ...this.opts.trailingInputs ?? []]));
    if (focused instanceof HTMLElement && focused !== document.activeElement && this.root.contains(focused))
      focused.focus();
  }
}
