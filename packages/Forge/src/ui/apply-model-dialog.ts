import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {APPLY_COLUMNS, ApplyModelRow, ApplyRequest, applyWithProgress, DEFAULT_BATCH_SIZE, featureSchemasOf,
  LoadedModel, loadedModelOf} from '../apply/apply-model';
import {ColumnMapping, compatibility, isSuggested, kindText, MappingProblem, mappingProblems, suggestMapping}
  from '../apply/column-matching';
import {applicableTables} from '../catalog/applicable-tables';
import {MINUTE_FORMAT, SECOND_FORMAT} from '../constants';
import {Engine} from '../engines/engine';
import {EngineRegistry} from '../engines/engine-registry';
import {ForgeError} from '../forge-error';
import {forgeDb, ModelRow} from '../generated/db';
import {missingColumnsOf, missingValuesProblems} from '../preparation/missing-values';
import {ColumnSchema} from '../training/train-model';
import {ButtonGate} from './button-gate';
import {CollapsibleGroup} from './collapsible-group';
import {MissingValuesInputs} from './missing-values-inputs';
import {reportError} from './report-error';

const MAX_MODELS = 10000;
const MAX_BATCH_SIZE = 100000;
const ROWS_HEIGHT_SHARE = 0.4;
const MIN_ROWS_HEIGHT = 160;

export interface ApplyDialogOptions {
  table: DG.DataFrame;
  modelId?: string;
  /** Show the table's view after applying (the catalog opens the dialog from another view). */
  switchToTable?: boolean;
  /** With [modelId]: the table to open on while it is open (the catalog's **Applicable to**). */
  preferredTable?: DG.DataFrame | null;
}

export async function applyModelDialog(options: ApplyDialogOptions): Promise<DG.Dialog> {
  const rows: ApplyModelRow[] = await forgeDb.models.query().select(...APPLY_COLUMNS).top(MAX_MODELS);
  if (rows.length === 0)
    throw new ForgeError('No Forge models yet. Train and save one with ML | Forge | Train... first.');
  const preset = rows.find((row) => row.id === options.modelId);
  const table = preset === undefined ? options.table : presetTable(options, preset);
  return new ApplyForm(rows, EngineRegistry.discover(), {...options, table}).dialog;
}

/** The table the dialog opens on for a preset [model]: the preferred table while it is open, else the current table
 * if the model fits it, else the first open table it fits, else the caller's. */
function presetTable(options: ApplyDialogOptions, model: ApplyModelRow): DG.DataFrame {
  const preferred = options.preferredTable;
  if (preferred && isOpen(preferred))
    return preferred;
  const fits = (t: DG.DataFrame | null): t is DG.DataFrame => t !== null && applicableTables(model, [t]).length > 0;
  const current = grok.shell.currentTable;
  return fits(current) ? current : grok.shell.tables.find(fits) ?? options.table;
}

/** Whether [table] is still open: a table input or a caller may hold a closed one. */
export function isOpen(table: DG.DataFrame): boolean {
  return grok.shell.tables.some((t) => t.dart === table.dart);
}

/** Unique list labels of [rows]: the name; with the creation time when names repeat, to the minute, or to the second
 * when the minute repeats too; then ` #2`, ` #3`, ... in creation order. */
export function modelLabels<T extends Pick<ModelRow, 'id' | 'name' | 'created_on'>>(rows: T[]): Map<string, T> {
  const countsOf = (labels: string[]) => {
    const counts = new Map<string, number>();
    for (const label of labels)
      counts.set(label, (counts.get(label) ?? 0) + 1);
    return counts;
  };
  const timed = (row: T, format: string) => `${row.name} (${row.created_on.format(format)})`;
  const names = countsOf(rows.map((row) => row.name));
  const isNameRepeated = (row: T) => (names.get(row.name) ?? 0) > 1;
  const entries = rows.map((row) => ({row, byMinute: isNameRepeated(row) ? timed(row, MINUTE_FORMAT) : row.name}));
  const minutes = countsOf(entries.map((e) => e.byMinute));
  entries.sort((a, b) => a.row.created_on.valueOf() - b.row.created_on.valueOf() || a.row.id.localeCompare(b.row.id));
  const labels = new Map<string, T>();
  for (const {row, byMinute} of entries) {
    const base = isNameRepeated(row) && (minutes.get(byMinute) ?? 0) > 1 ? timed(row, SECOND_FORMAT) : byMinute;
    let label = base;
    for (let i = 2; labels.has(label); i++)
      label = `${base} #${i}`;
    labels.set(label, row);
  }
  return labels;
}

/** Shows the Apply dialog on [table], or "Open a table first." without one: the menu's and the catalog's boundary. */
export async function openApplyDialog(table: DG.DataFrame | null,
  options: Omit<ApplyDialogOptions, 'table'> = {}): Promise<void> {
  try {
    if (table === null)
      grok.shell.warning('Open a table first.');
    else
      (await applyModelDialog({...options, table})).show();
  } catch (e) {
    reportError(e);
  }
}

/** Runs [request] under the task-bar progress, then names the new column in a balloon and, with [switchToTable], shows
 * the table's view: the Apply dialog's OK. Outside the dialog, so the running application holds the request only, not
 * a closed dialog's form and models. */
async function applyAndReport(request: ApplyRequest, switchToTable: boolean): Promise<void> {
  const table = request.table;
  try {
    const {column, skippedRows} = await applyWithProgress(request, 'ui');
    const skipped = skippedRows === 0 ? '' :
      `; ${skippedRows} ${skippedRows === 1 ? 'row' : 'rows'} skipped (missing values)`;
    grok.shell.info(`Added the column "${column.name}" to ${table.name}${skipped}.`);
    const view = switchToTable ? grok.shell.getTableView(table.name) : null;
    if (view)
      grok.shell.v = view;
  } catch (e) {
    reportError(e);
  }
}

class ApplyForm {
  readonly dialog: DG.Dialog;
  /** The models by their list labels, which are unique. */
  private readonly models: Map<string, ApplyModelRow>;
  private readonly engines: Engine[];
  private readonly switchToTable: boolean;
  private readonly tableInput: DG.InputBase<DG.DataFrame | null>;
  private readonly modelInput: DG.ChoiceInput<string | null>;
  private readonly missingValues: MissingValuesInputs;
  private readonly batchInput: DG.InputBase<number | null>;
  private readonly host = ui.div([]);
  private readonly rowsBlock = ui.div([], 'forge-apply-rows');
  private readonly columns = new CollapsibleGroup('Columns', [this.rowsBlock], true, () => {
    this.rowsBlock.scrollTop = 0;
  });
  private readonly moreOptions = new CollapsibleGroup('More options', [], false);
  // One input per feature name, reused across models: the dialog keeps every input it was given.
  private readonly rowInputs = new Map<string, DG.InputBase<DG.Column | null>>();
  private readonly subs: rxjs.Subscription[] = [];
  private readonly okGate: ButtonGate;
  private table: DG.DataFrame;
  private tableSubs: rxjs.Subscription[] = [];
  private model: LoadedModel | null = null;
  private modelError: string | null = null;
  private problems: MappingProblem[] = [];
  private blocker: string | null = null;
  private isUpdating = false;

  constructor(rows: ApplyModelRow[], engines: Engine[], options: ApplyDialogOptions) {
    this.engines = engines;
    this.table = options.table;
    this.switchToTable = options.switchToTable ?? false;
    this.models = modelLabels(rows);
    const order = this.order();
    const preset = [...this.models].find(([, row]) => row.id === options.modelId)?.[0];

    this.dialog = ui.dialog('Apply predictive model');
    this.tableInput = ui.input.table('Table', {items: grok.shell.tables, value: this.table, nullable: false,
      tooltipText: 'Table to add the prediction to.', onValueChanged: (t) => this.changeTable(t)});
    this.modelInput = ui.input.choice('Model', {items: order, value: preset ?? order[0], nullable: false,
      tooltipText: 'Saved model to apply.', onValueChanged: () => this.render()});
    this.modelInput.addValidator(() => this.modelError);
    this.missingValues = new MissingValuesInputs(() => this.update(), () => this.missingValuesProblem());
    this.batchInput = ui.input.int('Batch size', {value: DEFAULT_BATCH_SIZE, min: 1, max: MAX_BATCH_SIZE,
      nullable: false, tooltipText: 'Rows predicted in one step. Lower it for heavy methods.',
      onValueChanged: () => this.update()});
    const rowsHeight = Math.max(MIN_ROWS_HEIGHT, Math.round(window.innerHeight * ROWS_HEIGHT_SHARE));
    this.rowsBlock.style.maxHeight = `${rowsHeight}px`;

    this.dialog.onOK(() => this.apply());
    this.okGate = new ButtonGate(this.dialog.getButton('OK'), () => this.blocker);
    this.dialog.onClose.subscribe(() => {
      this.unsubscribeTable();
      for (const sub of this.subs)
        sub.unsubscribe();
    });
    this.subscribeTable();
    // The rows are given to the dialog before the other inputs: on show() the dialog focuses the input it was given
    // last, and a focused row would scroll the rows block down. The last one is Batch size, folded and unfocusable.
    this.render();
    this.dialog.add(this.tableInput).add(this.modelInput).add(this.host);
    for (const input of this.missingValues.inputs)
      this.dialog.add(input);
    // Added to the dialog first, so the dialog validates it, then moved into More options.
    this.dialog.add(this.batchInput).add(this.moreOptions.root);
    this.moreOptions.body.append(this.batchInput.root);
    this.subs.push(...this.moreOptions.expandOnError([this.batchInput]));
  }

  /** Labels with the models that fit the table first, newest first in each group. */
  private order(): string[] {
    // Retrained models share their feature lists: each list is matched against the table once.
    const fitting = new Map<string, boolean>();
    const entries = [...this.models].map(([label, row]) => {
      const features = featureSchemasOf(row.features);
      const key = JSON.stringify(features);
      let fits = fitting.get(key);
      if (fits === undefined) {
        fits = features !== null && isSuggested(features, this.table);
        fitting.set(key, fits);
      }
      return {label, created: row.created_on.valueOf(), isSuggested: fits};
    });
    entries.sort((a, b) => a.isSuggested === b.isSuggested ? b.created - a.created : a.isSuggested ? -1 : 1);
    return entries.map((e) => e.label);
  }

  private changeTable(table: DG.DataFrame | null): void {
    if (this.isUpdating || table === null)
      return;
    this.table = table;
    this.subscribeTable();
    const model = this.modelInput.value;
    this.isUpdating = true;
    try {
      this.modelInput.items = this.order();
      this.modelInput.value = model;
    } finally {
      this.isUpdating = false;
    }
    this.render();
  }

  /** Rebuilds the feature rows for the chosen model, prefilled with the suggested columns; **Columns** starts
   * collapsed when every feature has a valid column, expanded otherwise. */
  private render(): void {
    if (this.isUpdating)
      return;
    this.isUpdating = true;
    try {
      const row = this.models.get(this.modelInput.value ?? '');
      this.model = null;
      this.modelError = null;
      try {
        this.model = row === undefined ? null : loadedModelOf(row, this.engines);
      } catch (e) {
        if (!(e instanceof ForgeError))
          throw e;
        this.modelError = e.message;
      }
      for (const input of this.rowInputs.values()) {
        input.nullable = true;
        input.value = null;
      }
      ui.empty(this.host);
      ui.empty(this.rowsBlock);
      const model = this.model;
      if (model === null)
        this.host.append(ui.divText(this.modelError ?? ''));
      else {
        const prefill = suggestMapping(model.features, this.table);
        for (const feature of model.features)
          this.rowsBlock.append(this.rowInput(model, feature, prefill.get(feature.name)).root);
        this.host.append(this.columns.root);
      }
    } finally {
      this.isUpdating = false;
    }
    this.update();
    this.columns.setExpanded(this.problems.length > 0);
    this.rowsBlock.scrollTop = 0;
  }

  private rowInput(model: LoadedModel, feature: ColumnSchema, prefill: string | undefined):
    DG.InputBase<DG.Column | null> {
    const filter = (c: DG.Column) => compatibility(feature, c, model.engine.name).kind !== 'error';
    let input = this.rowInputs.get(feature.name);
    if (input !== undefined) {
      input.nullable = false;
      ui.input.setColumnInputTable(input, this.table, filter);
    } else {
      input = ui.input.column(feature.name, {table: this.table, filter, nullable: false,
        onValueChanged: () => this.update()});
      input.addValidator(() => this.problems.find((p) => p.feature === feature.name)?.message ?? null);
      this.dialog.add(input);
      this.rowInputs.set(feature.name, input);
      this.subs.push(...this.columns.expandOnError([input]));
    }
    input.value = prefill === undefined ? null : this.table.getCol(prefill);
    return input;
  }

  /** After any change: problems, tooltips, the Missing values block, validation marks and the OK button. */
  private update(): void {
    if (this.isUpdating)
      return;
    this.isUpdating = true;
    try {
      const model = this.model;
      const mapping = this.mapping();
      this.problems = model === null ? [] : mappingProblems(model.features, mapping, this.table, model.engine.name);
      const featureCount = model?.features.length ?? 0;
      const unmatched = new Set(this.problems.map((p) => p.feature)).size;
      this.columns.setSummary(`${featureCount - unmatched} of ${featureCount} matched`, unmatched > 0);
      const missing = missingColumnsOf(this.mappedColumns(mapping));
      this.missingValues.update(missing);
      const gaps = new Map(missing.map((m) => [m.name, m.count]));
      for (const feature of model?.features ?? [])
        this.rowInputs.get(feature.name)?.setTooltip(this.rowTooltip(feature, mapping, gaps));
      const inputs = this.visibleInputs();
      for (const input of inputs)
        input.validate();
      this.blocker = model === null ? this.modelError ?? 'Choose a model.' : this.problems[0]?.message ??
        inputs.map((input) => input.validity).find((v): v is string => v !== null) ?? null;
    } finally {
      this.isUpdating = false;
    }
    this.okGate.update();
  }

  private mapping(): ColumnMapping {
    const mapping: ColumnMapping = new Map();
    for (const feature of this.model?.features ?? []) {
      const col = this.rowInputs.get(feature.name)?.value;
      if (col !== null && col !== undefined)
        mapping.set(feature.name, col.name);
    }
    return mapping;
  }

  private mappedColumns(mapping: ColumnMapping): DG.Column[] {
    return [...mapping.values()].map((name) => this.table.col(name)).filter((c): c is DG.Column => c !== null);
  }

  private visibleInputs(): DG.InputBase[] {
    const rows = (this.model?.features ?? []).map((f) => this.rowInputs.get(f.name))
      .filter((input): input is DG.InputBase<DG.Column | null> => input !== undefined);
    return [this.tableInput, this.modelInput, ...rows, ...this.missingValues.visibleInputs, this.batchInput];
  }

  /** [gaps] holds the missing-value counts of the mapped columns that have any. */
  private rowTooltip(feature: ColumnSchema, mapping: ColumnMapping, gaps: Map<string, number>): string {
    const parts = [`Table column used as the feature '${feature.name}' (${kindText(feature.type)}).`];
    const name = mapping.get(feature.name);
    const col = name === undefined ? null : this.table.col(name);
    if (col !== null) {
      const count = gaps.get(col.name) ?? 0;
      if (count > 0)
        parts.push(`'${col.name}' has ${count} missing ${count === 1 ? 'value' : 'values'}.`);
      const fit = compatibility(feature, col, this.model?.engine.name ?? '');
      if (fit.kind === 'hint')
        parts.push(fit.message);
    }
    return parts.join(' ');
  }

  private missingValuesProblem(): string | null {
    if (!this.missingValues.isImpute)
      return null;
    const problems = missingValuesProblems(this.mappedColumns(this.mapping()), this.missingValues.settings());
    return problems.length > 0 ? problems.join(' ') : null;
  }

  /** The OK handler: returns at once, so the dialog closes, and leaves the application to the task bar. */
  private apply(): void {
    const model = this.model;
    if (model === null || this.blocker !== null) {
      grok.shell.warning(this.blocker ?? 'Choose a model.');
      return;
    }
    void applyAndReport({model, table: this.table, mapping: this.mapping(),
      batchSize: this.batchInput.value ?? DEFAULT_BATCH_SIZE, missingValues: this.missingValues.settings()},
    this.switchToTable);
  }

  private subscribeTable(): void {
    this.unsubscribeTable();
    // Through an input's change event (its handler is update()), so the dialog also re-checks the inputs Enter needs.
    const update = () => this.batchInput.fireChanged();
    this.tableSubs = [this.table.onColumnsRemoved.subscribe(update), this.table.onColumnsAdded.subscribe(update)];
  }

  private unsubscribeTable(): void {
    for (const sub of this.tableSubs)
      sub.unsubscribe();
    this.tableSubs = [];
  }
}
