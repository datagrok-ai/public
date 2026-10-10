import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {featureFit} from '../catalog/applicable-tables';
import {ModelActivity, modelActivity} from '../catalog/model-activity';
import {MINUTE_FORMAT, MODEL_TYPE, SECOND_FORMAT} from '../constants';
import {ApplicationRow, ApplicationSource, ApplicationStatus, forgeDb, ModelRow} from '../generated/db';
import {METRIC_DESCRIPTIONS, METRIC_IDS, METRIC_LABELS, MetricId} from '../metrics/metrics';
import {preparationOptionsOf} from '../preparation/preparation-options';
import {datasetRefCaption, storedDatasetRef} from '../storage/dataset-ref';
import {tagsOf, tagsText} from '../storage/model-fields';
import {MetricsRecord, metricsRecordOf} from '../training/train-model';
import {readOnlyGrid, textColumn} from './data-grid';
import {reportError} from './report-error';
import {STORAGE_CAPTIONS, writeModelInfo} from './save-model-dialog';
import {tagsInput, tagsOfInput} from './tags-input';

const NO_SHARE = 'Sharing a model needs the Share permission on this model; ask an administrator.';
const METRIC_HEADERS: {[column: string]: string} = {
  'Metric': 'Quality measure of the model.',
  'Train': 'Value on the training rows, from the model trained on all of them.',
  'Validation': 'Value on rows the model did not see: the pooled out-of-fold predictions of 5-fold cross-validation.',
};
const ACTIVITY_HEADERS: {[column: string]: string} = {
  'When': 'Time the model was applied.',
  'Who': 'User who applied the model.',
  'Table': 'Table the model was applied to.',
  'Rows': 'Number of rows in that table.',
  'Prediction column': 'Column the application added; empty when it added none.',
  'Status': 'Outcome of the application.',
  'Source': 'Where the application was started.',
  'Duration (ms)': 'Time the application took, in milliseconds.',
};
const STATUS_TOOLTIPS: Record<ApplicationStatus, string> = {
  completed: 'Completed: the column was added.',
  failed: 'Failed: ',
  cancelled: 'Cancelled: stopped between batches, no column added.',
};
const SOURCE_TOOLTIPS: Record<ApplicationSource, string> = {
  ui: 'The Apply dialog or the catalog.',
  api: 'A script, through Forge:applyModel.',
};

let sharingSub: rxjs.Subscription | undefined;

type Shares = {view?: DG.Group[]; edit?: DG.Group[]};

export interface MetricsSummary {
  metrics: MetricsRecord;
  rowCount: number;
  skippedRows: number;
  folds?: number;
  seed?: number;
}

type Split = 'train' | 'validation';
const SPLIT_COLUMNS: [string, Split][] = [['Train', 'train'], ['Validation', 'validation']];

/** The [split]'s value of each of [ids], null where it has none. */
function splitValues(metrics: MetricsRecord, split: Split, ids: MetricId[]): (number | null)[] {
  return ids.map((id) => metrics[split][id] ?? null);
}

/** The metrics [ids], one row each: **Metric** (the label), **Train**, **Validation**. */
function metricsFrame(metrics: MetricsRecord, ids: MetricId[]): DG.DataFrame {
  const values = SPLIT_COLUMNS.map(([name, split]) => {
    const col = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, splitValues(metrics, split, ids));
    col.meta.format = '0.000';
    return col;
  });
  return DG.DataFrame.fromColumns([textColumn('Metric', ids.map((id) => METRIC_LABELS[id])), ...values]);
}

/** The metrics grid (hover a header or a metric for its meaning), then a bullet list of the rows, the validation, the
 * seed and the positive class: the Train view's Results and the model's Performance pane. {@link update} shows
 * another summary of the same metric rows in the same grid, which keeps its column widths and scroll. */
export class MetricsTable {
  readonly root: HTMLElement;
  private readonly ids: MetricId[];
  private readonly grid: DG.Grid;
  private readonly list = ui.element('ul');
  private summary: MetricsSummary;

  constructor(summary: MetricsSummary) {
    this.summary = summary;
    this.ids = MetricsTable.idsOf(summary);
    // A value in full precision, with what the column measures; read from the current summary.
    const valueTooltip = (column: string, split: Split) => (row: number) => {
      const value = this.summary.metrics[split][this.ids[row]];
      return value === undefined ? '' : `${value}. ${METRIC_HEADERS[column]}`;
    };
    this.grid = readOnlyGrid(metricsFrame(summary.metrics, this.ids), {headers: METRIC_HEADERS, cells: {
      'Metric': (row) => METRIC_DESCRIPTIONS[this.ids[row]],
      'Train': valueTooltip('Train', 'train'),
      'Validation': valueTooltip('Validation', 'validation'),
    }});
    this.root = ui.divV([this.grid.root, this.list]);
    this.fillList();
  }

  /** Shows [summary] and returns true when it has the metric rows of the grid; false (nothing changed) when not. */
  update(summary: MetricsSummary): boolean {
    const ids = MetricsTable.idsOf(summary);
    if (ids.length !== this.ids.length || ids.some((id, i) => id !== this.ids[i]))
      return false;
    this.summary = summary;
    for (const [name, split] of SPLIT_COLUMNS) {
      const values = splitValues(summary.metrics, split, ids);
      this.grid.dataFrame.getCol(name).init((i) => values[i]);
    }
    this.fillList();
    return true;
  }

  private static idsOf({metrics}: MetricsSummary): MetricId[] {
    return METRIC_IDS.filter((id) => metrics.validation[id] !== undefined);
  }

  private fillList(): void {
    const {metrics, rowCount, skippedRows, folds, seed} = this.summary;
    const items: (string | HTMLElement)[][] = [];
    if (skippedRows > 0)
      items.push([`Rows: ${rowCount} used, ${skippedRows} skipped (missing values)`]);
    if (folds !== undefined && seed !== undefined)
      items.push([`Validation: ${folds}-fold cross-validation on ${rowCount} rows`], seedLine(seed));
    if (metrics.positiveClass !== undefined)
      items.push([`Positive class: ${metrics.positiveClass}`]);
    ui.empty(this.list);
    for (const item of items) {
      const li = ui.element('li');
      li.append(...item);
      this.list.append(li);
    }
    ui.setDisplay(this.list, items.length > 0);
  }
}

/** `Seed: <seed>` with the seed selectable and a copy icon right after it. */
function seedLine(seed: number): (string | HTMLElement)[] {
  const copy = ui.icons.copy(() => {
    navigator.clipboard.writeText(`${seed}`)
      .then(() => grok.shell.info(`Seed ${seed} copied.`))
      .catch((e) => reportError(e));
  }, 'Copy the seed');
  return ['Seed: ', ui.span([`${seed}`], 'forge-seed'), ' ', copy];
}

/** The model icon of the built-in tool. */
export function modelIcon(): HTMLElement {
  return ui.iconSvg('model');
}

/** The model's context panel: the title the platform gives a row (icon, favorites star, name, context actions), then
 * **Details** (open), **Performance**, **Activity**, **Sharing** and **History**, each built when first opened. [row]
 * is the shown row, [model] its full values. */
export function modelAccordion(row: DG.DomainRow, model: ModelRow): DG.Accordion {
  let activity: Promise<ModelActivity> | undefined;
  const activityOf = () => activity ??= modelActivity(model.id);
  const accordion = ui.accordion(MODEL_TYPE);
  accordion.addTitle(ui.span([modelIcon(), ui.star(row.id), ui.label(model.name), ui.contextActions(row)]));
  accordion.addPane('Details', () => detailsPane(model, activityOf), true);
  accordion.addPane('Performance', () => performancePane(model));
  accordion.addPane('Activity', () => activityPane(activityOf));
  accordion.addPane('Sharing', () => sharingPane(row));
  accordion.addPane('History', () => DG.DomainObjectHandler.auditPane(row));
  return accordion;
}

function detailsPane(model: ModelRow, activity: () => Promise<ModelActivity>): HTMLElement {
  const host = ui.divV([]);
  fillDetails(host, model, activity);
  return host;
}

function fillDetails(host: HTMLElement, model: ModelRow, activity: () => Promise<ModelActivity>): void {
  const {names: features, tables} = featureFit(model, grok.shell.tables);
  const rows = model.row_count === undefined ? '' : ` (${model.row_count} rows)`;
  const details: {[caption: string]: string | HTMLElement} = {
    'Author': ui.wait(async () => ui.render(await grok.dapi.users.find(model.author_id))),
    'Created': model.created_on.format(MINUTE_FORMAT),
    'Updated': model.updated_on.format(MINUTE_FORMAT),
    'Table': `${model.dataset_name ?? ''}${rows}`,
    ...storageDetails(model),
    'Last run': ui.wait(async () => ui.divText(lastRunText((await activity()).lastRun))),
    'Applications': ui.wait(async () => ui.divText(`${(await activity()).count}`)),
    'Features': features.join(', '),
    'Target': model.target_name,
    'Method': model.engine_name,
    'Task': model.task,
  };
  if (tables.length > 0)
    details['Applicable to'] = tables.map((t) => t.name).join(', ');
  ui.empty(host);
  host.append(...(model.description ? [ui.divText(model.description)] : []), ui.tableFromMap(details),
    ui.form([detailsTags(host, model, activity)]));
}

/** **Data storage**, then **Data source** of a reference or **Data copy** (the uploaded table, or `missing`). */
function storageDetails(model: ModelRow): {[caption: string]: string | HTMLElement} {
  const details: {[caption: string]: string | HTMLElement} = {'Data storage': STORAGE_CAPTIONS[model.storage_mode]};
  const ref = storedDatasetRef(model.dataset_ref);
  if (ref !== null)
    details['Data source'] = datasetRefCaption(ref);
  const tableId = model.dataset_table_id;
  if (tableId) {
    details['Data copy'] = ui.wait(async () => {
      const info: DG.TableInfo | undefined = await grok.dapi.tables.find(tableId);
      return info ? ui.render(info) : ui.divText('missing');
    });
  }
  return details;
}

/** The editable **Tags** of Details: every change is written; a model changed elsewhere meanwhile asks to reload or
 * overwrite. */
function detailsTags(host: HTMLElement, model: ModelRow, activity: () => Promise<ModelActivity>): DG.InputBase {
  const input = tagsInput(tagsOf(model.tags));
  let version = model.version;
  // The input reports a change while it fills itself too: only tags that differ from the stored ones are written.
  let stored = tagsText(tagsOf(model.tags));
  const reload = async (): Promise<void> => {
    const fresh: ModelRow | null = await forgeDb.models.get(model.id);
    if (fresh !== null)
      fillDetails(host, fresh, activity);
  };
  const write = async (): Promise<void> => {
    const tags = tagsOfInput(input);
    if (tagsText(tags) === stored)
      return;
    const written = await writeModelInfo(model, version, {tags}, reload);
    if (written !== null) {
      version = written;
      stored = tagsText(tags);
    }
  };
  // One write at a time, each with the version the previous one returned: quick chip changes are not a conflict.
  let writes = Promise.resolve();
  input.onChanged.subscribe(() => {
    writes = writes.then(write).catch((e) => reportError(e));
  });
  return input;
}

function lastRunText(lastRun: ModelActivity['lastRun']): string {
  if (lastRun === undefined)
    return 'Never';
  const when = lastRun.when.format(MINUTE_FORMAT);
  return lastRun.status === 'completed' ? when : `${when} (${lastRun.status})`;
}

function performancePane(model: ModelRow): HTMLElement {
  const metrics = metricsRecordOf(model.metrics);
  if (metrics === null)
    return ui.divText('No metrics were recorded for this model.');
  const folds: unknown = model.splitting?.folds;
  return new MetricsTable({metrics, rowCount: model.row_count ?? 0, seed: model.seed,
    skippedRows: preparationOptionsOf(model.options).missingValues?.skippedRows ?? 0,
    folds: typeof folds === 'number' ? folds : undefined}).root;
}

function activityPane(activity: () => Promise<ModelActivity>): HTMLElement {
  return ui.wait(async () => {
    const {applications, count} = await activity();
    if (count === 0)
      return ui.divText('Not applied yet.');
    // One read of the authors; a deleted user is not found and shows no login.
    const users = await grok.dapi.getEntities([...new Set(applications.map((a) => a.author_id))]);
    const logins = new Map(users.filter((u): u is DG.User => u instanceof DG.User).map((u) => [u.id, u.login]));
    const at = (row: number) => applications[row];
    const grid = readOnlyGrid(activityFrame(applications, logins), {headers: ACTIVITY_HEADERS, cells: {
      'When': (row) => at(row).created_on.format(SECOND_FORMAT),
      'Status': (row) => `${STATUS_TOOLTIPS[at(row).status]}${at(row).status === 'failed' ? at(row).error ?? '' : ''}`,
      'Source': (row) => SOURCE_TOOLTIPS[at(row).source],
    }});
    return ui.divV([ui.divText(`${count} ${count === 1 ? 'application' : 'applications'}`), grid.root]);
  });
}

/** One row per application: **When**, **Who**, **Table**, **Rows**, **Prediction column**, **Status**, **Source**,
 * **Duration (ms)**. */
function activityFrame(applications: ApplicationRow[], logins: Map<string, string>): DG.DataFrame {
  return DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 'When', applications.map((a) => a.created_on)),
    textColumn('Who', applications.map((a) => logins.get(a.author_id))),
    textColumn('Table', applications.map((a) => a.table_name)),
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'Rows', applications.map((a) => a.row_count)),
    textColumn('Prediction column', applications.map((a) => a.column_name)),
    textColumn('Status', applications.map((a) => a.status)),
    textColumn('Source', applications.map((a) => a.source)),
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'Duration (ms)', applications.map((a) => a.duration_ms ?? null)),
  ]);
}

/** **Sharing**: the groups the model was shared with and **Share...**, read again whenever the platform reports the
 * model shared (the sharing dialog's OK). */
function sharingPane(row: DG.DomainRow): HTMLElement {
  const host = ui.divV([ui.loader()]);
  const entity = grok.dapi.getEntities([row.id]).then((found) => found[0] ?? null);
  refreshSharing(host, row, entity).catch((e) => reportError(e));
  // The context panel shows one model: the newest Sharing pane is the one to refresh.
  sharingSub?.unsubscribe();
  sharingSub = grok.events.onEntityShared.subscribe((shared: DG.Entity | null) => {
    if (shared?.id === row.id && host.isConnected)
      refreshSharing(host, row, entity).catch((e) => reportError(e));
  });
  return host;
}

/** Fills [host] with the Sharing content of the model and its [entity]. The platform reads the shares of an entity
 * through its wrapping project, where Share... writes them; the author's own grants sit on the row and are not
 * listed. */
export async function refreshSharing(host: HTMLElement, row: DG.DomainRow, entity: Promise<DG.Entity | null>):
  Promise<void> {
  const target = await entity;
  const shares: Promise<{shares: Shares} | {error: unknown}> = target === null ? Promise.resolve({shares: {}}) :
    grok.dapi.permissions.get(target).then((found: Shares) => ({shares: found}), (error: unknown) => ({error}));
  // One column: the system columns are not selectable and come along anyway.
  const [[access], read] = await Promise.all([forgeDb.models.query({filter: {property: 'id', operator: '=',
    value: row.id}, columns: ['name'], withAccess: true, limit: 1}), shares]);
  const canShare = access?.['~can_share'] === true;
  // Only a user without Share is expected to be refused the shares: he sees the limitation line instead of them.
  if ('error' in read && canShare)
    throw read.error;
  const lines: HTMLElement[] = [];
  if ('shares' in read) {
    const view = read.shares.view ?? [];
    const edit = read.shares.edit ?? [];
    const groups = (list: DG.Group[]) => ui.divV(list.map((g) => ui.render(g)));
    lines.push(view.length + edit.length === 0 ?
      ui.divText('Not shared yet. Only its author and administrators can see this model.') :
      ui.tableFromMap({'Can view': groups(view), 'Can edit': groups(edit)}));
  }
  // shareRow can refuse before any dialog opens (the promotion check), which the platform does not report.
  lines.push(ui.button('Share...', async () => {
    try {
      await DG.DomainObjectHandler.shareRow(row);
    } catch (e) {
      reportError(e);
    }
  }));
  if (!canShare)
    lines.push(ui.divText(NO_SHARE));
  ui.empty(host);
  host.append(ui.divV(lines));
}
