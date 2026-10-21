import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {applicationRecorded} from '../apply/application-store';
import {featureFit} from '../catalog/applicable-tables';
import {loginsOf} from '../catalog/list-models';
import {ModelActivity, modelActivity} from '../catalog/model-activity';
import {textColumn} from '../columns';
import {MINUTE_FORMAT, MODEL_TYPE, SECOND_FORMAT} from '../constants';
import {ApplicationRow, ApplicationSource, ApplicationStatus, forgeDb, ModelRow} from '../generated/db';
import {datasetRefCaption, storedDatasetRef} from '../storage/dataset-ref';
import {tagsOf, tagsText} from '../storage/model-fields';
import {ownBlob} from '../storage/model-store';
import {readOnlyGrid} from './data-grid';
import {runEvaluation} from './evaluation-view';
import {MetricsTable, storedSummary} from './metrics-table';
import {reportError} from './report-error';
import {STORAGE_CAPTIONS, writeModelInfo} from './save-model-dialog';
import {tagsInput, tagsOfInput} from './tags-input';

const NO_SHARE = 'Sharing a model needs the Share permission on this model; ask an administrator.';
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
let activitySub: rxjs.Subscription | undefined;
/** The model the newest accordion shows, at the version it shows (its own Tags writes included). */
let shownModel: Pick<ModelRow, 'id' | 'version'> | undefined;
type Refill = () => void;
/** By Details host: shows again its cells that follow the model's applications (Last run, Applications). */
const ACTIVITY_CELLS = new WeakMap<HTMLElement, Refill>();
/** The authors by user id, read once per session for the tooltips and Details; a failed read is tried again. */
const USERS = new Map<string, Promise<DG.User>>();

type Shares = {view?: DG.Group[]; edit?: DG.Group[]};
type FactsRow = Pick<ModelRow, 'author_id' | 'created_on' | 'updated_on' | 'dataset_name' | 'row_count' | 'features' |
  'options'>;

/** What both the model's Details and the table pane's tooltip show: the author (a link), the creation and update
 * times, the training table with its rows, the last run, the feature names and the open tables the model fits ('' for
 * none). */
export interface ModelFacts {
  author: HTMLElement;
  created: string;
  updated: string;
  table: string;
  lastRun: HTMLElement;
  features: string[];
  applicableTo: string;
}

export function modelFacts(model: FactsRow, lastRun: () => Promise<ModelActivity['lastRun']>): ModelFacts {
  const {names, tables} = featureFit(model, grok.shell.tables);
  return {
    author: ui.wait(async () => ui.render(await userOf(model.author_id))),
    created: model.created_on.format(MINUTE_FORMAT),
    updated: model.updated_on.format(MINUTE_FORMAT),
    table: `${model.dataset_name ?? ''}${model.row_count === undefined ? '' : ` (${model.row_count} rows)`}`,
    lastRun: lastRunCell(lastRun),
    features: names,
    applicableTo: tables.map((t) => t.name).join(', '),
  };
}

/** The user [id], read once per session; a failed read is not kept, so the next call tries again. */
export function userOf(id: string): Promise<DG.User> {
  let user = USERS.get(id);
  if (user === undefined) {
    user = grok.dapi.users.find(id);
    USERS.set(id, user);
    user.catch(() => USERS.delete(id));
  }
  return user;
}

/** Whether the newest accordion shows [model] at its version. */
export function isShownVersion(model: Pick<ModelRow, 'id' | 'version'>): boolean {
  return shownModel?.id === model.id && shownModel.version === model.version;
}

/** The model icon of the built-in tool. */
export function modelIcon(): HTMLElement {
  return ui.iconSvg('model');
}

/** The model's context panel: the title the platform gives a row (icon, favorites star, name, context actions), then
 * **Details** (open), **Performance**, **Activity**, **Sharing** and **History**, each built when first opened. [row]
 * is the shown row, [model] its full values. An application of the model recorded in this session while the panel is
 * shown refreshes Details' Last run and Applications and the Activity pane in place. */
export function modelAccordion(row: DG.DomainRow, model: ModelRow): DG.Accordion {
  let activity: Promise<ModelActivity> | undefined;
  const activityOf = () => activity ??= modelActivity(model.id);
  const details = ui.divV([]);
  const activityHost = ui.divV([], 'forge-pane-fill');
  const accordion = ui.accordion(MODEL_TYPE);
  // Bound for the tooltip and the drag (onto a Space); the commands are the context actions'.
  const name = ui.bind(row, ui.label(model.name, 'forge-model-name'), {contextMenu: false});
  accordion.addTitle(ui.span([modelIcon(), ui.star(row.id), name, ui.contextActions(row)]));
  accordion.addPane('Details', () => fillDetails(details, model, activityOf), true);
  accordion.addPane('Performance', () => performancePane(model));
  accordion.addPane('Activity', () => fillActivity(activityHost, activityOf));
  accordion.addPane('Sharing', () => sharingPane(row));
  accordion.addPane('History', () => DG.DomainObjectHandler.auditPane(row));
  shownModel = {id: model.id, version: model.version};
  // The context panel shows one model: the newest accordion follows its model's applications.
  activitySub?.unsubscribe();
  activitySub = applicationRecorded.subscribe((id) => {
    if (id !== model.id || !accordion.root.isConnected)
      return;
    activity = undefined;
    // Details keeps its Tags and their version; a pane not built yet reads the activity when it is opened.
    ACTIVITY_CELLS.get(details)?.();
    if (activityHost.hasChildNodes())
      fillActivity(activityHost, activityOf);
  });
  return accordion;
}

/** Fills [host] with `N application(s)` and the grid of the newest applications, or "Not applied yet.". */
function fillActivity(host: HTMLElement, activity: () => Promise<ModelActivity>): HTMLElement {
  ui.empty(host);
  host.append(ui.wait(async () => {
    const {applications, count} = await activity();
    if (count === 0)
      return ui.divText('Not applied yet.');
    const logins = await loginsOf(applications.map((a) => a.author_id));
    const at = (row: number) => applications[row];
    const grid = readOnlyGrid(activityFrame(applications, logins), {headers: ACTIVITY_HEADERS, cells: {
      'When': (row) => at(row).created_on.format(SECOND_FORMAT),
      'Status': (row) => `${STATUS_TOOLTIPS[at(row).status]}${at(row).status === 'failed' ? at(row).error ?? '' : ''}`,
      'Source': (row) => SOURCE_TOOLTIPS[at(row).source],
    }});
    return ui.divV([ui.divText(`${count} ${count === 1 ? 'application' : 'applications'}`), grid.root]);
  }));
  return host;
}

function fillDetails(host: HTMLElement, model: ModelRow, activity: () => Promise<ModelActivity>): HTMLElement {
  const lastRun = async () => (await activity()).lastRun;
  const facts = modelFacts(model, lastRun);
  const lastRunHost = ui.div([facts.lastRun]);
  const countHost = ui.div([]);
  const showCount = () => countHost.append(ui.wait(async () => ui.divText(`${(await activity()).count}`)));
  showCount();
  ACTIVITY_CELLS.set(host, () => {
    ui.empty(lastRunHost);
    ui.empty(countHost);
    lastRunHost.append(lastRunCell(lastRun));
    showCount();
  });
  const details: {[caption: string]: string | HTMLElement} = {
    'Author': facts.author,
    'Created': facts.created,
    'Updated': facts.updated,
    'Table': facts.table,
    ...storageDetails(model),
    'Last run': lastRunHost,
    'Applications': countHost,
    'Features': facts.features.join(', '),
    'Target': model.target_name,
    'Method': model.engine_name,
    'Task': model.task,
  };
  if (facts.applicableTo !== '')
    details['Applicable to'] = facts.applicableTo;
  ui.empty(host);
  host.append(...(model.description ? [ui.divText(model.description)] : []), ui.tableFromMap(details),
    ui.form([detailsTags(host, model, activity)]));
  return host;
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
      if (shownModel?.id === model.id)
        shownModel = {id: model.id, version};
    }
  };
  // One write at a time, each with the version the previous one returned: quick chip changes are not a conflict.
  let writes = Promise.resolve();
  input.onChanged.subscribe(() => {
    writes = writes.then(write).catch((e) => reportError(e));
  });
  return input;
}

function lastRunCell(lastRun: () => Promise<ModelActivity['lastRun']>): HTMLElement {
  return ui.wait(async () => ui.divText(lastRunText(await lastRun())));
}

function lastRunText(lastRun: ModelActivity['lastRun']): string {
  if (lastRun === undefined)
    return 'Never';
  const when = lastRun.when.format(MINUTE_FORMAT);
  return lastRun.status === 'completed' ? when : `${when} (${lastRun.status})`;
}

/** The stored metrics, then **Run Evaluation** for a model with Forge's own model file. */
function performancePane(model: ModelRow): HTMLElement {
  const summary = storedSummary(model);
  const metrics = summary === null ? ui.divText('No metrics were recorded for this model.') :
    new MetricsTable(summary).root;
  const run = ownBlob(model.blob) === null ? [] : [ui.button('Run Evaluation', () => void runEvaluation(model),
    'Apply the model to its training data, or to an open table, and show how well it predicts.')];
  return ui.divV([metrics, ...run], 'forge-pane-fill');
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
