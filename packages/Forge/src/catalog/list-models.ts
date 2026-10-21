import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {featureSchemasOf} from '../apply/apply-model';
import {textColumn} from '../columns';
import {MAX_MODELS} from '../constants';
import {ForgeError} from '../forge-error';
import {forgeDb, ModelRow, ModelTask} from '../generated/db';
import {tagsOf, tagsText} from '../storage/model-fields';
import {applicableTables} from './applicable-tables';
import {metricColumns} from './compare-models';

const DEFAULT_LIMIT = 100;
const TASKS: readonly ModelTask[] = ['regression', 'classification'];

/** The `model` columns `Forge:listModels` reads (the system columns come with any query). */
const LIST_COLUMNS = ['name', 'description', 'engine_name', 'task', 'target_name', 'features', 'options',
  'row_count', 'storage_mode', 'tags', 'metrics'] as const;
type ListModelRow = Pick<ModelRow, typeof LIST_COLUMNS[number] | 'id' | 'author_id' | 'created_on'>;

/** The inputs of `Forge:listModels`; a missing optional input is null or undefined. */
export interface ListModelsInput {
  name?: string | null;
  target?: string | null;
  task?: string | null;
  method?: string | null;
  tags?: string | null;
  table?: DG.DataFrame | null;
  limit?: number | null;
}

/** The body of `Forge:listModels`: the visible models, newest first, at most [ListModelsInput.limit]; `task` and
 * `method` filter on the server, the name part, the target, the tags and the table in the browser over at most
 * {@link MAX_MODELS} rows. Nothing is uploaded. */
export async function runListModels(input: ListModelsInput): Promise<DG.DataFrame> {
  const task = input.task || null;
  if (task !== null && !TASKS.some((t) => t === task))
    throw new ForgeError('The task must be regression or classification.');
  const limit = input.limit ?? DEFAULT_LIMIT;
  if (!Number.isInteger(limit) || limit < 1)
    throw new ForgeError('The limit must be a whole number of at least 1.');
  const query = forgeDb.models.query();
  if (task !== null)
    query.where('task', '=', task);
  if (input.method)
    query.where('engine_name', '=', input.method);
  const name = input.name?.toLowerCase();
  const target = input.target?.toLowerCase();
  const tags = tagsOf(input.tags).map((t) => t.toLowerCase());
  const table = input.table || null;
  const filtered = !!name || !!target || tags.length > 0 || table !== null;
  const rows: ListModelRow[] = await query.select(...LIST_COLUMNS).orderBy('created_on', true)
    .top(filtered ? MAX_MODELS : Math.min(limit, MAX_MODELS));
  const fits = (row: ListModelRow): boolean => {
    if (name && !row.name.toLowerCase().includes(name) || target && row.target_name.toLowerCase() !== target)
      return false;
    if (tags.length > 0) {
      const rowTags = tagsOf(row.tags).map((t) => t.toLowerCase());
      if (!tags.every((tag) => rowTags.includes(tag)))
        return false;
    }
    return table === null || applicableTables(row, [table]).length > 0;
  };
  const listed: ListModelRow[] = [];
  for (let i = 0; i < rows.length && listed.length < limit; i++) {
    if (fits(rows[i]))
      listed.push(rows[i]);
  }
  return modelsFrame(listed, await loginsOf(listed.map((row) => row.author_id)));
}

/** The users' logins by user id, from one read of the unique [ids]; a deleted or unknown user is not found. */
export async function loginsOf(ids: string[]): Promise<Map<string, string>> {
  const unique = [...new Set(ids)];
  const users = unique.length === 0 ? [] : await grok.dapi.getEntities(unique);
  return new Map(users.filter((u): u is DG.User => u instanceof DG.User).map((u) => [u.id, u.login]));
}

/** **Forge models**: one row per model with **Id**, **Name**, **Description**, **Method**, **Task**, **Target**,
 * **Features**, **Training rows**, **Data storage**, **Tags**, **Author**, **Created**, then `<metric> (validation)`
 * for every metric a listed model has. */
function modelsFrame(rows: ListModelRow[], logins: Map<string, string>): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    textColumn('Id', rows.map((r) => r.id)),
    textColumn('Name', rows.map((r) => r.name)),
    textColumn('Description', rows.map((r) => r.description)),
    textColumn('Method', rows.map((r) => r.engine_name)),
    textColumn('Task', rows.map((r) => r.task)),
    textColumn('Target', rows.map((r) => r.target_name)),
    textColumn('Features', rows.map((r) => (featureSchemasOf(r.features) ?? []).map((f) => f.name).join(', '))),
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'Training rows', rows.map((r) => r.row_count ?? null)),
    textColumn('Data storage', rows.map((r) => r.storage_mode)),
    textColumn('Tags', rows.map((r) => tagsText(tagsOf(r.tags)))),
    textColumn('Author', rows.map((r) => logins.get(r.author_id))),
    DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 'Created', rows.map((r) => r.created_on)),
    ...metricColumns(rows, ['validation']),
  ]);
  df.name = 'Forge models';
  return df;
}
