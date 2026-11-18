import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ForgeError} from '../forge-error';
import {isRecord} from '../preparation/preparation-options';

export const DATASET_REF_KINDS = ['file', 'query', 'script'] as const;
export type DatasetRefKind = typeof DATASET_REF_KINDS[number];

/** Where a table came from, in reference mode: `script` is its creation script, one call per line, the first
 * assigning the table's variable; `path` for a file, `id` and `name` for a query. */
export interface DatasetRef { kind: DatasetRefKind; script: string; path?: string; id?: string; name?: string }

// The platform appends `//{"timestamp": <ms>}` to every line of a creation script (data_history.dart).
const TIMESTAMP_COMMENT = /\s*\/\/\{"timestamp":\s*\d+\}\s*$/;
const ASSIGNMENT = /^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.+)$/;
const OPEN_FILE = /^(?:OpenFile|OpenServerFile)\("([^"]+)"/;
const CALL_LINE = /^[A-Za-z_]\w*\s*=\s*[A-Za-z_]\w*(?::[A-Za-z_]\w*)?\((.*)\)(?:\.[A-Za-z_]\w*)?$/;
const QUOTED = /"(?:[^"\\]|\\.)*"/g;
const SOURCE_VARIABLE = 'data';

/** The origin of [table] from the tags the platform writes: a query, a server file or another creation script;
 * null when the platform recorded none (a table built in code, a local file, data history turned off). */
export function datasetRefOf(table: DG.DataFrame): DatasetRef | null {
  const lines = (table.getTag(DG.Tags.CreationScript) ?? '').split('\n')
    .map((line) => line.replace(TIMESTAMP_COMMENT, '').trim())
    .filter((line) => line !== '');
  const first = ASSIGNMENT.exec(lines[0] ?? '');
  if (first === null)
    return fileRefOf(table.getTag(DG.Tags.SourceFile));
  const script = lines.join('\n');
  const queryId = table.getTag(DG.Tags.DataQueryId);
  const queryName = table.getTag(DG.Tags.DataQueryName);
  if (queryId)
    return queryName ? {kind: 'query', script, id: queryId, name: queryName} : {kind: 'query', script, id: queryId};
  const path = OPEN_FILE.exec(first[2])?.[1];
  return path === undefined ? {kind: 'script', script} : {kind: 'file', script, path};
}

/** What [ref] points at, for the user: the file path, the query name or `a script`. */
export function datasetRefCaption(ref: DatasetRef): string {
  return ref.path ?? ref.name ?? 'a script';
}

/** Opens the referenced data; the table is not added to the workspace. A file is opened by its path; a query, once
 * found, and any other script replay the script as the platform's data sync does: the lines run one by one in a
 * fresh context, and the table is the first line's variable. A stored reference is editable, so only lines that
 * assign a single call with plain arguments are run. */
export async function openDatasetRef(ref: DatasetRef): Promise<DG.DataFrame> {
  if (ref.kind === 'file' && ref.path !== undefined)
    return await grok.data.files.openTable(ref.path);
  if (ref.kind === 'query' && (ref.id === undefined || !(await grok.dapi.queries.find(ref.id))))
    throw new ForgeError('The query the model was trained on has been moved or deleted.');
  const lines = ref.script.split('\n');
  const variable = ASSIGNMENT.exec(lines[0])?.[1];
  if (variable === undefined || !lines.every(isCallLine))
    throw new ForgeError('The data source of the model is not a script Forge can run. Choose another table.');
  const context = DG.Context.create();
  for (const line of lines)
    await grok.functions.eval(line, context);
  const table: unknown = context.getVariable(variable);
  if (!(table instanceof DG.DataFrame))
    throw new ForgeError('The data source of the model did not give a table. It may have been moved or deleted.');
  return table;
}

/** `<variable> = <call>(<arguments>)`, optionally with an output accessor; the arguments, quoted strings aside, hold
 * no call and no statement separator. */
function isCallLine(line: string): boolean {
  const args = CALL_LINE.exec(line.trim())?.[1];
  return args !== undefined && !/[();]/.test(args.replace(QUOTED, '""'));
}

/** A stored `dataset_ref` value, or null when it is not a reference Forge can open. */
export function storedDatasetRef(value: unknown): DatasetRef | null {
  if (!isRecord(value))
    return null;
  const {kind, script, path, id, name} = value;
  if (!isDatasetRefKind(kind) || typeof script !== 'string' || script === '')
    return null;
  const ref: DatasetRef = {kind, script};
  if (typeof path === 'string')
    ref.path = path;
  if (typeof id === 'string')
    ref.id = id;
  if (typeof name === 'string')
    ref.name = name;
  return ref;
}

function isDatasetRefKind(value: unknown): value is DatasetRefKind {
  return DATASET_REF_KINDS.some((kind) => kind === value);
}

/** A server file path recorded by a file handler; a local file's tag holds only its name, which has no `:`. */
function fileRefOf(sourceFile: string | null): DatasetRef | null {
  if (!sourceFile?.includes(':'))
    return null;
  return {kind: 'file', script: `${SOURCE_VARIABLE} = OpenFile(${JSON.stringify(sourceFile)})`, path: sourceFile};
}
