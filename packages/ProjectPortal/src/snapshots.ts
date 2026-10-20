/** Dated program snapshots against the platform: the program's context frozen into a
 * `program_snapshot` row, and its program dataset written as a CSV file to the program's workspace
 * folder, with the file's SHA-256 on the row. Written once; nothing here updates or deletes a
 * snapshot. The rules live in ./model/snapshots. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {textColumn} from './columns';
import dayjs from 'dayjs';
import {portalDb, ProgramRow, ProgramSnapshotRow} from './generated/db';
import {programContext} from './portal';
import {programAccess} from './program-roles';
import {selectedRunsOf} from './evidence';
import {datasetsOf, programTable} from './program-data';
import {existingSpaces} from './workspace';
import {MAIN_DATASET} from './model/program-data';
import {
  checkTable, defaultSnapshotLabel, describeSnapshot, snapshotArtifacts, snapshotCriteria, SnapshotCriterion, snapshotDate,
  snapshotFileName, SnapshotSummary, TableCheck,
} from './model/snapshots';
import {newestFirstBy} from './model/time';
import {filterValue, ROW_LIMIT} from './rows';

/** SHA-256 of a text, as hex. */
export async function sha256(text: string): Promise<string> {
  const digest = await crypto.subtle.digest('SHA-256', new TextEncoder().encode(text));
  return Array.from(new Uint8Array(digest)).map((b) => b.toString(16).padStart(2, '0')).join('');
}

/** The frozen lists, as written: json columns hold objects, so each list sits under one key. */
export const frozenCriteria = (s: ProgramSnapshotRow): SnapshotCriterion[] => (s.criteria as {items?: SnapshotCriterion[]} | undefined)?.items ?? [];
export const frozenRuns = (s: ProgramSnapshotRow): ReturnType<typeof snapshotArtifacts> =>
  (s.artifacts as {runs?: ReturnType<typeof snapshotArtifacts>} | undefined)?.runs ?? [];
export const frozenColumns = (s: ProgramSnapshotRow): string[] => (s.table_columns as {names?: string[]} | undefined)?.names ?? [];

export const summaryOf = (s: ProgramSnapshotRow): SnapshotSummary => ({
  takenOn: s.taken_on ? dayjs(s.taken_on as any).toISOString() : '', label: s.label,
  stageName: s.stage_name ?? null, criteria: frozenCriteria(s),
  tableRows: s.table_rows ?? null,
});

export interface TakenSnapshot {
  snapshot: ProgramSnapshotRow;
  /** Why the table is not kept, or other things worth saying. */
  notes: string[];
}

/** Freezes [program] now. Contributors and approvers of the program may; the table is kept when the
 * program has a workspace folder and a program dataset, else the snapshot says why not. */
export async function takeSnapshot(given: ProgramRow, options: {label?: string, note?: string, reason?: string} = {}):
  Promise<TakenSnapshot> {
  const ctx = await programContext(given.id);
  // the row just read: the caller's may predate a decision (stage, status, folder)
  const program = ctx.program;
  const access = await programAccess(ctx.portfolio, program.id);
  if (!access.approver && !access.contributor)
    throw new Error('Only the portfolio\'s or the program\'s contributors and approvers take snapshots');
  const takenOn = new Date().toISOString();
  const notes: string[] = [];
  const stage = ctx.stages.find((s) => s.id === program.stage_id) ?? null;
  const criteria = snapshotCriteria(ctx.criteria, ctx.stages, ctx.recorded, ctx.readiness);
  const runs = await selectedRunsOf(program.id).catch(() => []);

  let table: {file?: string, folder?: string, rows?: number, columns?: string[], sha?: string, dataset?: string} = {};
  const main = datasetsOf(ctx.portfolio).find((d) => d.code === MAIN_DATASET);
  if (!main)
    notes.push(`${ctx.portfolio.code} declares no program dataset: no table`);
  else {
    try {
      const t = await programTable(ctx.portfolio.code, program.code);
      const csv = t.df.toCsv();
      table = {rows: t.df.rowCount, columns: t.df.columns.names(), sha: await sha256(csv), dataset: t.dataset.name};
      const live = program.space_id && (await existingSpaces([program.space_id]).catch(() => new Set<string>())).has(program.space_id);
      if (!live)
        notes.push(`${program.code} has no workspace folder: the table (${t.df.rowCount} rows) is not kept, only its fingerprint`);
      else {
        const files = grok.dapi.spaces.id(program.space_id!).files;
        let file = snapshotFileName(program.code, takenOn);
        // two snapshots in the same second: the second must not overwrite the first one's table
        if (await files.exists(file).catch(() => false))
          file = file.replace(/Z\.csv$/, `.${takenOn.slice(20, 23)}Z.csv`);
        await files.writeString(file, csv);
        table = {...table, file, folder: program.space_id!};
      }
    } catch (e: any) {
      notes.push(`the table was not kept: ${e?.message ?? e}`);
    }
  }
  const label = options.label?.trim() || defaultSnapshotLabel(program.code, takenOn, options.reason);
  const insert = () => portalDb.programSnapshots.insert({program_id: program.id, taken_on: takenOn as any,
    taken_by: DG.User.current().id, label, note: options.note || notes.join('; ') || undefined,
    stage_code: stage?.code, stage_name: stage?.name, status: program.status ?? 'active',
    // a json column holds an object, not a list
    criteria: {items: criteria}, artifacts: {runs: snapshotArtifacts(runs)},
    dataset: table.dataset, table_file: table.file, table_folder: table.folder, table_rows: table.rows,
    table_columns: table.columns ? {names: table.columns} : undefined, table_sha256: table.sha});
  let row: {id: string};
  try {
    [row] = await insert();
  } catch (e) {
    // no record, no file: a table nobody points at would only confuse the folder's readers
    if (table.file)
      await grok.dapi.spaces.id(table.folder!).files.delete(table.file).catch(() => undefined);
    throw e;
  }
  return {snapshot: await portalDb.programSnapshots.get(row.id), notes};
}

/** The program's snapshots, newest first. */
export async function snapshotsOf(programId: string): Promise<ProgramSnapshotRow[]> {
  return (await portalDb.programSnapshots.query({filter: `program_id = ${filterValue(programId)}`, limit: ROW_LIMIT}))
    .sort(newestFirstBy((x) => x.taken_on ?? x.created_on));
}

/** The snapshot's table as it was frozen, and whether the file still holds it. */
export async function snapshotTable(s: ProgramSnapshotRow): Promise<{df: DG.DataFrame | null, check: TableCheck}> {
  if (!s.table_file || !s.table_folder)
    return {df: null, check: 'not kept'};
  const csv = await grok.dapi.spaces.id(s.table_folder).files.readAsString(s.table_file).catch(() => null);
  const check = checkTable(s.table_sha256, csv === null ? null : await sha256(csv), true);
  if (csv === null)
    return {df: null, check};
  const df = DG.DataFrame.fromCsv(csv);
  const program = await portalDb.programs.get(s.program_id).catch(() => null);
  df.name = `${program?.code ?? 'Program'} · snapshot ${snapshotDate(summaryOf(s).takenOn)}`;
  df.setTag('portal.snapshot', s.id);
  return {df, check};
}

/** One line per snapshot, for reports. */
export const describeSnapshotRow = (s: ProgramSnapshotRow) => describeSnapshot(summaryOf(s));

export function snapshotsFrame(rows: ProgramSnapshotRow[], name: string): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    textColumn('id', rows.map((r) => r.id)),
    DG.Column.fromStrings('taken_on', rows.map((r) => summaryOf(r).takenOn)),
    textColumn('label', rows.map((r) => r.label)),
    textColumn('stage', rows.map((r) => r.stage_code ?? '')),
    textColumn('status', rows.map((r) => r.status ?? '')),
    textColumn('summary', rows.map(describeSnapshotRow)),
    textColumn('table_file', rows.map((r) => r.table_file ?? '')),
    textColumn('note', rows.map((r) => r.note ?? '')),
  ]);
  df.name = name;
  return df;
}
