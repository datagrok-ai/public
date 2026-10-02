/** Parking orphan folders against the platform (FR-22): the folders the orphan report lists are
 * moved into `_Parked/<date>` under the portfolio's root with `SpaceClient.addEntity` (a move: the
 * folder keeps its id, files, subfolders and grants), each move recorded in `parked_folder`, and
 * `unparkFolder` moves one back. Nothing is deleted, except an empty `_Parked/<date>` or `_Parked`
 * folder that unparking leaves behind. The rules live in ./model/parking. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {ParkedFolderRow, portalDb, PortfolioRow} from './generated/db';
import {isApprover, portfolioByCode} from './portal';
import {findRoot, subspaces, withRootLock, workspaceOrphans} from './workspace';
import {
  checkParking, checkUnpark, isParkingFolder, parkDay, ParkedRef, ParkPlan, planParking, PARKED_FOLDER, summarizeParking,
} from './model/parking';
import {isSpaceNamed} from './model/workspace-plan';

const asRef = (e: DG.Entity) => ({id: e.id, name: e.name, friendlyName: (e as any).friendlyName ?? null});

export const parkedRef = (r: ParkedFolderRow): ParkedRef => ({id: r.id, name: r.name, status: r.status as ParkedRef['status'],
  fromPath: r.from_path, parkedPath: r.parked_path});

/** The portfolio's parked folders, newest first; only those still parked unless [all]. */
export async function parkedFoldersOf(portfolio: PortfolioRow, all = false): Promise<ParkedFolderRow[]> {
  return (await portalDb.parkedFolders.query({filter: `portfolio_id = "${portfolio.id}"`}))
    .filter((r) => all || r.status === 'parked')
    .sort((a, b) => String(b.created_on).localeCompare(String(a.created_on)));
}

/** The child of [parentId] called [name], made when missing. */
async function ensureChild(parentId: string, name: string): Promise<DG.Project> {
  const found = (await subspaces(parentId)).find((c) => isSpaceNamed(asRef(c), name));
  return found ?? await grok.dapi.spaces.id(parentId).addSubspace(name);
}

export interface ParkResult {
  portfolio: string;
  dryRun: boolean;
  plan: ParkPlan;
  /** The parked_folder rows written. */
  parked: ParkedFolderRow[];
  errors: string[];
}

/** Parks the folders the orphan report lists (or, with [only], those of them named by path or id)
 * into `_Parked/<date>` under the root. A dry run (the default) only plans. Applying needs one of
 * the portfolio's approvers. Each move is recorded, so `unparkFolder` can undo it. */
export async function parkOrphans(code: string, dryRun = true, only?: string[]): Promise<ParkResult> {
  const portfolio = await portfolioByCode(code);
  const refused = checkParking({isApprover: await isApprover(portfolio), dryRun});
  if (refused.length)
    throw new Error(refused.join('; '));
  const rootName = portfolio.workspace_root;
  if (!rootName)
    throw new Error(`portfolio ${code} has no workspaceRoot`);
  const report = await workspaceOrphans(code);
  const parkedIds = (await parkedFoldersOf(portfolio)).map((r) => r.space_id);
  const day = parkDay(new Date());
  const plan = planParking(rootName, report.orphans, day, parkedIds, only);
  const result: ParkResult = {portfolio: code, dryRun, plan, parked: [], errors: []};
  if (dryRun || plan.moves.length === 0)
    return result;
  await withRootLock(rootName, async () => {
    const root = await findRoot(rootName);
    if (!root)
      throw new Error(`the workspace of ${code} is not built`);
    const parking = await ensureChild(root.id, PARKED_FOLDER);
    const dayFolder = await ensureChild(parking.id, day);
    for (const m of plan.moves) {
      try {
        await grok.dapi.spaces.id(dayFolder.id).addEntity(m.id, false);
        // the move is checked, not assumed: the folder must now be under the day folder
        if (!(await subspaces(dayFolder.id)).some((c) => c.id === m.id))
          throw new Error('the server did not move it');
        const [r] = await portalDb.parkedFolders.insert({portfolio_id: portfolio.id, space_id: m.id, name: m.name, kind: m.kind,
          reason: m.reason, from_parent_id: root.id, from_path: m.from, parked_in_id: dayFolder.id, parked_path: m.to,
          parked_by: DG.User.current().id, parked_on: dayjs(), status: 'parked'});
        result.parked.push(await portalDb.parkedFolders.get(r.id));
      } catch (e: any) {
        result.errors.push(`${m.from}: not parked (${e?.message ?? e})`);
      }
    }
  });
  return result;
}

export const summarizeParkResult = (r: ParkResult): string => summarizeParking(r.portfolio, r.plan, r.dryRun, r.errors);

/** Moves a parked folder back to where it was, by its parked_folder id or its space id. Approvers
 * of its portfolio only. An empty day folder, and then an empty `_Parked`, is removed. */
export async function unparkFolder(id: string): Promise<{row: ParkedFolderRow, line: string, notes: string[]}> {
  // by the record's id, else the folder's: the latest record of it still parked
  const row = (await portalDb.parkedFolders.query({filter: `id = "${id}"`}))[0] ??
    (await portalDb.parkedFolders.query({filter: `space_id = "${id}"`}))
      .sort((a, b) => String(b.created_on).localeCompare(String(a.created_on))).find((r) => r.status === 'parked');
  if (!row)
    throw new Error(`no parked folder ${id}`);
  const portfolio = await portalDb.portfolios.get(row.portfolio_id);
  const refused = checkParking({isApprover: await isApprover(portfolio), dryRun: false});
  const notes: string[] = [];
  await withRootLock(portfolio.workspace_root ?? portfolio.code, async () => {
    refused.push(...checkUnpark(parkedRef(row), (await subspaces(row.from_parent_id)).map(asRef)));
    if (refused.length)
      throw new Error(refused.join('; '));
    await grok.dapi.spaces.id(row.from_parent_id).addEntity(row.space_id, false);
    if (!(await subspaces(row.from_parent_id)).some((c) => c.id === row.space_id))
      throw new Error(`${row.name}: the server did not move it back`);
    await portalDb.parkedFolders.update(row.id, {status: 'restored', restored_by: DG.User.current().id, restored_on: dayjs()});
    notes.push(...await tidyParking(row.parked_in_id, row.from_parent_id));
  });
  return {row: await portalDb.parkedFolders.get(row.id), line: `${row.parked_path} → ${row.from_path}: unparked`, notes};
}

/** Whether a folder holds nothing: no subfolder, no file, no other entity. Every space lists its
 * own storage connection, `Files`, among its children; that one does not count. */
async function isEmpty(spaceId: string): Promise<boolean> {
  const children = await grok.dapi.spaces.id(spaceId).children.list({pageSize: 20}).catch(() => null);
  return Array.isArray(children) && children.every((c) => c instanceof DG.DataConnection && c.name === 'Files');
}

/** Removes the day folder [dayId] when unparking emptied it, then `_Parked` under [rootId] when
 * that is empty too. Only folders named as the parking folder are ever removed. */
async function tidyParking(dayId: string, rootId: string): Promise<string[]> {
  const notes: string[] = [];
  const parking = (await subspaces(rootId)).find((c) => isParkingFolder(asRef(c)));
  if (!parking)
    return notes;
  const day = (await subspaces(parking.id)).find((c) => c.id === dayId);
  if (day && await isEmpty(day.id)) {
    await grok.dapi.spaces.delete(day);
    notes.push(`${PARKED_FOLDER}/${day.friendlyName || day.name} was empty: removed`);
  }
  if (await isEmpty(parking.id)) {
    await grok.dapi.spaces.delete(parking);
    notes.push(`${PARKED_FOLDER} was empty: removed`);
  }
  return notes;
}

/** One row per parked folder: name, why, from, to, who and when, status. */
export function parkedFrame(rows: ParkedFolderRow[], name: string): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('id', rows.map((r) => r.id)),
    DG.Column.fromStrings('folder', rows.map((r) => r.name)),
    DG.Column.fromStrings('kind', rows.map((r) => r.kind ?? '')),
    DG.Column.fromStrings('reason', rows.map((r) => r.reason ?? '')),
    DG.Column.fromStrings('from', rows.map((r) => r.from_path)),
    DG.Column.fromStrings('parked in', rows.map((r) => r.parked_path)),
    DG.Column.fromStrings('parked on', rows.map((r) => r.parked_on ? dayjs(r.parked_on as any).toISOString() : '')),
    DG.Column.fromStrings('status', rows.map((r) => r.status)),
    DG.Column.fromStrings('restored on', rows.map((r) => r.restored_on ? dayjs(r.restored_on as any).toISOString() : '')),
  ]);
  df.name = name;
  return df;
}
