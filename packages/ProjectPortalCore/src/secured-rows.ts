/** Records of one portfolio placed where only its people read them, and the move off the tables
 * every user could read. The rules and why are in ./model/secured-rows.
 *
 * A new request or sync log row is linked into its space as it is written ([placeApproval],
 * [placeSyncLog]). Linking needs Edit on the space and on the row: the row's author holds the
 * second (the server grants it on promotion), and the portfolio's contributors hold Edit on its
 * workspace. A row that cannot be linked stays with its author and admins; [moveToSecuredTables]
 * links what is left. A sync run, and a request with no program folder, go into a record folder under
 * the root (`_Sync history`, `_Approvals`), so the root lists the program folders only;
 * [moveRecordsOutOfRoot] moves those linked into the root before 0.3.0. */
import * as grok from 'datagrok-api/grok';
import {ApprovalRow, portalDb, PortfolioRow, SyncLogRow} from './generated/db';
import {isApprover, portfolioByCode} from './portal';
import {folderLinks} from './program-rows';
import {filterValue, queryIn, ROW_LIMIT} from './rows';
import {subspaces} from './workspace';
import {APPROVALS_FOLDER, isSpaceNamed, SYNC_HISTORY_FOLDER} from './model/workspace-plan';
import {cardPortfolios} from './model/catalogue';
import {
  copyValues, describeMove, describeRecordMoves, hasContent, inPortfolios, linkValuesOf, MoveReport, Placement, planLinks,
  planRecordMoves, portfoliosOfGroups, requestSpace, toCopy, COMPOUND_CONTENT,
} from './model/secured-rows';

/** The portfolio's workspace root space id, or null when it has none (or the user cannot see it). */
export async function rootIdOf(portfolio: PortfolioRow): Promise<string | null> {
  if (!portfolio.workspace_root)
    return null;
  const found = await grok.dapi.spaces.filter(`name = "${portfolio.workspace_root}"`).list({pageSize: 100}).catch(() => []);
  const root = found.find((s: any) => isSpaceNamed({id: s.id, name: s.name, friendlyName: s.friendlyName ?? null},
    portfolio.workspace_root!));
  return root?.id ?? null;
}

/** The record folder [name] under the portfolio's root, made when missing (unless not [create]).
 * The root itself when the folder is not there and cannot or may not be made: the root reaches the
 * same groups. Null when the portfolio has no root. */
export async function recordFolderId(portfolio: PortfolioRow, name: string, create = true): Promise<string | null> {
  const rootId = await rootIdOf(portfolio);
  if (!rootId)
    return null;
  const found = (await subspaces(rootId).catch(() => [])).find((c) => isSpaceNamed({id: c.id, name: c.name,
    friendlyName: (c as any).friendlyName ?? null}, name));
  if (found || !create)
    return found?.id ?? rootId;
  return (await grok.dapi.spaces.id(rootId).addSubspace(name).catch(() => null))?.id ?? rootId;
}

/** Links [rowId] into [folderId], or into the root when the folder cannot take it. Null when placed. */
async function linkIntoRecords(rowId: string, portfolio: PortfolioRow, folderId: string | null): Promise<string | null> {
  const why = await linkInto(rowId, folderId);
  const rootId = why && folderId ? await rootIdOf(portfolio) : null;
  return rootId && rootId !== folderId ? linkInto(rowId, rootId) : why;
}

/** Links [rowId] into [spaceId] unless it is there. Null when it is linked; else why not. */
export async function linkInto(rowId: string, spaceId: string | null): Promise<string | null> {
  if (!spaceId)
    return 'no workspace to link it into: its author and admins read it';
  try {
    const links = await folderLinks([spaceId]);
    if (!links.has(spaceId))
      return 'its folder cannot be read';
    if (!links.get(spaceId)!.has(rowId))
      await grok.dapi.spaces.id(spaceId).addEntity(rowId, true);
    return null;
  } catch (e: any) {
    return `not linked into its folder: ${e?.message ?? e}`;
  }
}

/** Links a request into its program's folder, or the portfolio's `_Approvals`. Null when placed. */
export async function placeApproval(row: ApprovalRow, portfolio?: PortfolioRow): Promise<string | null> {
  const program = row.program_id ? await portalDb.programs.get(row.program_id).catch(() => null) : null;
  const pf = portfolio ?? await portalDb.portfolios.get(row.portfolio_id);
  const folder = program?.space_id ?? null;
  if (!folder)
    return linkIntoRecords(row.id, pf, await recordFolderId(pf, APPROVALS_FOLDER));
  const first = await linkInto(row.id, requestSpace(folder, null));
  // a program folder that cannot take it: the portfolio's root still reaches the portfolio's groups
  if (first)
    return linkInto(row.id, await rootIdOf(pf));
  return first;
}

/** Links a sync log row into the portfolio's `_Sync history`. Null when placed. */
export async function placeSyncLog(row: SyncLogRow | string, portfolio: PortfolioRow): Promise<string | null> {
  return linkIntoRecords(typeof row === 'string' ? row : row.id, portfolio, await recordFolderId(portfolio, SYNC_HISTORY_FOLDER));
}

/** Moves the sync runs and requests linked into portfolio [code]'s root into its record folders:
 * each is linked into the folder first, then its link in the root removed (the record itself stays).
 * A dry run (the default) only counts. Applying needs one of the portfolio's approvers (admins are).
 * Idempotent: a second run finds nothing in the root. */
export async function moveRecordsOutOfRoot(code: string, dryRun = true): Promise<string[]> {
  const portfolio = await portfolioByCode(code);
  if (!dryRun && !await isApprover(portfolio))
    throw new Error(`Only the approvers of ${code} move its records`);
  const rootId = await rootIdOf(portfolio);
  if (!rootId)
    return [`${code}: no workspace root`];
  // only what is linked into the root: a link's removal leaves the record, an owned child's would not
  const rootLinks = (await folderLinks([rootId])).get(rootId) ?? new Set<string>();
  const children = (await grok.dapi.spaces.id(rootId).children.ofTypes('EntityRecord', true).list({pageSize: 100000}))
    .filter((c) => rootLinks.has(c.id));
  const folders = new Map<string, string>();
  for (const name of [SYNC_HISTORY_FOLDER, APPROVALS_FOLDER]) {
    const id = await recordFolderId(portfolio, name, false);
    if (id && id !== rootId)
      folders.set(name, id);
  }
  const links = await folderLinks([...folders.values()]);
  const inFolder = new Map([...folders].map(([name, id]) => [name, links.get(id) ?? new Set<string>()]));
  const moves = planRecordMoves(children.map((c) => ({id: c.id, name: c.name})), inFolder);
  const errors: string[] = [];
  if (!dryRun) {
    for (const m of moves) {
      try {
        if (!folders.has(m.folder)) {
          const id = await recordFolderId(portfolio, m.folder);
          if (!id || id === rootId)
            throw new Error(`${m.folder} could not be made`);
          folders.set(m.folder, id);
        }
        if (!m.linked)
          await grok.dapi.spaces.id(folders.get(m.folder)!).addEntity(m.id, true);
        // linked into the folder first, so the portfolio's groups never lose it
        await grok.dapi.spaces.id(rootId).removeEntity(m.id);
      } catch (e: any) {
        errors.push(`${m.id} not moved: ${e?.message ?? e}`);
      }
    }
  }
  return describeRecordMoves(code, moves, dryRun, errors);
}

/** The ids of the portfolios the current user belongs to, or null for an admin (all of them). */
export async function myPortfolioIds(): Promise<Set<string> | null> {
  const mine = await grok.dapi.groups.currentUserGroups();
  if (mine.some((g) => g.name === 'Administrators' || g.friendlyName === 'Administrators'))
    return null;
  const portfolios = await portalDb.portfolios.query({limit: ROW_LIMIT});
  return portfoliosOfGroups(portfolios, new Set(mine.map((g) => g.id)));
}

/** The codes of the current user's portfolios, or null for an admin. */
export async function myPortfolioCodes(): Promise<Set<string> | null> {
  const ids = await myPortfolioIds();
  if (ids === null)
    return null;
  return new Set((await portalDb.portfolios.query({limit: ROW_LIMIT})).filter((p) => ids.has(p.id)).map((p) => p.code));
}

/** Whether a catalogued model is the current user's to see: one of their portfolios' catalogues
 * lists it, or none does (a model loaded without a portfolio belongs to every portfolio). */
export async function inMyCatalogues(portfolios: string[]): Promise<boolean> {
  return cardPortfolios(portfolios, await myPortfolioCodes()) !== null;
}

/** [rows] of the current user's portfolios only. */
export async function ofMyPortfolios<T extends {portfolio_id: string}>(rows: T[]): Promise<T[]> {
  return inPortfolios(rows, await myPortfolioIds());
}

/** Moves portfolio [code]'s records off the tables every user could read: its approval requests
 * into `approval` and its sync runs into `sync_log` (same ids), each linked into its space; its
 * compound links given the compound's values, and compounds whose every link has them emptied.
 * Idempotent; with [dryRun] only reports. Needs an admin: it reads the closed tables. */
export async function moveToSecuredTables(code: string, dryRun = false): Promise<string[]> {
  if (await myPortfolioIds() !== null)
    throw new Error('Only admins move records: the old tables are closed to everyone else');
  const portfolio = await portfolioByCode(code);
  const by = `portfolio_id = ${filterValue(portfolio.id)}`;
  const report: MoveReport = {portfolio: code, copied: {}, already: {}, linked: {}, unplaced: [], linkValues: 0, emptied: 0, errors: []};
  const programs = await portalDb.programs.query({filter: by, limit: ROW_LIMIT});
  const folderOf = new Map(programs.map((p) => [p.id, p.space_id ?? null]));
  const codeOf = new Map(programs.map((p) => [p.id, p.code]));

  const move = async (name: string, oldRows: any[], table: any, placement: (r: any) => Placement) => {
    const copied = new Set<string>((await table.query({filter: by, limit: ROW_LIMIT})).map((r: any) => r.id as string));
    const pending = toCopy(oldRows, copied);
    report.copied[name] = pending.length;
    report.already[name] = oldRows.length - pending.length;
    if (!dryRun) {
      for (const r of pending)
        await table.insert(copyValues(r)).catch((e: any) => report.errors.push(`${name} ${r.id} not copied: ${e?.message ?? e}`));
    }
    const rows = dryRun ? oldRows : await table.query({filter: by, limit: ROW_LIMIT});
    const places = rows.map(placement);
    const spaces = [...new Set<string>(places.map((p: Placement) => p.spaceId).filter((x: string | null): x is string => !!x))];
    const plan = planLinks(places, await folderLinks(spaces));
    report.linked[name] = plan.link.length;
    report.unplaced.push(...plan.unplaced.map((p) => `${name} "${p.label}"`));
    if (!dryRun) {
      for (const p of plan.link) {
        const why = await linkInto(p.id, p.spaceId);
        if (why)
          report.errors.push(`${name} "${p.label}" ${why}`);
      }
    }
  };

  try {
    const oldRequests = await grok.dapi.domains.table('portal.approval_request').query({filter: by, limit: ROW_LIMIT});
    // a record folder is made only for records that go there, and not by a dry run (it counts against the root)
    const approvalsId = await recordFolderId(portfolio, APPROVALS_FOLDER,
      !dryRun && oldRequests.some((r: any) => !(r.program_id && folderOf.get(r.program_id))));
    await move('approval', oldRequests, portalDb.approvals, (r) => ({id: r.id, label: r.title,
      spaceId: requestSpace(r.program_id ? folderOf.get(r.program_id) : null, approvalsId)}));
    const oldRuns = await grok.dapi.domains.table('portal.sync_run').query({filter: by, limit: ROW_LIMIT});
    const syncId = await recordFolderId(portfolio, SYNC_HISTORY_FOLDER, !dryRun && oldRuns.length > 0);
    await move('sync_log', oldRuns, portalDb.syncLogs, (r) => ({id: r.id, label: `sync ${String(r.started_on)}`, spaceId: syncId}));
  } catch (e: any) {
    report.errors.push(`not moved: ${e?.message ?? e}`);
  }

  // compounds: each link of this portfolio carries its compound's values
  try {
    const links = await queryIn(portalDb.programCompounds, 'program_id', programs.map((p) => p.id));
    const compounds = await queryIn(portalDb.compounds, 'id', links.map((l) => l.compound_id));
    const byId = new Map(compounds.map((c) => [c.id, c]));
    const bare = links.filter((l) => !l.compound_code && byId.has(l.compound_id));
    report.linkValues = bare.length;
    if (!dryRun) {
      for (const l of bare) {
        await portalDb.programCompounds.update(l.id, linkValuesOf(byId.get(l.compound_id)!) as any)
          .catch((e: any) => report.errors.push(`${codeOf.get(l.program_id)}/${byId.get(l.compound_id)!.code}: values not copied: ${e?.message ?? e}`));
      }
    }
    // a compound is emptied once every link to it, in any portfolio, carries the values
    const allLinks = await queryIn(portalDb.programCompounds, 'compound_id', compounds.map((c) => c.id));
    const done = new Set(bare.map((l) => l.id));
    const full = compounds.filter((c) => hasContent(c) &&
      allLinks.filter((l) => l.compound_id === c.id).every((l) => l.compound_code || done.has(l.id)));
    report.emptied = full.length;
    if (!dryRun) {
      const empty = Object.fromEntries(COMPOUND_CONTENT.map((k) => [k, null]));
      for (const c of full)
        await portalDb.compounds.update(c.id, empty as any).catch((e: any) => report.errors.push(`${c.code} not emptied: ${e?.message ?? e}`));
    }
  } catch (e: any) {
    report.errors.push(`compounds not moved: ${e?.message ?? e}`);
  }
  return describeMove(report, dryRun);
}
