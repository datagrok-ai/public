/** Records of one portfolio placed where only its people read them, and the move off the tables
 * every user could read. The rules and why are in ./model/secured-rows.
 *
 * A new request or sync log row is linked into its space as it is written ([placeApproval],
 * [placeSyncLog]). Linking needs Edit on the space and on the row: the row's author holds the
 * second (the server grants it on promotion), and the portfolio's contributors hold Edit on its
 * workspace. A row that cannot be linked stays with its author and admins; [moveToSecuredTables]
 * links what is left. */
import * as grok from 'datagrok-api/grok';
import {ApprovalRow, portalDb, PortfolioRow, SyncLogRow} from './generated/db';
import {portfolioByCode} from './portal';
import {folderLinks} from './program-rows';
import {filterValue, queryIn, ROW_LIMIT} from './rows';
import {isSpaceNamed} from './model/workspace-plan';
import {
  copyValues, describeMove, hasContent, inPortfolios, linkValuesOf, MoveReport, Placement, planLinks,
  portfoliosOfGroups, requestSpace, toCopy, COMPOUND_CONTENT,
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

/** Links a request into its program's folder, or the portfolio's root. Null when placed. */
export async function placeApproval(row: ApprovalRow, portfolio?: PortfolioRow): Promise<string | null> {
  const program = row.program_id ? await portalDb.programs.get(row.program_id).catch(() => null) : null;
  const pf = portfolio ?? await portalDb.portfolios.get(row.portfolio_id);
  const folder = program?.space_id ?? null;
  const first = await linkInto(row.id, requestSpace(folder, folder ? null : await rootIdOf(pf)));
  // a program folder that cannot take it: the portfolio's root still reaches the portfolio's groups
  if (first && folder)
    return linkInto(row.id, await rootIdOf(pf));
  return first;
}

/** Links a sync log row into the portfolio's root. Null when placed. */
export async function placeSyncLog(row: SyncLogRow | string, portfolio: PortfolioRow): Promise<string | null> {
  return linkInto(typeof row === 'string' ? row : row.id, await rootIdOf(portfolio));
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
  const mine = await myPortfolioCodes();
  return mine === null || !portfolios.length || portfolios.some((p) => mine.has(p));
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
  const portfolio = await portfolioByCode(code);
  const by = `portfolio_id = ${filterValue(portfolio.id)}`;
  const report: MoveReport = {portfolio: code, copied: {}, already: {}, linked: {}, unplaced: [], linkValues: 0, emptied: 0, errors: []};
  const rootId = await rootIdOf(portfolio);
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
    await move('approval', oldRequests, portalDb.approvals, (r) => ({id: r.id, label: r.title,
      spaceId: requestSpace(r.program_id ? folderOf.get(r.program_id) : null, rootId)}));
    const oldRuns = await grok.dapi.domains.table('portal.sync_run').query({filter: by, limit: ROW_LIMIT});
    await move('sync_log', oldRuns, portalDb.syncLogs, (r) => ({id: r.id, label: `sync ${String(r.started_on)}`, spaceId: rootId}));
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
