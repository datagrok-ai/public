/** Whether each portfolio's groups can query the data connections its datasets and feeds use,
 * read from the server and never changed: the rules are in ./model/data-access.
 *
 * The js-api's `grok.dapi.permissions.get` lists only the groups with View and with Edit, so it
 * cannot tell whether a group holds `DataConnection.Query`. The server's privileges route, the
 * one the sharing dialog reads, lists every permission by name:
 * `GET /privileges/permissions?entityId=<connection>&groupId=<group>&all=true` returns the
 * group's effective grants on the connection, its parent groups' (All users among them) and
 * those inherited through projects included. So a group covered through its parent group counts
 * as covered. The caller needs View on the connection; without it the answer is "unknown". */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {portalDb, PortfolioRow} from './generated/db';
// portal.ts imports from here too; the cycle is safe because neither uses the other at load
import {COMPOUND_FEED, findConnection, portfolioByCode} from './portal';
import {datasetsOf} from './program-data';
import {PortfolioConfig} from './model/config';
import {programDatasets} from './model/program-data';
import {
  Access, accessFrom, AccessCheck, AccessRole, ACCESS_ROLES, accessStatus, accessFix, ConnectionUse, connectionUses,
  describeAccessProblems, missingPermissions,
} from './model/data-access';
import {filterValue, ROW_LIMIT} from './rows';

/** The effective permission names [groupId] holds on [connectionId]; null when the caller may not
 * read them. */
export async function effectivePermissions(connectionId: string, groupId: string): Promise<string[] | null> {
  try {
    const r = await fetch(`${grok.dapi.root}/privileges/permissions?entityId=${encodeURIComponent(connectionId)}` +
      `&groupId=${encodeURIComponent(groupId)}&all=true`,
    {headers: {Authorization: grok.dapi.impersonationToken ?? grok.dapi.token}});
    if (!r.ok)
      return null;
    const rows = await r.json();
    return Array.isArray(rows) ? rows.map((p: any) => String(p?.permission?.name ?? '')).filter((n) => n) : null;
  } catch (_) {
    return null;
  }
}

const groupName = async (id: string): Promise<string> => {
  const g = await grok.dapi.groups.find(id).catch(() => null);
  return g ? g.friendlyName || g.name : id;
};

/** Checks [uses] against the groups of portfolio [code] ([groups]: role → group id). Read-only. */
export async function checkConnections(code: string, groups: Partial<Record<AccessRole, string | null | undefined>>,
  uses: ConnectionUse[]): Promise<AccessCheck[]> {
  const roles = ACCESS_ROLES.filter((r) => groups[r]);
  const names = new Map<string, string>();
  for (const r of roles)
    names.set(groups[r]!, await groupName(groups[r]!));
  const checks: AccessCheck[] = [];
  for (const u of uses) {
    const conn = await findConnection(u.connection).catch(() => null);
    for (const role of roles) {
      const id = groups[role]!;
      let access: Access | null = null;
      let reason: string | undefined;
      if (!conn)
        reason = 'the connection is not found or not shared with you';
      else {
        const perms = await effectivePermissions(conn.id, id);
        if (perms === null)
          reason = 'you may not read the connection\'s sharing';
        else
          access = accessFrom(perms);
      }
      checks.push({portfolio: code, connection: u.connection, uses: u.uses, needed: u.needed, role,
        group: names.get(id) ?? id, access, ...(reason ? {reason} : {})});
    }
  }
  return checks;
}

/** The connections a loaded portfolio uses: its datasets and its enabled feeds. */
export async function connectionsOf(p: PortfolioRow): Promise<ConnectionUse[]> {
  const feeds = await portalDb.registryFeeds.query({filter: `portfolio_id = ${filterValue(p.id)}`, limit: ROW_LIMIT});
  return connectionUses(datasetsOf(p), feeds.filter((f) => f.enabled !== false)
    .map((f) => ({code: f.code, connection: f.connection, target: f.target ?? null})));
}

/** The checks of one loaded portfolio. */
export async function portfolioAccess(p: PortfolioRow): Promise<AccessCheck[]> {
  return checkConnections(p.code, {viewers: p.viewers, contributors: p.contributors, approvers: p.approvers},
    await connectionsOf(p));
}

/** The checks of portfolio [code], or of every portfolio the caller can see. */
export async function dataAccessChecks(code?: string): Promise<AccessCheck[]> {
  const portfolios = code ? [await portfolioByCode(code)] :
    (await portalDb.portfolios.query({limit: ROW_LIMIT})).sort((a, b) => a.code.localeCompare(b.code));
  const out: AccessCheck[] = [];
  for (const p of portfolios)
    out.push(...await portfolioAccess(p));
  return out;
}

/** The checks as a table: one row per portfolio, connection and group. */
export function accessFrame(checks: AccessCheck[]): DG.DataFrame {
  const yn = (c: AccessCheck, k: keyof Access) => c.access === null ? 'unknown' : c.access[k] ? 'yes' : 'no';
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('portfolio', checks.map((c) => c.portfolio)),
    DG.Column.fromStrings('connection', checks.map((c) => c.connection)),
    DG.Column.fromStrings('used by', checks.map((c) => c.uses.join(', '))),
    DG.Column.fromStrings('role', checks.map((c) => c.role)),
    DG.Column.fromStrings('group', checks.map((c) => c.group)),
    DG.Column.fromStrings('view', checks.map((c) => yn(c, 'view'))),
    DG.Column.fromStrings('query', checks.map((c) => yn(c, 'query'))),
    DG.Column.fromStrings('status', checks.map((c) => accessStatus(c))),
    DG.Column.fromStrings('fix', checks.map((c) => accessStatus(c) === 'missing' ?
      accessFix(c.connection, c.group) : c.access === null ? c.reason ?? '' : missingPermissions(c.access).length ?
        'feeds run as whoever syncs; groups need not query it' : '')),
  ]);
  df.name = 'Data access';
  return df;
}

/** The loader's warnings for a portfolio file just loaded: one per group that cannot query a
 * dataset's connection, with the fix. Never throws: a failed check is one line. */
export async function dataAccessWarnings(c: PortfolioConfig, groups: Partial<Record<AccessRole, string | null>>): Promise<string[]> {
  try {
    const feeds = [
      ...c.feeds.map((f) => ({code: f.code, connection: f.connection, target: 'program'})),
      ...(c.compoundFeed ? [{code: COMPOUND_FEED, connection: c.compoundFeed.connection, target: 'compound'}] : []),
      ...(c.workUnitFeeds ?? []).map((f) => ({code: f.code, connection: f.connection, target: f.kind})),
    ];
    const uses = connectionUses(programDatasets(c.programData, c.programDatasets), feeds);
    return describeAccessProblems(await checkConnections(c.code, groups, uses)).map((l) => `warning: ${l}`);
  } catch (e: any) {
    return [`warning: ${c.code}: data access not checked: ${e?.message ?? e}`];
  }
}
