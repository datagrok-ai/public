/** Whether a portfolio's groups can query the data connections its datasets and feeds use. A
 * program dataset runs its SQL in the viewer's own session, so the server asks for View and
 * `DataConnection.Query` on the connection: without Query a plain viewer gets 'Connection "Query"
 * permission required'. The portal shares its own tables and folders when a portfolio loads, but a
 * connection belongs to whoever set it up, so this only reports: it never grants.
 *
 * Feeds run in the session of whoever syncs (an admin or the scheduler), not the groups', so a
 * connection only feeds use is reported but not required of the groups.
 *
 * Pure: the caller resolves the connections and reads each group's effective grants. */

export type AccessRole = 'viewers' | 'contributors' | 'approvers';
export const ACCESS_ROLES: AccessRole[] = ['viewers', 'contributors', 'approvers'];

/** The platform's names of the two permissions a dataset query needs. */
export const VIEW = 'View';
export const QUERY = 'DataConnection.Query';

/** One connection a portfolio uses, as its file spells it, and what uses it. */
export interface ConnectionUse {
  connection: string;
  /** `dataset program`, `dataset assays`, `feed REG`, `compound feed`, ... in declaration order. */
  uses: string[];
  /** Whether the groups run it: true when a dataset uses it, false when only feeds do. */
  needed: boolean;
}

/** What a group may do on a connection, from its effective grants. */
export interface Access {
  view: boolean;
  query: boolean;
}

export type AccessStatus = 'ok' | 'missing' | 'unknown' | 'not needed';

/** One group of one portfolio against one connection. [access] null: the grants could not be
 * read (the connection is not visible to the caller, or the caller may not read its sharing);
 * [reason] says why. */
export interface AccessCheck {
  portfolio: string;
  connection: string;
  uses: string[];
  needed: boolean;
  role: AccessRole;
  group: string;
  access: Access | null;
  reason?: string;
}

/** The connections of a portfolio's datasets and feeds, one entry per connection, in order. */
export function connectionUses(datasets: {code: string, connection: string}[],
  feeds: {code: string, connection: string, target?: string | null}[]): ConnectionUse[] {
  const out = new Map<string, ConnectionUse>();
  const add = (connection: string, use: string, needed: boolean) => {
    if (!connection)
      return;
    const u = out.get(connection) ?? {connection, uses: [], needed: false};
    if (!u.uses.includes(use))
      u.uses.push(use);
    u.needed = u.needed || needed;
    out.set(connection, u);
  };
  for (const d of datasets)
    add(d.connection, `dataset ${d.code}`, true);
  for (const f of feeds)
    add(f.connection, !f.target || f.target === 'program' ? `feed ${f.code}` : `${f.target} feed ${f.code}`, false);
  return [...out.values()];
}

/** View and Query from a group's effective permission names on the connection: its own grants,
 * its parent groups' (All users among them) and the ones inherited through projects, as the
 * server lists them. Edit does not stand in for View: the server checks each by name. */
export function accessFrom(permissions: string[]): Access {
  return {view: permissions.includes(VIEW), query: permissions.includes(QUERY)};
}

export function accessStatus(c: AccessCheck): AccessStatus {
  if (c.access === null)
    return 'unknown';
  if (c.access.view && c.access.query)
    return 'ok';
  return c.needed ? 'missing' : 'not needed';
}

/** The permissions the group lacks, as the sharing dialog names them: View, Query. */
export function missingPermissions(a: Access | null): string[] {
  if (a === null)
    return [];
  return [...(a.view ? [] : ['View']), ...(a.query ? [] : ['Query'])];
}

/** The one-line fix for a group that cannot query a connection. */
export function accessFix(connection: string, group: string): string {
  return `share ${connection} with ${group}: View + Query`;
}

/** One warning per group that cannot query a connection a dataset uses, and per check that could
 * not be made. Groups that are fine, and feed-only connections, say nothing. */
export function describeAccessProblems(checks: AccessCheck[]): string[] {
  const lines: string[] = [];
  for (const c of checks) {
    const status = accessStatus(c);
    const where = `${c.portfolio}: connection ${c.connection} (${c.uses.join(', ')})`;
    if (status === 'missing')
      lines.push(`${where}: ${c.group} (${c.role}) lacks ${missingPermissions(c.access).join(' and ')}; ` +
        `its members cannot open the data. Fix: ${accessFix(c.connection, c.group)}`);
    else if (status === 'unknown' && c.needed)
      lines.push(`${where}: access of ${c.group} (${c.role}) unknown: ${c.reason ?? 'its sharing could not be read'}`);
  }
  return lines;
}

/** The board's notice: null when every group can query every dataset connection. */
export function accessNotice(checks: AccessCheck[]): {missing: number, unknown: number, text: string} | null {
  const needed = checks.filter((c) => c.needed);
  const missing = needed.filter((c) => accessStatus(c) === 'missing');
  const unknown = needed.filter((c) => accessStatus(c) === 'unknown');
  if (!missing.length && !unknown.length)
    return null;
  const conns = (cs: AccessCheck[]) => [...new Set(cs.map((c) => c.connection))].join(', ');
  const parts: string[] = [];
  if (missing.length) {
    const groups = [...new Set(missing.map((c) => c.group))].join(', ');
    parts.push(`${groups} cannot query ${conns(missing)}, so the program data does not open for them`);
  }
  if (unknown.length)
    parts.push(`access to ${conns(unknown)} could not be checked`);
  return {missing: missing.length, unknown: unknown.length, text: parts.join('; ')};
}

/** The whole audit as text: per portfolio and connection, each group's status, then the fixes. */
export function summarizeAccess(checks: AccessCheck[]): string {
  if (!checks.length)
    return 'no data connections to check';
  const lines: string[] = [];
  let last = '';
  for (const c of checks) {
    const key = `${c.portfolio}|${c.connection}`;
    if (key !== last)
      lines.push(`${c.portfolio}: ${c.connection} (${c.uses.join(', ')})${c.needed ? '' : ' — feeds only, run as whoever syncs'}`);
    last = key;
    const status = accessStatus(c);
    const detail = status === 'unknown' ? `: ${c.reason ?? 'its sharing could not be read'}` :
      status === 'ok' ? '' : `: lacks ${missingPermissions(c.access).join(' and ')}`;
    lines.push(`  ${c.group} (${c.role}): ${status}${detail}`);
  }
  const fixes = checks.filter((c) => accessStatus(c) === 'missing').map((c) => accessFix(c.connection, c.group));
  if (fixes.length)
    lines.push('fixes:', ...[...new Set(fixes)].map((f) => `  ${f}`));
  return lines.join('\n');
}
