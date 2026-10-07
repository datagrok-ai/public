/** Who the server lets write and read the portal's tables, as far as table grants decide it. The
 * app checks who may change a role, review a request or park a folder, but package code runs in
 * the user's browser, so that check is not enforcement: anyone can call the domain API from the
 * console. The server enforces table grants: Edit on the table is needed to insert, update or
 * delete (Delete for deletes), View to read a table-secured table at all.
 *
 * A table grant covers every row, so it cannot keep one portfolio's rows from another's. Records
 * of one portfolio or one user therefore live in row-secured tables (approval, sync_log,
 * following) or under the program (program_compound), and the old table-secured tables are
 * closed to All users ([CLOSED]).
 *
 * Tables secured through the program row (`master`) need Edit on the program row instead, which
 * only the row's author and admins hold today (GROK-21041); they are not managed here.
 *
 * Pure: the caller reads the current grants and the portfolios' groups. */

/** The group id the platform's "All users" group has on every server. */
export const ALL_USERS = 'a4b45840-9a50-11e6-9cc9-8546b8bf62e6';

export type PortfolioRole = 'approvers' | 'contributors';

/** Which portfolio group writes each table. Approvers are nested in contributors (ensureGroups),
 * so a contributors grant covers the approvers too. `approval` and `sync_log` are row-secured with
 * `defaultRowVisibility: table`: Edit on the table lets the group write any row, while reads go
 * per row (the folder a row is linked into), since no group holds View on the table. */
export const WRITERS: Record<string, PortfolioRole> = {
  program_role: 'approvers',
  parked_folder: 'approvers',
  approval: 'contributors',
  sync_log: 'approvers',
};

/** Which portfolio group reads each table that All users no longer read. Table security is one
 * grant for the whole table, so the group reads every portfolio's rows. */
export const READERS: Record<string, PortfolioRole> = {
  parked_folder: 'approvers',
};

/** Tables All users must not read or write, and why. A table here loses All users' View, Edit and
 * Delete; a redeploy never revokes a grant, so this runs on load. */
export const CLOSED: Record<string, string> = {
  approval_request: 'superseded by approval: every user read every portfolio\'s requests and rationales',
  sync_run: 'superseded by sync_log: every user read every portfolio\'s sync history',
  follow: 'superseded by following: every user read everyone\'s follows',
  compound: 'structures moved to program_compound, under each program\'s security',
  parked_folder: 'one portfolio\'s folder names and reasons; its approvers read it',
  approval: 'a View grant on the table would show every row (defaultRowVisibility: table)',
  sync_log: 'a View grant on the table would show every row (defaultRowVisibility: table)',
};

/** Tables every user still reads or writes, and why they are not narrowed. */
export const OPEN_TABLES: Record<string, string> = {
  following: 'each user writes their own follows; row security with no row grants keeps everyone ' +
    'else\'s row out of reach',
  program_role: 'every user reads which groups hold a role in a program to know their own access; ' +
    'it holds group ids and roles, no program content',
};

/** Every table whose grants are managed here. */
export const MANAGED_TABLES = [...new Set([...Object.keys(WRITERS), ...Object.keys(READERS), ...Object.keys(CLOSED)])];

export interface PortfolioGroups {
  code: string;
  approvers?: string | null;
  contributors?: string | null;
}

/** One permission row on a table: [group id, permission]. */
export type Grant = [string, string];

export interface GrantChange {
  table: string;
  group: string;
  permission: 'View' | 'Edit' | 'Delete';
  /** The portfolio the grant is for, or null for All users. */
  portfolio: string | null;
}

export interface WriteGrantPlan {
  grant: GrantChange[];
  revoke: GrantChange[];
}

/** The group of [p] that writes [table]: its contributors, or its approvers when it names none. */
export function writerGroup(p: PortfolioGroups, table: string): string | null {
  const role = WRITERS[table];
  if (!role)
    return null;
  return (role === 'contributors' ? p.contributors || p.approvers : p.approvers) || null;
}

/** The group of [p] that reads [table], when the table has one ([READERS]). */
export function readerGroup(p: PortfolioGroups, table: string): string | null {
  const role = READERS[table];
  if (!role)
    return null;
  return (role === 'contributors' ? p.contributors || p.approvers : p.approvers) || null;
}

/** What makes the grants right for [portfolios]: Edit to each portfolio's writer group and View to
 * its reader group where missing; All users' Edit and Delete taken from [WRITERS]' tables and
 * View, Edit and Delete from [CLOSED]'s. Grants to other groups are left alone: an admin may have
 * given them on purpose. [current]: each table's grants. */
export function planWriteGrants(portfolios: PortfolioGroups[], current: Record<string, Grant[]>): WriteGrantPlan {
  const plan: WriteGrantPlan = {grant: [], revoke: []};
  for (const table of MANAGED_TABLES) {
    const has = new Set((current[table] ?? []).map(([g, p]) => `${g}:${p}`));
    for (const [permission, groupOf] of [['Edit', writerGroup], ['View', readerGroup]] as const) {
      const planned = new Set<string>();
      for (const p of portfolios) {
        const group = groupOf(p, table);
        if (!group || has.has(`${group}:${permission}`) || planned.has(group))
          continue;
        planned.add(group);
        plan.grant.push({table, group, permission, portfolio: p.code});
      }
    }
    const taken = CLOSED[table] ? ['View', 'Edit', 'Delete'] as const : WRITERS[table] ? ['Edit', 'Delete'] as const : [];
    for (const permission of taken) {
      if (has.has(`${ALL_USERS}:${permission}`))
        plan.revoke.push({table, group: ALL_USERS, permission, portfolio: null});
    }
  }
  return plan;
}

/** One line per change, for the report. */
export function describeWriteGrants(plan: WriteGrantPlan, groupName: (id: string) => string = (id) => id): string[] {
  const lines = plan.grant.map((g) => `${g.table}: ${g.permission} to ${groupName(g.group)} (${g.portfolio})`);
  lines.push(...plan.revoke.map((r) => `${r.table}: ${r.permission} taken from All users`));
  return lines;
}
