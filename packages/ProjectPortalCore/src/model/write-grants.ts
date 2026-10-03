/** Who the server lets write the portal's table-secured workflow tables. The app checks who may
 * change a role, review a request or park a folder, but package code runs in the user's browser,
 * so that check is not enforcement: anyone can call the domain API from the console. The server
 * enforces table grants: Edit on a table-secured table is needed to insert, update or delete
 * (Delete for deletes). So these tables are View to All users and Edit to the portfolios' groups
 * that the app lets write them.
 *
 * Tables secured through the program row (`master`) need Edit on the program row instead, which
 * only the row's author and admins hold today (GROK-21041); they are not managed here.
 *
 * Pure: the caller reads the current grants and the portfolios' groups. */

/** The group id the platform's "All users" group has on every server. */
export const ALL_USERS = 'a4b45840-9a50-11e6-9cc9-8546b8bf62e6';

export type PortfolioRole = 'approvers' | 'contributors';

/** Which portfolio group writes each table. Approvers are nested in contributors (ensureGroups),
 * so a contributors grant covers the approvers too. */
export const WRITERS: Record<string, PortfolioRole> = {
  program_role: 'approvers',
  parked_folder: 'approvers',
  approval_request: 'contributors',
};

/** Tables every user still writes, and why they are not narrowed. */
export const OPEN_TABLES: Record<string, string> = {
  follow: 'each user writes their own follow; the sync updates everyone\'s. Narrowing it to the ' +
    'author needs row security, and the server refuses a securityMode change on a deployed table',
};

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
  permission: 'Edit' | 'Delete';
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

/** What makes the grants of [WRITERS]' tables right for [portfolios]: Edit to each portfolio's
 * writer group where missing, and All users' Edit and Delete taken away. Grants to other groups
 * are left alone: an admin may have given them on purpose. [current]: each table's grants. */
export function planWriteGrants(portfolios: PortfolioGroups[], current: Record<string, Grant[]>): WriteGrantPlan {
  const plan: WriteGrantPlan = {grant: [], revoke: []};
  for (const table of Object.keys(WRITERS)) {
    const has = new Set((current[table] ?? []).map(([g, p]) => `${g}:${p}`));
    const planned = new Set<string>();
    for (const p of portfolios) {
      const group = writerGroup(p, table);
      if (!group || has.has(`${group}:Edit`) || planned.has(group))
        continue;
      planned.add(group);
      plan.grant.push({table, group, permission: 'Edit', portfolio: p.code});
    }
    for (const permission of ['Edit', 'Delete'] as const) {
      if (has.has(`${ALL_USERS}:${permission}`))
        plan.revoke.push({table, group: ALL_USERS, permission, portfolio: null});
    }
  }
  return plan;
}

/** One line per change, for the report. */
export function describeWriteGrants(plan: WriteGrantPlan, groupName: (id: string) => string = (id) => id): string[] {
  const lines = plan.grant.map((g) => `${g.table}: Edit to ${groupName(g.group)} (${g.portfolio})`);
  lines.push(...plan.revoke.map((r) => `${r.table}: ${r.permission} taken from All users`));
  return lines;
}
