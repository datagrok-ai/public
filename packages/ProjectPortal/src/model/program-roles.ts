/** Program roles (F-16, FR-90): viewers, contributors and approvers of one program, on top of the
 * portfolio's three groups. A role is held by a group or by a user (through the user's personal
 * group); roles nest, so an approver may do what a contributor may, and a contributor what a
 * viewer may. A program role only adds: it never takes away what the portfolio's groups give.
 *
 * Where each role takes effect today:
 * - in the app, as the portfolio's groups do: approvers record gate decisions, contributors and
 *   approvers select runs;
 * - in the program's workspace folder, which accepts grants: View to viewers, Edit to contributors
 *   and approvers.
 * Program rows do not take grants yet (GROK-21041). The row is linked into its folder instead
 * (./program-rows), so a role's View reaches the row through the folder; Edit on the row does not.
 * `grantsFor` says what the row grants will be.
 * Pure: the caller reads the roles and the user's groups. */

export type ProgramRole = 'viewer' | 'contributor' | 'approver';
export const PROGRAM_ROLES: ProgramRole[] = ['viewer', 'contributor', 'approver'];
export const PROGRAM_ROLE_LABELS: Record<ProgramRole, string> = {viewer: 'Viewer', contributor: 'Contributor', approver: 'Approver'};

const RANK: Record<ProgramRole, number> = {viewer: 1, contributor: 2, approver: 3};

/** One role of one principal (a group, or a user's personal group) on a program. */
export interface RoleGrant {
  groupId: string;
  role: ProgramRole;
}

/** The highest of [roles] held by any of [myGroupIds], or null. */
export function programRoleOf(roles: RoleGrant[], myGroupIds: Iterable<string>): ProgramRole | null {
  const mine = new Set(myGroupIds);
  let best: ProgramRole | null = null;
  for (const r of roles) {
    if (mine.has(r.groupId) && (!best || RANK[r.role] > RANK[best]))
      best = r.role;
  }
  return best;
}

/** What the current user may do on a program: what the portfolio's groups give, plus the program's
 * roles. A program role never removes a right. */
export interface ProgramAccess {
  approver: boolean;
  contributor: boolean;
  /** The program role that added to the portfolio's, if any; for saying why. */
  programRole: ProgramRole | null;
}

export function accessOf(portfolio: {approver: boolean, contributor: boolean}, roles: RoleGrant[],
  myGroupIds: Iterable<string>): ProgramAccess {
  const role = programRoleOf(roles, myGroupIds);
  return {
    approver: portfolio.approver || role === 'approver',
    contributor: portfolio.contributor || role === 'approver' || role === 'contributor',
    programRole: role,
  };
}

/** Why the current user may not change a program's roles, empty when they may. The portfolio's
 * approvers manage them; a program approver does not, so a role cannot grant itself more. */
export function checkRoleChange(i: {isPortfolioApprover: boolean, role: string | null}): string[] {
  const reasons: string[] = [];
  if (!i.isPortfolioApprover)
    reasons.push('Only the portfolio\'s approvers manage program roles');
  if (i.role !== null && !PROGRAM_ROLES.includes(i.role as ProgramRole))
    reasons.push(`role "${i.role}" is not one of ${PROGRAM_ROLES.join(', ')}, or none to remove it`);
  return reasons;
}

/** The grants a program's roles imply on an entity of the program (its workspace folder; its row
 * once GROK-21041 is fixed): View to viewers, Edit to contributors and approvers. One per group,
 * the strongest. */
export function grantsFor(roles: RoleGrant[]): {groupId: string, edit: boolean}[] {
  const byGroup = new Map<string, ProgramRole>();
  for (const r of roles) {
    const was = byGroup.get(r.groupId);
    if (!was || RANK[r.role] > RANK[was])
      byGroup.set(r.groupId, r.role);
  }
  return [...byGroup.entries()].map(([groupId, role]) => ({groupId, edit: role !== 'viewer'}));
}

/** What to change on an entity whose current grants are [current] (group id → edit) so that it
 * carries [wanted] for the groups [managed] covers: grants to add or change, and grants to
 * remove. A group [managed] does not cover (the portfolio's groups, the author) is left alone. */
export function grantChanges(current: Map<string, boolean>, wanted: {groupId: string, edit: boolean}[],
  managed: Iterable<string>): {grant: {groupId: string, edit: boolean}[], revoke: string[]} {
  const want = new Set(wanted.map((w) => w.groupId));
  const grant = wanted.filter((w) => current.get(w.groupId) !== w.edit);
  const revoke = [...new Set(managed)].filter((g) => !want.has(g) && current.has(g));
  return {grant, revoke};
}

/** "Approver: ann; Viewer: Bio, Chemists" */
export function describeRoles(roles: {name: string, role: ProgramRole}[]): string {
  if (!roles.length)
    return 'No program roles: the portfolio\'s groups decide';
  return PROGRAM_ROLES.slice().reverse().map((role) => {
    const names = roles.filter((r) => r.role === role).map((r) => r.name).sort();
    return names.length ? `${PROGRAM_ROLE_LABELS[role]}: ${names.join(', ')}` : null;
  }).filter((x) => x).join('; ');
}
