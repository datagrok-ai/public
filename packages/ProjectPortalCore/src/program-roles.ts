/** Program roles against the platform (F-16): the `program_role` rows, the current user's access
 * to a program, and the program's workspace folder shared to match. The rules are in
 * ./model/program-roles.
 *
 * Row grants wait on GROK-21041: the server refuses `permissions.grant` on a promoted domain row,
 * so a program role cannot make a program row visible yet. `syncRowGrants` applies the same grants
 * to the row as to the folder; it runs only when ROW_GRANTS is true. Once the platform accepts the
 * grant, set ROW_GRANTS to true: the `program_role` table then drives the row's grants too, with no
 * change to the table or to who manages it. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {portalDb, PortfolioRow, ProgramRoleRow, ProgramRow} from './generated/db';
import {isApprover, isContributor} from './portal';
import {accessOf, checkRoleChange, grantChanges, grantsFor, ProgramAccess, ProgramRole, RoleGrant} from './model/program-roles';
import {filterValue, ROW_LIMIT} from './rows';

/** Whether program roles also grant on the program row. Off until GROK-21041 is fixed. */
export const ROW_GRANTS = false;

/** The program's roles in force (not removed). */
export async function rolesOf(programId: string): Promise<ProgramRoleRow[]> {
  return (await portalDb.programRoles.query({filter: `program_id = ${filterValue(programId)}`, limit: ROW_LIMIT})).filter((r) => r.active !== false);
}

const asGrants = (rows: ProgramRoleRow[]): RoleGrant[] => rows.map((r) => ({groupId: r.group_id, role: r.role as ProgramRole}));

/** The ids of the current user's groups, their personal group included. */
export async function myGroupIds(): Promise<Set<string>> {
  const ids = new Set((await grok.dapi.groups.currentUserGroups()).map((g) => g.id));
  const personal = DG.User.current().group?.id;
  if (personal)
    ids.add(personal);
  return ids;
}

/** What the current user may do on a program: the portfolio's groups, plus the program's roles. */
export async function programAccess(portfolio: PortfolioRow, programId: string): Promise<ProgramAccess> {
  const [approver, contributor, roles, groups] = await Promise.all([isApprover(portfolio), isContributor(portfolio),
    rolesOf(programId).catch(() => [] as ProgramRoleRow[]), myGroupIds()]);
  return accessOf({approver, contributor}, asGrants(roles), groups);
}

/** A group by its name, or a user by login (their personal group). */
export async function resolvePrincipal(name: string): Promise<DG.Group> {
  const n = name.trim();
  if (!n)
    throw new Error('Name a group or a user');
  const group = await grok.dapi.groups.filter(`friendlyName = "${n}"`).first().catch(() => null) ??
    await grok.dapi.groups.filter(`name = "${n}"`).first().catch(() => null);
  if (group)
    return group;
  const user = await grok.dapi.users.filter(`login = "${n}"`).include('group').first().catch(() => null);
  if (user?.group)
    return user.group;
  throw new Error(`No group or user "${n}"`);
}

/** The program's roles with the names of their groups, strongest first. */
export async function namedRolesOf(programId: string): Promise<(ProgramRoleRow & {name: string})[]> {
  const rows = await rolesOf(programId);
  const names = new Map<string, string>();
  for (const id of new Set(rows.map((r) => r.group_id))) {
    const g = await grok.dapi.groups.find(id).catch(() => null);
    names.set(id, g ? g.friendlyName || g.name : id);
  }
  const rank = {approver: 0, contributor: 1, viewer: 2} as Record<string, number>;
  return rows.map((r) => ({...r, name: names.get(r.group_id) ?? r.group_id}))
    .sort((a, b) => rank[a.role] - rank[b.role] || a.name.localeCompare(b.name));
}

export interface RoleChange {
  line: string;
  /** What happened to the workspace folder's sharing. */
  notes: string[];
}

/** Gives [groupId] [role] on [program], or with null removes its role. Only the portfolio's
 * approvers may (asked here, not only by the page). The program's workspace folder is shared to
 * match; a failure there is reported, not thrown: the role stands. */
export async function setProgramRole(program: ProgramRow, groupId: string, role: ProgramRole | null): Promise<RoleChange> {
  const portfolio = await portalDb.portfolios.get(program.portfolio_id);
  const refused = checkRoleChange({isPortfolioApprover: await isApprover(portfolio), role});
  if (refused.length)
    throw new Error(refused.join('; '));
  const existing = (await portalDb.programRoles.query({filter: `program_id = ${filterValue(program.id)} and group_id = ${filterValue(groupId)}`, limit: ROW_LIMIT}))[0];
  const stamp = {granted_by: DG.User.current().id, granted_on: dayjs()};
  if (role === null) {
    if (existing && existing.active !== false)
      await portalDb.programRoles.update(existing.id, {active: false, ...stamp});
  } else if (existing)
    await portalDb.programRoles.update(existing.id, {role, active: true, ...stamp});
  else
    await portalDb.programRoles.insert({program_id: program.id, group_id: groupId, role, active: true, ...stamp});
  const group = await grok.dapi.groups.find(groupId).catch(() => null);
  const who = group?.friendlyName || group?.name || groupId;
  const notes: string[] = [];
  try {
    notes.push(...await syncFolderGrants(program, portfolio));
  } catch (e: any) {
    notes.push(`workspace folder not shared: ${e?.message ?? e}`);
  }
  if (ROW_GRANTS) {
    try {
      notes.push(...await syncRowGrants(program, portfolio));
    } catch (e: any) {
      notes.push(`program row not shared: ${e?.message ?? e}`);
    }
  }
  return {line: `${program.code}: ${who} ${role ? `is ${role}` : 'has no program role'}`, notes};
}

/** The groups program roles manage on an entity: every group that holds or held a role on the
 * program, never the portfolio's own groups, whose grants the sync and the builder make. */
async function managedGroups(program: ProgramRow, portfolio: PortfolioRow): Promise<string[]> {
  const all = await portalDb.programRoles.query({filter: `program_id = ${filterValue(program.id)}`, limit: ROW_LIMIT});
  const own = new Set([portfolio.viewers, portfolio.contributors, portfolio.approvers].filter((x) => x));
  return [...new Set(all.map((r) => r.group_id))].filter((g) => !own.has(g));
}

/** Makes [entity]'s grants match the program's roles: View to viewers, Edit to contributors and
 * approvers, nothing to a group whose role was removed. Other grants are left alone. */
async function syncGrants(entity: DG.Entity, program: ProgramRow, portfolio: PortfolioRow, what: string): Promise<string[]> {
  const managed = await managedGroups(program, portfolio);
  if (!managed.length)
    return [];
  const own = new Set([portfolio.viewers, portfolio.contributors, portfolio.approvers].filter((x) => x));
  const wanted = grantsFor(asGrants(await rolesOf(program.id))).filter((g) => !own.has(g.groupId));
  // a space nobody has been granted comes back without the view/edit keys
  const perms = await grok.dapi.permissions.get(entity) as {view?: DG.Group[], edit?: DG.Group[]};
  const current = new Map<string, boolean>();
  for (const g of perms.view ?? [])
    current.set(g.id, false);
  for (const g of perms.edit ?? [])
    current.set(g.id, true);
  const changes = grantChanges(current, wanted, managed);
  const notes: string[] = [];
  for (const g of changes.grant) {
    const group = await grok.dapi.groups.find(g.groupId);
    // a grant only adds; going from edit to view needs the edit taken away first
    if (current.has(g.groupId))
      await grok.dapi.permissions.revoke(entity, group);
    await grok.dapi.permissions.grant(entity, group, g.edit);
    notes.push(`${what} shared with ${group.friendlyName || group.name} (${g.edit ? 'edit' : 'view'})`);
  }
  for (const id of changes.revoke) {
    const group = await grok.dapi.groups.find(id);
    await grok.dapi.permissions.revoke(entity, group);
    notes.push(`${what} no longer shared with ${group.friendlyName || group.name}`);
  }
  return notes;
}

/** Shares the program's workspace folder to match its roles. Without a folder there is nothing to
 * share; the builder applies the roles when it creates one. */
export async function syncFolderGrants(program: ProgramRow, portfolio?: PortfolioRow, folder?: DG.Entity): Promise<string[]> {
  const pf = portfolio ?? await portalDb.portfolios.get(program.portfolio_id);
  const space = folder ?? (program.space_id ? (await grok.dapi.getEntities([program.space_id]))[0] : null);
  if (!space)
    return ['no workspace folder yet: building it applies the program roles'];
  return syncGrants(space, program, pf, 'workspace folder');
}

/** The program row's grants from its roles. The server refuses them today (GROK-21041); runs only
 * when ROW_GRANTS is true. */
export async function syncRowGrants(program: ProgramRow, portfolio?: PortfolioRow): Promise<string[]> {
  const pf = portfolio ?? await portalDb.portfolios.get(program.portfolio_id);
  await portalDb.programs.promote(program.id);
  const [entity] = await grok.dapi.getEntities([program.id]);
  if (!entity)
    throw new Error(`program ${program.code} is not an entity after promotion`);
  return syncGrants(entity, program, pf, 'program row');
}

/** One row per role: group or user, role, who granted it, when. */
export function programRolesFrame(rows: (ProgramRoleRow & {name: string})[], name: string): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('group', rows.map((r) => r.name)),
    DG.Column.fromStrings('role', rows.map((r) => r.role)),
    DG.Column.fromStrings('granted on', rows.map((r) => r.granted_on ? dayjs(r.granted_on as any).format('YYYY-MM-DD HH:mm') : '')),
  ]);
  df.name = name;
  return df;
}
