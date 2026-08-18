import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {
  APPROVAL_COLUMNS, CURATION_COLUMNS, programGroupNames,
  T_ALIGNMENT, T_COMPOUND, T_PROGRAM, T_TAG,
  UMBRELLA_APPROVERS, UMBRELLA_CONTRIBUTORS,
} from './constants';
import {shareEntity} from './rest';

export async function findGroup(name: string): Promise<DG.Group | null> {
  const groups = await grok.dapi.groups.filter(`friendlyName = "${name}"`).list();
  return groups.find((g) => g.friendlyName === name) ?? null;
}

export async function ensureGroup(name: string): Promise<DG.Group> {
  return await findGroup(name) ?? await grok.dapi.groups.createNew(name);
}

async function allUsers(): Promise<DG.Group> {
  return grok.dapi.groups.find(DG.Group.defaultGroupsIds['All users']);
}

/** One-time (idempotent) security setup for the schema: umbrella groups, registry
 * visibility, and per-column write gating on the approval/curation columns.
 *
 * Column security is table-wide, so Edit on the approval columns goes to an umbrella
 * group containing every program's approvers. Per-program authority (an Alpha approver
 * must not approve Beta rows) additionally requires Edit on the row's delegate program
 * row, and the service re-checks membership in the row's own approvers group before
 * flipping status. */
export async function setupSchemaSecurity(): Promise<{approvers: DG.Group, contributors: DG.Group}> {
  const approvers = await ensureGroup(UMBRELLA_APPROVERS);
  const contributors = await ensureGroup(UMBRELLA_CONTRIBUTORS);
  const all = await allUsers();

  const compound = grok.dapi.domains.table(T_COMPOUND);
  const tag = grok.dapi.domains.table(T_TAG);
  await compound.grant(all.id, 'View');
  await tag.grant(all.id, 'View');
  await tag.grant(contributors.id, 'Edit');

  const alignment = grok.dapi.domains.table(T_ALIGNMENT);
  for (const column of APPROVAL_COLUMNS) {
    await alignment.shareColumn(column, all.id, 'View');
    await alignment.shareColumn(column, approvers.id, 'Edit');
  }
  for (const column of CURATION_COLUMNS) {
    await alignment.shareColumn(column, all.id, 'View');
    await alignment.shareColumn(column, contributors.id, 'Edit');
  }
  return {approvers, contributors};
}

let _setupDone: Promise<{approvers: DG.Group, contributors: DG.Group}> | null = null;
export function setupSchemaSecurityOnce(): Promise<{approvers: DG.Group, contributors: DG.Group}> {
  return _setupDone ??= setupSchemaSecurity();
}

export interface ProgramSpec {
  code: string;
  name: string;
  description?: string;
  externalRef?: string;
  stage?: string;
  status?: string;
  indication?: string;
}

export interface ProgramInfo {
  id: string;
  code: string;
  viewers: DG.Group;
  contributors: DG.Group;
  approvers: DG.Group;
}

/** Idempotently creates a program with its three audience groups and row grants:
 * viewers get View on the program row (alignment rows inherit it through the master
 * delegate); contributors and approvers get Edit (write-coupling accepted per the
 * design — inserting a version row routes through the delegate). */
export async function ensureProgram(spec: ProgramSpec): Promise<ProgramInfo> {
  const names = programGroupNames(spec.code);
  const viewers = await ensureGroup(names.viewers);
  const contributors = await ensureGroup(names.contributors);
  const approvers = await ensureGroup(names.approvers);

  const umbrellas = await setupSchemaSecurityOnce();
  await grok.dapi.groups.includeTo(contributors, umbrellas.contributors);
  await grok.dapi.groups.includeTo(approvers, umbrellas.approvers);

  const programs = grok.dapi.domains.table(T_PROGRAM);
  const existing = await programs.getByKey({code: spec.code});
  const values = {
    code: spec.code, name: spec.name, description: spec.description,
    external_ref: spec.externalRef, stage: spec.stage, status: spec.status,
    indication: spec.indication,
    viewers_group: viewers.id, contributors_group: contributors.id, approvers_group: approvers.id,
  };
  let programId: string;
  if (existing == null)
    programId = (await programs.insert(values))[0].id;
  else {
    await programs.update(existing.id, values);
    programId = existing.id;
  }

  await programs.promote(programId);
  await shareEntity(programId, names.viewers, 'View');
  await shareEntity(programId, names.contributors, 'Edit');
  await shareEntity(programId, names.approvers, 'Edit');
  return {id: programId, code: spec.code, viewers, contributors, approvers};
}

export async function getProgramGroups(programId: string):
    Promise<{viewers: DG.Group | null, approvers: DG.Group | null, contributors: DG.Group | null}> {
  const program = await grok.dapi.domains.table(T_PROGRAM).get(programId);
  const load = async (id?: string | null) => id == null ? null : grok.dapi.groups.find(id);
  return {
    viewers: await load(program?.viewers_group),
    approvers: await load(program?.approvers_group),
    contributors: await load(program?.contributors_group),
  };
}

/** Transitive membership check of the current user in [group] (walks membership
 * parents breadth-first, capped at 6 hops). */
export async function currentUserIn(group: DG.Group | null): Promise<boolean> {
  if (group == null)
    return false;
  const user = await grok.dapi.users.current();
  const seen = new Set<string>();
  let frontier: string[] = user.group?.id != null ? [user.group.id] : [];
  for (let depth = 0; depth < 6 && frontier.length > 0; depth++) {
    if (frontier.includes(group.id))
      return true;
    const next: string[] = [];
    for (const id of frontier) {
      if (seen.has(id))
        continue;
      seen.add(id);
      const g = await grok.dapi.groups.find(id);
      for (const parent of [...g?.memberships ?? [], ...g?.adminMemberships ?? []]) {
        if (!seen.has(parent.id))
          next.push(parent.id);
      }
    }
    frontier = next;
  }
  return frontier.includes(group.id);
}
