import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {
  APPROVAL_COLUMNS, CURATION_COLUMNS, programGroupNames,
  T_ALIGNMENT, T_ALIGNMENT_COMPOUND, T_ALIGNMENT_TAG, T_COMPOUND, T_HISTORY,
  T_PROGRAM, T_PROGRAM_COMPOUND, T_STUDY, T_TAG,
  UMBRELLA_APPROVERS, UMBRELLA_CONTRIBUTORS, UMBRELLA_VIEWERS,
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

/** One-time (idempotent) security setup: umbrella groups, registry visibility, column
 * gating — docs/artifact-alignment.html § Permissions (as implemented). */
export async function setupSchemaSecurity():
    Promise<{viewers: DG.Group, approvers: DG.Group, contributors: DG.Group}> {
  const viewers = await ensureGroup(UMBRELLA_VIEWERS);
  const approvers = await ensureGroup(UMBRELLA_APPROVERS);
  const contributors = await ensureGroup(UMBRELLA_CONTRIBUTORS);
  const all = await allUsers();

  // Grants are table-wide by platform necessity; per-program isolation is
  // service-side — doc § Permissions and § Known design gaps.
  for (const table of [T_PROGRAM, T_STUDY, T_ALIGNMENT, T_HISTORY,
    T_PROGRAM_COMPOUND, T_ALIGNMENT_COMPOUND, T_ALIGNMENT_TAG]) {
    const client = grok.dapi.domains.table(table);
    for (const group of [viewers, contributors, approvers])
      await client.grant(group.id, 'View');
    for (const group of [contributors, approvers])
      await client.grant(group.id, 'Edit');
    try {
      await client.revoke(all.id, 'View');
    } catch (_) {/* nothing to revoke */}
  }

  const compound = grok.dapi.domains.table(T_COMPOUND);
  const tag = grok.dapi.domains.table(T_TAG);
  await compound.grant(all.id, 'View');
  await tag.grant(all.id, 'View');
  await tag.grant(contributors.id, 'Edit');
  // publish creates missing compounds by registration code, same as tags
  await compound.grant(contributors.id, 'Edit');

  // The archive delete authorizes against the SECURING program table, not the
  // alignment table (verified) — doc § Permissions.
  const programs = grok.dapi.domains.table(T_PROGRAM);
  for (const group of [contributors, approvers])
    await programs.grant(group.id, 'Delete');

  const alignment = grok.dapi.domains.table(T_ALIGNMENT);
  for (const column of APPROVAL_COLUMNS) {
    await alignment.shareColumn(column, all.id, 'View');
    await alignment.shareColumn(column, approvers.id, 'Edit');
  }
  for (const column of CURATION_COLUMNS) {
    await alignment.shareColumn(column, all.id, 'View');
    await alignment.shareColumn(column, contributors.id, 'Edit');
  }
  return {viewers, approvers, contributors};
}

let _setupDone: Promise<{viewers: DG.Group, approvers: DG.Group, contributors: DG.Group}> | null = null;
export function setupSchemaSecurityOnce():
    Promise<{viewers: DG.Group, approvers: DG.Group, contributors: DG.Group}> {
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

/** Idempotently creates a program with its audience groups, umbrella nesting, and
 * row grants — doc § Permissions (as implemented). */
export async function ensureProgram(spec: ProgramSpec): Promise<ProgramInfo> {
  const names = programGroupNames(spec.code);
  const viewers = await ensureGroup(names.viewers);
  const contributors = await ensureGroup(names.contributors);
  const approvers = await ensureGroup(names.approvers);

  const umbrellas = await setupSchemaSecurityOnce();
  await grok.dapi.groups.includeTo(viewers, umbrellas.viewers);
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
