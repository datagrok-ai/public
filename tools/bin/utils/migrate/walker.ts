/// Docs: [Entity export / import](/docs/features/grok-tool/export-import/DESIGN.md)
import {NodeDapi} from '../node-dapi';
import {BytesKind, TYPES, Ref, TypeOptions, UUID_RE, inNamespace, isBuiltinGroup, nqNameOf, untransferableReason} from './registry';
import {BundleEntity} from './bundle';

export interface Selection {
  types: string[];
  names: string[];
  name?: string;
  namespace?: string;
  space?: string;
  author?: string;
  tag?: string;
  since?: string;
  filter?: string;
  /** Per-type listing rules of the `--type` aliases (`space`, `dashboard`). */
  typeOptions?: Record<string, TypeOptions>;
}

export type Note = (row: {name: string; entityType: string; action: 'warn' | 'info'; reason: string; detail?: string}) => void;

const ALL_USERS = 'a4b45840-9a50-11e6-9cc9-8546b8bf62e6';

const quoted = (v: string): string => `"${v.replace(/"/g, '\\"')}"`;

/**
 * Compiles the selection flags into one smart-filter expression. A structured `name`
 * clause is unreliable — the server rewrites `name` to `friendlyName`
 * (`repository_query.dart` whereSmart) and `/projects` never resolves it — so `--name`
 * becomes a free-text search, narrowed client-side by `matchesName`. Free text is not
 * part of the grammar, so it is sent only when it is the whole filter.
 */
export function compileFilter(sel: Selection, type?: string): string {
  const clauses: string[] = [];
  const typeClause = type ? sel.typeOptions?.[type]?.clause : undefined;
  const namespace = sel.namespace ?? sel.space;
  if (namespace)
    clauses.push(`namespace starts ${quoted(`${namespace.replace(/:$/, '')}:`)}`);
  if (sel.author)
    clauses.push(`author = ${quoted(sel.author)}`);
  if (sel.since)
    clauses.push(`updatedOn > ${sel.since}`);
  if (sel.filter)
    clauses.push(`(${sel.filter})`);
  if (sel.name && !clauses.length && !typeClause)
    return sel.name.replace(/[*?]/g, ' ').trim();
  if (typeClause)
    clauses.push(typeClause);
  return clauses.join(' and ');
}

/** `--name` glob, matched against both the internal name and the user-visible one. */
export function matchesName(glob: string, json: any): boolean {
  const pattern = glob.replace(/[.+^${}()|[\]\\]/g, '\\$&').replace(/\*/g, '.*').replace(/\?/g, '.');
  const re = new RegExp(`^${pattern}$`, 'i');
  return re.test(json?.name ?? '') || re.test(json?.friendlyName ?? '');
}

/**
 * minimist swallows `--since -2w` (the value parses as a flag), so `2w` means `-2w`.
 * An absolute date has to reach the filter quoted — `updatedOn > 2020-01-01` parses
 * but matches nothing.
 */
export function normalizeSince(value: any): string | undefined {
  if (value === undefined || value === null) return undefined;
  if (value === true || value === '')
    throw new Error('`--since` needs a value: `--since=-2w` or `--since 2w`');
  const text = String(value).trim();
  if (/^-?\d+[dwmy]$/.test(text))
    return text.startsWith('-') ? text : `-${text}`;
  if (/^\d{4}-\d{2}-\d{2}([T ][\d:.]+Z?)?$/.test(text))
    return quoted(text);
  throw new Error(`\`--since\` takes a timespan ('2w', '-30d') or an ISO date ('2026-08-01'), got '${text}'`);
}

/** `GET /entities/{id}` answers with an array holding one typed entity. */
export async function findEntity(dapi: NodeDapi, id: string): Promise<any> {
  const found = await dapi.internal('/entities').find(id);
  return (Array.isArray(found) ? found[0] : found) ?? null;
}

/** One unreachable entity is a warning, not the end of the run. */
async function tryFind(dapi: NodeDapi, type: string, id: string, name: string, note: Note): Promise<any> {
  try {
    return await dapi.internal(TYPES[type].route).find(id);
  } catch (err: any) {
    note({name, entityType: type, action: 'warn', reason: 'fetch_failed', detail: err?.message ?? String(err)});
    return null;
  }
}

/** UUID, or `namespace:name` (`name` alone for entities in the root namespace). */
export async function resolveEntity(dapi: NodeDapi, token: string, type?: string): Promise<any> {
  if (UUID_RE.test(token)) {
    const one = await findEntity(dapi, token);
    if (!one) throw new Error(`No entity with id '${token}'`);
    return one;
  }
  const cut = token.lastIndexOf(':');
  const namespace = cut === -1 ? '' : token.slice(0, cut + 1);
  const name = token.slice(cut + 1);
  // The namespace is sent even when empty — omitting it matches the name in every namespace.
  const matches = (await dapi.internal('/entities').list({namespace, name}))
    .filter((m: any) => inNamespace(m, namespace) && (!type || m['#type'] === type));
  if (!matches.length)
    throw new Error(`No ${type ?? 'entity'} named '${token}'` +
      (cut === -1 ? ` in the root namespace — qualify it, e.g. Admin:${name}` : ''));
  if (matches.length > 1)
    throw new Error(`'${token}' matches ${matches.length} entities: ${matches.map((m: any) => `${m['#type']} ${m.id}`).join(', ')}`);
  return matches[0];
}

/**
 * What never travels, however the entity was chosen — `--no-deps` skips the walk, not these
 * rules. Returns true when the entity was dropped (and noted).
 */
async function untransferable(dapi: NodeDapi, type: string, json: any, note: Note,
                              packageNames: Map<string, string>): Promise<boolean> {
  const drop = (action: 'warn' | 'info', reason: string, detail?: string) => {
    note({name: nqNameOf(json), entityType: type, action, reason, detail});
    return true;
  };
  if (json.package)
    return drop('warn', 'package_entity', await packageName(dapi, json.package.id, packageNames));
  const refuse = untransferableReason(type, json);
  if (refuse)
    return drop(refuse.action, refuse.reason);
  // Only a blob of its own migrates: a file inside a share is a row the target's own
  // share re-creates, and its connection never travels (same rule as the 1.28 executor).
  if (type === 'FileInfo' && json.connection?.id)
    return drop('info', 'file_in_share_not_migratable');
  return false;
}

export async function select(dapi: NodeDapi, sel: Selection, note: Note): Promise<Map<string, BundleEntity>> {
  const picked = new Map<string, BundleEntity>();
  const packageNames = new Map<string, string>();
  const add = async (lite: any) => {
    const type = lite?.['#type'];
    if (!TYPES[type]) {
      note({name: lite?.name ?? lite?.id ?? '', entityType: type ?? '?', action: 'warn', reason: 'type_unsupported'});
      return;
    }
    if (picked.has(lite.id)) return;
    const json = await tryFind(dapi, type, lite.id, nqNameOf(lite), note);
    if (json && !await untransferable(dapi, type, json, note, packageNames))
      picked.set(lite.id, {type, json});
  };

  if (sel.space)
    await add(await resolveEntity(dapi, sel.space));
  for (const n of sel.names)
    await add(await resolveEntity(dapi, n));

  if (!sel.names.length || compileFilter(sel) || sel.tag) {
    for (const type of sel.types) {
      const spec = TYPES[type];
      if (sel.tag && !spec.tags) {
        note({name: '', entityType: type, action: 'warn', reason: 'tag_unsupported'});
        continue;
      }
      // Types without a list route of their own are listed polymorphically by type id.
      const route = spec.listVia === 'entities' ? '/entities' : spec.route;
      const params = {
        text: compileFilter(sel, type) || undefined,
        tags: sel.tag,
        typeId: spec.typeId,
        ...(sel.typeOptions?.[type]?.params ?? {}),
      };
      for (const lite of await dapi.internal(route).listAll(params))
        if (!sel.name || matchesName(sel.name, lite))
          await add(lite);
    }
  }
  return picked;
}

export async function expand(dapi: NodeDapi, selected: Map<string, BundleEntity>, note: Note): Promise<Map<string, BundleEntity>> {
  const out = new Map<string, BundleEntity>();
  const seen = new Set<string>();
  const queue: {type: string; id: string; json?: any}[] = [];
  const packageNames = new Map<string, string>();
  // Datasync scripts name the same query over and over; one lookup per nqName is enough.
  const resolved = new Map<string, any>();

  const enqueue = (type: string, id: string, json?: any) => {
    if (!id || seen.has(id) || !TYPES[type]) return;
    seen.add(id);
    queue.push({type, id, json});
  };

  for (const [id, e] of selected)
    enqueue(e.type, id, e.json);

  while (queue.length) {
    const {type, id, json: known} = queue.shift()!;
    const json = known ?? await tryFind(dapi, type, id, id, note);
    if (!json) continue;

    if (await untransferable(dapi, type, json, note, packageNames)) continue;
    out.set(id, {type, json});

    if (type === 'Project') {
      for (const r of await projectRelations(dapi, id, json, note))
        enqueue(r?.entity?.['#type'], r?.entity?.id);
      for (const route of ['/views', '/layouts'])
        for (const v of await dapi.internal(route).listAll({projectId: id}))
          enqueue(v['#type'], v.id);
    }

    for (const dep of TYPES[type].deps?.(json) ?? []) {
      if (dep.id || !dep.nqName) {
        enqueue(dep.type ?? '', dep.id ?? '');
        continue;
      }
      const key = `${dep.type ?? ''}|${dep.nqName}`;
      if (!resolved.has(key))
        resolved.set(key, await resolveDep(dapi, dep, note));
      const found = resolved.get(key);
      if (found)
        enqueue(found['#type'], found.id);
    }
  }
  await expandGrants(dapi, out, note);
  return out;
}

/**
 * `GET /projects/relations?include=entity` fails on a project that links domain-table rows
 * (`projects_service.dart:425`), so the project's own `relations[]` stands in — typed as
 * well, just fetched together with the project.
 */
async function projectRelations(dapi: NodeDapi, id: string, json: any, note: Note): Promise<any[]> {
  try {
    return await dapi.internal('/projects/relations').listAll({projectId: id, include: 'entity'});
  } catch (err: any) {
    note({name: nqNameOf(json), entityType: 'Project', action: 'warn', reason: 'relations_degraded',
      detail: err?.message ?? String(err)});
    return json.relations ?? [];
  }
}

async function resolveDep(dapi: NodeDapi, dep: Ref, note: Note): Promise<any> {
  try {
    return await resolveEntity(dapi, dep.nqName!, dep.type);
  } catch (err: any) {
    note({name: dep.nqName ?? '', entityType: dep.type ?? 'Entity', action: 'warn', reason: 'dependency_not_found', detail: err?.message});
    return null;
  }
}

/** Entities carry the id of the PUBLISHED package, which only `/packages/published` resolves. */
async function packageName(dapi: NodeDapi, id: string, cache: Map<string, string>): Promise<string> {
  if (!cache.has(id))
    cache.set(id, (await dapi.internal('/packages/published').find(id))?.name ?? '');
  return cache.get(id)!;
}

export const visibleGroup = (g: any): boolean => !!g?.id && g.personal !== true && g.hidden !== true;

/** One `GET /groups/{id}` per group per run, shared by everything that resolves grants. */
export function groupCache(dapi: NodeDapi): (id: string) => Promise<any> {
  const fetched = new Map<string, any>();
  return async (id: string) => {
    if (!fetched.has(id))
      fetched.set(id, await dapi.internal('/groups').find(id));
    return fetched.get(id);
  };
}

/**
 * Grants held on the entity itself by groups a human can belong to. `all=true` is the
 * only listing that hydrates `permission`; it also returns rows inherited from containing
 * entities, which are granted on those entities, not here.
 */
export async function grantsOf(dapi: NodeDapi, id: string, group: (gid: string) => Promise<any> = groupCache(dapi),
): Promise<{groupId: string; group: string; permission: string}[]> {
  const out: {groupId: string; group: string; permission: string}[] = [];
  for (const perm of await dapi.internal('/privileges/permissions').list({entityId: id, all: 'true'})) {
    if (perm.entityId !== id || !perm.permission?.name || !perm.userGroup?.id) continue;
    const g = await group(perm.userGroup.id);
    if (!visibleGroup(g)) continue;
    out.push({groupId: g.id, group: g.friendlyName ?? g.name, permission: perm.permission.name});
  }
  return out;
}

/**
 * Grant-holder groups, their parents and their members. Group and user ids are not
 * portable, so grants and memberships are recorded under bundle-only `_grants` /
 * `_members` keys and replayed by name and login on the target.
 */
async function expandGrants(dapi: NodeDapi, out: Map<string, BundleEntity>, note: Note): Promise<void> {
  const group = groupCache(dapi);
  const queue: any[] = [];
  for (const [id, e] of out) {
    if (e.type === 'UserGroup') { queue.push(e.json); continue; }
    const grants = await grantsOf(dapi, id, group);
    if (!grants.length) continue;
    e.json._grants = grants.map((g) => ({group: g.group, permission: g.permission}));
    for (const g of grants)
      queue.push(await group(g.groupId));
  }

  const done = new Set<string>();
  while (queue.length) {
    const g = queue.shift()!;
    if (done.has(g.id)) continue;
    done.add(g.id);
    // The grant still travels by name; the group itself belongs to the target instance.
    if (isBuiltinGroup(g)) {
      note({name: g.friendlyName ?? g.name, entityType: 'UserGroup', action: 'info', reason: 'platform_group'});
      continue;
    }
    if (!out.has(g.id))
      out.set(g.id, {type: 'UserGroup', json: g});
    const json = out.get(g.id)!.json;
    const members = await groupMembers(dapi, json, group, queue, note);
    if (members.length)
      json._members = members;
    for (const rel of json.parents ?? []) {
      const parent = rel.parent?.id ? await group(rel.parent.id) : null;
      if (visibleGroup(parent) && parent.id !== ALL_USERS)
        queue.push(parent);
    }
  }
}

/** Member groups are bundled too — replaying a membership by name needs the group on the target. */
async function groupMembers(dapi: NodeDapi, g: any, group: (id: string) => Promise<any>, queue: any[], note: Note): Promise<any[]> {
  const members: any[] = [];
  for (const rel of g.children ?? []) {
    const child = rel.child?.id ? await group(rel.child.id) : null;
    if (!child) continue;
    const isAdmin = rel.isAdmin === true;
    if (child.personal !== true) {
      if (!visibleGroup(child)) continue;
      members.push({kind: 'group', name: child.friendlyName ?? child.name, isAdmin});
      queue.push(child);
      continue;
    }
    const user = await dapi.internal(`/groups/${child.id}`).find('user');
    if (user?.login)
      members.push({kind: 'user', login: user.login, email: user.email, isAdmin});
    else
      note({name: g.friendlyName ?? g.name, entityType: 'UserGroup', action: 'warn',
        reason: 'member_unresolved', detail: child.friendlyName ?? child.id});
  }
  return members;
}

export async function pullBytes(dapi: NodeDapi, entities: Map<string, BundleEntity>, note: Note,
                                kinds: BytesKind[]): Promise<Map<string, Buffer>> {
  const bytes = new Map<string, Buffer>();
  for (const [id, {type, json}] of entities) {
    const spec = TYPES[type];
    if (!spec.bytes || !kinds.includes(spec.bytes.kind)) continue;
    try {
      bytes.set(id, await dapi.client.getBytes(spec.bytes.get(id)));
    } catch (err: any) {
      note({name: json.name ?? id, entityType: type, action: 'warn', reason: 'no_data', detail: err?.message});
    }
  }
  return bytes;
}
