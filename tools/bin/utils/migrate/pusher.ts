/// Docs: [Entity export / import](/docs/features/grok-tool/export-import/DESIGN.md)
import * as fs from 'fs';
import {randomUUID} from 'crypto';
import {NodeDapi} from '../node-dapi';
import {TYPES, inNamespace, nqNameOf, rankOf, untransferableReason} from './registry';
import {Bundle, BundleEntity, bytesPath, hashOf, hashView, stripPrivate, writeIdmap} from './bundle';
import {nestedIds, rewrite} from './rewriter';
import {findEntity, grantsOf, groupCache} from './walker';

export type Action = 'create' | 'update' | 'identical' | 'skip' | 'failed' | 'warn' | 'info' | 'needs-credentials';
export type ConflictPolicy = 'fail' | 'skip' | 'duplicate' | 'adopt';

export interface Row {name: string; entityType: string; action: Action; reason: string; detail?: string}

export interface Op {id: string; type: string; json: any; row: Row; creds?: Record<string, any>; expectedNamespace?: string}

export interface PushOptions {
  dryRun?: boolean;
  onConflict?: ConflictPolicy;
  concurrency?: number;
  creds?: Record<string, any>;
}

export interface PushResult {items: Row[]; counts: Record<string, number>; status: string; remoteUrl: string}

export async function plan(dapi: NodeDapi, bundle: Bundle,
                           opts: {onConflict: ConflictPolicy; creds?: Record<string, any>; idmap?: Record<string, string>},
): Promise<{rows: Row[]; ops: Op[]; planned: Map<string, Row>; effective: Map<string, BundleEntity>}> {
  const rows: Row[] = [];
  const ops: Op[] = [];
  const planned = new Map<string, Row>();
  const effective = new Map<string, BundleEntity>();
  const conflicts: string[] = [];
  const idmap = opts.idmap ?? {};
  const orphans = new Set<string>();
  const referrers = new Map<string, string[]>();
  const onTarget = new Set<string>();
  const pusherNamespace = await currentNamespace(dapi);

  // Every twin is resolved before the first payload is rewritten: an adoption discovered
  // halfway through would leave the references of everything rewritten before it stale.
  const resolved = new Map<string, {target: any; twin: any}>();
  for (const entry of bundle.manifest.order) {
    const {type, json: source} = bundle.entities.get(entry.id)!;
    const target = await dapi.internal(TYPES[type].route).find(idmap[entry.id] ?? source.id);
    const twin = target ? null : await findByNqName(dapi, bundle, type, source, pusherNamespace);
    resolved.set(entry.id, {target, twin});
    if (!twin || opts.onConflict !== 'adopt') continue;
    // The twin's own nested rows must not be hijacked by the bundle's row ids; the fresh
    // ones go into the idmap so the next push produces the same payload.
    idmap[entry.id] = twin.id;
    for (const nested of nestedIds(source))
      idmap[nested] = randomUUID();
  }

  for (const entry of bundle.manifest.order) {
    const {type, json: source} = bundle.entities.get(entry.id)!;
    const mine = new Set<string>();
    const json = rewrite(source, idmap, mine);
    for (const orphan of mine) {
      orphans.add(orphan);
      referrers.set(orphan, [...(referrers.get(orphan) ?? []), nqNameOf(json)]);
    }
    const row: Row = {name: nqNameOf(json), entityType: type, action: 'create', reason: ''};
    rows.push(row);
    planned.set(entry.id, row);

    const {target, twin} = resolved.get(entry.id)!;
    // A hand-edited bundle must not be able to overwrite what the target owns itself.
    const refuse = untransferableReason(type, json);
    if (refuse) {
      row.action = 'skip';
      row.reason = refuse.reason;
      if (target) onTarget.add(json.id);
      effective.set(entry.id, {type, json});
      continue;
    }
    if (target) {
      onTarget.add(json.id);
      compare(rows, row, type, target, json, '');
    } else if (twin) {
      row.reason = `name taken by ${twin.id}`;
      if (opts.onConflict === 'fail') {
        conflicts.push(`${type} ${row.name} → ${twin.id}`);
        row.action = 'failed';
      } else if (opts.onConflict === 'skip')
        row.action = 'skip';
      else if (opts.onConflict === 'adopt')
        compare(rows, row, type, twin, json, `adopted ${twin.id}`);
    }
    effective.set(entry.id, {type, json});

    const creds = type === 'DataConnection' ? credentialsFor(opts.creds, json, pusherNamespace) : undefined;
    // A new secret is invisible in the payload, so a covered connection is written anyway.
    if (creds && row.action === 'identical') {
      row.action = 'update';
      row.reason = 'credentials';
    }
    if (!['create', 'update'].includes(row.action)) continue;
    if (type === 'DataConnection' && !planCredentials(json, row, rows, creds)) continue;
    ops.push({id: entry.id, type, json, row, creds, expectedNamespace: expectedNamespace(bundle, json)});
  }

  if (conflicts.length)
    throw new Error(`Name conflicts on the target (use --on-conflict skip|duplicate|adopt):\n  ${conflicts.join('\n  ')}`);

  failDependants(effective, planned, onTarget);

  // Nested row ids belong to the bundle as much as entity ids do — only a reference to
  // something that was never pulled is an orphan.
  const known = new Set<string>(Object.keys(idmap));
  for (const [id, e] of bundle.entities) {
    known.add(id);
    for (const nested of nestedIds(e.json))
      known.add(nested);
  }
  // An id the bundle points at but does not carry is only harmless if the target happens to hold
  // it too (a platform row keeps its id everywhere). One that exists on neither side is dangling
  // on the source, and whatever needs it will be refused there with the server's own error.
  for (const orphan of orphans) {
    if (known.has(orphan)) continue;
    const by = [...new Set(referrers.get(orphan) ?? [])].join(', ');
    rows.push(await findEntity(dapi, orphan)
      ? {name: orphan, entityType: 'Entity', action: 'info', reason: 'orphan_ref',
        detail: `not in the bundle, but the target has this id — referenced by ${by}`}
      : {name: orphan, entityType: 'Entity', action: 'warn', reason: 'dependency_missing',
        detail: `on neither the source nor the target; the save of ${by} may be refused`});
  }
  return {rows, ops: ops.filter((o) => ['create', 'update'].includes(o.row.action)), planned, effective};
}

/** Relations are outside the hash: what matters is that the target links everything the bundle does. */
function relationsCovered(target: any, json: any): boolean {
  const have = new Set((target.relations ?? []).map((r: any) => r?.entity?.id));
  return (json.relations ?? []).every((r: any) => have.has(r.entity?.id));
}

function compare(rows: Row[], row: Row, type: string, target: any, json: any, reason: string): void {
  const covered = relationsCovered(target, json);
  const sameHash = hashOf(type, target) === hashOf(type, json);
  const same = covered && sameHash;
  row.action = same ? 'identical' : 'update';
  row.reason = same ? reason : reason || (sameHash ? 'relations missing on target' : 'hash differs');
  row.detail = same ? undefined : [...changedKeys(type, target, json), ...(covered ? [] : ['relations'])].join(', ');

  // Coverage semantics: relations are only ever added, so a link dropped from the bundle stays.
  const wanted = new Set((json.relations ?? []).map((r: any) => r?.entity?.id));
  const kept = (target.relations ?? [])
    .filter((r: any) => r?.entity?.id && !wanted.has(r.entity.id) && r.entity.parameters?.isProject !== true)
    .map((r: any) => r.entity.id);
  if (kept.length)
    rows.push({name: row.name, entityType: type, action: 'info', reason: 'relation_not_removed', detail: kept.join(', ')});
}

/** What has to be on the target before this entity can land: bundle-internal references. */
function referencedIds(json: any, groupIds: Map<string, string>, type: string): string[] {
  const ids: string[] = (TYPES[type].deps?.(json) ?? []).map((r) => r.id).filter(Boolean) as string[];
  for (const rel of json.relations ?? [])
    if (rel?.entity?.id)
      ids.push(rel.entity.id);
  for (const grant of json._grants ?? [])
    if (groupIds.has(grant.group))
      ids.push(groupIds.get(grant.group)!);
  return ids;
}

/**
 * A skipped entity is not on the target, so anything pointing at it would be rejected
 * there — report that here instead of letting the server fail the save.
 */
function failDependants(effective: Map<string, BundleEntity>, planned: Map<string, Row>, onTarget: Set<string>): void {
  const blocked = new Map<string, Row>();
  const groupIds = new Map<string, string>();
  for (const [id, {type, json}] of effective) {
    if (type === 'UserGroup')
      groupIds.set(json.friendlyName ?? json.name, json.id);
    // A connection skipped because every parameter was masked is still on the target. A platform
    // group is too, under whatever id that instance gave it — grants name it, so it never blocks.
    if (planned.get(id)?.action === 'skip' && !onTarget.has(json.id) && planned.get(id)!.reason !== 'platform_group')
      blocked.set(json.id, planned.get(id)!);
  }
  for (let changed = true; changed;) {
    changed = false;
    for (const [id, {type, json}] of effective) {
      const row = planned.get(id)!;
      if (!['create', 'update'].includes(row.action)) continue;
      const on = referencedIds(json, groupIds, type).find((ref) => blocked.has(ref));
      if (!on) continue;
      row.action = 'failed';
      row.reason = 'dependency_skipped';
      row.detail = blocked.get(on)!.name;
      blocked.set(json.id, row);
      changed = true;
    }
  }
}

/**
 * Where the entity should end up, or undefined where the target legitimately decides: a personal
 * namespace becomes the pusher's, and a root space has none to keep.
 */
function expectedNamespace(bundle: Bundle, json: any): string | undefined {
  const source: string = json.namespace ?? '';
  return !source || source === bundle.manifest.source.userNamespace ? undefined : source;
}

/**
 * The `--creds` file is authored for the target, so its keys are the connection nqName as
 * the bundle spells it, or as the target does once a personal namespace is remapped (R11).
 */
function credentialsFor(creds: Record<string, any> | undefined, json: any, pusherNamespace: string): Record<string, any> | undefined {
  return creds?.[nqNameOf(json)] ?? creds?.[`${pusherNamespace}${json.name ?? ''}`];
}

/**
 * Secrets never travel in a bundle: `_credentials` lists the parameters the source masked.
 * When they were the only parameters there is nothing left worth pushing.
 */
function planCredentials(json: any, row: Row, rows: Row[], creds?: Record<string, any>): boolean {
  const masked: string[] = json._credentials ?? [];
  if (masked.length && !creds && !Object.keys(json.parameters ?? {}).length) {
    row.action = 'skip';
    row.reason = 'every parameter was masked on the source';
    return false;
  }
  if (!creds)
    rows.push({name: row.name, entityType: 'DataConnection', action: 'needs-credentials', reason: '',
      detail: masked.length ? masked.join(', ') : 'credentials are not migrated'});
  return true;
}

export async function push(dapi: NodeDapi, bundle: Bundle, opts: PushOptions, log: (rows: Row[]) => void): Promise<PushResult> {
  const onConflict = opts.onConflict ?? 'fail';
  const extra: Row[] = [];

  const target = await dapi.serverInfo();
  const source = bundle.manifest.source;
  if (minor(source.version) !== minor(target.version))
    extra.push({
      name: dapi.client.baseUrl, entityType: 'Server', action: 'warn', reason: 'version_mismatch',
      detail: `${source.version} → ${target.version}; unknown fields are ignored by the server, failures below may be version-related`,
    });

  const installed = new Set((await dapi.packages.listFull()).map((p: any) => String(p.name).toLowerCase()));
  for (const name of bundle.manifest.packages)
    if (!installed.has(name.toLowerCase()))
      extra.push({name, entityType: 'Package', action: 'warn', reason: 'package_not_installed'});

  const idmap = {...bundle.idmap};
  const {rows, ops, planned, effective} = await plan(dapi, bundle, {onConflict, creds: opts.creds, idmap});
  const items = [...extra, ...rows];
  log(items);
  if (opts.dryRun)
    return summarize(items, dapi, 'dry-run');

  for (const rank of [...new Set(ops.map((o) => rankOf(o.type)))].sort((a, b) => a - b))
    await pool(ops.filter((o) => rankOf(o.type) === rank), opts.concurrency ?? 6,
      (op) => saveOne(dapi, bundle, op, items, idmap));
  // A FileInfo save can answer with an existing row's id: recording it keeps the next push
  // stable, and the passes below have to point at the id the entity actually landed under.
  if (Object.keys(idmap).length > Object.keys(bundle.idmap).length) {
    writeIdmap(bundle.dir, idmap);
    for (const [id, e] of effective)
      effective.set(id, {type: e.type, json: rewrite(e.json, idmap)});
  }

  await pushRelations(dapi, effective, planned, items);
  await pushTags(dapi, effective, planned, items);
  await pushMemberships(dapi, effective, planned, items);
  const {visible, projects} = await pushGrants(dapi, effective, planned, items);
  const written = items.filter((r) => r.action === 'create' || r.action === 'update').length;
  if (!visible && projects && written)
    items.push({name: '', entityType: 'Project', action: 'warn', reason: 'not_visible',
      detail: 'Pushed but not visible to any user on the remote — share it with a group first'});

  return summarize(items, dapi, items.some((r) => r.action === 'failed') ? 'failed' : 'ok');
}

async function saveOne(dapi: NodeDapi, bundle: Bundle, op: Op, rows: Row[], idmap: Record<string, string>): Promise<void> {
  const spec = TYPES[op.type];
  const payload = stripPrivate(JSON.parse(JSON.stringify(op.json)));
  const targetId = payload.id ?? op.id;
  if (payload.metaParams && typeof payload.metaParams === 'object')
    payload.metaParams.sync_id = op.id;

  // The server encrypts and masks password-class parameters itself, so target-side
  // secrets go in as plain parameters and never touch the bundle.
  if (op.creds)
    payload.parameters = {...payload.parameters, ...op.creds};

  try {
    if (spec.bytes && !spec.bytes.afterSave)
      await pushBytes(dapi, bundle, op.type, op.id, targetId);

    const saved = await dapi.internal(spec.saveRoute ?? spec.route).save(payload);
    const peerId = saved?.id ?? targetId;
    if (peerId !== targetId)
      idmap[op.id] = peerId;
    if (spec.bytes && spec.bytes.afterSave)
      await pushBytes(dapi, bundle, op.type, op.id, peerId);

    const verified = await dapi.internal(spec.route).find(peerId);
    if (!verified) {
      op.row.action = 'failed';
      op.row.reason = 'Save reported success but not on target';
      return;
    }
    if (verified.name !== op.json.name)
      rows.push({name: op.row.name, entityType: op.type, action: 'warn', reason: 'renamed', detail: `${op.json.name} → ${verified.name}`});
    // A namespace is a label; placement comes from the owning space's relations, so an entity
    // whose space is absent (or that is not among its relations) lands under the pusher instead.
    if (op.expectedNamespace !== undefined && (verified.namespace ?? '') !== op.expectedNamespace)
      rows.push({name: op.row.name, entityType: op.type, action: 'warn', reason: 'namespace_not_preserved',
        detail: `${op.expectedNamespace} → ${verified.namespace ?? ''}; pull the owning space so the entity travels in its relations`});
    if (op.type === 'PredictiveModelInfo')
      rows.push({name: op.row.name, entityType: op.type, action: 'info', reason: 'model_blob_skipped',
        detail: 'the trained model itself stays on the source — retrain or copy it separately'});
  } catch (err: any) {
    op.row.action = 'failed';
    op.row.reason = err?.message ?? String(err);
  }
}

/** Bytes live under the bundle id and land under the id the target ended up using. */
async function pushBytes(dapi: NodeDapi, bundle: Bundle, type: string, bundleId: string, targetId: string): Promise<void> {
  const bytes = TYPES[type].bytes!;
  const file = bytesPath(bundle.dir, bytes.kind, bundleId);
  if (fs.existsSync(file))
    await dapi.client.putBytes(bytes.put(targetId), fs.readFileSync(file));
}

/**
 * Relations need every member to exist, so they are re-attached once everything is saved.
 * `saveRelations=true` replaces the project's whole relation set (`projects_repository.dart`
 * `_deleteRelations`), so the payload starts from what the target already links and the
 * bundle only ever adds to it.
 */
async function pushRelations(dapi: NodeDapi, effective: Map<string, BundleEntity>, planned: Map<string, Row>, rows: Row[]): Promise<void> {
  const projects = dapi.internal('/projects');
  const inBundle = new Set([...effective.values()].map((e) => e.json.id));
  for (const [id, {type, json}] of effective) {
    if (type !== 'Project' || !json.relations?.length) continue;
    if (!['create', 'update'].includes(planned.get(id)?.action ?? 'create')) continue;

    const target = await projects.find(json.id);
    if (!target) continue;
    const wanted: any[] = (target.relations ?? [])
      .filter((r: any) => r?.entity?.id)
      .map((r: any) => ({id: r.id, entity: {'#type': 'EntityRecord', id: r.entity.id}, isLink: r.isLink ?? false}));
    const linked = new Set<string>(wanted.map((r) => r.entity.id));
    let added = 0;
    for (const rel of json.relations) {
      if (linked.has(rel.entity.id)) continue;
      // The target stamps its own Files connection onto a space; re-attaching the source's
      // would leave the space with two of them.
      if (!inBundle.has(rel.entity.id)) {
        const existing = await findEntity(dapi, rel.entity.id);
        if (!existing || existing.parameters?.isProject === true) continue;
      }
      wanted.push(rel);
      linked.add(rel.entity.id);
      added++;
    }
    if (!added) continue;
    delete target.storage;
    target.relations = wanted;

    for (let attempt = 0; ; attempt++) {
      try {
        await projects.save(target, 'saveRelations=true');
        break;
      } catch (err: any) {
        const text = String(err?.message ?? err);
        if (attempt >= 2 || !(text.includes('40P01') || text.toLowerCase().includes('deadlock'))) {
          rows.push({name: nqNameOf(json), entityType: 'Project', action: 'failed', reason: 'relations', detail: text});
          break;
        }
        await new Promise((r) => setTimeout(r, 250 * (attempt + 1)));
      }
    }
  }
}

/**
 * Tag rows are server-managed and duplicate on every POST, so only the tags the
 * target is missing are replayed.
 */
async function pushTags(dapi: NodeDapi, effective: Map<string, BundleEntity>, planned: Map<string, Row>, rows: Row[]): Promise<void> {
  const missing = new Map<string, string[]>();
  try {
    for (const [id, {type, json}] of effective) {
      const tags: string[] = json._tags ?? [];
      if (!tags.length || ['skip', 'failed'].includes(planned.get(id)?.action ?? '')) continue;
      const target = await dapi.internal(TYPES[type].route).find(json.id);
      if (!target) continue;
      const have = new Set((target.entityTags ?? []).map((t: any) => t?.tag));
      for (const tag of tags)
        if (!have.has(tag))
          missing.set(tag, [...(missing.get(tag) ?? []), json.id]);
    }
    for (const [tag, ids] of missing)
      await dapi.client.post(`/entities/tag?tag=${encodeURIComponent(tag)}`, ids);
  } catch (err: any) {
    rows.push({name: [...missing.keys()].join(', '), entityType: 'Entity', action: 'warn', reason: 'tags',
      detail: err?.message ?? String(err)});
  }
}

/** Groups travel bare, so their members are re-attached by login and by group name. */
async function pushMemberships(dapi: NodeDapi, effective: Map<string, BundleEntity>, planned: Map<string, Row>, rows: Row[]): Promise<void> {
  for (const [id, {type, json}] of effective) {
    const members: any[] = json._members ?? [];
    if (type !== 'UserGroup' || !members.length) continue;
    if (['skip', 'failed'].includes(planned.get(id)?.action ?? '')) continue;

    let matched = 0;
    let missing = 0;
    try {
      for (const kind of ['user', 'group'])
        for (const isAdmin of [false, true]) {
          const wanted = members.filter((m) => m.kind === kind && (m.isAdmin === true) === isAdmin);
          if (!wanted.length) continue;
          const names = wanted.map((m) => (kind === 'user' ? m.login : m.name));
          for (const result of await dapi.groups.addMembers(json.id, names, isAdmin, kind === 'user')) {
            if (result.status !== 'error') { matched++; continue; }
            missing++;
            rows.push({name: nqNameOf(json), entityType: 'UserGroup', action: 'warn', reason: 'member_not_found',
              detail: `${result.member}: ${result.error}`});
          }
        }
    } catch (err: any) {
      rows.push({name: nqNameOf(json), entityType: 'UserGroup', action: 'warn', reason: 'members',
        detail: err?.message ?? String(err)});
      continue;
    }
    const row = planned.get(id);
    if (row)
      row.detail = [row.detail, `members: ${matched} matched${missing ? `, ${missing} not on remote` : ''}`].filter(Boolean).join('; ');
  }
}

/**
 * Without a grant to a real group the pushed content is invisible on the target, so the
 * source grants are replayed by group name. Returns whether anything is visible at all.
 */
async function pushGrants(dapi: NodeDapi, effective: Map<string, BundleEntity>, planned: Map<string, Row>, rows: Row[],
): Promise<{visible: boolean; projects: number}> {
  let anyVisible = false;
  let projects = 0;
  const group = groupCache(dapi);
  for (const [id, {type, json}] of effective) {
    if (type !== 'Project' || ['skip', 'failed'].includes(planned.get(id)?.action ?? '')) continue;
    projects++;
    try {
      if (!await dapi.internal('/projects').find(json.id)) continue;
      const have = new Set<string>();
      for (const perm of await grantsOf(dapi, json.id, group)) {
        have.add(`${perm.group.toLowerCase()}|${perm.permission}`);
        anyVisible = true;
      }
      for (const grant of json._grants ?? []) {
        if (have.has(`${grant.group.toLowerCase()}|${grant.permission}`)) continue;
        if (!['View', 'Edit'].includes(grant.permission)) {
          rows.push({name: nqNameOf(json), entityType: type, action: 'info', reason: 'unsupported_grant',
            detail: `${grant.group}: ${grant.permission}`});
          continue;
        }
        const peer = (await dapi.groups.lookup(grant.group))
          .find((g: any) => g?.personal !== true && (g.friendlyName ?? g.name) === grant.group);
        if (!peer) {
          rows.push({name: nqNameOf(json), entityType: type, action: 'warn', reason: 'group_not_found', detail: grant.group});
          continue;
        }
        await dapi.shares.share(json.id, grant.group, grant.permission);
        anyVisible = true;
      }
    } catch (err: any) {
      rows.push({name: nqNameOf(json), entityType: type, action: 'warn', reason: 'grants',
        detail: err?.message ?? String(err)});
    }
  }
  return {visible: anyVisible, projects};
}

export async function pool<T>(items: T[], concurrency: number, work: (item: T) => Promise<void>): Promise<void> {
  const queue = items.slice();
  const workers: Promise<void>[] = [];
  for (let i = 0; i < Math.min(concurrency, queue.length); i++)
    workers.push((async () => {
      while (queue.length)
        await work(queue.shift()!);
    })());
  await Promise.all(workers);
}

export function summarize(items: Row[], dapi: NodeDapi, status: string): PushResult {
  const counts: Record<string, number> = {};
  for (const r of items)
    counts[r.action] = (counts[r.action] ?? 0) + 1;
  return {items, counts, status, remoteUrl: dapi.client.baseUrl};
}

/** Top-level keys whose normalized values differ — the same view of the JSON `hashOf` takes. */
export function changedKeys(type: string, a: any, b: any): string[] {
  const na = hashView(type, a), nb = hashView(type, b);
  return [...new Set([...Object.keys(na), ...Object.keys(nb)])]
    .filter((k) => JSON.stringify(na[k]) !== JSON.stringify(nb[k]));
}

function minor(version: string): string {
  return String(version ?? '').split('.').slice(0, 2).join('.');
}

async function currentNamespace(dapi: NodeDapi): Promise<string> {
  const user = await dapi.client.get('/users/current');
  return user?.project?.name ? `${user.project.name}:` : '';
}

/**
 * An entity pulled from the source author's personal namespace lands under the pusher's
 * own namespace on the target, so that is where a same-name twin would be.
 */
async function findByNqName(dapi: NodeDapi, bundle: Bundle, type: string, json: any, pusherNamespace: string): Promise<any> {
  // Not `expectedNamespace`: a twin for a namespace-less entity is looked up in the root.
  const sourceNamespace: string = json.namespace ?? '';
  const namespace = sourceNamespace === bundle.manifest.source.userNamespace ? pusherNamespace : sourceNamespace;
  const matches = await dapi.internal('/entities').list({namespace, name: json.name});
  return matches.find((m: any) => m['#type'] === type && inNamespace(m, namespace)) ?? null;
}
