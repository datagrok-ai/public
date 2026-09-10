/// Docs: [Entity export / import](/docs/features/grok-tool/export-import/DESIGN.md)
import * as fs from 'fs';
import {randomUUID} from 'crypto';
import {NodeDapi} from '../node-dapi';
import {TYPES, inNamespace, nqNameOf, rankOf, untransferableReason} from './registry';
import {Bundle, BundleEntity, bytesPath, hashOf, hashView, listShares, sharePath, stripPrivate, writeIdmap} from './bundle';
import {pool} from './pool';
import {nestedIds, rewrite} from './rewriter';
import {Progress, findEntity, grantsOf, groupCache} from './walker';

/** A relation set is written whole; a project holding thousands of them is slow, not stuck. */
const BULK_WRITE_MS = Number(process.env['GROK_HTTP_BULK_TIMEOUT'] ?? 600000);
const RELATION_REFUSED_RE = /entity\s+([0-9a-f-]{36})/i;
const WRITE_CONCURRENCY = 6;

export type Action = 'create' | 'update' | 'identical' | 'skip' | 'failed' | 'warn' | 'info' | 'needs-credentials';
export type ConflictPolicy = 'fail' | 'skip' | 'duplicate' | 'adopt';

export interface Row {name: string; entityType: string; action: Action; reason: string; detail?: string}

export interface Op {id: string; type: string; json: any; row: Row; creds?: Record<string, any>;
  expectedNamespace?: string; dangling?: string[]}

/** What a save could not judge on its own, because placement had not happened yet. */
interface Deferred {renamed: {op: Op; row: Row}[]; misplaced: Op[]}

export interface PushOptions {
  dryRun?: boolean;
  onConflict?: ConflictPolicy;
  concurrency?: number;
  creds?: Record<string, any>;
  progress?: Progress;
}

export interface PushResult {items: Row[]; counts: Record<string, number>; status: string; remoteUrl: string}

export async function plan(dapi: NodeDapi, bundle: Bundle,
                           opts: {onConflict: ConflictPolicy; creds?: Record<string, any>; idmap?: Record<string, string>},
): Promise<{rows: Row[]; ops: Op[]; planned: Map<string, Row>; effective: Map<string, BundleEntity>; missing: Set<string>}> {
  const rows: Row[] = [];
  const ops: Op[] = [];
  const planned = new Map<string, Row>();
  const effective = new Map<string, BundleEntity>();
  const conflicts: string[] = [];
  const idmap = opts.idmap ?? {};
  const orphans = new Set<string>();
  const missing = new Set<string>();
  const referrers = new Map<string, string[]>();
  const onTarget = new Set<string>();
  const pusherNamespace = await currentNamespace(dapi);

  // A personal root project exists on the target already, under the target's own id, because it
  // was created with the user. Pointing at that one before anything is rewritten is what keeps
  // the content under it in its owner's namespace instead of the pusher's.
  const owners = bundle.manifest.order
    .map((entry) => bundle.entities.get(entry.id)!.json._personalOf).filter(Boolean);
  if (owners.length) {
    const personal = await personalProjects(dapi);
    for (const entry of bundle.manifest.order) {
      const owner = bundle.entities.get(entry.id)!.json._personalOf;
      if (owner && personal.has(owner))
        idmap[entry.id] = personal.get(owner)!;
    }
    const absent = [...new Set(owners)].filter((login) => !personal.has(login)).sort();
    if (absent.length)
      rows.push({name: absent.join(', '), entityType: 'User', action: 'warn', reason: 'user_missing',
        detail: `create ${absent.length === 1 ? 'this user' : 'these users'} on the target first, ` +
          'or their content lands under the pushing account'});
  }

  // A package entity gets a fresh id on every instance, so it only resolves by name.
  await pool(bundle.manifest.externals ?? [], WRITE_CONCURRENCY, async (external) => {
    // A mapping the bundle carries was learned from whichever target it was last pushed at, so it
    // is only reusable while this target still answers to it.
    const mapped = idmap[external.id];
    const spec = TYPES[external.type];
    if (mapped && spec && await dapi.internal(spec.route).find(mapped).catch(() => null)) return;
    const twin = await resolveByNqName(dapi, external.type, external.nqName);
    if (!twin) {
      delete idmap[external.id];
      rows.push({name: external.nqName, entityType: external.type, action: 'warn', reason: 'external_missing',
        detail: 'stayed on the source and the target has nothing by that name — install or publish what owns it'});
      return;
    }
    if (twin.id !== external.id)
      idmap[external.id] = twin.id;
    else
      delete idmap[external.id];
  });

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
    if (json._personalOf) {
      row.action = 'skip';
      row.reason = 'personal_project';
      row.detail = `the target keeps ${json._personalOf}'s own`;
      onTarget.add(json.id);
      effective.set(entry.id, {type, json});
      continue;
    }
    // A hand-edited bundle must not be able to overwrite what the target owns itself.
    const refuse = untransferableReason(type, json) ?? bytesMissing(bundle, type, json);
    if (refuse) {
      row.action = 'skip';
      row.reason = refuse.reason;
      row.detail = refuse.detail;
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
    ops.push({id: entry.id, type, json, row, creds, expectedNamespace: expectedNamespace(bundle, json),
      dangling: [...mine]});
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
    if (await findEntity(dapi, orphan)) {
      rows.push({name: orphan, entityType: 'Entity', action: 'info', reason: 'orphan_ref',
        detail: `not in the bundle, but the target has this id — referenced by ${by}`});
      continue;
    }
    missing.add(orphan);
    rows.push(bundle.manifest.dangling?.includes(orphan)
      ? {name: orphan, entityType: 'Entity', action: 'info', reason: 'dead_on_source',
        detail: `nothing on the source answers to it either — ${by} cannot be migrated intact`}
      : {name: orphan, entityType: 'Entity', action: 'warn', reason: 'dependency_missing',
        detail: `not in the bundle and not on the target — install the package that owns it, ` +
          `or pull it; the save of ${by} may be refused`});
  }
  return {rows, ops: ops.filter((o) => ['create', 'update'].includes(o.row.action)), planned, effective, missing};
}

/**
 * A file with no share is stored as a blob under its own id (`files_service.dart`:
 * `addToUserProject: f.connection == null`), so without its bytes there is nothing to create.
 */
function bytesMissing(bundle: Bundle, type: string, json: any): {reason: string; detail: string} | null {
  const bytes = TYPES[type].bytes;
  if (type !== 'FileInfo' || json?.connection?.id || !bytes) return null;
  return fs.existsSync(bytesPath(bundle.dir, bytes.kind, json.id)) ? null
    : {reason: 'file_bytes_missing', detail: 'a file with no share is its bytes — pull with --include-files'};
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
  const {rows, ops, planned, effective, missing} = await plan(dapi, bundle, {onConflict, creds: opts.creds, idmap});
  const items = [...extra, ...rows];
  log(items);
  if (opts.dryRun)
    return summarize(items, dapi, 'dry-run');

  const deferred: Deferred = {renamed: [], misplaced: []};
  const progress = opts.progress ?? (() => {});
  let saved = 0;
  for (const rank of [...new Set(ops.map((o) => rankOf(o.type)))].sort((a, b) => a - b))
    await pool(ops.filter((o) => rankOf(o.type) === rank), opts.concurrency ?? WRITE_CONCURRENCY,
      async (op) => {
        await saveOne(dapi, bundle, op, items, idmap, deferred, missing);
        progress('saving', ++saved, ops.length);
      });
  // A FileInfo save can answer with an existing row's id: recording it keeps the next push
  // stable, and the passes below have to point at the id the entity actually landed under.
  // Comparing sizes would miss a run that dropped one stale mapping and learned another.
  if (JSON.stringify(idmap) !== JSON.stringify(bundle.idmap)) {
    writeIdmap(bundle.dir, idmap);
    for (const [id, e] of effective)
      effective.set(id, {type: e.type, json: rewrite(e.json, idmap)});
  }

  await pushShares(dapi, bundle, items, onConflict);
  await pushRelations(dapi, effective, planned, items, progress);
  progress('restoring names');
  const renamed = await restoreNames(dapi, deferred.renamed, idmap);
  // A save re-homes the entity into the pusher's own root (`repository_query.dart`:
  // `addToUserProject`), undoing the placement the namespace is derived from.
  if (renamed)
    await pushRelations(dapi, effective, planned, items, progress);
  progress('checking placement');
  await reportPlacement(dapi, deferred.misplaced, items, idmap);
  progress('tags and memberships');
  await pushTags(dapi, effective, planned, items);
  await pushMemberships(dapi, effective, planned, items);
  progress('sharing');
  const {visible, projects} = await pushGrants(dapi, effective, planned, items);
  const written = items.filter((r) => r.action === 'create' || r.action === 'update').length;
  if (!visible && projects && written)
    items.push({name: '', entityType: 'Project', action: 'warn', reason: 'not_visible',
      detail: 'Pushed but not visible to any user on the remote — share it with a group first'});

  return summarize(items, dapi, items.some((r) => r.action === 'failed') ? 'failed' : 'ok');
}

async function saveOne(dapi: NodeDapi, bundle: Bundle, op: Op, rows: Row[], idmap: Record<string, string>,
                       deferred: Deferred, missing: Set<string>): Promise<void> {
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
    // Names and namespaces are settled by placement, so a clash here is only a candidate.
    if (verified.name !== op.json.name) {
      const row: Row = {name: op.row.name, entityType: op.type, action: 'warn', reason: 'renamed',
        detail: `${op.json.name} → ${verified.name}`};
      rows.push(row);
      deferred.renamed.push({op, row});
    }
    if (op.expectedNamespace !== undefined && (verified.namespace ?? '') !== op.expectedNamespace)
      deferred.misplaced.push(op);
    if (op.type === 'PredictiveModelInfo')
      rows.push({name: op.row.name, entityType: op.type, action: 'info', reason: 'model_blob_skipped',
        detail: 'the trained model itself stays on the source — retrain or copy it separately'});
  } catch (err: any) {
    // The server reports its own refusal, not what caused it; a reference that exists on
    // neither instance is the cause often enough to be worth naming here.
    const cause = (op.dangling ?? []).filter((id) => missing.has(id));
    const dead = cause.filter((id) => bundle.manifest.dangling?.includes(id));
    const refusal = err?.message ?? String(err);
    // The server names its own refusal, not what caused it, so a dead reference is inferred — but
    // only from a refusal. A request that timed out or lost its connection says nothing about the
    // entity, and calling that "nothing to migrate" would hide a transport problem as clean data.
    const transport = /no answer in|deadlock|40P01|ECONN|socket|fetch failed/i.test(refusal);
    const named = transport ? [] : dead;
    op.row.action = named.length ? 'skip' : 'failed';
    op.row.reason = named.length ? 'dead_on_source' : refusal;
    if (named.length)
      op.row.detail = `references ${named.join(', ')}, which the source no longer has either — ` +
        `nothing to migrate, delete it there or ignore (${refusal})`;
    else if (cause.length)
      op.row.detail = `references ${cause.join(', ')}, which the target does not have — ` +
        'a package connection is the usual cause, install the package first';
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
async function pushRelations(dapi: NodeDapi, effective: Map<string, BundleEntity>, planned: Map<string, Row>, rows: Row[],
                             progress: Progress = () => {}): Promise<void> {
  const projects = dapi.internal('/projects');
  const relations = dapi.internal('/projects/relations');
  // Anything this push did not write has to prove it exists before being offered as a relation.
  const landed = new Set<string>();
  for (const [id, {json}] of effective)
    if (['create', 'update', 'identical'].includes(planned.get(id)?.action ?? '')) landed.add(json.id);
  // Containment is a tree: each entity goes to the deepest project claiming it, so a space and a
  // dashboard listing the same table do not erase each other.
  const claimants = new Map<string, string[]>();
  for (const [, {type, json}] of effective) {
    if (type !== 'Project') continue;
    for (const rel of json.relations ?? [])
      if (rel?.entity?.id) claimants.set(rel.entity.id, [...(claimants.get(rel.entity.id) ?? []), json.id]);
  }
  const depthOf = (id: string): number => {
    let depth = 0;
    for (let at = (claimants.get(id) ?? [])[0]; at && depth < 64; at = (claimants.get(at) ?? [])[0]) depth++;
    return depth;
  };
  const ownerOf = (entityId: string): string | undefined => {
    const holders = claimants.get(entityId) ?? [];
    return holders.length < 2 ? holders[0]
      : holders.reduce((deepest, id) => depthOf(id) > depthOf(deepest) ? id : deepest, holders[0]);
  };
  const order = [...effective].filter(([, e]) => e.type === 'Project' && e.json.relations?.length)
    .sort((a, b) => depthOf(a[1].json.id) - depthOf(b[1].json.id));

  const tally = {wanted: 0, notWritable: 0, absent: 0, nothingToAdd: 0, written: 0, reasserted: 0};
  const written: typeof order = [];
  for (const [id, {type, json}] of order) {
    progress('placing', ++tally.wanted, order.length);
    // Every project the bundle places is re-asserted, not only the ones whose own row was
    // written: one project's write takes contained entities from every other project holding them.
    const row = planned.get(id);
    if (['skip', 'failed'].includes(row?.action ?? '') && row?.reason !== 'personal_project') {
      tally.notWritable++;
      continue;
    }
    // Reading the relation rows on their own, rather than the whole project, is what keeps this
    // stage usable: nearly every project already holds what the bundle wants, and fetching each
    // one in full to discover that costs hours on a stand-sized push.
    // `GET /projects/relations` fails on a project that links domain-table rows, the same way it
    // does on the pull side — and reading "holds nothing" would strip every relation the target
    // has, because the write replaces the whole set. Fall back to the project's own copy.
    const held = await relations.listAll({projectId: json.id, include: 'entity'}).catch(() => null)
      ?? (await projects.find(json.id).catch(() => null))?.relations;
    if (!held) {
      tally.absent++;
      continue;
    }
    // `projects_repository.dart`: a non-link relation deletes the entity's other non-link rows.
    const contains = (entityId: string): boolean => ownerOf(entityId) === json.id;
    const claimed = new Set<string>((json.relations ?? []).map((r: any) => r?.entity?.id).filter(Boolean));
    // The server derives `is_link` — it keeps one container per entity and marks every other holder
    // a link — so only a claim is worth writing for. Writing the release direction is ignored and
    // recomputed, which would rewrite the project on every push, and a relation write costs tens of
    // seconds on a stand-sized target.
    let corrected = 0;
    const wanted: any[] = [];
    for (const r of held) {
      if (!r?.entity?.id) continue;
      const isLink = claimed.has(r.entity.id) ? !contains(r.entity.id) : (r.isLink ?? false);
      if ((r.isLink ?? false) && !isLink) corrected++;
      wanted.push({id: r.id, entity: {'#type': 'EntityRecord', id: r.entity.id}, isLink});
    }
    const linked = new Set<string>(wanted.map((r) => r.entity.id));
    let added = 0;
    for (const rel of json.relations) {
      if (linked.has(rel.entity.id)) continue;
      // The target stamps its own Files connection onto a space; re-attaching the source's
      // would leave the space with two of them.
      if (!landed.has(rel.entity.id)) {
        const existing = await findEntity(dapi, rel.entity.id);
        if (!existing || existing.parameters?.isProject === true) continue;
      }
      wanted.push(contains(rel.entity.id) ? rel : {...rel, isLink: true});
      linked.add(rel.entity.id);
      added++;
    }
    if (!added && !corrected) {
      tally.nothingToAdd++;
      continue;
    }
    // Only now is the project itself needed: the write posts it back whole.
    const target = await projects.find(json.id).catch(() => null);
    if (!target) {
      tally.absent++;
      continue;
    }
    tally.written++;
    // Persistently non-zero means the server is not keeping the flag as sent.
    if (!added) tally.reasserted++;
    delete target.storage;
    target.relations = wanted;

    const refused: string[] = [];
    for (let attempt = 0; ;) {
      try {
          await dapi.client.post('/projects?saveRelations=true', stripPrivate(JSON.parse(JSON.stringify(target))), BULK_WRITE_MS);
        break;
      } catch (err: any) {
        const text = String(err?.message ?? err);
        // The write is all-or-nothing, so a single entity the target refuses to link would cost
        // the whole space its placement. Drop the one it named and write the rest.
        // `projects_repository.dart` refuses a relation for want of a permission in several
        // wordings; dropping the entity then reports "the target would not link it" for what is
        // really the pusher's own access, and quietly migrates less.
        const denied = /privileges|you do not have/i.test(text);
        const bad = denied ? undefined : RELATION_REFUSED_RE.exec(text)?.[1];
        if (bad && refused.length < 3 && target.relations.some((r: any) => r.entity?.id === bad)) {
          target.relations = target.relations.filter((r: any) => r.entity?.id !== bad);
          refused.push(bad);
          continue;
        }
        if (attempt++ >= 2 || !(text.includes('40P01') || text.toLowerCase().includes('deadlock'))) {
          rows.push({name: nqNameOf(json), entityType: 'Project', action: 'failed', reason: 'relations', detail: text});
          break;
        }
        await new Promise((r) => setTimeout(r, 250 * attempt));
      }
    }
    if (refused.length)
      rows.push({name: nqNameOf(json), entityType: 'Project', action: 'warn', reason: 'relations_refused',
        detail: `${refused.length} left out because the target would not link them: ${refused.slice(0, 3).join(', ')}`});
    written.push([id, {type, json}]);
  }
  const short: string[] = [];
  for (const [, entity] of written) {
    const back = await projects.find(entity.json.id).catch(() => null);
    const kept = new Set<string>((back?.relations ?? []).map((r: any) => r?.entity?.id));
    if ((entity.json.relations ?? []).some((r: any) => r?.entity?.id && !kept.has(r.entity.id)))
      short.push(nqNameOf(entity.json));
  }
  tally.written -= short.length;
  if (short.length)
    rows.push({name: short.slice(0, 3).join(', '), entityType: 'Project',
      action: 'warn', reason: 'relations_not_kept',
      detail: `${short.length} project(s) lost relations the target dropped after writing them`});
  // Placement is what makes a dashboard openable, and it is invisible in a per-entity report.
  if (tally.wanted)
    rows.push({name: '', entityType: 'Project', action: 'info', reason: 'relations_written',
      detail: `${tally.written} of ${tally.wanted} placed` +
        (tally.reasserted ? ` (${tally.reasserted} only re-stating ownership)` : '') +
        (tally.notWritable ? `; ${tally.notWritable} not written` : '') +
        (tally.absent ? `; ${tally.absent} missing on target` : '') +
        (tally.nothingToAdd ? `; ${tally.nothingToAdd} already linked or unresolvable` : '')});
}

/**
 * Files a datasync table reads, written back into the share of the same name. A personal
 * `Home` share resolves to the target's own connection for that user, so each owner's files
 * land in their own share rather than the pusher's.
 */
async function pushShares(dapi: NodeDapi, bundle: Bundle, rows: Row[], onConflict: ConflictPolicy): Promise<void> {
  const conflicts: string[] = [];
  for (const remote of listShares(bundle.dir)) {
    try {
      // The only write that replaces a file the target already has, so it answers to the same
      // policy as everything else rather than overwriting whatever is there.
      const already = await dapi.files.readBytes(remote).then(() => true).catch(() => false);
      if (already && onConflict === 'skip') {
        rows.push({name: remote, entityType: 'File', action: 'skip', reason: 'share_file_exists'});
        continue;
      }
      if (already && onConflict === 'fail')
        conflicts.push(remote);
      await dapi.files.writeBytes(remote, fs.readFileSync(sharePath(bundle.dir, remote)));
      rows.push({name: remote, entityType: 'File', action: already ? 'update' : 'create', reason: 'share_file'});
    } catch (err: any) {
      rows.push({name: remote, entityType: 'File', action: 'warn', reason: 'share_file_not_written',
        detail: err?.message ?? String(err)});
    }
  }
  if (conflicts.length)
    throw new Error('The target already has these share files (use --on-conflict skip|adopt): ' +
      conflicts.join(', '));
}

/**
 * The server keeps a name unique within the namespace the entity is in, and an entity created
 * by a push is in the pusher's namespace until its space's relations are written. A name taken
 * there is usually free once the entity is placed. Returns true when anything was written back,
 * since a save re-homes the entity and the placement has to be asserted again.
 */
async function restoreNames(dapi: NodeDapi, renamed: {op: Op; row: Row}[],
                            idmap: Record<string, string>): Promise<boolean> {
  let saved = false;
  for (const {op, row} of renamed) {
    const spec = TYPES[op.type];
    const current = await dapi.internal(spec.route).find(idmap[op.id] ?? op.json.id).catch(() => null);
    if (!current || current.name === op.json.name) continue;
    current.name = op.json.name;
    delete current.storage;
    const restored = await dapi.internal(spec.saveRoute ?? spec.route)
      .save(stripPrivate(JSON.parse(JSON.stringify(current)))).catch(() => null);
    if (!restored) continue;
    // The entity was written either way, so placement has to be re-asserted either way.
    saved = true;
    if (restored.name !== op.json.name) continue;
    row.action = 'info';
    row.reason = 'name_restored';
    row.detail = `${op.json.name} was taken until the entity was placed`;
  }
  return saved;
}

/**
 * A namespace is a label derived from the owning space, so it is only worth judging once
 * relations and names are settled.
 */
async function reportPlacement(dapi: NodeDapi, misplaced: Op[], rows: Row[],
                               idmap: Record<string, string>): Promise<void> {
  for (const op of misplaced) {
    const now = await dapi.internal(TYPES[op.type].route).find(idmap[op.id] ?? op.json.id).catch(() => null);
    if (!now || (now.namespace ?? '') === op.expectedNamespace) continue;
    rows.push({name: op.row.name, entityType: op.type, action: 'warn', reason: 'namespace_not_preserved',
      detail: `${op.expectedNamespace} → ${now.namespace ?? ''}; pull the owning space so the entity travels in its relations`});
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
      const target = await dapi.internal(TYPES[type].route).find(json.id).catch(() => null);
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

/** Login -> the id of that user's personal root project on this instance. */
async function personalProjects(dapi: NodeDapi): Promise<Map<string, string>> {
  const byLogin = new Map<string, string>();
  for (const user of await dapi.internal('/users').listAll({limit: 500}))
    if (user?.login && user?.project?.id)
      byLogin.set(user.login, user.project.id);
  return byLogin;
}

async function currentNamespace(dapi: NodeDapi): Promise<string> {
  const user = await dapi.client.get('/users/current');
  return user?.project?.name ? `${user.project.name}:` : '';
}

/**
 * An entity pulled from the source author's personal namespace lands under the pusher's
 * own namespace on the target, so that is where a same-name twin would be.
 */
/** The target's own entity of that qualified name, whatever id it gave it. */
async function resolveByNqName(dapi: NodeDapi, type: string, nqName: string): Promise<any> {
  const cut = nqName.lastIndexOf(':');
  const namespace = cut === -1 ? '' : nqName.slice(0, cut + 1);
  const name = nqName.slice(cut + 1);
  if (!name) return null;
  const matches = await dapi.internal('/entities').list({namespace, name}).catch(() => []);
  return matches.find((m: any) => m['#type'] === type && inNamespace(m, namespace)) ?? null;
}

async function findByNqName(dapi: NodeDapi, bundle: Bundle, type: string, json: any, pusherNamespace: string): Promise<any> {
  // Not `expectedNamespace`: a twin for a namespace-less entity is looked up in the root.
  const sourceNamespace: string = json.namespace ?? '';
  const namespace = sourceNamespace === bundle.manifest.source.userNamespace ? pusherNamespace : sourceNamespace;
  const matches = await dapi.internal('/entities').list({namespace, name: json.name});
  return matches.find((m: any) => m['#type'] === type && inNamespace(m, namespace)) ?? null;
}
