/// Docs: [Entity export / import](/docs/features/grok-tool/export-import/DESIGN.md)
import * as fs from 'fs';
import * as path from 'path';
import {createHash} from 'crypto';
import {BytesKind, TYPES, fileNameFor, nqNameOf, rankOf} from './registry';

export interface ManifestEntry {type: string; id: string; file: string}

export interface Manifest {
  formatVersion: 1;
  source: {url: string; version: string; commit?: string; userNamespace: string};
  pulls: {at: string; args: string[]; ids: string[]}[];
  order: ManifestEntry[];
  packages: string[];
}

export interface BundleEntity {type: string; json: any; file?: string}

export interface Bundle {
  dir: string;
  manifest: Manifest;
  entities: Map<string, BundleEntity>;
  idmap: Record<string, string>;
}

// `entityTags` rows carry their own primary keys and come back in an unstable order:
// pushing them makes the hash flap AND re-points the source's tag rows at the copy.
const VOLATILE = ['createdOn', 'updatedOn', 'author', 'pictureId', 'encryptedParametersId',
  'keyKid', 'isAvailable', 'entityTags'];

export function normalize(type: string, json: any): any {
  const copy = JSON.parse(JSON.stringify(json));
  // Tag rows are volatile, but the tag names are worth migrating: keep them as a
  // bundle-only list the pusher replays through `POST /entities/tag`.
  const tags = [...new Set((copy.entityTags ?? []).map((t: any) => t?.tag).filter(Boolean))].sort();
  if (tags.length)
    copy._tags = tags;
  for (const k of VOLATILE)
    delete copy[k];
  if (copy.package)
    copy.package = {id: copy.package.id};
  // The stamp is written by the pusher, so it must not make an unchanged entity look different.
  if (copy.metaParams)
    delete copy.metaParams.sync_id;
  TYPES[type]?.strip?.(copy);
  return sortKeys(copy);
}

function sortKeys(v: any): any {
  if (Array.isArray(v)) return v.map(sortKeys);
  if (v === null || typeof v !== 'object') return v;
  const out: any = {};
  for (const k of Object.keys(v).sort())
    out[k] = sortKeys(v[k]);
  return out;
}

/** Bundle-only keys (`_credentials`, `_grants`, `_members`) never travel to the server. */
export function stripPrivate(json: any): any {
  for (const k of Object.keys(json))
    if (k.startsWith('_'))
      delete json[k];
  return json;
}

/**
 * What a comparison sees: the payload without the bundle-only keys and without the
 * namespace, which the server computes from ownership (`Askalkin:` on the source,
 * `Admin:` on the target) and would otherwise make every cross-instance push a rewrite.
 * The bundle file keeps it — `findByNqName` resolves the twin by it.
 */
export function hashView(type: string, json: any): any {
  const view = stripPrivate(normalize(type, json));
  delete view.namespace;
  // Relations have a pass of their own, which only ever adds: their row ids belong to the
  // target, and a target that links more than the bundle is not stale. The pusher compares
  // them by coverage instead.
  delete view.relations;
  return view;
}

export function hashOf(type: string, json: any): string {
  return createHash('sha256').update(JSON.stringify(hashView(type, json))).digest('hex');
}

export function bytesPath(dir: string, kind: BytesKind, id: string): string {
  return path.join(dir, kind, kind === 'tables' ? `${id}.d42` : id);
}

const EMPTY: Manifest = {
  formatVersion: 1,
  source: {url: '', version: '', userNamespace: ''},
  pulls: [],
  order: [],
  packages: [],
};

function readManifest(dir: string): Manifest {
  const file = path.join(dir, 'manifest.json');
  if (!fs.existsSync(file)) return JSON.parse(JSON.stringify(EMPTY));
  return JSON.parse(fs.readFileSync(file, 'utf8'));
}

export function write(dir: string, entities: Map<string, BundleEntity>,
                      meta: {source: Manifest['source']; args: string[]; packages: string[]},
                      opts: {replace?: boolean},
                      bytes: Map<string, Buffer> = new Map()): Manifest {
  if (opts.replace)
    fs.rmSync(dir, {recursive: true, force: true});
  const manifest = readManifest(dir);
  if (!manifest.source.url)
    manifest.source = meta.source;
  else if (manifest.source.url !== meta.source.url)
    throw new Error(`${dir} was pulled from ${manifest.source.url}; ` +
      `use --replace to start a new bundle from ${meta.source.url}`);
  fs.mkdirSync(dir, {recursive: true});

  const byId = new Map<string, ManifestEntry>(manifest.order.map((e) => [e.id, e]));
  const takenBy = new Map<string, string>([...byId.values()].map((e) => [e.file, e.id]));
  const ids: string[] = [];
  for (const [id, {type, json}] of entities) {
    const base = `${type}/${fileNameFor(json)}`;
    // Names are not unique across a bundle (two files of the same name in different shares);
    // whoever claimed the plain name keeps it, the rest are suffixed by their id.
    const owner = takenBy.get(`${base}.json`);
    const file = !owner || owner === id ? `${base}.json` : `${base}-${id.slice(0, 8)}.json`;
    const previous = byId.get(id);
    if (previous && previous.file !== file) {
      fs.rmSync(path.join(dir, previous.file), {force: true});
      takenBy.delete(previous.file);
    }
    takenBy.set(file, id);
    fs.mkdirSync(path.join(dir, type), {recursive: true});
    fs.writeFileSync(path.join(dir, file), JSON.stringify(normalize(type, json), null, 2));
    byId.set(id, {type, id, file});
    ids.push(id);
  }

  for (const [id, buf] of bytes) {
    const kind = TYPES[entities.get(id)!.type].bytes!.kind;
    fs.mkdirSync(path.join(dir, kind), {recursive: true});
    fs.writeFileSync(bytesPath(dir, kind, id), buf);
  }

  manifest.formatVersion = 1;
  manifest.pulls.push({at: new Date().toISOString(), args: meta.args, ids});
  manifest.packages = [...new Set([...manifest.packages, ...meta.packages])].sort();
  const pullOf = new Map<string, number>();
  manifest.pulls.forEach((p, i) => p.ids.forEach((id) => pullOf.has(id) || pullOf.set(id, i)));
  manifest.order = [...byId.values()].sort((a, b) =>
    rankOf(a.type) - rankOf(b.type) || pullOf.get(a.id)! - pullOf.get(b.id)! || a.file.localeCompare(b.file));
  fs.writeFileSync(path.join(dir, 'manifest.json'), JSON.stringify(manifest, null, 2));
  return manifest;
}

export function read(dir: string): Bundle {
  if (!fs.existsSync(path.join(dir, 'manifest.json')))
    throw new Error(`Not a bundle directory (no manifest.json): ${dir}`);
  const manifest = readManifest(dir);
  const entities = new Map<string, BundleEntity>();
  for (const e of manifest.order) {
    const file = path.join(dir, e.file);
    if (!fs.existsSync(file))
      throw new Error(`Bundle is missing ${e.file} listed in manifest.json`);
    entities.set(e.id, {type: e.type, json: JSON.parse(fs.readFileSync(file, 'utf8')), file: e.file});
  }
  const idmapFile = path.join(dir, 'idmap.json');
  const idmap = fs.existsSync(idmapFile) ? JSON.parse(fs.readFileSync(idmapFile, 'utf8')) : {};
  return {dir, manifest, entities, idmap};
}

/** Adopted `sourceId → targetId` pairs, so the next push of this bundle is stable. */
export function writeIdmap(dir: string, idmap: Record<string, string>): void {
  fs.writeFileSync(path.join(dir, 'idmap.json'), JSON.stringify(idmap, null, 2));
}

export function list(dir: string): any[] {
  const {manifest, entities} = read(dir);
  return manifest.order.map((e) => {
    const pull = manifest.pulls.findIndex((p) => p.ids.includes(e.id));
    const json = entities.get(e.id)!.json;
    return {
      type: e.type,
      nqName: nqNameOf(json),
      file: e.file,
      pulledOn: manifest.pulls[pull]?.at ?? '',
      pull: pull + 1,
    };
  });
}
