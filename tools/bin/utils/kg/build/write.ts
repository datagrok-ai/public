/// The build outputs (build-plan.md WO-1, Decisions): `data/nodes/<type>.jsonl`, `data/edges/<type|refname>.jsonl`,
/// `manifest.json` and `reports/`, deterministic line by line; the public projection of a graph.
import * as fs from 'fs';
import * as path from 'path';
import {createHash} from 'crypto';
import {spawnSync} from 'child_process';
import {TypeSystem, isSubtype} from '../types';
import {Graph} from './emitter';
import {Row} from './normalize';

export interface Manifest {
  built_at: string;
  schema_version: number;
  builder: string;
  batch: string;
  mode: string;
  revisions: Record<string, string>;
  sources: Record<string, string>;
  counts: {nodes: Record<string, number>, edges: Record<string, number>};
  problems: Record<string, number>;
}

export interface BuildInfo {
  mode: string;
  batch: string;
  builder: string;
  schemaVersion: number;
  revisions: Record<string, string>;
}

/** Node types the public snapshot carries (Decisions "Public mode"); a doc-page only when it is a help page. */
const PUBLIC_TYPES = ['feature', 'concept', 'package', 'library', 'doc-page', 'doc-anchor', 'sample', 'scenario'];
const HEAD_KEYS = ['id', 'type', 'name', 'from', 'to'];
const INVALID_CAP = 500;

export function toolsVersion(): string {
  return JSON.parse(fs.readFileSync(path.join(__dirname, '..', '..', '..', '..', 'package.json'), 'utf8')).version;
}

/** HEAD of the monorepo, HEAD of the public checkout and the submodule pin; `unknown` where git cannot answer. */
export function gitRevisions(repoRoot: string): Record<string, string> {
  const sha = (cwd: string, ref: string) => {
    const r = spawnSync('git', ['rev-parse', '--verify', ref], {cwd, encoding: 'utf8'});
    return r.status === 0 ? r.stdout.trim() : 'unknown';
  };
  return {reddata: sha(repoRoot, 'HEAD'), public: sha(path.join(repoRoot, 'public'), 'HEAD'), public_pin: sha(repoRoot, 'HEAD:public')};
}

/** Content-addressed: two builds of the same revisions with the same builder and schema share a batch. */
export function batchId(revisions: Record<string, string>, schemaVersion: number, builder: string): string {
  const hash = createHash('sha1').update([revisions.reddata, revisions.public, schemaVersion, builder].join('\n')).digest('hex');
  return `b-${hash.slice(0, 12)}`;
}

/** Public logical nodes with their descriptions; the home, the owner, non-public paths and every edge with a non-public end stay behind. */
export function projectPublic(graph: Graph, system: TypeSystem): Graph {
  const keep = (row: Row) => {
    const type = String(row.type);
    if (row.visibility !== 'public' || !PUBLIC_TYPES.some((t) => isSubtype(system, type, t))) return false;
    return !isSubtype(system, type, 'doc-page') || row.kind === 'help';
  };
  const nodes = graph.nodes.filter(keep);
  const ids = new Set(nodes.map((n) => n.id as string));
  const isPublicPath = (p: unknown) => typeof p === 'string' && (p.startsWith('public/') || p.startsWith('landing/'));
  const projected = nodes.map((row) => {
    const out: Row = {};
    const members = system.nodes.get(row.type as string)!.members;
    for (const [k, v] of Object.entries(row)) {
      if (k === 'home' || k === 'owner') continue;
      const m = members[k];
      if (m?.kind === 'ref' && !(m.list ? (v as string[]).every((id) => ids.has(id)) : ids.has(v as string))) continue;
      if (m?.scalar === 'Path' && !(m.list ? (v as string[]).every(isPublicPath) : isPublicPath(v))) continue;
      out[k] = v;
    }
    return out;
  });
  const edges = graph.edges.filter((e) => ids.has(e.from as string) && ids.has(e.to as string) && !(e.type === 'ref' && e.name === 'owner')).map((e) => {
    const evidence = (e.evidence as string[] | undefined)?.filter(isPublicPath);
    const out = {...e};
    delete out.evidence;
    return evidence?.length ? {...out, evidence} : out;
  });
  return {...graph, nodes: projected, edges};
}

/** Writes everything under [outRoot] and returns the manifest. `data/` and `reports/` are replaced whole. */
export function writeBuild(graph: Graph, outRoot: string, info: BuildInfo): Manifest {
  const dataDir = path.join(outRoot, 'data');
  const reportsDir = path.join(outRoot, 'reports');
  for (const dir of [dataDir, reportsDir]) fs.rmSync(dir, {recursive: true, force: true});
  fs.mkdirSync(path.join(dataDir, 'nodes'), {recursive: true});
  fs.mkdirSync(path.join(dataDir, 'edges'), {recursive: true});
  fs.mkdirSync(reportsDir, {recursive: true});

  const nodes = groupBy(graph.nodes, (n) => String(n.type));
  const edges = groupBy(graph.edges, (e) => e.type === 'ref' ? String(e.name) : String(e.type));
  const counts: Manifest['counts'] = {nodes: {}, edges: {}};
  for (const [type, rows] of nodes) {
    rows.sort((a, b) => compare(String(a.id), String(b.id)));
    writeJsonl(path.join(dataDir, 'nodes', `${type}.jsonl`), rows);
    counts.nodes[type] = rows.length;
  }
  for (const [name, rows] of edges) {
    rows.sort((a, b) => compare(String(a.type), String(b.type)) || compare(String(a.from), String(b.from)) || compare(String(a.to), String(b.to)) || compare(String(a.name ?? ''), String(b.name ?? '')));
    writeJsonl(path.join(dataDir, 'edges', `${name}.jsonl`), rows);
    counts.edges[name] = rows.length;
  }
  const claims = [...graph.claims].sort((a, b) => compare(a.file, b.file) || compare(a.feature, b.feature) || a.rung - b.rung);
  writeJsonl(path.join(reportsDir, 'claims.jsonl'), claims as unknown as Row[]);
  if (graph.invalid.length) writeJsonl(path.join(reportsDir, 'invalid.jsonl'), graph.invalid.slice(0, INVALID_CAP));
  if (Object.keys(graph.details).length) fs.writeFileSync(path.join(reportsDir, 'problems.json'), `${JSON.stringify(graph.details, null, 2)}\n`);

  const manifest: Manifest = {
    built_at: new Date().toISOString(),
    schema_version: info.schemaVersion,
    builder: info.builder,
    batch: info.batch,
    mode: info.mode,
    revisions: info.revisions,
    sources: sortKeys(graph.sources),
    counts: {nodes: sortKeys(counts.nodes), edges: sortKeys(counts.edges)},
    problems: sortKeys(graph.problems),
  };
  fs.writeFileSync(path.join(outRoot, 'manifest.json'), `${JSON.stringify(manifest, null, 2)}\n`);
  return manifest;
}

function groupBy(rows: Row[], key: (row: Row) => string): Map<string, Row[]> {
  const out = new Map<string, Row[]>();
  for (const row of rows) {
    let list = out.get(key(row));
    if (!list) out.set(key(row), list = []);
    list.push(row);
  }
  return new Map([...out].sort(([a], [b]) => compare(a, b)));
}

/** One object per line, keys in a fixed order: id, type, name, from, to, then the rest alphabetically. */
function writeJsonl(file: string, rows: Row[]): void {
  const lines = rows.map((row) => {
    const keys = [...HEAD_KEYS.filter((k) => row[k] !== undefined), ...Object.keys(row).filter((k) => !HEAD_KEYS.includes(k)).sort(compare)];
    return JSON.stringify(Object.fromEntries(keys.map((k) => [k, row[k]])));
  });
  fs.writeFileSync(file, lines.length ? `${lines.join('\n')}\n` : '');
}

function sortKeys<T>(record: Record<string, T>): Record<string, T> {
  return Object.fromEntries(Object.entries(record).sort(([a], [b]) => compare(a, b)));
}

/** Code-point order, the same on every platform and locale. */
function compare(a: string, b: string): number {
  return a < b ? -1 : a > b ? 1 : 0;
}
