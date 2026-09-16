/// The build outputs (build-plan.md WO-1, Decisions): `data/nodes/<type>.jsonl`, `data/edges/<type|refname>.jsonl`,
/// `manifest.json` and `reports/`, deterministic line by line; the public projection of a graph.
import * as fs from 'fs';
import * as path from 'path';
import {createHash} from 'crypto';
import {spawnSync} from 'child_process';
import {TypeSystem, isSubtype} from '../types';
import {Graph} from './emitter';
import {Row, isOrdered, compare} from '../normalize';
import {Manifest, dataDir, manifestFile} from '../generation';

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

/** What a build actually reads, so that two graphs cannot share a batch (kg-codex-review-3 #12). */
export interface BuildInputs {
  repoRoot: string;
  mode: string;
  schemaVersion: number;
  builder: string;
  revisions: Record<string, string>;
  /** The extractors this build runs; order does not matter. */
  extractors: string[];
  backlogDir?: string;
  /** Where this build writes: its own output is not one of its inputs. */
  outRoot?: string;
}

/** Content-addressed: the revisions, the dirty tree of both repositories, the schema and builder versions, the
 * mode, the extractor selection and the backlog snapshot. */
export function buildInputs(inputs: BuildInputs): Record<string, string> {
  return {
    reddata: inputs.revisions.reddata,
    public: inputs.revisions.public,
    dirty: dirtyDigest(inputs.repoRoot, inputs.outRoot),
    schema_version: String(inputs.schemaVersion),
    builder: inputs.builder,
    mode: inputs.mode,
    extractors: [...inputs.extractors].sort(compare).join(','),
    backlog: backlogInput(inputs.backlogDir),
  };
}

export function batchId(inputs: Record<string, string>): string {
  const text = Object.entries(inputs).sort(([a], [b]) => compare(a, b)).map(([k, v]) => `${k}=${v}`).join('\n');
  return `b-${createHash('sha1').update(text).digest('hex').slice(0, 12)}`;
}

/** `git status --porcelain` of the monorepo and of the public checkout, with the content of every modified
 * tracked file: an uncommitted edit changes the graph, so it changes the batch. What the build itself writes
 * under [outRoot] is left out — its own output is not an input, gitignored or not. */
function dirtyDigest(repoRoot: string, outRoot?: string): string {
  const hash = createHash('sha1');
  const out = outRoot === undefined ? undefined : path.resolve(outRoot);
  const seen = new Set<string>();
  for (const dir of [repoRoot, path.join(repoRoot, 'public')]) {
    // porcelain paths are relative to the repository root, and `public/` is one only when it is the submodule
    const top = spawnSync('git', ['rev-parse', '--show-toplevel'], {cwd: dir, encoding: 'utf8'});
    if (top.status !== 0) {
      hash.update(`${dir}: unknown\n`);
      continue;
    }
    const cwd = path.resolve(top.stdout.trim());
    if (seen.has(cwd)) continue;
    seen.add(cwd);
    const r = spawnSync('git', ['status', '--porcelain', '-z'], {cwd, encoding: 'utf8', maxBuffer: 256 * 1024 * 1024});
    if (r.status !== 0) {
      hash.update(`${cwd}: unknown\n`);
      continue;
    }
    const entries = r.stdout.split('\0');
    for (let i = 0; i < entries.length; i++) {
      const entry = entries[i];
      if (!entry) continue;
      const status = entry.slice(0, 2);
      const file = path.resolve(cwd, entry.slice(3));
      const renamed = (status.startsWith('R') || status.startsWith('C')) ? entries[++i] ?? '' : '';
      if (out !== undefined && inside(out, file)) continue;
      hash.update(`${entry}\n${renamed}`);
      if (status !== '??') hash.update(digestOf(file));
    }
  }
  return hash.digest('hex');
}

/** The backlog watermark and the digest of the snapshot it came from; `missing` when there is none. */
function backlogInput(dir?: string): string {
  if (!dir || !fs.existsSync(path.join(dir, 'index.jsonl'))) return 'missing';
  const file = path.join(dir, 'index.jsonl');
  const text = fs.readFileSync(file, 'utf8');
  let watermark = '';
  for (const line of text.split('\n')) {
    if (!line.trim()) continue;
    try {
      const row = JSON.parse(line);
      if (typeof row.updated === 'string' && row.updated > watermark) watermark = row.updated;
    }
    catch {
      // a line the process layer will count as invalid; the digest below covers it
    }
  }
  return `${watermark || 'none'}:${createHash('sha1').update(text).digest('hex')}`;
}

function inside(dir: string, file: string): boolean {
  const r = path.relative(dir, file);
  return !r.startsWith('..') && !path.isAbsolute(r);
}

/** The sha1 of a file, or `absent`; a directory and an unreadable path are `absent` too. */
function digestOf(file: string): string {
  try {
    const stat = fs.statSync(file);
    return stat.isFile() ? createHash('sha1').update(fs.readFileSync(file)).digest('hex') : 'absent';
  }
  catch {
    return 'absent';
  }
}

/** Public logical nodes with their descriptions; the home, the owner, non-public paths, every edge with a non-public end, the
 * claims and the full-graph problem records stay behind. */
export function projectPublic(graph: Graph, system: TypeSystem): Graph {
  const pages = new Map<string, Row>();
  for (const row of graph.nodes)
    if (isSubtype(system, String(row.type), 'doc-page')) pages.set(String(row.id), row);
  /** A heading is no more public than the page that carries it, and never outlives it (kg-codex-review-3 #1). */
  const keep = (row: Row): boolean => {
    const type = String(row.type);
    if (row.visibility !== 'public' || !PUBLIC_TYPES.some((t) => isSubtype(system, type, t))) return false;
    if (isSubtype(system, type, 'doc-page')) return row.kind === 'help';
    if (!isSubtype(system, type, 'doc-anchor')) return true;
    const page = pages.get(String(row.page));
    return page !== undefined && keep(page);
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
  const stubs = graph.stubs.filter((id) => ids.has(id));
  const problems = Object.fromEntries(Object.keys(graph.problems).map((k) => [k, k === 'partial_stubs' ? stubs.length : 0]));
  problems.projection_dropped = closure(projected, edges, ids, system);
  return {nodes: projected, edges, stubs, claims: [], sources: graph.sources, problems, details: {}, invalid: [], reports: {}, manifest: {}};
}

/** The projection is closed or it is not published: every edge end and every reference must name a node that
 * survived. Anything left over is removed here and counted, so a leak is visible in the manifest. */
function closure(nodes: Row[], edges: Row[], ids: Set<string>, system: TypeSystem): number {
  let dropped = 0;
  for (let i = edges.length - 1; i >= 0; i--)
    if (!ids.has(edges[i].from as string) || !ids.has(edges[i].to as string)) {
      edges.splice(i, 1);
      dropped++;
    }
  for (const row of nodes) {
    const members = system.nodes.get(row.type as string)!.members;
    for (const [k, v] of Object.entries(row)) {
      const m = members[k];
      if (m?.kind !== 'ref') continue;
      const ok = m.list ? (v as string[]).every((id) => ids.has(id)) : ids.has(v as string);
      if (ok) continue;
      delete row[k];
      dropped++;
    }
  }
  return dropped;
}

/** Writes the data and the reports of one generation into [genDir] and returns its manifest, which the caller
 * writes last (after the index, if any): a generation without a manifest is an interrupted build. The public
 * snapshot gets no reports and only the public revision. */
export function writeBuild(graph: Graph, genDir: string, info: BuildInfo): Manifest {
  const isPublic = info.mode === 'public';
  const nodesDir = dataDir(genDir, 'nodes'), edgesDir = dataDir(genDir, 'edges');
  const reportsDir = path.join(genDir, 'reports');
  for (const dir of [path.join(genDir, 'data'), reportsDir]) fs.rmSync(dir, {recursive: true, force: true});
  fs.mkdirSync(nodesDir, {recursive: true});
  fs.mkdirSync(edgesDir, {recursive: true});

  const nodes = groupBy(graph.nodes, (n) => String(n.type));
  const edges = groupBy(graph.edges, (e) => e.type === 'ref' ? String(e.name) : String(e.type));
  const counts: Manifest['counts'] = {nodes: {}, edges: {}};
  for (const [type, rows] of nodes) {
    rows.sort((a, b) => compare(String(a.id), String(b.id)));
    writeJsonl(path.join(nodesDir, `${type}.jsonl`), rows);
    counts.nodes[type] = rows.length;
  }
  for (const [name, rows] of edges) {
    rows.sort((a, b) => compare(String(a.type), String(b.type)) || compare(String(a.from), String(b.from)) || compare(String(a.to), String(b.to)) || compare(String(a.name ?? ''), String(b.name ?? '')) || compare(JSON.stringify(a), JSON.stringify(b)));
    writeJsonl(path.join(edgesDir, `${name}.jsonl`), rows);
    counts.edges[name] = rows.length;
  }
  if (!isPublic) {
    fs.mkdirSync(reportsDir, {recursive: true});
    const claims = [...graph.claims].sort((a, b) => compare(a.file, b.file) || compare(a.feature, b.feature) || a.rung - b.rung);
    writeJsonl(path.join(reportsDir, 'claims.jsonl'), claims as unknown as Row[]);
    if (graph.invalid.length) writeJsonl(path.join(reportsDir, 'invalid.jsonl'), graph.invalid.slice(0, INVALID_CAP));
    if (Object.keys(graph.details).length) fs.writeFileSync(path.join(reportsDir, 'problems.json'), `${JSON.stringify(graph.details, null, 2)}\n`);
    for (const [name, data] of Object.entries(graph.reports).sort(([a], [b]) => compare(a, b)))
      fs.writeFileSync(path.join(reportsDir, `${name}.json`), `${JSON.stringify(data, null, 2)}\n`);
  }

  const manifest: Manifest = {
    built_at: new Date().toISOString(),
    schema_version: info.schemaVersion,
    builder: info.builder,
    batch: info.batch,
    mode: info.mode,
    revisions: isPublic ? {public: info.revisions.public} : info.revisions,
    sources: sortKeys(graph.sources),
    counts: {nodes: sortKeys(counts.nodes), edges: sortKeys(counts.edges)},
    problems: sortKeys(graph.problems),
    ...graph.manifest,
  };
  return manifest;
}

/** The last write of a generation: until this file exists the generation is a build in flight. */
export function writeManifest(genDir: string, manifest: Manifest): void {
  fs.writeFileSync(manifestFile(genDir), `${JSON.stringify(manifest, null, 2)}\n`);
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

/** One object per line, keys in a fixed order: id, type, name, from, to, then the rest alphabetically; every list
 * that is not an ordered member (normalize.ts `isOrdered`) is written sorted, so that the order its rows arrived
 * in never reaches the bytes. */
function writeJsonl(file: string, rows: Row[]): void {
  const lines = rows.map((row) => {
    const keys = [...HEAD_KEYS.filter((k) => row[k] !== undefined), ...Object.keys(row).filter((k) => !HEAD_KEYS.includes(k)).sort(compare)];
    return JSON.stringify(Object.fromEntries(keys.map((k) => [k, sorted(k, row[k])])));
  });
  fs.writeFileSync(file, lines.length ? `${lines.join('\n')}\n` : '');
}

function sorted(key: string, value: unknown): unknown {
  return Array.isArray(value) && !isOrdered(key) ? [...value].sort((a, b) => compare(String(a), String(b))) : value;
}

function sortKeys<T>(record: Record<string, T>): Record<string, T> {
  return Object.fromEntries(Object.entries(record).sort(([a], [b]) => compare(a, b)));
}

