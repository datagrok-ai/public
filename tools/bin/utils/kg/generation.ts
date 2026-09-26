/// A generation is one build's output directory (build-plan.md WO-11): `<out>/gen/<batch>-<suffix>/` holding
/// `manifest.json`, `data/nodes/<type>.jsonl`, `data/edges/<name>.jsonl` and `reports/`, and `<out>/current`
/// naming the one readers open. The layout is spelled here and nowhere else; `build/write.ts` fills it, the
/// index, the reports and the browser read it.
import * as fs from 'fs';
import * as path from 'path';
import {compare} from './normalize';

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
  dart_packages?: Record<string, number>;
  /** How deep the Dart pass reads: `lexical` while it is a regex pass over the sources. */
  dart_depth?: string;
  inventory?: Record<string, number>;
  /** The labels of the concrete edge types under each `edges/` folder, in schema order (`types.ts edgeGroups`). */
  edge_groups?: Record<string, string[]>;
  /** What each source contributes, by source name, from the extractors that ran (`Extractor.describes`). */
  provides?: Record<string, string>;
  /** The generation the index in this directory was loaded from; absent when it has none (`--no-db`). */
  indexed_batch?: string;
  /** Peak resident memory of the load, in MB, as telemetry; readers size their pool from `--memory`,
   * `KG_KUZU_MEMORY` or the 512 MB default. */
  index_memory_mb?: number;
  index_platform?: string;
}

/** One generation per build under `<out>/gen/`, and the one-line file naming the generation to read. */
export const GENERATIONS = 'gen';
export const CURRENT = 'current';
/** How long a directory under `gen/` without a manifest may be an interrupted build before `gc` removes it. */
const INTERRUPTED_MS = 3600_000;

export type DataKind = 'nodes' | 'edges';

export function dataDir(genDir: string, kind: DataKind): string {
  return path.join(genDir, 'data', kind);
}

/** `data/nodes/<type>.jsonl` or `data/edges/<edge name>.jsonl`. */
export function dataFile(genDir: string, kind: DataKind, name: string): string {
  return path.join(dataDir(genDir, kind), `${name}.jsonl`);
}

export function manifestFile(genDir: string): string {
  return path.join(genDir, 'manifest.json');
}

export function* readJsonl(file: string): Generator<Record<string, unknown>> {
  if (!fs.existsSync(file)) return;
  for (const line of fs.readFileSync(file, 'utf8').split('\n'))
    if (line) yield JSON.parse(line);
}

export function readManifest(genDir: string): Manifest | undefined {
  const file = manifestFile(genDir);
  if (!fs.existsSync(file)) return undefined;
  try {
    return JSON.parse(fs.readFileSync(file, 'utf8')) as Manifest;
  }
  catch {
    return undefined;
  }
}

/** One directory per build, named by its batch and a suffix: the batch is the content identity, the directory the
 * publication, so a rebuild never overwrites what a reader holds open. */
export function generationDir(outRoot: string, name: string): string {
  return path.join(outRoot, GENERATIONS, name);
}

/** A fresh generation directory, `<batch>-<six characters>`, created atomically. */
export function newGeneration(outRoot: string, batch: string): string {
  fs.mkdirSync(path.join(outRoot, GENERATIONS), {recursive: true});
  return fs.mkdtempSync(path.join(outRoot, GENERATIONS, `${batch}-`));
}

/** The generation `<out>/current` names, when it is complete. */
export function readCurrent(outRoot: string): string | undefined {
  const file = path.join(outRoot, CURRENT);
  if (!fs.existsSync(file)) return undefined;
  const name = fs.readFileSync(file, 'utf8').trim();
  return name && fs.existsSync(manifestFile(generationDir(outRoot, name))) ? name : undefined;
}

/** The directory to read a graph from: the current generation, or [outRoot] itself when it holds a build from
 * before generations existed (a committed public snapshot, an `--out` folder written by an older builder). */
export function currentDir(outRoot: string): string | undefined {
  const name = readCurrent(outRoot);
  if (name) return generationDir(outRoot, name);
  return fs.existsSync(manifestFile(outRoot)) ? outRoot : undefined;
}

/** Switches the pointer to the generation [name]: the file is replaced by a rename. */
export function publish(outRoot: string, name: string): void {
  const tmp = path.join(outRoot, `${CURRENT}.tmp`);
  fs.writeFileSync(tmp, `${name}\n`);
  fs.renameSync(tmp, path.join(outRoot, CURRENT));
}

export interface Generation {
  /** The directory name, `<batch>-<suffix>`. */
  name: string;
  batch: string;
  dir: string;
  built_at: string;
  current: boolean;
  indexed: boolean;
}

/** Every complete generation under `<out>/gen/`, newest first. */
export function generations(outRoot: string): Generation[] {
  const dir = path.join(outRoot, GENERATIONS);
  if (!fs.existsSync(dir)) return [];
  const current = readCurrent(outRoot);
  const found: Generation[] = [];
  for (const name of fs.readdirSync(dir)) {
    if (name === CURRENT) continue;
    const manifest = readManifest(path.join(dir, name));
    if (!manifest) continue;
    found.push({name, batch: manifest.batch, dir: path.join(dir, name), built_at: manifest.built_at, current: name === current,
      indexed: manifest.indexed_batch === manifest.batch});
  }
  return found.sort((a, b) => compare(b.built_at, a.built_at));
}

export interface GcResult {
  removed: string[];
  kept: string[];
  locked: string[];
}

/** Keeps the [keep] newest generations and the current one, removes the rest and every interrupted build; one
 * whose index a reader holds open cannot be removed on every platform, and is reported instead. */
export function gc(outRoot: string, keep: number): GcResult {
  const all = generations(outRoot);
  const result: GcResult = {removed: [], kept: [], locked: []};
  const remove = (name: string, dir: string) => {
    try {
      fs.rmSync(dir, {recursive: true, force: true});
      result.removed.push(name);
    }
    catch {
      result.locked.push(name);
    }
  };
  for (const [i, gen] of all.entries()) {
    if (i < keep || gen.current) result.kept.push(gen.name);
    else remove(gen.name, gen.dir);
  }
  for (const name of interrupted(outRoot, new Set(all.map((g) => g.name))))
    remove(name, path.join(outRoot, GENERATIONS, name));
  return result;
}

/** Directories under `gen/` with no manifest and no write for an hour: a build that died before it finished. */
function interrupted(outRoot: string, complete: Set<string>): string[] {
  const dir = path.join(outRoot, GENERATIONS);
  if (!fs.existsSync(dir)) return [];
  const out: string[] = [];
  for (const name of fs.readdirSync(dir)) {
    if (name === CURRENT || complete.has(name)) continue;
    const stat = fs.statSync(path.join(dir, name), {throwIfNoEntry: false});
    if (stat?.isDirectory() && Date.now() - stat.mtimeMs > INTERRUPTED_MS) out.push(name);
  }
  return out;
}
