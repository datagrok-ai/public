/// Membership resolution (conventions.md §8, build-plan.md WO-4): the one feature that owns each
/// file, the features that only participate in it, the tests that follow an owned file, and the
/// report of what stayed ambiguous or orphaned. Runs after every extractor has made its claims.
import * as fs from 'fs';
import * as path from 'path';
import {Emitter, Claim} from '../emitter';
import {Row} from '../../normalize';
import {BuildContext, Extractor} from '../context';
import {REPO_PREFIX, GLOB_MAGIC} from '../../homes';
import {fileId, posix, docId, CODE_ROOTS, countLines, sourceFileRow, unitOf} from '../../ids';
import {homesOf} from './markers';

/** Claim properties an ownership edge carries; anything else the `code:` item said stays in the claim. */
const OWNER_PROPS = ['role', 'layer'];
const RUNGS: (1 | 2 | 3)[] = [1, 2, 3];
/** A TypeScript file no mirror may point at: in a tests folder, or named as a test. */
const TS_TEST_PATH = /(?:^|\/)(?:tests|__tests__)\/|\.(?:test|spec)\.(?:tsx?|[cm]?js)$/;
const MIRROR_NAMES: [RegExp, string][] = [[/^(.+)_test\.dart$/, '_'], [/^(.+)-tests\.ts$/, '-'], [/^(.+)(?:\.test|-test)\.ts$/, '']];
/** The Test Track specs, whose folder names the feature. */
const TEST_TRACK = 'public/packages/UsageAnalysis/files/TestTrack/';
const TEST_WORDS = ['tests', 'test', 'spec', 'ui'];
const NAME_CONFIDENCE = 0.7;

export const membershipExtractor: Extractor = {
  name: 'membership',
  describes: {membership: 'file ownership'},
  modes: ['full'],
  run(ctx: BuildContext, emitter: Emitter): void {
    new Membership(ctx, emitter).run();
  },
};

interface Owner {
  feature: string;
  rung: 1 | 2 | 3 | 4;
}

class Membership {
  /** Home file per feature, for the evidence of a claim and for the folders of rung 4. */
  private homeFile = new Map<string, string>();
  private featureName = new Map<string, string>();
  /** Folder -> the features rung 4 may inherit from it; a folder several features claim is no candidate. */
  private folders = new Map<string, Set<string>>();
  private owners = new Map<string, string>();
  private ambiguous: {file: string, features: string[], rung: number}[] = [];
  private chained: {file: string, owner: string, over: string[]}[] = [];
  private orphans: {file: string, loc: number}[] = [];
  /** The denominator every ownership number is a fraction of: what this build actually saw. */
  private inventory = {observed_files: 0, observed_loc: 0, owned_files: 0, owned_loc: 0, participating_files: 0};

  constructor(private ctx: BuildContext, private emitter: Emitter) {}

  run(): void {
    this.readHomes();
    const claims = new Map<string, Claim[]>();
    for (const claim of this.emitter.claimed) {
      const list = claims.get(claim.file);
      if (list) list.push(claim);
      else claims.set(claim.file, [claim]);
    }
    const files = new Map<string, Row>();
    for (const row of this.emitter.rowsOf('source-file')) files.set(String(row.path), row);
    for (const file of claims.keys())
      if (!files.has(file)) {
        const row = this.fileNode(file);
        if (row) files.set(file, row);
      }
    for (const file of [...files.keys()].sort()) this.resolve(file, files.get(file)!, claims.get(file) ?? []);
    this.testEdges();
    this.nameEdges();
    this.mirrorEdges(files);
    this.helpEdges();
    this.citationEdges(files);
    this.orphans.sort((a, b) => b.loc - a.loc || (a.file < b.file ? -1 : 1));
    this.emitter.manifest('inventory', this.inventory);
    this.emitter.report('ownership', {inventory: this.inventory, ambiguous: this.ambiguous, orphans: this.orphans, resolved_by_chain: this.chained});
  }

  /** Feature homes: where each one lives, and the folder a file with no claim of its own inherits from (rung 4). */
  private readHomes(): void {
    for (const home of homesOf(this.ctx).homes) {
      if (home.type.root !== 'feature') continue;
      this.homeFile.set(home.id, home.file);
      this.featureName.set(home.id, home.name);
      const roots = codeRoots(home.data.code);
      // every root is a candidate, however many the feature has: rung 4 asks whether the folder is unambiguous, not the feature
      for (const folder of roots.map((r) => this.rootFolder(r)))
        if (folder !== undefined) this.claimFolder(folder, home.id);
      if (!roots.some(under)) this.claimFolder(path.posix.dirname(home.file), home.id);
    }
  }

  private claimFolder(folder: string, feature: string): void {
    const features = this.folders.get(folder);
    if (features) features.add(feature);
    else this.folders.set(folder, new Set([feature]));
  }

  /** The folder a `code:` root stands for: its part before any glob, when that is a directory on disk. */
  private rootFolder(target: string): string | undefined {
    const segments = target.split('/');
    const magic = segments.findIndex((s) => GLOB_MAGIC.test(s));
    const clean = (magic < 0 ? segments : segments.slice(0, magic)).join('/');
    const full = path.join(this.ctx.repoRoot, clean);
    return clean && fs.existsSync(full) && fs.statSync(full).isDirectory() ? clean : undefined;
  }

  /** The node a claimed path needs when no extractor emitted one; a cited folder gets none. */
  private fileNode(file: string): Row | undefined {
    const full = path.join(this.ctx.repoRoot, file);
    if (!fs.existsSync(full) || !fs.statSync(full).isFile()) return undefined;
    const row = sourceFileRow(file, {loc: countLines(fs.readFileSync(full))});
    return this.emitter.node(row).accepted ? row : undefined;
  }

  /**
   * §8 for one file: the first rung yielding exactly one owner wins, rung 2 breaking a nested claim by the part-of
   * chain, rung 3 owning only what it is also the nearest home of, rung 4 inheriting from the folder. Every other
   * feature that claimed the file participates in it, ancestors of the owner excepted: part-of already says that.
   */
  private resolve(file: string, row: Row, claims: Claim[]): void {
    const loc = Number(row.loc ?? 0);
    this.inventory.observed_files++;
    this.inventory.observed_loc += loc;
    const owning = claims.filter((c) => c.mode !== 'participates');
    let owner: Owner | undefined;
    let ambiguous = false;
    for (const rung of RUNGS) {
      const features = [...new Set(owning.filter((c) => c.rung === rung).map((c) => c.feature))].sort();
      if (!features.length) continue;
      let picked = features;
      if (rung === 3) picked = features.filter((f) => this.nearestHomes(file).includes(f));
      else if (rung === 2 && features.length > 1) {
        const leaf = chainLeaf(features);
        if (leaf) {
          picked = [leaf];
          this.chained.push({file, owner: leaf, over: features.filter((f) => f !== leaf)});
        }
      }
      if (picked.length === 1) {
        owner = {feature: picked[0], rung};
        break;
      }
      if (rung === 3 && !picked.length) continue;
      this.ambiguous.push({file, features, rung});
      this.emitter.problem('ambiguous_owners', `${file}: ${features.join(' and ')} claim it at rung ${rung}`);
      ambiguous = true;
      break;
    }
    if (!owner && !ambiguous) {
      const inherited = this.folderOwner(file);
      if (inherited) owner = {feature: inherited, rung: 4};
    }
    const id = fileId(file);
    if (owner) {
      this.owners.set(file, owner.feature);
      this.inventory.owned_files++;
      this.inventory.owned_loc += loc;
      const claim = owning.find((c) => c.rung === owner!.rung && c.feature === owner!.feature);
      this.emitter.edge({type: 'is-implemented-in', from: owner.feature, to: id, derived_by: owner.rung === 4 ? 'filesystem' : 'annotation',
        confidence: owner.rung === 4 ? 0.9 : 1, evidence: this.evidence(owner.feature, file, owner.rung), ...pick(claim?.props, OWNER_PROPS)});
    }
    else if (CODE_ROOTS.some((root) => file.startsWith(root))) {
      this.orphans.push({file, loc});
      this.emitter.problem('orphans');
    }
    let participates = false;
    for (const claim of claims) {
      if (claim.feature === owner?.feature || (owner && owner.feature.startsWith(`${claim.feature}/`))) continue;
      participates = true;
      this.emitter.edge({type: 'participates-in', from: id, to: claim.feature, derived_by: 'annotation', confidence: 1,
        evidence: this.evidence(claim.feature, file, claim.rung), ...pick(claim.props, ['role'])});
    }
    if (participates) this.inventory.participating_files++;
  }

  /** Where the relation is visible: the home that claimed the file, or the file itself for a marker in it. */
  private evidence(feature: string, file: string, rung: number): string[] | undefined {
    const where = rung === 1 ? file : this.homeFile.get(feature);
    return where ? [where] : undefined;
  }

  /** The feature homes whose folder is the longest ancestor of [file]'s; more than one when they share it. */
  private nearestHomes(file: string): string[] {
    const dir = path.posix.dirname(file);
    let longest = -1;
    let winners: string[] = [];
    for (const [feature, home] of this.homeFile) {
      const folder = path.posix.dirname(home);
      if (dir !== folder && !dir.startsWith(`${folder}/`)) continue;
      if (folder.length > longest) {
        longest = folder.length;
        winners = [feature];
      }
      else if (folder.length === longest) winners.push(feature);
    }
    return winners;
  }

  /** Rung 4: the nearest ancestor folder that exactly one feature's `code:` root or home folder names. */
  private folderOwner(file: string): string | undefined {
    for (let dir = path.posix.dirname(file); dir && dir !== '.'; dir = path.posix.dirname(dir)) {
      const features = this.folders.get(dir);
      if (features?.size === 1) return [...features][0];
    }
    return undefined;
  }

  /** A help page a file names documents the feature that owns the file (conventions.md §4, the Dart pass). */
  private helpEdges(): void {
    for (const {file, page} of this.emitter.helpRefs) {
      const feature = this.owners.get(file);
      if (!feature) continue;
      this.emitter.edge({type: 'documents', from: docId(page), to: feature, derived_by: 'ast', confidence: 0.8, evidence: [file]});
    }
  }

  /** A page that cites a source file by its repo path mentions it; a cited path the graph holds no file for draws nothing. */
  private citationEdges(files: Map<string, Row>): void {
    for (const {page, file} of this.emitter.citations)
      if (files.has(file)) this.emitter.edge({type: 'mentions', from: docId(page), to: fileId(file), derived_by: 'annotation', confidence: 0.8, evidence: [page]});
  }

  private testEdges(): void {
    for (const test of this.emitter.rowsOf('test')) {
      const feature = this.owners.get(String(test.path));
      if (!feature) continue;
      this.emitter.edge({type: 'tests', from: String(test.id), to: feature, derived_by: 'filesystem', confidence: 0.9,
        evidence: [String(test.path)]});
    }
  }

  /** §8.1 name: a test whose file's base name, a segment of its category or, under Test Track, its folder spells a feature's
   * leaf id or name in the same words; a spelling several features share draws nothing and is reported once. */
  private nameEdges(): void {
    const byWords = new Map<string, string[]>();
    for (const [feature, name] of this.featureName)
      for (const key of new Set([nameWords(feature.slice(feature.lastIndexOf('/') + 1)), nameWords(name)]))
        if (key) byWords.set(key, [...byWords.get(key) ?? [], feature]);
    const reported = new Set<string>();
    for (const test of this.emitter.rowsOf('test')) {
      const file = String(test.path);
      const base = path.posix.basename(file).replace(/\..*$/, '');
      const spellings = new Set([nameWords(base, true), ...String(test.category ?? '').split(/[:|]/).map((s) => nameWords(s))]);
      if (file.startsWith(TEST_TRACK)) spellings.add(nameWords(path.posix.basename(path.posix.dirname(file))));
      for (const key of spellings) {
        const features = byWords.get(key);
        if (!features) continue;
        if (features.length > 1) {
          const detail = `${key}: ${[...features].sort().join(' and ')}`;
          if (!reported.has(detail)) this.emitter.problem('ambiguous_test_names', detail);
          reported.add(detail);
        }
        else if (!this.emitter.hasEdge('tests', String(test.id), features[0]))
          this.emitter.edge({type: 'tests', from: String(test.id), to: features[0], derived_by: 'name', confidence: NAME_CONFIDENCE, evidence: [file]});
      }
    }
  }

  /** §8.1 mirror: a test file named after a source file of its unit, by the whole base name or the part before a suffix,
   * when exactly one file of the unit bears that name (Dart under `lib/`, TypeScript outside the tests folders). */
  private mirrorEdges(files: Map<string, Row>): void {
    const tested = new Set(this.emitter.rowsOf('test').map((t) => String(t.path)));
    const byName = new Map<string, string[]>();
    for (const file of files.keys()) {
      const unit = unitOf(file);
      if (!unit || tested.has(file) || (file.endsWith('.dart') ? !file.startsWith(`${unit}/lib/`) : TS_TEST_PATH.test(file))) continue;
      const key = `${unit} ${path.posix.basename(file)}`;
      byName.set(key, [...byName.get(key) ?? [], file]);
    }
    for (const test of [...tested].sort()) {
      const unit = unitOf(test);
      const form = MIRROR_NAMES.map(([re, sep]) => ({m: re.exec(path.posix.basename(test)), sep})).find((f) => f.m);
      if (!unit || !form) continue;
      const parts = form.sep ? form.m![1].split(form.sep) : [form.m![1]];
      for (let n = parts.length; n >= 1; n--) {
        const hits = byName.get(`${unit} ${parts.slice(0, n).join(form.sep)}${path.posix.extname(test)}`);
        if (!hits) continue;
        if (hits.length === 1)
          this.emitter.edge({type: 'mirrors', from: fileId(test), to: fileId(hits[0]), derived_by: 'filesystem', confidence: n === parts.length ? 0.9 : 0.8, evidence: [test]});
        else this.emitter.problem('ambiguous_mirrors', `${test}: ${hits.join(' and ')} share the name`);
        break;
      }
    }
  }
}

/** The `code:` items of a home as repo paths, without the `path#Anchor` form, which names a declaration and claims nothing. */
function codeRoots(code: unknown): string[] {
  const out: string[] = [];
  for (const item of Array.isArray(code) ? code : []) {
    const target = typeof item === 'string' ? item : item && typeof item === 'object' ? (item as Record<string, unknown>).path : undefined;
    if (typeof target !== 'string' || target.includes('#')) continue;
    const repo = REPO_PREFIX.exec(target.trim());
    out.push(posix(repo ? `${repo[1]}/${repo[2]}` : target.trim()).replace(/\/+$/, ''));
  }
  return out;
}

/** `ScatterPlot`, `scatter-plot`, `scatter_plot tests` and `Scatter plot` are all `scatterplot`: lower-case words joined,
 * a trailing test word dropped from a file's base name. */
function nameWords(text: string, fileName = false): string {
  const words = text.replace(/([a-z0-9])([A-Z])/g, '$1 $2').toLowerCase().split(/[^a-z0-9]+/).filter(Boolean);
  if (fileName && words.length > 1 && TEST_WORDS.includes(words[words.length - 1])) words.pop();
  return words.join('');
}

/** The deepest of [features] when they form one part-of chain (`a`, `a/b`, `a/b/c`), otherwise nothing. */
function chainLeaf(features: string[]): string | undefined {
  const sorted = [...features].sort((a, b) => a.length - b.length);
  for (let i = 1; i < sorted.length; i++)
    if (!sorted[i].startsWith(`${sorted[i - 1]}/`)) return undefined;
  return sorted[sorted.length - 1];
}

function under(root: string): boolean {
  return CODE_ROOTS.some((r) => root.startsWith(r));
}

function pick(props: Record<string, unknown> | undefined, keys: string[]): Record<string, unknown> {
  return Object.fromEntries(keys.filter((k) => props?.[k] !== undefined).map((k) => [k, props![k]]));
}
