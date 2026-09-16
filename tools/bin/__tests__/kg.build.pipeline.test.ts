/// `grok kg build` end to end (build-plan.md WO-10): the whole pipeline over the mini monorepo under
/// fixtures/kg/build — every extractor, the row counts each of them produces, byte-identical repeats,
/// the Dart pass over the fixture's sources, the public projection and the five reports; then, where the
/// optional `kuzu` binding is installed, the index, `query`, the four operations and a CSV round-trip
/// of values the loader cannot put through a CSV.
/// Named `.pipeline.` and not `.integration.`, which vitest.config.mts reserves for the suites that need
/// a live Datagrok server: this one needs nothing but the fixture, and runs with `npm test`.
import {describe, it, expect, vi, beforeAll, afterAll} from 'vitest';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {fileURLToPath} from 'url';
import {spawnSync} from 'child_process';
import {loadTypeSystem, TypeSystem} from '../utils/kg/types';
import {loadKuzu, load as loadIndex, open, run} from '../utils/kg/kuzu';
import {find, explain, impact, testsFor, resolveTarget, printOps, OpsResult} from '../utils/kg/ops';
import {OutputFormat} from '../utils/server-output';
import {readGraph, makeReport, REPORT_NAMES, ReportName} from '../utils/kg/report';
import {currentDir} from '../utils/kg/build/write';
import {kg} from '../commands/kg';

const fixture = path.join(path.dirname(fileURLToPath(import.meta.url)), 'fixtures', 'kg', 'build');
const KG_DIR = path.join('core', 'docs', 'knowledge-graph');
const CACHING = 'core/docs/CACHING.md';
const BIO_HOME = 'public/help/domains/bio/bio.md';
const SEQUENCES = 'public/help/domains/bio/sequences.md';
const RENDERER = 'core/client/d4/lib/src/legends/legend_renderer.dart';
const LIMIT = {limit: 50};

interface Built {
  repo: string;
  out: string;
  manifest: any;
  rows: (file: string) => any[];
  report: (name: string) => any;
}

function git(repo: string, ...args: string[]): string {
  return spawnSync('git', ['-c', 'user.email=kg@test', '-c', 'user.name=kg', '-C', repo, ...args], {encoding: 'utf8'}).stdout.trim();
}

/** The fixture as a git repository with one commit. */
function makeRepo(): string {
  const repo = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-kg-int-'));
  fs.cpSync(fixture, repo, {recursive: true});
  git(repo, 'init', '-q');
  git(repo, 'add', '-A');
  git(repo, 'commit', '-q', '-m', 'fixture');
  stampRelease(repo);
  return repo;
}

/** The same fixture with `public/` as a repository of its own, which the monorepo commit records as a gitlink. */
function makeSubmoduleRepo(): string {
  const repo = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-kg-sub-'));
  fs.cpSync(fixture, repo, {recursive: true});
  const publicDir = path.join(repo, 'public');
  git(publicDir, 'init', '-q');
  git(publicDir, 'add', '-A');
  git(publicDir, 'commit', '-q', '-m', 'public');
  git(repo, 'init', '-q');
  git(repo, 'add', '-A');
  git(repo, 'commit', '-q', '-m', 'fixture');
  stampRelease(repo);
  return repo;
}

/** The picked commits of the release record are the fixture's own commit: git has to find them for a pick to exist. */
function stampRelease(repo: string): void {
  const file = path.join(repo, 'core', 'docs', 'release', '1.0.1.yaml');
  const sha = git(repo, 'rev-parse', '--short=10', 'HEAD');
  fs.writeFileSync(file, fs.readFileSync(file, 'utf8').replace(/\b(?:aaaaaaaaaa|bbbbbbbbbb)\b/g, sha));
}

async function build(repo: string, extra: Record<string, unknown> = {}): Promise<Built> {
  const log = vi.spyOn(console, 'log').mockImplementation(() => {});
  const error = vi.spyOn(console, 'error').mockImplementation(() => {});
  const before = process.exitCode;
  try {
    await kg({_: ['kg', 'build'], kg: path.join(repo, KG_DIR), backlog: path.join(repo, 'backlog'), db: false, output: 'json', ...extra});
    expect(error.mock.calls).toEqual([]);
    const out = currentDir(path.join(repo, ...(extra.public ? ['public', '.kg'] : ['.kg'])))!;
    const rows = (file: string) => {
      const p = path.join(out, file.startsWith('reports/') ? file : `data/${file}.jsonl`);
      return fs.existsSync(p) ? fs.readFileSync(p, 'utf8').split('\n').filter(Boolean).map((l) => JSON.parse(l)) : [];
    };
    const report = (name: string) => JSON.parse(fs.readFileSync(path.join(out, 'reports', `${name}.json`), 'utf8'));
    return {repo, out, manifest: JSON.parse(String(log.mock.calls[0][0])), rows, report};
  }
  finally {
    process.exitCode = before;
    log.mockRestore();
    error.mockRestore();
  }
}

/** Every JSONL file under `data/`, by name, so two builds can be compared byte for byte. */
function dataFiles(out: string): Record<string, string> {
  const files: Record<string, string> = {};
  for (const kind of ['nodes', 'edges'])
    for (const name of fs.readdirSync(path.join(out, 'data', kind)).sort())
      files[`${kind}/${name}`] = fs.readFileSync(path.join(out, 'data', kind, name), 'utf8');
  return files;
}

const built = makeRepo();
const graph = build(built);

describe('grok kg build over the fixture monorepo (build-plan.md WO-10)', () => {
  it('runs every extractor and writes the rows the fixture calls for', async () => {
    const {manifest} = await graph;
    expect(manifest.mode).toBe('full');
    expect(manifest.sources).toEqual({
      backlog: 'ok@2026-01-12T07:00:00Z', dart: 'ok', docs: 'partial', git: 'ok', homes: 'ok', membership: 'ok',
      people: 'partial', process: 'ok', releases: 'ok', 'ts-changelog': 'partial', 'ts-declarations': 'ok',
      'ts-functions': 'partial(2 rejected)', 'ts-imports': 'ok', 'ts-markers': 'ok', 'ts-packages': 'ok', 'ts-samples': 'partial',
      'ts-tests': 'ok', 'ts-uses': 'ok',
    });
    expect(manifest.counts.nodes).toEqual({
      app: 2, 'cell-renderer': 1, 'changelog-entry': 5, commit: 2, connection: 2, container: 3, customer: 2,
      declaration: 59, 'doc-anchor': 21, 'doc-page': 14, editor: 1, feature: 9, 'file-handler': 1,
      'file-viewer': 1, filter: 1, function: 13, library: 2, 'lifecycle-hook': 2, package: 7, panel: 2, person: 3,
      query: 3, release: 2, sample: 4, scenario: 7, script: 3, 'script-environment': 1, 'script-handler': 1,
      'sem-type-detector': 6, 'semantic-type': 7, 'source-file': 44, test: 16, 'test-suite': 7, ticket: 13, tutorial: 1, viewer: 2,
    });
    expect(manifest.counts.edges).toEqual({
      affects: 2, assignee: 3, automates: 3, base: 1, calls: 6, changes: 2, connection: 3, covers: 1,
      declares: 170, demonstrates: 1, 'depends-on': 6, documents: 5, environment: 1, extends: 4,
      implements: 1, imports: 26, includes: 2, 'is-implemented-in': 19, mentions: 16, owner: 9,
      package: 72, page: 21, 'part-of': 11, 'participates-in': 6, reporter: 4, 'requested-by': 2, resolves: 1,
      suite: 16, 'targets-release': 5, 'targets-semtype': 13, tests: 13, 'tracked-in': 2, user_help: 1, uses: 18,
    });
    expect(manifest.problems).toMatchObject({dangling_edges: 0, ambiguous_owners: 1, orphans: 28, partial_stubs: 23});
  }, 120_000);

  it('writes the same bytes twice, with the same content-addressed batch and a later built_at', async () => {
    const {out, manifest} = await graph;
    const first = dataFiles(out);
    const again = await build(built);
    expect(again.out).not.toBe(out);
    expect(dataFiles(again.out)).toEqual(first);
    expect(again.manifest.batch).toBe(manifest.batch);
    expect(again.manifest.built_at >= manifest.built_at).toBe(true);
    expect(Object.values(first).every((text) => !text.includes('built_at'))).toBe(true);
  }, 120_000);

  it('reads the Dart sources of the fixture: files, declarations, tests and their suite', async () => {
    const {manifest, rows} = await graph;
    expect(manifest.sources.dart).toBe('ok');
    expect(manifest.dart_depth).toBe('lexical');
    expect(manifest.dart_packages).toEqual({d4: 8, grok_shared: 1});
    const files = rows('nodes/source-file').filter((f) => f.language === 'dart');
    expect(files).toHaveLength(9);
    expect(files.filter((f) => f.generated).map((f) => f.path)).toEqual(['core/client/d4/lib/src/viewers/viewer.g.dart']);
    expect(rows('nodes/declaration').filter((d) => d.language === 'dart').map((d) => d.name).sort())
      .toEqual(['HelpUrl', 'Histogram', 'Legend', 'LegendCache', 'LegendRenderer', 'ScatterPlot', 'Viewer', 'ViewerProps']);
    expect(rows('nodes/test').filter((t) => t.framework === 'dart').map((t) => t.category)).toEqual([undefined, 'placement']);
    expect(rows('nodes/test-suite').filter((s) => s.framework === 'dart')).toHaveLength(1);
  }, 120_000);

  it('projects the public layer: public nodes only, no home, no owner, no source files, no reports', async () => {
    const {manifest, rows, out} = await build(makeRepo(), {public: true});
    expect(manifest.mode).toBe('public');
    expect(Object.keys(manifest.revisions)).toEqual(['public']);
    expect(Object.keys(manifest.counts.nodes).sort()).toEqual(['doc-anchor', 'doc-page', 'feature', 'library', 'package', 'sample', 'scenario']);
    const features = rows('nodes/feature');
    expect(features).toHaveLength(9);
    expect(features.every((f) => f.visibility === 'public' && f.home === undefined && f.owner === undefined)).toBe(true);
    expect(features.find((f) => f.id === 'domains/bio').description).toContain('Sequence analysis');
    expect(rows('nodes/source-file')).toEqual([]);
    expect(fs.existsSync(path.join(out, 'data', 'nodes', 'person.jsonl'))).toBe(false);
    expect(fs.existsSync(path.join(out, 'reports'))).toBe(false);
    // an edge survives only when both its ends did
    expect(rows('edges/owner')).toEqual([]);
    expect(rows('edges/covers')).toMatchObject([{from: 'TS:viewers/scatter-plot/ui', to: 'domains/bio'}]);
  }, 120_000);
});

describe('the reports build writes and the report verb prints (build-plan.md WO-8, WO-10)', () => {
  it('writes json and md for the four reports a build can answer on its own', async () => {
    const {out} = await graph;
    for (const name of REPORT_NAMES.filter((n) => n !== 'diff'))
      for (const ext of ['json', 'md']) expect(fs.existsSync(path.join(out, 'reports', `${name}.${ext}`)), `${name}.${ext}`).toBe(true);
    expect(fs.existsSync(path.join(out, 'reports', 'diff.json'))).toBe(false);
    expect(fs.readFileSync(path.join(out, 'reports', 'coverage.md'), 'utf8')).toContain('| visualize/viewers/scatter-plot |');
  });

  it('groups orphan files by package and core sub-project, largest first', async () => {
    const {report} = await graph;
    const orphans = report('orphans');
    expect(orphans.summary).toBe('28 files of 44 observed (16 owned, 5 participating) have no owner, 405 of 716 lines, in 8 groups; the 8 largest groups below.');
    expect(orphans.sections[0].rows[0]).toEqual({group: 'public/js-api', owner: '', files: 11, owned: 0, participating: 0, orphans: 11, loc: 155,
      largest: expect.stringContaining('src/dataframe.ts (45)')});
    expect(orphans.sections[0].rows.map((r: any) => r.group)).toContain('core/client/d4');
    expect(orphans.sections[0].rows.map((r: any) => r.loc)).toEqual([...orphans.sections[0].rows.map((r: any) => r.loc)].sort((a: number, b: number) => b - a));
  });

  it('names the citation targets, help pages and specs the graph can no longer reach', async () => {
    const {report} = await graph;
    const rows = report('stale').sections[0].rows;
    expect(rows).toContainEqual({kind: 'help-url', source: 'public/packages/ApiSamples/scripts/misc/missing.js',
      target: 'https://datagrok.ai/help/nowhere/at-all', reason: 'names no page under public/help'});
    expect(rows).toContainEqual({kind: 'realized-as', source: 'public/packages/UsageAnalysis/files/TestTrack/Viewers/ScatterPlot/scatterplot-legend.md',
      target: 'missing-spec.ts', reason: 'is not beside the scenario'});
    // the batch is ok here, so defined_by is checked and the note about it is absent
    expect(report('stale').notes).toEqual([]);
  });

  it('takes the citations check reports as stale references, with the code and the path that is gone', async () => {
    const repo = makeRepo();
    fs.rmSync(path.join(repo, 'core', 'client', 'd4', 'lib', 'src', 'legends', 'legend.dart'));
    const {report} = await build(repo);
    expect(report('stale').sections[0].rows).toContainEqual({kind: 'citation', source: 'core/client/d4/lib/src/legends/README.md:9',
      target: 'core/client/d4/lib/src/legends/legend.dart', reason: "missing-cited-path: cited path 'core/client/d4/lib/src/legends/legend.dart' does not exist"});
  }, 120_000);

  it('gives every feature a coverage row, counting its tests and the scenarios that cover it', async () => {
    const {report} = await graph;
    const rows = report('coverage').sections[0].rows;
    expect(rows.map((r: any) => r.feature)).toEqual(['domains', 'domains/bio', 'platform', 'platform/caching', 'visualize',
      'visualize/legends', 'visualize/viewers', 'visualize/viewers/histogram', 'visualize/viewers/scatter-plot']);
    expect(rows.find((r: any) => r.feature === 'domains/bio')).toEqual({feature: 'domains/bio', owner: 'P:jane', status: 'active', stub: false,
      tests_runnable: 4, tests_skipped: 1, tests_dynamic: 0, inherited: 0, scenarios: 1, scenario_automated: 0,
      user_help: '', developer_help: '', no_tests: false, no_docs: false, no_description: false});
    // the page named by user_help alone, with no documents edge, is documentation all the same
    expect(rows.find((r: any) => r.feature === 'visualize/viewers/histogram'))
      .toMatchObject({user_help: 'doc:public/help/visualize/viewers/histogram.md', developer_help: '', no_docs: false});
    expect(rows.find((r: any) => r.feature === 'platform/caching')).toMatchObject({no_tests: true, no_docs: true, no_description: false});
    expect(rows.find((r: any) => r.feature === 'domains')).toMatchObject({owner: '', status: 'proposed', no_description: true});
  });

  it('proposes a feature id for every unowned package folder, and none for a library', async () => {
    const {report} = await graph;
    const rows = report('proposed').sections[0].rows;
    expect(rows).toEqual([
      {path: 'public/packages/Demo', files: 12, loc: 370, unowned_files: 7, unowned_loc: 153, suggested: 'domains/demo'},
      {path: 'public/packages/Plain', files: 3, loc: 35, unowned_files: 3, unowned_loc: 35, suggested: 'domains/plain'},
      {path: 'public/packages/Tutorials', files: 3, loc: 29, unowned_files: 3, unowned_loc: 29, suggested: 'domains/tutorials'},
      {path: 'public/libraries/utils', files: 1, loc: 16, unowned_files: 1, unowned_loc: 16, suggested: ''},
      {path: 'public/packages/ApiTests', files: 1, loc: 7, unowned_files: 1, unowned_loc: 7, suggested: 'domains/api-tests'},
    ]);
    // the Tested package is owned file for file through the viewers home, and only that drops out
    expect(rows.some((r: any) => r.path === 'public/packages/Tested')).toBe(false);
  });

  it('names the features a branch touches through the home it changed, with the tests that cover them', async () => {
    const repo = makeRepo();
    const base = git(repo, 'rev-parse', 'HEAD');
    fs.appendFileSync(path.join(repo, ...CACHING.split('/')), '\nOne more line.\n');
    // the batch envelope was stamped after the first commit, so only the home goes into the second
    git(repo, 'commit', '-q', '-m', 'touch caching', '--', CACHING);
    const {out} = await build(repo);
    const system: TypeSystem = loadTypeSystem(path.join(repo, KG_DIR));
    const data = readGraph(out, repo, system, 'diff');
    const report = makeReport('diff', data, {system, repoRoot: repo, base});
    expect(report.title).toBe(`Features touched by ${base}...HEAD`);
    expect(report.sections[0].rows).toEqual([{feature: 'platform/caching', name: 'Caching', owner: 'P:jane', relation: 'home', confidence: 1}]);
    expect(report.sections[1].rows).toEqual([{file: 'core/docs/CACHING.md', feature: 'platform/caching', relation: 'home', confidence: 1, change: 'modified'}]);
    expect(report.summary).toContain('1 files changed (0 deleted), 0 of them in no feature');
  }, 120_000);

  it('separates owned, participating and orphan files per group, over the inventory it observed', async () => {
    const {report} = await graph;
    const orphans = report('orphans');
    // one owned file used to make the whole Demo package look accounted for; it is 5 owned, 1 participating, 7 orphans
    // and, whatever the features say, the package.json author owns all twelve
    expect(orphans.sections[0].rows).toContainEqual({group: 'public/packages/Demo', owner: 'P:jane', files: 12, owned: 5, participating: 1, orphans: 7,
      loc: 153, largest: expect.stringContaining('detectors.js (69)')});
    expect(orphans.notes).toContain('the inventory covers the files the extractors observed, not every file in the repositories');
  });

  it('reports a tracked-in ticket as unknown, not as absent, when the backlog snapshot was not read', async () => {
    const repo = makeRepo();
    const {report} = await build(repo, {backlog: path.join(repo, 'no-such-backlog')});
    const rows = report('stale').sections[0].rows.filter((r: any) => r.kind === 'ticket');
    expect(rows.length).toBeGreaterThan(0);
    expect(rows.every((r: any) => r.reason === 'unknown: the backlog snapshot is missing')).toBe(true);
    expect(report('stale').notes).toContain('the backlog snapshot was not read, so a tracked-in ticket is reported as unknown, not as absent');
  }, 120_000);

  it('splits runnable, skipped and dynamic tests, counts what a feature inherits, and marks a stub', async () => {
    const {report} = await graph;
    const rows = report('coverage').sections[0].rows;
    // 6 tests carry ~visualize/viewers: four that run, one skipped, one whose name is a template
    expect(rows.find((r: any) => r.feature === 'visualize/viewers')).toMatchObject({stub: false,
      tests_runnable: 4, tests_skipped: 1, tests_dynamic: 1, inherited: 0});
    // a stub has no home of its own; what it "has" is what the features under it have
    expect(rows.find((r: any) => r.feature === 'visualize')).toMatchObject({stub: true, tests_runnable: 0, inherited: 8});
    expect(rows.find((r: any) => r.feature === 'domains')).toMatchObject({stub: true, inherited: 5});
    expect(report('coverage').summary).toBe('9 features, 3 of them stubs; of the 6 with a home, 3 have no test or scenario, ' +
      '3 no document beside the home, 0 no first paragraph.');
  });

  it('ranks a folder by the code no feature owns, not by whether anything in it is owned', async () => {
    const {report} = await graph;
    const rows = report('proposed').sections[0].rows;
    // the Demo package has an owned file and 7 unowned ones; the unowned remainder is what puts it first
    expect(rows[0]).toEqual({path: 'public/packages/Demo', files: 12, loc: 370, unowned_files: 7, unowned_loc: 153, suggested: 'domains/demo'});
    expect(report('proposed').summary).toContain('1 of them is partly owned already');
  });

  it('keeps a deleted file with the owner the graph still has, and warns when the graph is behind the tree', async () => {
    const repo = makeRepo();
    const {out} = await build(repo);
    const base = git(repo, 'rev-parse', 'HEAD');
    fs.rmSync(path.join(repo, ...CACHING.split('/')));
    fs.rmSync(path.join(repo, 'public', 'packages', 'Plain', 'src', 'package.js'));
    git(repo, 'commit', '-q', '-m', 'remove', '--', CACHING, 'public/packages/Plain/src/package.js');
    const system: TypeSystem = loadTypeSystem(path.join(repo, KG_DIR));
    const report = makeReport('diff', readGraph(out, repo, system, 'diff'), {system, repoRoot: repo, base});
    expect(report.summary).toContain('2 files changed (2 deleted)');
    // the file is gone from the tree, so the graph built before the deletion is the only place its owner survives
    expect(report.sections[1].rows).toContainEqual({file: CACHING, feature: 'platform/caching', relation: 'home', confidence: 1, change: 'deleted'});
    expect(report.sections[1].rows).toContainEqual({file: 'public/packages/Plain/src/package.js', feature: '', relation: 'deleted', confidence: '', change: 'deleted'});
    expect(report.summary).toContain('0 of them in no feature');
    expect(report.notes).toContain('1 deleted file had no owner in this graph either; they are listed as deleted');
    expect(report.notes).toContainEqual(expect.stringContaining('the graph was built at'));
  }, 120_000);

  it('takes the public baseline from the gitlink the core base recorded', async () => {
    const repo = makeSubmoduleRepo();
    const publicDir = path.join(repo, 'public');
    const base = git(repo, 'rev-parse', 'HEAD');
    fs.appendFileSync(path.join(publicDir, 'packages', 'Plain', 'src', 'package.ts'), '\n// one more line\n');
    git(publicDir, 'commit', '-q', '-m', 'touch plain', '--', 'packages/Plain/src/package.ts');
    git(repo, 'commit', '-q', '-m', 'bump submodule', '--', 'public');
    const {out} = await build(repo);
    const system: TypeSystem = loadTypeSystem(path.join(repo, KG_DIR));
    const report = makeReport('diff', readGraph(out, repo, system, 'diff'), {system, repoRoot: repo, base});
    // the monorepo sha names nothing inside the submodule: without the gitlink this diff failed and found no public file
    expect(report.notes.some((n: string) => n.startsWith('public/:'))).toBe(false);
    expect(report.sections[1].rows).toContainEqual({file: 'public/packages/Plain/src/package.ts', feature: '', relation: '', confidence: '', change: 'modified'});
  }, 120_000);

  it('refuses a report it does not have, and diff without a base', async () => {
    const {repo} = await graph;
    const error = vi.spyOn(console, 'error').mockImplementation(() => {});
    const before = process.exitCode;
    try {
      await kg({_: ['kg', 'report', 'nonsense'], kg: path.join(repo, KG_DIR)});
      await kg({_: ['kg', 'report', 'diff'], kg: path.join(repo, KG_DIR)});
      await kg({_: ['kg', 'report', 'coverage'], kg: path.join(repo, KG_DIR), output: 'csv'});
      expect(error.mock.calls.map((c) => String(c[0]))).toEqual([
        `unknown report 'nonsense': ${REPORT_NAMES.join(', ')}`,
        'grok kg report diff needs the revision to compare against: --diff <ref>',
        "--output must be table, json or md, got 'csv'",
      ]);
    }
    finally {
      process.exitCode = before;
      error.mockRestore();
    }
  });
});

/** What an op writes to the console, line by line. */
function render(result: OpsResult, output: OutputFormat = 'table'): string[] {
  const log = vi.spyOn(console, 'log').mockImplementation(() => {});
  try {
    printOps(result, output);
    return log.mock.calls.map((c) => String(c[0]));
  }
  finally {
    log.mockRestore();
  }
}

function csv(result: OpsResult): string {
  return render(result, 'csv').join('\n');
}

describe('what an op prints (slice-results.md §10)', () => {
  it('prints a list cell as its items and says how many it left out, and marks a cut object', () => {
    const long = 'x'.repeat(100);
    const lines = render({op: 'explain', target: {id: 'pkg:Chem'}, sections: [{title: 'properties', total: 2, rows: [
      {property: 'settings', value: ['Sketcher', 'TemplatesPath', 'BuildingBlocksPath', 'ReagentsPath', 'MolecularFingerprints', 'SubstructureSearch']},
      {property: 'meta', value: {a: long}},
    ]}]}).join('\n');
    expect(lines).toContain('Sketcher, TemplatesPath, BuildingBlocksPath, ReagentsPath, MolecularFingerprints … (+1 more)');
    expect(lines).toContain(`{"a":"${'x'.repeat(73)}…`);
    expect(lines).not.toContain(long);
  });

  it('gives the edges section a sub-header per group in a table, and a group column in json and csv', () => {
    const edges = {title: 'edges', total: 3, rows: [
      {group: 'ownership', edge: 'IS_IMPLEMENTED_IN', direction: 'out', count: 2},
      {group: 'evidence', edge: 'TESTS', direction: 'in', count: 4},
      {group: 'reference', edge: 'owner', direction: 'out', count: 1},
    ]};
    const lines = render({op: 'explain', target: {id: 'visualize/viewers'}, sections: [edges]});
    expect(lines.filter((l) => l.startsWith('\n  '))).toEqual(['\n  ownership', '\n  evidence', '\n  reference']);
    expect(lines.join('\n')).not.toContain('group');
    expect(csv({op: 'explain', target: {id: 'visualize/viewers'}, sections: [edges]})).toContain('group,edge,direction,count');
  });
});

/** The fixture graph with one feature carrying values a CSV cannot express (build-plan.md WO-7 "List cells"). */
async function unloadable(): Promise<{out: string, feature: Record<string, unknown>}> {
  const {out} = await graph;
  const file = path.join(out, 'data', 'nodes', 'feature.jsonl');
  const rows = fs.readFileSync(file, 'utf8').split('\n').filter(Boolean).map((l) => JSON.parse(l));
  const feature = rows.find((r) => r.id === 'platform/caching')!;
  feature.aliases = ['cache, the', 'a "quoted" one'];
  feature.description = 'One line, with a comma.\nAnd a "second" one.';
  fs.writeFileSync(file, `${rows.map((r) => JSON.stringify(r)).join('\n')}\n`);
  return {out, feature};
}

const withKuzu = loadKuzu() ? it : it.skip;

describe('the index over the fixture graph (build-plan.md WO-7, WO-10)', () => {
  /** One index for the whole block: a kuzu database reserves its buffer pool, and two at once exhaust a worker. */
  let index: {out: string, feature: Record<string, unknown>, system: TypeSystem, opened: Awaited<ReturnType<typeof open>>};

  beforeAll(async () => {
    if (!loadKuzu()) return;
    const {out, feature} = await unloadable();
    const system = loadTypeSystem(path.join((await graph).repo, KG_DIR));
    await loadIndex(out, system);
    index = {out, feature, system, opened: await open(out, true)};
  }, 300_000);

  afterAll(async () => {
    if (!index?.opened) return;
    await index.opened.conn.close();
    await index.opened.db.close();
  });

  withKuzu('round-trips a list item with a comma and a quote, and a description with a newline', async () => {
    const {conn} = index.opened!;
    const {rows} = await run(conn, 'MATCH (n:Feature) WHERE n.`id` = $id RETURN n.`aliases` AS aliases, n.`description` AS description', {id: index.feature.id});
    expect(rows[0].aliases).toEqual(index.feature.aliases);
    expect(rows[0].description).toBe(index.feature.description);
    expect((await run(conn, 'MATCH (n:Feature) RETURN count(n) AS n')).rows[0].n).toBe(9);
    expect(fs.existsSync(path.join(index.out, 'tmp'))).toBe(false);
  }, 120_000);

  withKuzu('answers find, explain, impact and tests-for over the whole fixture', async () => {
    const {conn} = index.opened!;
    expect((await find(conn, 'caching', LIMIT)).sections[0].rows[0]).toMatchObject({id: 'platform/caching', type: 'feature'});
    expect((await find(conn, 'cache, the', LIMIT)).sections[0].rows[0]).toMatchObject({id: 'platform/caching'});

    const bio = (await resolveTarget(conn, '~domains/bio'))!;
    const explained = await explain(conn, bio, LIMIT);
    expect(explained.sections[0].rows).toContainEqual({property: 'home', value: 'public/help/domains/bio/bio.md'});
    expect(explained.sections[1].rows).toContainEqual(expect.objectContaining({edge: 'owner', direction: 'out', targets: 'P:jane'}));

    const tests = await testsFor(conn, bio, LIMIT);
    expect(tests.sections.find((s) => s.title === 'tests')!.rows.length).toBe(5);
    expect(tests.sections.find((s) => s.title === 'scenarios')!.rows).toMatchObject([{scenario: 'TS:viewers/scatter-plot/ui', manual_only: true}]);

    const renderer = (await resolveTarget(conn, RENDERER))!;
    expect(renderer).toMatchObject({id: `file:${RENDERER}`, root: 'Component'});
    const reached = await impact(conn, renderer, LIMIT);
    expect(reached.sections[0].rows).toContainEqual(expect.objectContaining({feature: 'visualize/viewers/scatter-plot', relation: 'owns'}));
    expect(reached.sections[1].rows).toContainEqual({feature: 'visualize/viewers/scatter-plot', owner: 'P:jane', name: 'Jane Dev'});
  }, 120_000);

  withKuzu('counts a section whole and pages it after, so a header tells a page from the total', async () => {
    const {conn} = index.opened!;
    const bio = (await resolveTarget(conn, '~domains/bio'))!;
    const paged = (await testsFor(conn, bio, {limit: 3})).sections.find((s) => s.title === 'tests')!;
    expect(paged.rows.length).toBe(3);
    expect(paged.total).toBe(5);
    expect(render({op: 'tests-for', target: bio, sections: [paged]})).toContain('\ntests (3 of 5; --limit to see more)');
    const whole = (await testsFor(conn, bio, LIMIT)).sections.find((s) => s.title === 'tests')!;
    expect(render({op: 'tests-for', target: bio, sections: [whole]})).toContain('\ntests (5)');
  }, 120_000);

  withKuzu('answers for a home document and for a page that documents a feature, which have no file: node of their own', async () => {
    const {conn} = index.opened!;
    const home = (await resolveTarget(conn, BIO_HOME))!;
    expect(home).toMatchObject({id: `doc:${BIO_HOME}`, root: 'Artifact', type: 'doc-page'});
    expect((await impact(conn, home, LIMIT)).sections[0].rows).toMatchObject([{feature: 'domains/bio', relation: 'home', name: 'Bioinformatics', status: 'active',
      via: `doc:${BIO_HOME} → home ← domains/bio`}]);
    expect((await testsFor(conn, home, LIMIT)).sections.find((s) => s.title === 'tests')!.rows.length).toBe(5);

    const page = (await resolveTarget(conn, SEQUENCES))!;
    expect((await impact(conn, page, LIMIT)).sections[0].rows).toMatchObject([
      {feature: 'domains/bio', relation: 'documents', name: 'Bioinformatics', status: 'active'},
      // a doc comment in viewer.dart names the page too, so it documents the feature that owns that file
      {feature: 'visualize/viewers', relation: 'documents', name: 'Viewers'},
    ]);
  }, 120_000);

  withKuzu('answers for a package through what it declares, and for a file through what it declares', async () => {
    const {conn} = index.opened!;
    // the package holds the code the viewers home claims, so it must answer with the same tests the feature does
    const pkg = (await resolveTarget(conn, 'pkg:Tested'))!;
    const viewers = (await resolveTarget(conn, '~visualize/viewers'))!;
    const byPackage = await testsFor(conn, pkg, LIMIT);
    const byFeature = await testsFor(conn, viewers, LIMIT);
    const total = (r: OpsResult, title: string) => r.sections.find((s) => s.title === title)!.total;
    const features = (r: OpsResult) => r.sections[0].rows.map((x) => String(x.feature)).sort();
    expect(features(byPackage)).toEqual(features(byFeature));
    expect(features(byPackage)).toEqual(['visualize/viewers', 'visualize/viewers/histogram', 'visualize/viewers/scatter-plot']);
    expect(total(byPackage, 'tests')).toBe(total(byFeature, 'tests'));
    expect(total(byPackage, 'tests')).toBe(6);
    expect(byPackage.sections[0].rows.find((r) => r.feature === 'visualize/viewers')!.via)
      .toBe('pkg:Tested → declares → file:public/packages/Tested/src/tests/demo-tests.ts → is-implemented-in ← visualize/viewers');

    // js-api/src/viewer.ts belongs to no feature; the declaration inside it is the scatter plot's api
    const file = (await resolveTarget(conn, 'public/js-api/src/viewer.ts'))!;
    const reached = await impact(conn, file, LIMIT);
    expect(reached.sections[0].rows).toMatchObject([{feature: 'visualize/viewers/scatter-plot', relation: 'owns',
      via: 'file:public/js-api/src/viewer.ts → declares → decl:public/js-api/src/viewer.ts#JsViewer → is-implemented-in ← visualize/viewers/scatter-plot'}]);
    expect(reached.sections[0].rows[0].path).toEqual(['file:public/js-api/src/viewer.ts', '→ declares →',
      'decl:public/js-api/src/viewer.ts#JsViewer', '→ is-implemented-in ←', 'visualize/viewers/scatter-plot']);
  }, 120_000);

  withKuzu('gives a file no feature owns the owner of the package that declares it', async () => {
    const {conn} = index.opened!;
    const file = (await resolveTarget(conn, 'public/packages/Demo/detectors.js'))!;
    const reached = await impact(conn, file, LIMIT);
    expect(reached.sections[0].rows).toEqual([]);
    expect(reached.sections[0].empty).toContain('no home document owns this file');
    expect(reached.sections[1].rows).toEqual([{package: 'pkg:Demo', owner: 'P:jane', name: 'Jane Dev',
      via: 'file:public/packages/Demo/detectors.js → declares ← pkg:Demo → owner → P:jane'}]);
    expect((await testsFor(conn, file, LIMIT)).sections[1].rows).toEqual(reached.sections[1].rows);
  }, 120_000);

  withKuzu('puts the evidence path behind an edge group', async () => {
    const {conn} = index.opened!;
    const viewers = (await resolveTarget(conn, '~visualize/viewers'))!;
    const edges = (await explain(conn, viewers, LIMIT)).sections.find((s) => s.title === 'edges')!;
    expect(edges.rows).toContainEqual(expect.objectContaining({edge: 'IS_IMPLEMENTED_IN', direction: 'out',
      derived_by: 'annotation', confidence: '1', evidence: 'core/docs/VIEWERS.md'}));
  }, 120_000);

  withKuzu('separates what a release targeted, what it includes and what it shipped', async () => {
    const {conn} = index.opened!;
    const release = (await resolveTarget(conn, '~Rel:1.0.1'))!;
    const sections = (await explain(conn, release, LIMIT)).sections;
    const titled = (title: string) => sections.find((s) => s.title === title)!;
    // the record's own `features:` and the backlog's fix versions are what it aimed at; the picks are what it carries
    expect(titled('targeted').rows.map((r) => r.ticket)).toEqual(['GROK-100', 'GROK-101', 'GROK-102']);
    expect(titled('targeted').rows.find((r) => r.ticket === 'GROK-101')!.confidence).toBe(0.7);
    expect(titled('included').rows.filter((r) => r.relation === 'includes').map((r) => String(r.item).split(':')[1]).sort()).toEqual(['public', 'reddata']);
    expect(titled('included').rows.filter((r) => r.relation === 'picked').map((r) => r.item)).toEqual(['GROK-100', 'GROK-101']);
    // a dry run in testing cannot have shipped anything, and the record says so instead of showing an empty list
    expect(titled('shipped').rows).toEqual([]);
    expect(titled('shipped').empty).toBe('shipped (0): not derivable (the record is a dry run; release state testing)');
  }, 120_000);

  withKuzu('keeps an exact match the substring scan would have cut, and still ranks the vocabulary first', async () => {
    const {conn} = index.opened!;
    // one row per table is all the substring scan may contribute, and `demo` matches a dozen components
    const rows = (await find(conn, 'demo', {limit: 50, scan: 1})).sections[0].rows;
    expect(rows.find((r) => r.id === 'pkg:Demo')).toMatchObject({type: 'package', match: 0});
    expect(rows.findIndex((r) => r.id === 'pkg:Demo')).toBeLessThan(rows.findIndex((r) => r.match === 2));
    // exactness is what the exact query matched, not what the id and the name happen to say: this is an alias
    const alias = await find(conn, 'cache, the', LIMIT);
    expect(alias.sections[0].rows[0]).toMatchObject({id: 'platform/caching', match: 0});
  }, 120_000);

  withKuzu('leaves the reports to the JSONL: every one of them answers with the index open', async () => {
    const {repo} = await graph;
    const system = index.system;
    for (const name of REPORT_NAMES) {
      const data = readGraph(index.out, repo, system, name as ReportName);
      const report = makeReport(name as ReportName, data, {system, repoRoot: repo, base: 'HEAD'});
      expect(report.name, name).toBe(name);
      expect(report.sections[0].title, name).toBeTruthy();
    }
  }, 120_000);
});
