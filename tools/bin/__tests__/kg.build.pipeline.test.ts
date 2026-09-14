/// `grok kg build` end to end (build-plan.md WO-10): the whole pipeline over the mini monorepo under
/// fixtures/kg/build — every extractor, the row counts each of them produces, byte-identical repeats,
/// the three states of the Dart batch, the public projection and the five reports; then, where the
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
import {find, explain, impact, testsFor, resolveTarget} from '../utils/kg/ops';
import {readGraph, makeReport, REPORT_NAMES, ReportName} from '../utils/kg/report';
import {kg} from '../commands/kg';

const fixture = path.join(path.dirname(fileURLToPath(import.meta.url)), 'fixtures', 'kg', 'build');
const KG_DIR = path.join('core', 'docs', 'knowledge-graph');
const BATCH = path.join('.kg', 'batches', 'kg-dart.jsonl');
const CACHING = 'core/docs/CACHING.md';
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

/** The fixture as a git repository with one commit, its Dart batch stamped with that revision. */
function makeRepo(): string {
  const repo = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-kg-int-'));
  fs.cpSync(fixture, repo, {recursive: true});
  git(repo, 'init', '-q');
  git(repo, 'add', '-A');
  git(repo, 'commit', '-q', '-m', 'fixture');
  stampBatch(repo, {built_at: new Date().toISOString(), revision: git(repo, 'rev-parse', 'HEAD')});
  return repo;
}

/** Rewrites the batch envelope; `null` deletes the file. */
function stampBatch(repo: string, envelope: Record<string, unknown> | null): void {
  const file = path.join(repo, BATCH);
  if (!envelope) {
    fs.rmSync(file);
    return;
  }
  const lines = fs.readFileSync(file, 'utf8').split('\n').filter(Boolean);
  lines[0] = JSON.stringify({...JSON.parse(lines[0]), ...envelope});
  fs.writeFileSync(file, `${lines.join('\n')}\n`);
}

async function build(repo: string, extra: Record<string, unknown> = {}): Promise<Built> {
  const log = vi.spyOn(console, 'log').mockImplementation(() => {});
  const error = vi.spyOn(console, 'error').mockImplementation(() => {});
  const before = process.exitCode;
  try {
    await kg({_: ['kg', 'build'], kg: path.join(repo, KG_DIR), backlog: path.join(repo, 'backlog'), db: false, output: 'json', ...extra});
    expect(error.mock.calls).toEqual([]);
    const out = path.join(repo, ...(extra.public ? ['public', '.kg'] : ['.kg']));
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
      backlog: 'ok@2026-01-12T07:00:00Z', dart: 'ok', docs: 'partial', git: 'partial', homes: 'ok', membership: 'ok',
      people: 'partial', process: 'ok', releases: 'ok', 'ts-changelog': 'partial', 'ts-declarations': 'ok',
      'ts-functions': 'partial', 'ts-imports': 'ok', 'ts-packages': 'ok', 'ts-samples': 'partial', 'ts-tests': 'ok', 'ts-uses': 'ok',
    });
    expect(manifest.counts.nodes).toEqual({
      app: 2, 'cell-renderer': 1, 'changelog-entry': 5, connection: 2, container: 3, customer: 2, 'db-table': 1,
      declaration: 105, 'doc-anchor': 20, 'doc-page': 12, editor: 1, endpoint: 1, feature: 8, 'file-handler': 1,
      'file-viewer': 1, filter: 1, function: 13, library: 2, 'lifecycle-hook': 2, package: 7, panel: 2, person: 2,
      query: 3, release: 2, sample: 4, scenario: 7, script: 3, 'script-environment': 1, 'script-handler': 1,
      'sem-type-detector': 6, 'semantic-type': 7, 'source-file': 39, test: 14, 'test-suite': 6, ticket: 12, tutorial: 1, viewer: 2,
    });
    expect(manifest.counts.edges).toEqual({
      affects: 2, assignee: 3, automates: 3, base: 1, calls: 6, changes: 2, connection: 3, container: 33, covers: 1,
      declared_in: 1, declares: 215, demonstrates: 1, 'depends-on': 5, documents: 1, environment: 1, extends: 4,
      handler: 1, implements: 1, imports: 25, 'in-suite': 14, 'is-implemented-in': 14, mentions: 15, owner: 5,
      package: 72, page: 20, 'part-of': 10, 'participates-in': 5, reporter: 4, 'requested-by': 2, router: 1,
      semtype: 1, target_semtype: 1, 'targets-release': 3, 'targets-semtype': 12, tests: 11, uses: 17,
    });
    expect(manifest.problems).toMatchObject({dangling_edges: 0, ambiguous_owners: 1, orphans: 27, partial_stubs: 21});
  }, 120_000);

  it('writes the same bytes twice, with the same content-addressed batch and a later built_at', async () => {
    const {out, manifest} = await graph;
    const first = dataFiles(out);
    const again = await build(built);
    expect(dataFiles(out)).toEqual(first);
    expect(again.manifest.batch).toBe(manifest.batch);
    expect(again.manifest.built_at >= manifest.built_at).toBe(true);
    expect(Object.values(first).every((text) => !text.includes('built_at'))).toBe(true);
  }, 120_000);

  it('reports the Dart batch as ok, stale and missing as the batch itself changes', async () => {
    const {manifest, rows} = await graph;
    expect(manifest.sources.dart).toBe('ok');
    expect(rows('nodes/endpoint')).toHaveLength(1);

    const stale = makeRepo();
    stampBatch(stale, {revision: 'deadbeefdeadbeefdeadbeefdeadbeefdeadbeef'});
    const staleBuild = await build(stale);
    expect(staleBuild.manifest.sources.dart).toBe('stale');
    expect(staleBuild.rows('nodes/db-table')).toHaveLength(1);

    const gone = makeRepo();
    stampBatch(gone, null);
    const goneBuild = await build(gone);
    expect(goneBuild.manifest.sources.dart).toBe('missing');
    expect(goneBuild.manifest.counts.nodes.endpoint).toBeUndefined();
    expect(goneBuild.manifest.counts.nodes['db-table']).toBeUndefined();
  }, 180_000);

  it('projects the public layer: public nodes only, no home, no owner, no source files, no reports', async () => {
    const {manifest, rows, out} = await build(makeRepo(), {public: true});
    expect(manifest.mode).toBe('public');
    expect(Object.keys(manifest.revisions)).toEqual(['public']);
    expect(Object.keys(manifest.counts.nodes).sort()).toEqual(['doc-anchor', 'doc-page', 'feature', 'library', 'package', 'sample', 'scenario']);
    const features = rows('nodes/feature');
    expect(features).toHaveLength(8);
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
    expect(orphans.summary).toBe('27 files with no owner in 7 groups, 386 lines; the 7 largest groups below.');
    expect(orphans.sections[0].rows[0]).toEqual({group: 'public/js-api', files: 11, loc: 154, largest: expect.stringContaining('src/dataframe.ts (45)')});
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
      'visualize/legends', 'visualize/viewers', 'visualize/viewers/scatter-plot']);
    expect(rows.find((r: any) => r.feature === 'domains/bio')).toEqual({feature: 'domains/bio', owner: 'P:jane', status: 'active',
      tests: 5, scenarios: 1, automated: 0, no_tests: false, no_docs: false, no_description: false});
    expect(rows.find((r: any) => r.feature === 'platform/caching')).toMatchObject({no_tests: true, no_docs: true, no_description: false});
    expect(rows.find((r: any) => r.feature === 'domains')).toMatchObject({owner: '', status: 'proposed', no_description: true});
  });

  it('proposes a feature id for every unowned package folder, and none for a library', async () => {
    const {report} = await graph;
    const rows = report('proposed').sections[0].rows;
    expect(rows).toEqual([
      {path: 'public/packages/Plain', files: 3, loc: 35, suggested: 'domains/plain'},
      {path: 'public/packages/Tutorials', files: 3, loc: 29, suggested: 'domains/tutorials'},
      {path: 'public/packages/ApiTests', files: 1, loc: 7, suggested: 'domains/api-tests'},
      {path: 'public/libraries/utils', files: 1, loc: 5, suggested: ''},
    ]);
    // the folders the Demo package owns through its home are not proposed again
    expect(rows.some((r: any) => r.path === 'public/packages/Demo')).toBe(false);
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
    expect(report.sections[1].rows).toEqual([{file: 'core/docs/CACHING.md', feature: 'platform/caching', relation: 'home', confidence: 1}]);
    expect(report.summary).toContain('1 files changed, 0 of them in no feature');
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
    // kuzu 0.11.3 kills the process when the type system is loaded while a Database is open, so it is read first
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
    expect((await run(conn, 'MATCH (n:Feature) RETURN count(n) AS n')).rows[0].n).toBe(8);
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

    const service = (await resolveTarget(conn, 'core/server/datlas/lib/src/services/bio_service.dart'))!;
    expect(service).toMatchObject({id: 'file:core/server/datlas/lib/src/services/bio_service.dart', root: 'Component'});
    const reached = await impact(conn, service, LIMIT);
    expect(reached.sections[0].rows).toMatchObject([{feature: 'domains/bio', relation: 'owns'}]);
    expect(reached.sections[1].rows).toEqual([{feature: 'domains/bio', owner: 'P:jane', name: 'Jane Dev'}]);
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
