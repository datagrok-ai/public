/// The graph index (build-plan.md WO-7): the DDL the type system dictates, and — when the optional
/// `kuzu` binding is installed — a load of the fixture graph with `query` and the four operations over it.
import {describe, it, expect, vi, beforeAll, afterAll} from 'vitest';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {fileURLToPath} from 'url';
import {edgeGroups, TypeSystem} from '../utils/kg/types';
import {ddl, load, loadKuzu, open, run, Ddl, KuzuConnection, KuzuQueryResult} from '../utils/kg/kuzu';
import {find, explain, impact, testsFor, resolveTarget, groupOf, groupOrder} from '../utils/kg/ops';
import {readManifest} from '../utils/kg/generation';
import {EXTRACTORS, provides} from '../utils/kg/build/registry';
import {coverageNote, caveatNotes} from '../utils/kg/answer';
import {printAnswer} from '../utils/kg/print';
import {kg} from '../commands/kg';
import {copyFixture, buildFixture as buildGraph, fixtureTypes} from './kg-fixture';

const here = path.dirname(fileURLToPath(import.meta.url));
const fixture = path.join(here, 'fixtures', 'kg', 'build');
const types = (): TypeSystem => fixtureTypes('build');
const schema: Ddl = ddl(types());
const LIMIT = {limit: 50, groups: edgeGroups(fixtureTypes('build'))};

function table(name: string): string {
  const found = schema.statements.find((s) => s.startsWith(`CREATE NODE TABLE \`${name}\`(`) || s.startsWith(`CREATE REL TABLE \`${name}\`(`));
  expect(found, `no table ${name}`).toBeDefined();
  return found!;
}

describe('the index DDL (build-plan.md WO-7)', () => {
  it('gives every root one node table with the union of its subtree, typed lists and quoted identifiers', () => {
    expect(schema.nodes.map((t) => t.name)).toEqual(['Feature', 'Concept', 'Component', 'Artifact', 'Work', 'Actor', 'Infra', 'Type']);
    expect(table('Feature')).toBe('CREATE NODE TABLE `Feature`(`id` STRING PRIMARY KEY, `type` STRING, `types` STRING[], `name` STRING, ' +
      '`description` STRING, `status` STRING, `visibility` STRING, `owner` STRING, `aliases` STRING[], `source_layer` STRING, `home` STRING, ' +
      '`provenance` STRING, `batch` STRING, `developer_help` STRING, `user_help` STRING)');
    // the component subtree: a number is a DOUBLE, a boolean a BOOLEAN, a list of strings a STRING[]
    expect(table('Component')).toContain('`loc` DOUBLE');
    expect(table('Component')).toContain('`exported` BOOLEAN');
    expect(table('Component')).toContain('`path_params` STRING[]');
    // `key`, `type` and `from` are Cypher keywords, so every identifier is quoted
    expect(table('Work')).toContain('`key` STRING');
  });

  it('gives every concrete edge type a rel table named by its graph label, over every admissible root pair', () => {
    expect(table('COVERS')).toBe('CREATE REL TABLE `COVERS`(FROM `Artifact` TO `Feature`, `derived_by` STRING, `confidence` DOUBLE, ' +
      '`evidence` STRING[], `batch` STRING, `level` STRING)');
    const mentions = table('MENTIONS');
    expect(mentions.match(/FROM `\w+` TO `\w+`/g)!.length).toBe(16);
    expect(mentions).toContain('FROM `Artifact` TO `Feature`');
    expect(table('PART_OF')).toContain('FROM `Feature` TO `Feature`');
  });

  it('gives every reference property a rel table named after it, and merges one that collides with an edge label', () => {
    expect(table('owner')).toBe('CREATE REL TABLE `owner`(FROM `Feature` TO `Actor`, FROM `Concept` TO `Actor`, FROM `Component` TO `Actor`, ' +
      'FROM `Artifact` TO `Actor`, FROM `Work` TO `Actor`, FROM `Actor` TO `Actor`, FROM `Infra` TO `Actor`, FROM `Type` TO `Actor`, ' +
      '`derived_by` STRING, `confidence` DOUBLE, `evidence` STRING[], `batch` STRING)');
    expect(table('router')).toContain('FROM `Component` TO `Component`');
    expect(table('suite')).toContain('FROM `Artifact` TO `Artifact`');
    expect(table('base')).toContain('FROM `Infra` TO `Infra`');
    expect(schema.statements.some((s) => s.includes('`IN_SUITE`'))).toBe(false);
    // kuzu identifiers are case-insensitive: the `extends` edge and the `extends` reference of a type share one table
    expect(table('EXTENDS')).toBe('CREATE REL TABLE `EXTENDS`(FROM `Component` TO `Component`, FROM `Type` TO `Type`, ' +
      '`derived_by` STRING, `confidence` DOUBLE, `evidence` STRING[], `batch` STRING)');
  });

  it('widens two scalar declarations of one column and refuses a scalar against a list', () => {
    // doc-anchor.level is a number and test.level an enum: the artifact table carries both as text
    expect(table('Artifact')).toContain('`level` STRING');
    const conflicting = types();
    conflicting.nodes.get('scenario')!.members.keywords = {...conflicting.nodes.get('doc-page')!.members.keywords, list: false, spec: 'string'};
    expect(() => ddl(conflicting)).toThrow(/kuzu table Artifact: doc-page declares keywords as string\[], scenario as string/);
  });

  it('refuses two tables that differ only in case', () => {
    const colliding = types();
    const owner = colliding.nodes.get('node')!.own.owner;
    colliding.nodes.get('concept')!.own.feature = {...owner, name: 'feature'};
    expect(() => ddl(colliding)).toThrow(/tables 'Feature' and 'feature' differ only in case/);
  });
});

/** The fixture graph as JSONL, with one feature carrying values CSV cannot express. */
async function buildFixture(): Promise<{kgDir: string, feature: Record<string, unknown>}> {
  const {out: kgDir} = await buildGraph(copyFixture('build'), 'homes,dart');
  const file = path.join(kgDir, 'data', 'nodes', 'feature.jsonl');
  const rows = fs.readFileSync(file, 'utf8').split('\n').filter(Boolean).map((l) => JSON.parse(l));
  const feature = rows.find((r) => r.id === 'platform/caching')!;
  feature.aliases = ['cache, the', 'a "quoted" one'];
  feature.description = 'One line, with a comma.\nAnd a "second" one.';
  fs.writeFileSync(file, rows.map((r) => JSON.stringify(r)).join('\n') + '\n');
  return {kgDir, feature};
}

/** Remembers what a query handed back, so a test can prove `run` closed it. */
function keep(results: KuzuQueryResult[], result: KuzuQueryResult | KuzuQueryResult[]): KuzuQueryResult | KuzuQueryResult[] {
  results.push(...(Array.isArray(result) ? result : [result]));
  return result;
}

const withKuzu = loadKuzu() ? it : it.skip;

describe('the index itself (build-plan.md WO-7, WO-10)', () => {
  /** One index for both tests: a kuzu database reserves its whole buffer pool, and two at once exhaust a test worker. */
  let index: {kgDir: string, feature: Record<string, unknown>, result: Awaited<ReturnType<typeof load>>, opened: Awaited<ReturnType<typeof open>>};

  beforeAll(async () => {
    if (!loadKuzu()) return;
    const built = await buildFixture();
    const result = await load(built.kgDir, types());
    index = {...built, result, opened: await open(built.kgDir, true)};
  });

  const conn = () => index.opened!.conn;

  afterAll(async () => {
    if (!index?.opened) return;
    await index.opened.conn.close();
    await index.opened.db.close();
  });

  withKuzu('loads the fixture graph and round-trips a list and a text CSV cannot carry', async () => {
    const {kgDir, feature, result} = index;
    expect(result.nodes.find((t) => t.table === 'Feature')!.rows).toBe(9);
    expect(result.nodes.find((t) => t.table === 'Component')!.rows).toBe(26);
    expect(result.rels.find((t) => t.table === 'PART_OF')!.rows).toBe(8);
    expect(result.rels.find((t) => t.table === 'owner')!.rows).toBe(6);
    expect(result.parameterized).toBe(1);
    expect(fs.existsSync(path.join(kgDir, 'tmp'))).toBe(false);

    const {conn} = index.opened!;
    const {rows} = await run(conn, 'MATCH (n:Feature) WHERE n.`id` = $id RETURN n.`aliases` AS aliases, n.`description` AS description, n.`types` AS types', {id: feature.id});
    expect(rows[0].aliases).toEqual(feature.aliases);
    expect(rows[0].description).toBe(feature.description);
    expect(rows[0].types).toEqual(['feature', 'node']);
    const count = await run(conn, 'MATCH (n:Feature) RETURN count(n) AS n');
    expect(count.rows[0].n).toBe(9);
    // the Dart pass reached the index, and an integral DOUBLE comes back as a number, not 4.0
    const loc = await run(conn, 'MATCH (n:Component) WHERE n.`id` = $id RETURN n.`loc` AS loc', {id: 'file:core/client/d4/lib/src/legends/legend.dart'});
    expect(loc.rows[0].loc).toBe(4);
  });

  withKuzu('closes every query result, and reads the type system with the database open', async () => {
    const {conn} = index.opened!;
    const results: KuzuQueryResult[] = [];
    const watched: KuzuConnection = {
      query: async (statement) => keep(results, await conn.query(statement)),
      prepare: (statement) => conn.prepare(statement),
      execute: async (prepared, params) => keep(results, await conn.execute(prepared, params)) as KuzuQueryResult,
      close: () => conn.close(),
    };
    expect((await run(watched, 'MATCH (n:Feature) RETURN count(n) AS n')).rows[0].n).toBe(9);
    expect((await run(watched, 'MATCH (n:Feature) WHERE n.`id` = $id RETURN n.`id` AS id', {id: 'platform/caching'})).rows).toHaveLength(1);
    // a result left open when the database closes kills the process at exit (kuzu.ts): both of these are closed already
    expect(results).toHaveLength(2);
    for (const result of results) await expect(result.getAll()).rejects.toThrow(/closed/);
    // and the type system is read with the database open, which the same teardown used to take down with it
    expect(types().nodes.size).toBeGreaterThan(0);
  });

  withKuzu('answers find, explain, impact and tests-for', async () => {
    const {conn} = index.opened!;
    {
      const found = await find(conn, 'caching', LIMIT);
      expect(found.sections[0].rows[0]).toMatchObject({id: 'platform/caching', type: 'feature', name: 'Caching'});
      // the alias with a comma is searchable, the way it was written
      expect((await find(conn, 'cache, the', LIMIT)).sections[0].rows[0]).toMatchObject({id: 'platform/caching'});

      const caching = (await resolveTarget(conn, '~platform/caching'))!;
      expect(caching).toMatchObject({id: 'platform/caching', root: 'Feature', type: 'feature'});
      const explained = await explain(conn, caching, LIMIT);
      expect(explained.sections[0].rows).toContainEqual({property: 'home', value: 'core/docs/CACHING.md'});
      expect(explained.sections[1].rows).toContainEqual(expect.objectContaining({group: 'reference', edge: 'owner', direction: 'out', count: 1, targets: [{id: 'P:jane', name: 'Jane Dev', type: 'person'}]}));
      expect(explained.sections[1].rows).toContainEqual(expect.objectContaining({group: 'tree', edge: 'PART_OF', direction: 'out', targets: [expect.objectContaining({id: 'platform'})]}));
      // the edges section is ordered by group, folder order first and the reference predicates last
      const groups = explained.sections[1].rows.map((r) => String(r.group));
      expect(groups).toEqual([...groups].sort((a, b) => groupOrder(LIMIT.groups).indexOf(a) - groupOrder(LIMIT.groups).indexOf(b)));

      const reached = await impact(conn, caching, LIMIT);
      expect(reached.sections[0].rows).toEqual([{feature: 'platform/caching', relation: 'self', name: 'Caching', status: 'active',
        via: 'platform/caching', path: ['platform/caching']}]);
      expect(reached.sections[1].rows).toEqual([{feature: 'platform/caching', owner: 'P:jane', name: 'Jane Dev'}]);

      const bio = (await resolveTarget(conn, 'domains/bio'))!;
      const tests = await testsFor(conn, bio, LIMIT);
      expect(tests.sections.find((s) => s.title === 'scenarios')!.rows).toMatchObject([{scenario: 'TS:viewers/scatter-plot/ui', feature: 'domains/bio', manual_only: true}]);
      expect(tests.sections.find((s) => s.title === 'tests')!.rows).toEqual([]);
      expect(coverageNote({dart: 'partial'})).toBe('Dart coverage partial (some markers did not resolve)');
      expect(coverageNote(undefined)).toBe('Dart coverage unknown (the dart extractor did not run)');
      expect(coverageNote({dart: 'ok'})).toBeUndefined();
    }
  });

  withKuzu('puts the evidence behind every edge group, and caveats every source that is not ok', async () => {
    const {conn} = index.opened!;
    const caching = (await resolveTarget(conn, '~platform/caching'))!;
    const edges = (await explain(conn, caching, LIMIT)).sections.find((s) => s.title === 'edges')!;
    expect(edges.rows).toContainEqual(expect.objectContaining({edge: 'owner', direction: 'out', derived_by: ['annotation'], confidence: [1, 1]}));
    expect(edges.rows).toContainEqual(expect.objectContaining({edge: 'PART_OF', direction: 'out', derived_by: ['filesystem'], confidence: [1, 1]}));
    expect(edges.rows.every((r) => r.derived_by !== undefined && r.confidence !== undefined && r.evidence !== undefined)).toBe(true);
    // the Dart clause used to be the only one; a partial docs or people source was silently passed off as complete
    expect(caveatNotes({sources: {dart: 'ok', backlog: 'missing', docs: 'partial', git: 'ok@2026-01-12T07:00:00Z', people: 'partial'}, provides: provides(EXTRACTORS)})).toEqual([
      'backlog missing: no coverage of tickets, their state and who they are assigned to',
      'docs partial: incomplete coverage of documentation pages, their headings and the mentions in them',
      'people partial: incomplete coverage of people, teams and customers',
    ]);
    expect(caveatNotes({sources: {dart: 'partial', homes: 'ok'}})).toEqual(['Dart coverage partial (some markers did not resolve)']);
    // a generation that predates `provides` still gets a caveat, without the description
    expect(caveatNotes({sources: {dart: 'ok', backlog: 'missing'}})).toEqual(['backlog missing: no coverage of what backlog contributes']);
  });

  withKuzu('reads the edge groups from the manifest, and answers ungrouped over a generation that predates them', async () => {
    const manifest = readManifest(index.kgDir)!;
    expect(manifest.edge_groups).toEqual(edgeGroups(types()));
    expect(manifest.provides).toEqual(provides(EXTRACTORS.filter((e) => ['homes', 'dart'].includes(e.name))));
    expect(groupOf('PART_OF', manifest.edge_groups)).toBe('tree');
    expect(groupOf('owner', manifest.edge_groups)).toBe('reference');
    expect(groupOf('user_help', manifest.edge_groups)).toBe('reference');
    expect(groupOf('PART_OF')).toBe('');
    const caching = (await resolveTarget(conn(), '~platform/caching'))!;
    const old = (await explain(conn(), caching, {limit: 50})).sections.find((s) => s.title === 'edges')!;
    expect(old.rows.length).toBeGreaterThan(0);
    expect(old.rows.map((r) => r.group)).toEqual(old.rows.map((r) => String(r.edge) === String(r.edge).toLowerCase() ? 'reference' : ''));
  });

  withKuzu('says why a section came back empty instead of printing nothing', async () => {
    const {conn} = index.opened!;
    const reason = 'tests (0): no test carries ~domains/bio and no owned file contains tests';
    const tests = await testsFor(conn, (await resolveTarget(conn, 'domains/bio'))!, LIMIT);
    expect(tests.sections.find((s) => s.title === 'tests')!.empty).toBe(reason);
    const log = vi.spyOn(console, 'log').mockImplementation(() => {});
    let printed: string[];
    try {
      printAnswer(tests, 'table');
      printed = log.mock.calls.map((c) => String(c[0]));
    }
    finally {
      log.mockRestore();
    }
    expect(printed).toContain(`\n${reason}`);
  });
});
