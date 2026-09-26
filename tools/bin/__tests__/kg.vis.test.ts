/// The graph browser's render tier and server (core/docs/knowledge-graph/vis/): the blob exported from a
/// fixture generation, and — when the optional `kuzu` binding is installed — the routes over its index.
import {describe, it, expect, beforeAll, afterAll} from 'vitest';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {fileURLToPath} from 'url';
import {TypeSystem} from '../utils/kg/types';
import {load, loadKuzu, open} from '../utils/kg/kuzu';
import {readManifest} from '../utils/kg/generation';
import {exportVis, hasVis, readBlobHeader, visDir, BLOB, INDEX, SCHEMA} from '../utils/kg/vis';
import {serve, Served} from '../utils/kg/serve';
import {loadQuestions, resolveParams, isoDate, ask} from '../utils/kg/questions';
import {kg} from '../commands/kg';
import {copyFixture, buildFixture as buildGraph, fixtureTypes} from './kg-fixture';

const questionsRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..', '..', '..', '..', 'core', 'docs', 'knowledge-graph');
/** The question set lives in the private monorepo: a checkout of public/ alone (CI) has none. */
const withQuestions = fs.existsSync(questionsRoot) ? it : it.skip;

describe('the question set (questions/README.md)', () => {
  withQuestions('loads every file, each with its id, cypher and used parameters', () => {
    const {questions, errors} = loadQuestions(questionsRoot);
    expect(errors).toEqual([]);
    expect(questions.length).toBeGreaterThanOrEqual(30);
    for (const q of questions) {
      expect(q.cypher).toMatch(/MATCH/);
      for (const name of Object.keys(q.params)) expect(q.cypher, q.id).toContain(`$${name}`);
      if (q.status === 'blocked') expect(q.blocked_by, q.id).toBeTruthy();
    }
  });

  withQuestions('resolves parameters over the defaults and relative dates to midnight UTC', () => {
    const {questions} = loadQuestions(questionsRoot);
    const since = questions.find((q) => q.id === 'tickets-since')!;
    const now = new Date('2026-09-15T13:45:00Z');
    expect(resolveParams(since, {}, now)).toEqual({since: '2026-09-08T00:00:00.000Z'});
    expect(resolveParams(since, {since: '2026-01-02'}, now)).toEqual({since: '2026-01-02T00:00:00.000Z'});
    expect(() => resolveParams(since, {since: 'yesterday'}, now)).toThrow(/ISO date or a relative day count/);
    expect(() => resolveParams(since, {nope: 1}, now)).toThrow(/takes no parameter 'nope'/);
    expect(isoDate('+1d', now, 'x')).toBe('2026-09-16T00:00:00.000Z');
  });

  it('refuses a file whose id, params or status are wrong', () => {
    const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-kg-questions-'));
    fs.mkdirSync(path.join(dir, 'questions'));
    fs.writeFileSync(path.join(dir, 'questions', 'a.yaml'), 'id: b\nquestion: q\nwhy: w\ncypher: MATCH (n) RETURN n\n');
    fs.writeFileSync(path.join(dir, 'questions', 'c.yaml'), 'id: c\nquestion: q\nwhy: w\nparams: {x: {type: string}}\ncypher: MATCH (n) RETURN n\n');
    fs.writeFileSync(path.join(dir, 'questions', 'd.yaml'), 'id: d\nquestion: q\nwhy: w\ncypher: MATCH (n) RETURN n\nstatus: blocked\n');
    const {questions, errors} = loadQuestions(dir);
    expect(questions).toEqual([]);
    expect(errors.map((e) => e.split(':')[0])).toEqual(['a.yaml', 'c.yaml', 'd.yaml']);
    expect(errors[0]).toMatch(/does not match the file name/);
    expect(errors[1]).toMatch(/not used in the cypher/);
    expect(errors[2]).toMatch(/blocked_by/);
  });
});

const fixture = path.join(path.dirname(fileURLToPath(import.meta.url)), 'fixtures', 'kg', 'build');
const types = (): TypeSystem => fixtureTypes('build');

async function buildFixture(): Promise<{repo: string, kgDir: string}> {
  const {repo, out: kgDir} = await buildGraph(copyFixture('build'), 'homes,dart');
  return {repo, kgDir};
}

function lines(file: string): number {
  return fs.existsSync(file) ? fs.readFileSync(file, 'utf8').split('\n').filter(Boolean).length : 0;
}

describe('the render tier (vis/plan.md WO-2)', () => {
  let kgDir: string;
  beforeAll(async () => {
    kgDir = (await buildFixture()).kgDir;
  });

  it('writes the blob, the ids and the type tree once, byte for byte the same twice', () => {
    expect(hasVis(kgDir)).toBe(false);
    const first = exportVis(kgDir, types(), 'b-test');
    expect(hasVis(kgDir)).toBe(true);
    const bytes = fs.readFileSync(path.join(visDir(kgDir), BLOB));
    expect(first.bytes).toBe(bytes.length);
    exportVis(kgDir, types(), 'b-test');
    expect(fs.readFileSync(path.join(visDir(kgDir), BLOB)).equals(bytes)).toBe(true);
    const manifest = readManifest(kgDir)!;
    const nodes = Object.values(manifest.counts.nodes).reduce((a, b) => a + b, 0);
    const edges = Object.values(manifest.counts.edges).reduce((a, b) => a + b, 0);
    expect(first.nodes).toBe(nodes);
    expect(first.edges + first.dropped).toBe(edges);
    expect(first.dropped).toBe(0);
  });

  it('carries every node and edge as typed sections the page can view without copying', () => {
    const header = readBlobHeader(path.join(visDir(kgDir), BLOB));
    expect(header.version).toBe(1);
    expect(header.batch).toBe('b-test');
    const index = JSON.parse(fs.readFileSync(path.join(visDir(kgDir), INDEX), 'utf8'));
    expect(index.ids.length).toBe(header.nodes);
    expect(index.names.length).toBe(header.nodes);
    for (const name of ['type', 'layer', 'visibility', 'status', 'provenance', 'degree', 'parent'])
      expect(header.sections[name].length, name).toBe(header.nodes);
    for (const name of ['from', 'to', 'kind', 'confidence', 'derivedBy'])
      expect(header.sections[name].length, name).toBe(header.edges);
    for (const s of Object.values(header.sections)) expect(s.offset % 4).toBe(0);
    // the feature file feeds the first node type and part-of the parent column
    expect(header.nodeTypes[0]).toBe('feature');
    expect(header.nodeTypes.length).toBe(fs.readdirSync(path.join(kgDir, 'data', 'nodes')).length);
    expect(header.edgeKinds).toContain('part-of');
    expect(header.edgeKinds.length).toBe(fs.readdirSync(path.join(kgDir, 'data', 'edges')).length);
    expect(header.enums.provenance).toContain('annotation');
    const blob = fs.readFileSync(path.join(visDir(kgDir), BLOB));
    const base = (8 + blob.readUInt32LE(4) + 3) & ~3;
    const parent = new Int32Array(blob.buffer.slice(blob.byteOffset + base + header.sections.parent.offset,
      blob.byteOffset + base + header.sections.parent.offset + header.nodes * 4));
    const child = index.ids.indexOf('platform/caching');
    expect(child).toBeGreaterThanOrEqual(0);
    expect(index.ids[parent[child]]).toBe('platform');
    expect(lines(path.join(kgDir, 'data', 'edges', 'part-of.jsonl'))).toBeGreaterThan(0);
  });

  it('describes the type tree with roots, chains, labels and the reference properties', () => {
    const schema = JSON.parse(fs.readFileSync(path.join(visDir(kgDir), SCHEMA), 'utf8'));
    expect(schema.roots).toEqual(['feature', 'concept', 'component', 'artifact', 'work', 'actor', 'infra', 'type']);
    const feature = schema.nodeTypes.find((t: any) => t.name === 'feature');
    expect(feature.extends).toBe('node');
    expect(feature.hierarchical).toBe(true);
    const partOf = schema.edgeTypes.find((e: any) => e.name === 'part-of');
    expect(partOf.label).toBe('PART_OF');
    expect(schema.refs.find((r: any) => r.name === 'owner').to).toContain('person');
  });
});

const withKuzu = loadKuzu() ? it : it.skip;
const withKuzuAndQuestions = loadKuzu() && fs.existsSync(questionsRoot) ? it : it.skip;

describe('the server (vis/plan.md WO-3)', () => {
  let kgDir: string;
  let repo: string;
  let served: Served;
  let opened: Awaited<ReturnType<typeof open>>;

  beforeAll(async () => {
    if (!loadKuzu()) return;
    ({repo, kgDir} = await buildFixture());
    await load(kgDir, types());
    exportVis(kgDir, types(), 'b-test');
    opened = await open(kgDir, true);
    served = await serve({genDir: kgDir, repoRoot: repo, manifest: readManifest(kgDir)!, system: types(), questions: loadQuestions(questionsRoot).questions,
      db: opened!.db, conn: opened!.conn, port: 0});
  });

  afterAll(async () => {
    await served?.close();
    await opened?.conn.close();
    await opened?.db.close();
  });

  const get = async (route: string) => {
    const r = await fetch(served.url + route.replace(/^\//, ''));
    return {status: r.status, body: r.headers.get('content-type')?.includes('json') ? await r.json() : await r.arrayBuffer()};
  };

  withKuzu('serves the page, the manifest with its caveats and the render tier', async () => {
    const page = await fetch(served.url);
    expect(page.status).toBe(200);
    expect(await page.text()).toContain('<title>Knowledge graph</title>');
    const manifest = await get('/api/manifest');
    expect(manifest.status).toBe(200);
    expect(manifest.body.batch).toBeDefined();
    expect(manifest.body.repo_root).not.toContain('\\');
    expect(Array.isArray(manifest.body.notes)).toBe(true);
    const blob = await get('/api/graph.bin');
    expect((blob.body as ArrayBuffer).byteLength).toBe(fs.statSync(path.join(visDir(kgDir), BLOB)).size);
    expect((await get('/api/schema')).body.roots).toContain('feature');
    expect((await get('/api/index')).body.ids).toContain('platform/caching');
  });

  withKuzu('answers a node with its edge groups, an edge, a bounded query and an operation', async () => {
    const node = await get('/api/node?id=platform/caching');
    expect(node.status).toBe(200);
    expect(node.body.node.name).toBe('Caching');
    const partOf = node.body.edges.find((g: any) => g.edge === 'PART_OF' && g.direction === 'out');
    expect(partOf.count).toBe(1);
    expect(partOf.targets[0].id).toBe('platform');
    expect(partOf.derived_by).toEqual(['filesystem']);
    const edge = await get('/api/edge?from=platform/caching&to=platform&kind=PART_OF');
    expect(edge.body.edge.derived_by).toBe('filesystem');
    const r = await fetch(served.url + 'api/query', {method: 'POST', headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({cypher: 'MATCH (f:Feature) RETURN f.id AS id ORDER BY id', limit: 2})});
    const query = await r.json();
    expect(query.columns).toEqual(['id']);
    expect(query.rows.length).toBe(2);
    expect(query.truncated).toBe(true);
    const explain = await get('/api/op/explain?target=~platform/caching');
    expect(explain.body.op).toBe('explain');
    expect(explain.body.target.id).toBe('platform/caching');
    expect(explain.body.sections.map((s: any) => s.title)).toContain('edges');
    expect(explain.body.notes).toBeDefined();
    const find = await get('/api/op/find?target=caching');
    expect(find.body.sections[0].rows.some((row: any) => row.id === 'platform/caching')).toBe(true);
  });

  withKuzuAndQuestions('lists the questions and answers one with its parameters bound', async () => {
    const list = await get('/api/questions');
    expect(list.body.length).toBeGreaterThanOrEqual(30);
    expect(list.body[0].file).toBeUndefined();
    const tree = await get('/api/ask/feature-tree');
    expect(tree.status).toBe(200);
    expect(tree.body.rows.some((r: any) => r.id === 'platform/caching' && r.parent === 'platform')).toBe(true);
    expect(tree.body.highlight).toEqual(['id', 'parent']);
    const tests = await get('/api/ask/tests-for-feature?feature=platform');
    expect(tests.status).toBe(200);
    expect(tests.body.params).toEqual({feature: 'platform'});
    const blocked = await get('/api/ask/tickets-per-feature');
    expect(blocked.body.notes[0]).toMatch(/^blocked: /);
    expect((await get('/api/ask/nothing')).status).toBe(404);
    expect((await get('/api/ask/tickets-since?since=yesterday')).status).toBe(400);
    const direct = await ask(opened!.conn, loadQuestions(questionsRoot).questions.find((q) => q.id === 'feature-tree')!, {});
    expect(direct.rows.length).toBeGreaterThan(0);
    const bound = await fetch(served.url + 'api/query', {method: 'POST', headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({cypher: 'MATCH (f:Feature) WHERE f.id = $id RETURN f.name AS name', params: {id: 'platform/caching'}})});
    expect((await bound.json()).rows).toEqual([{name: 'Caching'}]);
  });

  withKuzu('refuses what it does not serve', async () => {
    expect((await get('/api/node?id=nope')).status).toBe(404);
    expect((await get('/api/node')).status).toBe(400);
    expect((await get('/api/nothing')).status).toBe(404);
    expect((await get('/api/op/drop?target=x')).status).toBe(404);
    expect((await fetch(served.url + '../package.json')).status).toBe(404);
    expect((await fetch(served.url + 'api/query')).status).toBe(405);
    const bad = await fetch(served.url + 'api/query', {method: 'POST', body: JSON.stringify({cypher: 'MATCH (n RETURN n'})});
    expect(bad.status).toBe(500);
    expect((await bad.json()).error).toBeTruthy();
  });
});
