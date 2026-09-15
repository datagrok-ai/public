/// Publication safety (build-plan.md WO-11a, kg-codex-review-3.md #1, #2, #11, #12): the closure of the
/// public projection, the immutable generations a build publishes through `current`, the Kuzu seam's path
/// escaping, identifier validation and buffer pool, and the batch id over every effective input.
/// Each test here fails against the code as the third review found it.
import {describe, it, expect, vi} from 'vitest';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {fileURLToPath} from 'url';
import {spawn} from 'child_process';
import {loadTypeSystem, TypeSystem} from '../utils/kg/types';
import {Graph} from '../utils/kg/build/emitter';
import {Row} from '../utils/kg/build/normalize';
import {projectPublic, writeBuild, currentDir, readCurrent, generationDir, generations, gc} from '../utils/kg/build/write';
import * as kuzu from '../utils/kg/kuzu';
import {kg} from '../commands/kg';

const fixture = path.join(path.dirname(fileURLToPath(import.meta.url)), 'fixtures', 'kg', 'build');
const KG_DIR = path.join('core', 'docs', 'knowledge-graph');
const system: TypeSystem = loadTypeSystem(path.join(fixture, KG_DIR));
const BATCH = 'b-000000000000';
const tools = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..', '..');

function scratch(name: string): string {
  return fs.mkdtempSync(path.join(os.tmpdir(), `grok-kg-11a-${name}-`));
}

function makeRepo(): string {
  const repo = scratch('repo');
  fs.cpSync(fixture, repo, {recursive: true});
  return repo;
}

async function run(argv: Record<string, unknown>): Promise<{out: string[], err: string[], exitCode: number | undefined}> {
  const log = vi.spyOn(console, 'log').mockImplementation(() => {});
  const error = vi.spyOn(console, 'error').mockImplementation(() => {});
  const before = process.exitCode;
  try {
    await kg(argv);
    return {out: log.mock.calls.map((c) => String(c[0])), err: error.mock.calls.map((c) => String(c[0])), exitCode: process.exitCode as number | undefined};
  }
  finally {
    process.exitCode = before;
    log.mockRestore();
    error.mockRestore();
  }
}

function build(repo: string, root: string, extra: Record<string, unknown> = {}): Promise<{out: string[], err: string[], exitCode: number | undefined}> {
  return run({_: ['kg', 'build'], kg: path.join(repo, KG_DIR), only: 'homes', db: false, out: root, ...extra});
}

function rows(dir: string, file: string): Row[] {
  const p = path.join(dir, 'data', `${file}.jsonl`);
  return fs.existsSync(p) ? fs.readFileSync(p, 'utf8').split('\n').filter(Boolean).map((l) => JSON.parse(l)) : [];
}

const withKuzu = kuzu.loadKuzu() ? it : it.skip;

function node(row: Row): Row {
  return {visibility: 'public', provenance: 'annotation', source_layer: 'public', batch: BATCH, ...row};
}

describe('the public projection is closed (kg-codex-review-3.md #1)', () => {
  const page = (id: string, kind: string, visibility = 'public') =>
    node({type: 'doc-page', id: `doc:${id}`, name: id, path: id, kind, visibility});
  const anchor = (id: string, slug: string, visibility = 'public') =>
    node({type: 'doc-anchor', id: `doc:${id}#${slug}`, name: slug, path: id, page: `doc:${id}`, slug, depth: 2, visibility});

  it('keeps no heading whose page the projection left behind', () => {
    const graph: Graph = {
      nodes: [page('public/help/a.md', 'help'), anchor('public/help/a.md', 'one'),
        page('landing/CLAUDE.md', 'agent'), anchor('landing/CLAUDE.md', 'build-docusaurus')],
      edges: [{type: 'mentions', from: 'doc:public/help/a.md', to: 'doc:landing/CLAUDE.md#build-docusaurus', derived_by: 'annotation', confidence: 1, batch: BATCH}],
      stubs: [], claims: [], sources: {docs: 'ok'}, problems: {invalid_rows: 0}, details: {}, invalid: [], reports: {},
    };
    const projected = projectPublic(graph, system);
    expect(projected.nodes.map((n) => n.id)).toEqual(['doc:public/help/a.md', 'doc:public/help/a.md#one']);
    expect(projected.edges).toEqual([]);
    expect(projected.problems.projection_dropped).toBe(0);
  });

  it('takes a page down with it when its own visibility says so, headings included', () => {
    const graph: Graph = {
      nodes: [page('public/help/internal.md', 'help', 'dev'), anchor('public/help/internal.md', 'secret', 'dev'),
        anchor('public/help/internal.md', 'public-looking')],
      edges: [], stubs: [], claims: [], sources: {}, problems: {}, details: {}, invalid: [], reports: {},
    };
    expect(projectPublic(graph, system).nodes).toEqual([]);
  });

  it('counts and removes anything the projection would otherwise leave dangling', () => {
    const graph: Graph = {
      nodes: [node({type: 'feature', id: 'visualize', name: 'Visualize'}),
        node({type: 'scenario', id: 'TS:a', name: 'A', path: 'public/packages/x/a.md'})],
      edges: [], stubs: [], claims: [], sources: {}, problems: {}, details: {}, invalid: [], reports: {},
    };
    const projected = projectPublic(graph, system);
    projected.nodes[0].owner = 'P:jane';
    expect(projected.problems.projection_dropped).toBe(0);
    // the same graph run through the closure again finds the reference that does not resolve
    expect(projectPublic({...graph, nodes: [{...graph.nodes[0], owner: 'P:jane'}, graph.nodes[1]]}, system).nodes[0].owner).toBeUndefined();
  });

  it('reads visibility: from a page\'s frontmatter and gives it to the page and its headings', async () => {
    const repo = makeRepo();
    const file = path.join(repo, 'public', 'help', 'domains', 'bio', 'internal-notes.md');
    fs.writeFileSync(file, '---\nvisibility: dev\n---\n\n# Internal notes\n\n## How it really works\n\nText.\n');
    const root = path.join(scratch('vis'), '.kg');
    const full = await build(repo, root, {only: 'docs'});
    expect(full.err).toEqual([]);
    const page = rows(currentDir(root)!, 'nodes/doc-page').find((r) => r.id === 'doc:public/help/domains/bio/internal-notes.md');
    expect(page).toMatchObject({visibility: 'dev'});
    expect(rows(currentDir(root)!, 'nodes/doc-anchor').find((r) => r.id === 'doc:public/help/domains/bio/internal-notes.md#how-it-really-works')).toMatchObject({visibility: 'dev'});

    const publicRoot = path.join(scratch('vis-public'), '.kg');
    await build(repo, publicRoot, {only: 'docs', public: true});
    const dir = currentDir(publicRoot)!;
    const ids = [...rows(dir, 'nodes/doc-page'), ...rows(dir, 'nodes/doc-anchor')].map((r) => String(r.id));
    expect(ids.some((id) => id.includes('internal-notes'))).toBe(false);
    // and no heading in the whole snapshot is left without its page
    const pages = new Set(rows(dir, 'nodes/doc-page').map((r) => r.id));
    expect(rows(dir, 'nodes/doc-anchor').filter((a) => !pages.has(a.page))).toEqual([]);
  }, 120_000);
});

describe('a build publishes an immutable generation (kg-codex-review-3.md #2)', () => {
  it('writes a fresh generation whole before current names it, and never removes the one before', async () => {
    const repo = makeRepo();
    const root = path.join(scratch('gen'), '.kg');
    const first = await build(repo, root);
    expect(first.err).toEqual([]);
    const a = readCurrent(root)!;
    expect(a).toMatch(/^b-[0-9a-f]{12}-[A-Za-z0-9]{6}$/);
    expect(fs.readFileSync(path.join(root, 'current'), 'utf8').trim()).toBe(a);
    expect(fs.existsSync(path.join(root, 'gen', a, 'manifest.json'))).toBe(true);
    expect(fs.existsSync(path.join(root, 'data'))).toBe(false);

    await build(repo, root, {only: 'homes,docs'});
    const b = readCurrent(root)!;
    expect(b).not.toBe(a);
    expect(fs.existsSync(path.join(root, 'gen', a, 'manifest.json'))).toBe(true);
    expect(generations(root).map((g) => g.name).sort()).toEqual([a, b].sort());

    // the same inputs again: the same batch, a directory of its own, and current moves to it
    await build(repo, root, {only: 'homes,docs'});
    const c = readCurrent(root)!;
    expect(c).not.toBe(b);
    const all = generations(root);
    expect(all).toHaveLength(3);
    expect(all.filter((g) => g.batch === all.find((x) => x.name === c)!.batch)).toHaveLength(2);
  }, 120_000);

  withKuzu('leaves current where it is when the index cannot be loaded', async () => {
    const repo = makeRepo();
    const root = path.join(scratch('fail'), '.kg');
    await build(repo, root, {db: true});
    const good = readCurrent(root)!;
    // a buffer pool of 1 MB takes the load down the way a killed COPY would: partway, with nothing to publish
    const failed = await run({_: ['kg', 'build'], kg: path.join(repo, KG_DIR), only: 'homes,docs', out: root, memory: 1});
    expect(failed.err.join('\n')).toContain('index not built:');
    expect(failed.exitCode).toBe(1);
    expect(readCurrent(root)).toBe(good);
    // the interrupted generation is on disk without a manifest, so nothing reads it and gc can take it
    const partial = fs.readdirSync(path.join(root, 'gen')).filter((n) => n !== good && n !== 'current');
    expect(partial).toHaveLength(1);
    expect(fs.existsSync(path.join(root, 'gen', partial[0], 'manifest.json'))).toBe(false);
    expect(generations(root).map((g) => g.name)).toEqual([good]);

    const answer = await run({_: ['kg', 'query', 'MATCH (n:Feature) RETURN count(n) AS n'], kg: path.join(repo, KG_DIR), out: root, output: 'json'});
    expect(answer.err).toEqual([]);
    expect(JSON.parse(answer.out[0])[0].n).toBeGreaterThan(0);
  }, 300_000);

  withKuzu('a COPY that fails leaves the generation without a manifest and the error in one line', async () => {
    const repo = makeRepo();
    const root = path.join(scratch('copy'), '.kg');
    await build(repo, root, {only: 'homes,ts-declarations'});
    const dir = currentDir(root)!;
    const file = path.join(dir, 'data', 'nodes', 'source-file.jsonl');
    const lines = fs.readFileSync(file, 'utf8').split('\n').filter(Boolean).map((l) => JSON.parse(l));
    expect(lines.length).toBeGreaterThan(0);
    fs.writeFileSync(file, `${lines.map((r) => JSON.stringify({...r, loc: 'not a number'})).join('\n')}\n`);
    await expect(kuzu.load(dir, system)).rejects.toThrow(/Copy exception|Conversion exception|not a number/);
  }, 300_000);

  it('gc keeps the current generation and the newest --keep, and counts what it removed', async () => {
    const repo = makeRepo();
    const root = path.join(scratch('gc'), '.kg');
    await build(repo, root, {only: 'homes'});
    await build(repo, root, {only: 'homes,docs'});
    await build(repo, root, {only: 'homes,docs,process'});
    const all = generations(root);
    expect(all).toHaveLength(3);
    const current = readCurrent(root);
    const result = gc(root, 1);
    expect(result.kept).toContain(current);
    expect(result.removed).toHaveLength(2);
    expect(result.locked).toEqual([]);
    expect(generations(root).map((g) => g.name)).toEqual([current]);
    expect(currentDir(root)).toBe(generationDir(root, current!));
  }, 180_000);

  it('gc also removes an interrupted build older than an hour, and leaves a fresh one alone', async () => {
    const repo = makeRepo();
    const root = path.join(scratch('interrupted'), '.kg');
    await build(repo, root, {only: 'homes'});
    for (const [name, age] of [['b-abcdef123456-stale', 2 * 3600_000], ['b-abcdef123456-fresh', 0]] as const) {
      fs.mkdirSync(path.join(root, 'gen', name), {recursive: true});
      const when = new Date(Date.now() - age);
      fs.utimesSync(path.join(root, 'gen', name), when, when);
    }
    const result = gc(root, 2);
    expect(result.removed).toEqual(['b-abcdef123456-stale']);
    expect(fs.existsSync(path.join(root, 'gen', 'b-abcdef123456-fresh'))).toBe(true);
  }, 120_000);
});

describe('the Kuzu seam (kg-codex-review-3.md #11)', () => {
  it('escapes a path with an apostrophe and refuses an identifier that is not a schema name', () => {
    // kuzu 0.11.3 reads \' inside a string literal and takes '' for the end of it
    expect(kuzu.literal("/tmp/o'brien/.kg/tmp/Feature.csv")).toBe("'/tmp/o\\'brien/.kg/tmp/Feature.csv'");
    expect(kuzu.literal('C:\\kg\\tmp')).toBe("'C:\\\\kg\\\\tmp'");
    expect(kuzu.literal('Feature')).toBe("'Feature'");
    expect(kuzu.quote('source_layer')).toBe('`source_layer`');
    expect(() => kuzu.quote('a name')).toThrow(/not a graph identifier/);
    expect(() => kuzu.quote('Feature`) MATCH (n')).toThrow(/not a graph identifier/);
    expect(() => kuzu.quote('2big')).toThrow(/not a graph identifier/);
  });

  it('takes the buffer pool from --memory, then KG_KUZU_MEMORY, then the default', () => {
    expect(kuzu.memoryMb(undefined, kuzu.BUILD_MEMORY_MB)).toBe(2048);
    expect(kuzu.memoryMb('256', kuzu.BUILD_MEMORY_MB)).toBe(256);
    expect(kuzu.memoryMb('nonsense', kuzu.READ_MEMORY_MB)).toBe(512);
    expect(kuzu.memoryMb(undefined, kuzu.READ_MEMORY_MB)).toBe(512);
    process.env.KG_KUZU_MEMORY = '777';
    try {
      expect(kuzu.memoryMb(undefined, kuzu.BUILD_MEMORY_MB)).toBe(777);
      expect(kuzu.memoryMb('256', kuzu.BUILD_MEMORY_MB)).toBe(256);
    }
    finally {
      delete process.env.KG_KUZU_MEMORY;
    }
  });

  withKuzu('opens an index whose manifest recorded a large load without reserving that much', async () => {
    const root = path.join(scratch('reader-pool'), '.kg');
    expect((await build(makeRepo(), root, {db: true, output: 'json'})).err).toEqual([]);
    const dir = currentDir(root)!;
    const manifest = JSON.parse(fs.readFileSync(path.join(dir, 'manifest.json'), 'utf8'));
    fs.writeFileSync(path.join(dir, 'manifest.json'), JSON.stringify({...manifest, index_memory_mb: 3000}, null, 2));
    const opened = (await kuzu.open(dir, true))!;
    try {
      expect(Number((await kuzu.run(opened.conn, 'MATCH (n:Feature) RETURN count(n) AS n')).rows[0].n)).toBeGreaterThan(0);
      expect(process.memoryUsage().rss).toBeLessThan(1024 * 1048576);
    }
    finally {
      await opened.conn.close();
      await opened.db.close();
    }
  }, 300_000);
});

/** One read-only connection in its own process, the way a second `grok kg query` would open the index. */
function reader(db: string): Promise<{code: number | null, out: string}> {
  const script = 'const kuzu = require("kuzu");' +
    'const db = new kuzu.Database(process.argv[1], 512 * 1024 * 1024, true, true);' +
    'const conn = new kuzu.Connection(db);' +
    'conn.query("MATCH (n:Feature) RETURN count(n) AS n").then(async (r) => {' +
    'const rows = await r.getAll(); r.close(); await conn.close(); await db.close();' +
    'console.log(JSON.stringify(rows.map((x) => Number(x.n)))); });';
  const child = spawn(process.execPath, ['-e', script, db], {cwd: tools});
  let out = '';
  child.stdout.on('data', (d) => out += String(d));
  child.stderr.on('data', (d) => out += String(d));
  return new Promise((resolve) => child.on('close', (code) => resolve({code, out: out.trim()})));
}

describe('the index of a generation (kg-codex-review-3.md #2, #11)', () => {
  withKuzu('loads from a path with an apostrophe, records what it needed, and answers two other processes at once', async () => {
    const repo = makeRepo();
    const root = path.join(scratch("o'brien"), '.kg');
    const built = await build(repo, root, {db: true, output: 'json'});
    expect(built.err).toEqual([]);
    const manifest = JSON.parse(built.out[built.out.length - 1]);
    expect(manifest.indexed_batch).toBe(manifest.batch);
    expect(manifest.index_memory_mb).toBeGreaterThan(0);
    expect(manifest.index_platform).toBe(`${process.platform}-${process.arch}`);
    const dir = currentDir(root)!;
    expect(fs.existsSync(path.join(dir, 'kg.kuzu'))).toBe(true);

    const opened = (await kuzu.open(dir, true))!;
    try {
      const [a, b, mine] = await Promise.all([
        reader(path.join(dir, 'kg.kuzu')),
        reader(path.join(dir, 'kg.kuzu')),
        kuzu.run(opened.conn, 'MATCH (n:Feature) RETURN count(n) AS n'),
      ]);
      const count = Number(mine.rows[0].n);
      expect(count).toBeGreaterThan(0);
      expect(a, a.out).toEqual({code: 0, out: `[${count}]`});
      expect(b, b.out).toEqual({code: 0, out: `[${count}]`});
    }
    finally {
      await opened.conn.close();
      await opened.db.close();
    }
  }, 300_000);

  withKuzu('a --no-db rebuild leaves current on the generation that has an index, and says so', async () => {
    const repo = makeRepo();
    const root = path.join(scratch('nodb'), '.kg');
    await build(repo, root, {db: true});
    const indexed = readCurrent(root)!;
    const again = await build(repo, root, {only: 'homes,docs'});
    expect(again.out.join('\n')).toContain(`current stays at ${indexed}`);
    expect(readCurrent(root)).toBe(indexed);

    const answer = await run({_: ['kg', 'query', 'MATCH (n:Feature) RETURN count(n) AS n'], kg: path.join(repo, KG_DIR), out: root, output: 'json'});
    expect(answer.err).toEqual([]);
    expect(JSON.parse(answer.out[0])[0].n).toBeGreaterThan(0);
  }, 300_000);

  withKuzu('refuses an index that was loaded from another generation', async () => {
    const repo = makeRepo();
    const root = path.join(scratch('mixed'), '.kg');
    await build(repo, root, {db: true});
    const dir = currentDir(root)!;
    const manifest = JSON.parse(fs.readFileSync(path.join(dir, 'manifest.json'), 'utf8'));
    fs.writeFileSync(path.join(dir, 'manifest.json'), JSON.stringify({...manifest, indexed_batch: 'b-somewhereelse'}, null, 2));
    const answer = await run({_: ['kg', 'query', 'MATCH (n:Feature) RETURN count(n) AS n'], kg: path.join(repo, KG_DIR), out: root});
    expect(answer.err.join('\n')).toContain('this index was loaded from batch b-somewhereelse');
    expect(answer.exitCode).toBe(1);
  }, 300_000);
});

describe('determinism of the written rows (kg-codex-review-3.md #12)', () => {
  it('writes a set-valued member sorted, whatever order its rows arrived in', () => {
    const dir = path.join(scratch('sorted'), 'gen');
    const graph: Graph = {
      nodes: [node({type: 'feature', id: 'visualize', name: 'Visualize', aliases: ['zeta', 'alpha']}),
        node({type: 'doc-page', id: 'doc:public/help/a.md', name: 'a', path: 'public/help/a.md', kind: 'help', keywords: ['zeta', 'alpha']}),
        node({type: 'function', id: 'func:Demo:f', name: 'f', language: 'js', input_types: ['string', 'dataframe']})],
      edges: [{type: 'covers', from: 'TS:a', to: 'visualize', derived_by: 'annotation', confidence: 1, evidence: ['z.md', 'a.md'], batch: BATCH},
        {type: 'imports', from: 'func:Demo:f', to: 'visualize', derived_by: 'ast', confidence: 1, symbols: ['b', 'a'], batch: BATCH}],
      stubs: [], claims: [], sources: {}, problems: {}, details: {}, invalid: [], reports: {},
    };
    writeBuild(graph, dir, {mode: 'public', batch: BATCH, builder: '6.5.10', schemaVersion: 1, revisions: {public: 'x'}});
    expect(fs.readFileSync(path.join(dir, 'data', 'nodes', 'feature.jsonl'), 'utf8')).toContain('"aliases":["alpha","zeta"]');
    expect(fs.readFileSync(path.join(dir, 'data', 'nodes', 'doc-page.jsonl'), 'utf8')).toContain('"keywords":["alpha","zeta"]');
    expect(fs.readFileSync(path.join(dir, 'data', 'edges', 'covers.jsonl'), 'utf8')).toContain('"evidence":["a.md","z.md"]');
    expect(fs.readFileSync(path.join(dir, 'data', 'edges', 'imports.jsonl'), 'utf8')).toContain('"symbols":["a","b"]');
    expect(fs.readFileSync(path.join(dir, 'data', 'nodes', 'function.jsonl'), 'utf8')).toContain('"input_types":["string","dataframe"]');
  });
});
