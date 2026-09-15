/// `grok kg build` (build-plan.md WO-1, WO-2): the emitter's merge rules and finalize passes, the
/// deterministic writer and manifest, the public projection, and the home-document layer's rows
/// against the mini monorepo under fixtures/kg/good.
import {describe, it, expect, vi} from 'vitest';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {fileURLToPath} from 'url';
import {spawnSync} from 'child_process';
import {loadTypeSystem, TypeSystem} from '../utils/kg/types';
import {Emitter, Graph} from '../utils/kg/build/emitter';
import {normalizeRow} from '../utils/kg/build/normalize';
import {batchId, buildInputs, currentDir, readCurrent} from '../utils/kg/build/write';
import {parseId, declId, docId, testId, epId, ticketId, stubName} from '../utils/kg/build/ids';
import {firstParagraph} from '../utils/kg/build/extract/homes';
import {createRequire} from 'module';
import {readGraph, makeReport} from '../utils/kg/report';
import {kg} from '../commands/kg';

const fixtures = path.join(path.dirname(fileURLToPath(import.meta.url)), 'fixtures', 'kg');
const KG_DIR = path.join('core', 'docs', 'knowledge-graph');
const system: TypeSystem = loadTypeSystem(path.join(fixtures, 'good', KG_DIR));
const BATCH = 'b-test';
const homesModule = createRequire(import.meta.url)('../utils/kg/homes.js');

function makeRepo(): string {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-kg-build-'));
  fs.cpSync(path.join(fixtures, 'good'), dir, {recursive: true});
  return dir;
}

function write(root: string, file: string, text: string): void {
  fs.mkdirSync(path.dirname(path.join(root, file)), {recursive: true});
  fs.writeFileSync(path.join(root, file), text);
}

async function run(argv: Record<string, unknown>): Promise<{ok: boolean, out: string[], err: string[], exitCode: number | undefined}> {
  const log = vi.spyOn(console, 'log').mockImplementation(() => {});
  const error = vi.spyOn(console, 'error').mockImplementation(() => {});
  const before = process.exitCode;
  try {
    const ok = await kg(argv);
    return {ok, out: log.mock.calls.map((c) => String(c[0])), err: error.mock.calls.map((c) => String(c[0])), exitCode: process.exitCode as number | undefined};
  } finally {
    process.exitCode = before;
    log.mockRestore();
    error.mockRestore();
  }
}

/** Builds the homes layer of a fixture copy and returns the repo, the manifest and a JSONL reader. */
async function build(repo = makeRepo(), extra: Record<string, unknown> = {}): Promise<{repo: string, manifest: any, rows: (file: string) => any[], out: string, root: string}> {
  const root = typeof extra.out === 'string' ? extra.out : path.join(repo, ...(extra.public ? ['public', '.kg'] : ['.kg']));
  const result = await run({_: ['kg', 'build'], kg: path.join(repo, KG_DIR), only: 'homes', db: false, output: 'json', ...extra});
  expect(result.err).toEqual([]);
  expect(result.exitCode).toBeUndefined();
  const out = currentDir(root)!;
  const rows = (file: string) => {
    const p = path.join(out, file.startsWith('reports/') ? file : `data/${file}`) + (file.endsWith('.jsonl') ? '' : '.jsonl');
    return fs.existsSync(p) ? fs.readFileSync(p, 'utf8').split('\n').filter(Boolean).map((l) => JSON.parse(l)) : [];
  };
  return {repo, manifest: JSON.parse(result.out[0]), rows, out, root};
}

function node(graph: Graph, id: string): any {
  return graph.nodes.find((n) => n.id === id);
}

function person(id: string, extra: Record<string, unknown> = {}): Record<string, unknown> {
  return {type: 'person', id, name: id, email: 'x@example.com', provenance: 'annotation', source_layer: 'core', ...extra};
}

describe('normalizeRow (shared by check and build)', () => {
  it('drops nulls, applies defaults, canonicalizes dates, deduplicates lists and coerces a scalar into a list', () => {
    const {row, problems, defaulted} = normalizeRow(system, {type: 'person', id: 'P:a', name: 'A', email: 'a@x', slack: null, aliases: 'P:old', company: 'Team:core'});
    expect(problems).toEqual([]);
    expect(row).toEqual({type: 'person', id: 'P:a', name: 'A', email: 'a@x', aliases: ['P:old'], company: 'Team:core', status: 'active', tier: 'staff'});
    expect(defaulted).toEqual(['status', 'tier']);
    const ticket = normalizeRow(system, {type: 'ticket', id: 'GROK-1', name: 'x', tracker: 'jira', key: 'GROK-1', created: new Date(Date.UTC(2026, 8, 13)), updated: new Date(Date.UTC(2026, 8, 13, 10, 30))});
    expect(ticket.row).toMatchObject({created: '2026-09-13', updated: '2026-09-13T10:30:00.000Z'});
    expect(normalizeRow(system, {type: 'person', id: 'P:a', name: 'A', email: 'e', aliases: ['a', 'a', 'b']}).row.aliases).toEqual(['a', 'b']);
  });

  it('trims string and Text members, so a folded YAML scalar and a markdown paragraph read the same', () => {
    const {row} = normalizeRow(system, {type: 'concept', id: 'C:x', name: ' X ', description: 'A folded\nscalar.\n', aliases: [' C:y ']});
    expect(row).toMatchObject({name: 'X', description: 'A folded\nscalar.', aliases: ['C:y']});
  });

  it('reports shape problems with the check vocabulary and never the missing required members', () => {
    const {problems} = normalizeRow(system, {type: 'person', id: 'P:a', tier: 'boss', colour: 'red'});
    expect(problems).toEqual([
      {key: 'tier', code: 'bad-value', message: expect.stringMatching(/^tier: "boss" is not one of/)},
      {key: 'colour', code: 'unknown-key', message: "unknown key 'colour' for type person"},
    ]);
    expect(normalizeRow(system, {type: 'actor', id: 'x'}).problems).toEqual([{key: 'type', code: 'bad-type', message: "type 'actor' is abstract"}]);
    expect(normalizeRow(system, {type: 'nope', id: 'x'}).problems[0].message).toBe("unknown type 'nope'");
  });
});

describe('emitter merge rules (build-plan.md WO-1)', () => {
  it('resolves conflicting scalars by provenance rank whatever the order, fills missing fields and unions lists', () => {
    for (const order of [['annotation', 'filesystem'], ['filesystem', 'annotation']]) {
      const e = new Emitter(system, BATCH);
      for (const provenance of order)
        e.node(person('P:a', {name: provenance === 'annotation' ? 'Authored' : 'Derived', provenance, aliases: [`P:${provenance}`], ...(provenance === 'filesystem' ? {slack: 'a'} : {})}));
      const a = node(e.finalize(), 'P:a');
      expect(a).toMatchObject({name: 'Authored', slack: 'a', provenance: 'annotation', batch: BATCH});
      expect([...a.aliases].sort()).toEqual(['P:annotation', 'P:filesystem']);
    }
  });

  it('keeps the subtype when one type is an ancestor of the other, keeps the first of two unrelated types', () => {
    const e = new Emitter(system, BATCH);
    e.node(person('P:a'));
    e.node(person('P:a', {type: 'developer', company: 'Team:core', bitbucket: 'a', provenance: 'filesystem'}));
    e.node({type: 'team', id: 'Team:core', name: 'Core', provenance: 'annotation', source_layer: 'core'});
    e.node({type: 'ticket', id: 'P:a', name: 'x', tracker: 'jira', key: 'x', provenance: 'external', source_layer: 'process'});
    const graph = e.finalize();
    expect(node(graph, 'P:a')).toMatchObject({type: 'developer', company: 'Team:core', bitbucket: 'a', provenance: 'annotation'});
    expect(graph.problems.invalid_rows).toBe(0);
    expect(graph.details.duplicate_ids).toEqual(['P:a: ticket at external ignored; developer at annotation kept']);
  });

  it('keeps the higher confidence of a duplicate edge and unions evidence up to 20 paths', () => {
    const e = new Emitter(system, BATCH);
    e.node(person('P:a'));
    e.node({type: 'team', id: 'Team:t', name: 'T', provenance: 'annotation', source_layer: 'core', lead: 'P:a'});
    e.node({type: 'feature', id: 'visualize', name: 'V', provenance: 'annotation', source_layer: 'core'});
    e.node({type: 'concept', id: 'C:c', name: 'c', provenance: 'annotation', source_layer: 'core'});
    e.edge({type: 'uses-concept', from: 'visualize', to: 'C:c', derived_by: 'annotation', confidence: 0.5, evidence: ['core/docs/zz.md'], role: 'supporting'});
    for (let i = 0; i < 25; i++)
      e.edge({type: 'uses-concept', from: 'visualize', to: 'C:c', derived_by: 'annotation', confidence: i === 3 ? 1 : 0.5, evidence: [`core/docs/${i}.md`], role: i === 3 ? 'central' : 'supporting'});
    const graph = e.finalize();
    const uses = graph.edges.filter((x) => x.type === 'uses-concept');
    expect(uses).toHaveLength(1);
    expect(uses[0]).toMatchObject({confidence: 1, role: 'central', batch: BATCH});
    expect(uses[0].evidence).toHaveLength(20);
    expect(uses[0].evidence[0]).toBe('core/docs/0.md');
    expect(uses[0].evidence).not.toContain('core/docs/zz.md');
  });

  it('exempts a partial stub from required members, drops an incomplete real row into invalid.jsonl', () => {
    const e = new Emitter(system, BATCH);
    e.stub('GROK-9', 'ticket', 'GROK-9', 'annotation', {tracker: 'jira', bogus: 1, kind: 'unknown'});
    e.node({type: 'person', id: 'P:noemail', name: 'No email', provenance: 'annotation', source_layer: 'core'});
    e.node(person('P:ok'));
    const graph = e.finalize();
    expect(node(graph, 'GROK-9')).toEqual({id: 'GROK-9', type: 'ticket', name: 'GROK-9', tracker: 'jira', status: 'proposed', provenance: 'annotation', source_layer: 'synthetic', batch: BATCH, visibility: 'dev'});
    expect(node(graph, 'P:noemail')).toBeUndefined();
    expect(node(graph, 'P:ok')).toBeDefined();
    expect(graph.problems).toMatchObject({invalid_rows: 1, partial_stubs: 1});
    expect(graph.invalid[0]).toMatchObject({id: 'P:noemail', problems: ["missing required member 'email' for type person"]});
  });

  it('gives a stub only what created it: no defaults, source_layer from its path when it has one', () => {
    const e = new Emitter(system, BATCH);
    e.stub('TS:x', 'scenario', 'X', 'filesystem');
    e.stub('decl:public/a.ts#A', 'declaration', 'A', 'annotation', {path: 'public/a.ts', language: 'ts'});
    e.stub('P:ghost', 'person', 'Ghost', 'annotation');
    const graph = e.finalize();
    expect(node(graph, 'TS:x')).toEqual({id: 'TS:x', type: 'scenario', name: 'X', status: 'proposed', provenance: 'filesystem', source_layer: 'synthetic', batch: BATCH, visibility: 'dev'});
    expect(node(graph, 'decl:public/a.ts#A')).toEqual({id: 'decl:public/a.ts#A', type: 'declaration', name: 'A', path: 'public/a.ts', language: 'ts', status: 'proposed', provenance: 'annotation', source_layer: 'public', batch: BATCH, visibility: 'public'});
    expect(node(graph, 'P:ghost')).not.toHaveProperty('tier');
    e.node(person('P:ghost'));
    expect(node(e.finalize(), 'P:ghost')).toMatchObject({tier: 'staff', status: 'active'});
  });

  it('turns a real row for a stub id into a full node, the stub filling only the gaps', () => {
    const e = new Emitter(system, BATCH);
    e.stub('GROK-9', 'ticket', 'GROK-9', 'annotation', {tracker: 'jira', key: 'GROK-9'});
    e.node({type: 'ticket', id: 'GROK-9', name: 'Real title', tracker: 'jira', key: 'GROK-9', provenance: 'external', source_layer: 'process'});
    const graph = e.finalize();
    expect(node(graph, 'GROK-9')).toMatchObject({name: 'Real title', status: 'active', provenance: 'external', source_layer: 'process'});
    expect(graph.problems.partial_stubs).toBe(0);
  });
});

describe('emitter finalize (build-plan.md WO-1)', () => {
  it('derives part-of from id paths with title-cased stub parents, inherits owner and status, writes reference members as ref lines', () => {
    const e = new Emitter(system, BATCH);
    e.node(person('P:a'));
    e.node({type: 'feature', id: 'visualize/viewers', name: 'Viewers', owner: 'P:a', status: 'removed', provenance: 'annotation', source_layer: 'core'});
    e.node({type: 'feature', id: 'visualize/viewers/box-plot', name: 'Box plot', provenance: 'annotation', source_layer: 'core'});
    e.node({type: 'feature', id: 'visualize/viewers/kept', name: 'Kept', status: 'active', provenance: 'annotation', source_layer: 'core'});
    const graph = e.finalize();
    expect(node(graph, 'visualize')).toMatchObject({type: 'feature', name: 'Visualize', status: 'proposed', provenance: 'filesystem', source_layer: 'synthetic', visibility: 'public'});
    expect(node(graph, 'visualize/viewers/box-plot')).toMatchObject({owner: 'P:a', status: 'removed'});
    expect(node(graph, 'visualize/viewers/kept')).toMatchObject({owner: 'P:a', status: 'active'});
    expect(graph.edges.filter((x) => x.type === 'part-of').map((x) => `${x.from} -> ${x.to}`).sort()).toEqual([
      'visualize/viewers -> visualize', 'visualize/viewers/box-plot -> visualize/viewers', 'visualize/viewers/kept -> visualize/viewers',
    ]);
    expect(graph.edges.filter((x) => x.type === 'part-of')[0]).toMatchObject({derived_by: 'filesystem', confidence: 1, batch: BATCH});
    expect(graph.edges.filter((x) => x.type === 'ref' && x.name === 'owner').map((x) => x.from).sort()).toEqual(['visualize/viewers', 'visualize/viewers/box-plot', 'visualize/viewers/kept']);
    expect(graph.edges.find((x) => x.name === 'owner')).toMatchObject({type: 'ref', to: 'P:a', derived_by: 'annotation', confidence: 1});
  });

  it('drops an edge whose endpoints violate the type, stubs a missing target of an authored edge and drops an extracted one', () => {
    const e = new Emitter(system, BATCH);
    e.node({type: 'feature', id: 'visualize', name: 'V', provenance: 'annotation', source_layer: 'core'});
    e.node({type: 'feature', id: 'govern', name: 'G', provenance: 'annotation', source_layer: 'core'});
    e.edge({type: 'covers', from: 'visualize', to: 'govern', derived_by: 'annotation', confidence: 1});
    e.edge({type: 'tracked-in', from: 'visualize', to: 'GROK-77', derived_by: 'annotation', confidence: 1});
    e.edge({type: 'uses-concept', from: 'visualize', to: 'C:nowhere', derived_by: 'ast', confidence: 0.9});
    e.edge({type: 'supersedes', from: 'visualize', to: 'govern/legacy', derived_by: 'annotation', confidence: 1});
    const graph = e.finalize();
    expect(graph.edges.map((x) => `${x.type} ${x.from} -> ${x.to}`).sort()).toEqual(['part-of govern/legacy -> govern', 'supersedes visualize -> govern/legacy', 'tracked-in visualize -> GROK-77']);
    expect(node(graph, 'GROK-77')).toMatchObject({type: 'ticket', status: 'proposed', provenance: 'annotation'});
    expect(node(graph, 'govern/legacy')).toMatchObject({type: 'feature', name: 'Legacy', status: 'proposed'});
    expect(graph.problems).toMatchObject({dangling_edges: 2, partial_stubs: 2});
  });

  it('settles visibility: the type default narrowed by the home, a path narrowed by its location, internal as the ceiling', () => {
    const e = new Emitter(system, BATCH);
    e.node({type: 'feature', id: 'visualize', name: 'V', provenance: 'annotation', source_layer: 'core', home: 'core/docs/v.md'});
    e.node({type: 'feature', id: 'govern', name: 'G', visibility: 'dev', provenance: 'annotation', source_layer: 'core'});
    e.node({type: 'concept', id: 'C:c', name: 'c', visibility: 'internal', provenance: 'annotation', source_layer: 'core'});
    e.node(person('P:a', {visibility: 'public'}));
    e.node({type: 'source-file', id: 'file:core/x.dart', name: 'x.dart', path: 'core/x.dart', language: 'dart', provenance: 'filesystem', source_layer: 'core'});
    e.node({type: 'source-file', id: 'file:public/x.ts', name: 'x.ts', path: 'public/x.ts', language: 'ts', provenance: 'filesystem', source_layer: 'public'});
    const graph = e.finalize();
    expect(Object.fromEntries(graph.nodes.map((n) => [n.id, n.visibility]))).toEqual({
      visualize: 'public', govern: 'dev', 'C:c': 'internal', 'P:a': 'internal', 'file:core/x.dart': 'dev', 'file:public/x.ts': 'public',
    });
  });
});

describe('ids (build-plan.md "Common contracts")', () => {
  it('constructs and parses the extracted id forms', () => {
    expect(declId('public\\js-api\\src\\viewer.ts', 'Viewer.root', 'get')).toBe('decl:public/js-api/src/viewer.ts#Viewer.root:get');
    expect(parseId('decl:public/js-api/src/viewer.ts#Viewer.root:get')).toMatchObject({form: 'scheme', scheme: 'decl', path: 'public/js-api/src/viewer.ts', anchor: 'Viewer.root:get'});
    expect(parseId(docId('public/help/a.md', 'usage'))).toMatchObject({scheme: 'doc', path: 'public/help/a.md', anchor: 'usage'});
    expect(parseId(testId('dg', 'public/packages/Chem/src/tests/a.ts', 'Chem', 'smiles'))).toMatchObject({scheme: 'test', parts: ['dg', 'public/packages/Chem/src/tests/a.ts'], anchor: 'Chem/smiles'});
    expect(epId('get', '/spaces/{id}')).toBe('ep:GET /spaces/{id}');
    expect(parseId('C:dataframe')).toEqual({form: 'prefix', prefix: 'C', local: 'dataframe'});
    expect(parseId('GROK-12')).toMatchObject({form: 'tracker', tracker: 'jira'});
    expect(parseId(ticketId('#4062'))).toMatchObject({form: 'tracker', tracker: 'github', anchor: '4062'});
    expect(parseId('tutorial:chem/activity-cliffs')).toMatchObject({form: 'scheme', scheme: 'tutorial', parts: ['chem/activity-cliffs']});
    expect(parseId('visualize/viewers')).toEqual({form: 'bare', local: 'visualize/viewers'});
    expect([stubName('visualize/scatter-plot'), stubName('decl:a/b.ts#X'), stubName('doc:a/b.md'), stubName('func:Chem:detect'), stubName('GROK-1')]).toEqual(['Scatter plot', 'X', 'b.md', 'detect', 'GROK-1']);
  });

  it('hashes the batch id from every effective input, and from nothing else', () => {
    const inputs = {reddata: '1', public: '2', dirty: 'd', schema_version: '1', builder: '6.5.10', mode: 'full',
      extractors: 'docs,homes', backlog: '2026-01-12T07:00:00Z:abc', dart: 'def'};
    const a = batchId(inputs);
    expect(a).toMatch(/^b-[0-9a-f]{12}$/);
    expect(batchId({...inputs})).toBe(a);
    // the same inputs in another key order are the same batch; a change in any one of them is not
    expect(batchId(Object.fromEntries(Object.entries(inputs).reverse()))).toBe(a);
    for (const key of Object.keys(inputs))
      expect(batchId({...inputs, [key]: 'other'}), key).not.toBe(a);
  });

  it('takes the extractor selection unordered and the backlog watermark into the inputs', () => {
    const repo = makeRepo();
    const of = (extra: Record<string, unknown>) => buildInputs({repoRoot: repo, mode: 'full', schemaVersion: 1, builder: '6.5.10',
      revisions: {reddata: '1', public: '2'}, extractors: ['homes', 'docs'], ...extra});
    expect(of({}).extractors).toBe('docs,homes');
    expect(of({extractors: ['docs', 'homes']})).toEqual(of({}));
    expect(of({extractors: ['docs']})).not.toEqual(of({}));
    const backlog = path.join(repo, 'backlog');
    fs.mkdirSync(backlog, {recursive: true});
    fs.writeFileSync(path.join(backlog, 'index.jsonl'), `${JSON.stringify({id: 'GROK-1', updated: '2026-01-01T00:00:00Z'})}\n`);
    const first = of({backlogDir: backlog});
    expect(first.backlog).toMatch(/^2026-01-01T00:00:00Z:[0-9a-f]{40}$/);
    fs.writeFileSync(path.join(backlog, 'index.jsonl'), `${JSON.stringify({id: 'GROK-1', updated: '2026-02-02T00:00:00Z'})}\n`);
    expect(of({backlogDir: backlog}).backlog).not.toBe(first.backlog);
    expect(batchId(of({backlogDir: backlog}))).not.toBe(batchId(first));
  });
});

describe('grok kg build: writer, manifest and public projection (build-plan.md WO-1)', () => {
  it('writes JSONL sorted with fixed key order, byte-identical across two builds; the manifest has the schema.yaml shape', async () => {
    const first = await build();
    const snapshot = (dir: string) => Object.fromEntries(['nodes', 'edges'].flatMap((d) => fs.readdirSync(path.join(dir, 'data', d)).map((f) => [`${d}/${f}`, fs.readFileSync(path.join(dir, 'data', d, f), 'utf8')])));
    const before = snapshot(first.out);
    const second = await build(first.repo);
    expect(second.out).not.toBe(first.out);
    expect(snapshot(second.out)).toEqual(before);
    expect(second.manifest.batch).toBe(first.manifest.batch);
    expect(first.manifest).toMatchObject({schema_version: 1, builder: expect.stringMatching(/^\d+\.\d+\.\d+/), batch: expect.stringMatching(/^b-[0-9a-f]{12}$/), mode: 'full',
      revisions: {reddata: expect.any(String), public: expect.any(String), public_pin: expect.any(String)}, sources: {homes: 'ok'},
      problems: {invalid_rows: 0, dangling_edges: 0, unresolved_ids: 0, ambiguous_owners: 0, orphans: 0}});
    expect(first.manifest.built_at).toMatch(/Z$/);
    const features = before['nodes/feature.jsonl'].split('\n').filter(Boolean);
    expect(features.map((l) => JSON.parse(l).id)).toEqual([...features.map((l) => JSON.parse(l).id)].sort());
    expect(features[0]).toMatch(/^\{"id":"platform","type":"feature","name":"Platform","batch":/);
    expect(Object.keys(JSON.parse(features[1])).slice(3)).toEqual([...Object.keys(JSON.parse(features[1])).slice(3)].sort());
    expect(before['edges/part-of.jsonl']).not.toContain('built_at');
  });

  it('--out redirects the output root, --only refuses an unknown extractor, --output table prints one line', async () => {
    const repo = makeRepo();
    const out = path.join(repo, 'elsewhere');
    const {manifest} = await build(repo, {out});
    expect(fs.existsSync(path.join(out, 'gen', readCurrent(out)!, 'manifest.json'))).toBe(true);
    expect(fs.readFileSync(path.join(out, 'current'), 'utf8').trim()).toMatch(new RegExp(`^${manifest.batch}-`));
    expect(fs.existsSync(path.join(repo, '.kg'))).toBe(false);
    const bad = await run({_: ['kg', 'build'], kg: path.join(repo, KG_DIR), only: 'homes,nope', db: false});
    expect(bad.exitCode).toBe(1);
    expect(bad.err).toEqual(['--only names unknown extractors: nope (known: homes, ts-packages, ts-functions, ts-declarations, ts-imports, ts-uses, ts-tests, ts-samples, ts-changelog, docs, ts-markers, dart, process, membership)']);
    const table = await run({_: ['kg', 'build'], kg: path.join(repo, KG_DIR), only: 'homes', db: false, out});
    expect(table.out).toHaveLength(1);
    expect(table.out[0]).toMatch(/^wrote .*elsewhere\/gen\/b-[0-9a-f]{12}-\w{6}: \d+ nodes \(concept 2, .*feature 7.*\), \d+ edges \(.*part-of 5.*\); sources: homes ok; problems: partial_stubs \d+; batch b-[0-9a-f]{12} \(full\)$/);
  });

  it('--public keeps public logical nodes without home or owner, no people, no files, and only edges with both ends public', async () => {
    const {rows, manifest} = await build(makeRepo(), {public: true});
    expect(manifest.mode).toBe('public');
    expect(Object.keys(manifest.counts.nodes)).toEqual(['concept', 'doc-page', 'feature', 'scenario']);
    const features = rows('nodes/feature');
    expect(features.map((f) => f.id)).toEqual(['platform', 'platform/caching', 'platform/caching/invalidation', 'visualize', 'visualize/viewers', 'visualize/viewers/old-scatter', 'visualize/viewers/scatter-plot']);
    expect(features.some((f) => 'home' in f || 'owner' in f)).toBe(false);
    expect(features.find((f) => f.id === 'visualize/viewers/scatter-plot')).toMatchObject({description: 'Points in two dimensions.', visibility: 'public'});
    expect(rows('nodes/person')).toEqual([]);
    expect(rows('nodes/source-file')).toEqual([]);
    expect(rows('nodes/doc-page').map((d) => d.id)).toEqual(['doc:public/help/visualize/viewers/scatter-plot-tips.md']);
    expect(rows('edges/owner')).toEqual([]);
    expect(rows('edges/is-implemented-in')).toEqual([]);
    expect(rows('edges/tracked-in')).toEqual([]);
    expect(rows('edges/covers')).toHaveLength(1);
    // the tips page declares `documents:` and the home cites it: one edge, both pages as evidence
    expect(rows('edges/documents')[0].evidence).toEqual(['public/help/visualize/viewers/scatter-plot-tips.md', 'public/help/visualize/viewers/scatter-plot.md']);
    expect(rows('edges/mentions')).toEqual([]);
    expect(rows('edges/uses-concept')).toHaveLength(2);
  });

  it('--public writes nothing private: no reports, only the public revision, no core/ path and no reddata sha anywhere', async () => {
    const repo = makeRepo();
    const git = (cwd: string, ...args: string[]) => spawnSync('git', ['-c', 'user.name=t', '-c', 'user.email=t@x', ...args], {cwd, encoding: 'utf8'});
    for (const dir of [path.join(repo, 'public'), repo]) {
      git(dir, 'init', '-q');
      git(dir, 'add', '-A');
      git(dir, 'commit', '-q', '-m', 'fixture');
    }
    const reddata = git(repo, 'rev-parse', 'HEAD').stdout.trim();
    const pub = git(path.join(repo, 'public'), 'rev-parse', 'HEAD').stdout.trim();
    expect(reddata).toMatch(/^[0-9a-f]{40}$/);
    expect(pub).not.toBe(reddata);
    const {manifest, out} = await build(repo, {public: true});
    expect(manifest.mode).toBe('public');
    expect(manifest.revisions).toEqual({public: pub});
    expect(manifest.problems.partial_stubs).toBe(2);
    expect(Object.entries(manifest.problems).filter(([k]) => k !== 'partial_stubs').every(([, n]) => n === 0)).toBe(true);
    expect(fs.existsSync(path.join(out, 'reports'))).toBe(false);
    const files = fs.readdirSync(out, {recursive: true, withFileTypes: true}).filter((d) => d.isFile()).map((d) => path.join(d.parentPath ?? d.path, d.name));
    expect(files.length).toBeGreaterThan(5);
    for (const file of files) {
      const text = fs.readFileSync(file, 'utf8');
      expect(text, file).not.toContain('core/');
      expect(text, file).not.toContain(reddata);
    }
  });

  it('writes the public projection under public/.kg by default', async () => {
    const {repo, out} = await build(makeRepo(), {public: true});
    expect(out.startsWith(path.join(repo, 'public', '.kg'))).toBe(true);
    expect(fs.existsSync(path.join(out, 'manifest.json'))).toBe(true);
    expect(fs.existsSync(path.join(repo, '.kg'))).toBe(false);
  });
});

describe('homes extractor (build-plan.md WO-2)', () => {
  it('emits one node per home with its members, the description, the home path and the layer', async () => {
    const {rows} = await build();
    const features = rows('nodes/feature');
    expect(features.find((f) => f.id === 'visualize/viewers/scatter-plot')).toEqual({
      id: 'visualize/viewers/scatter-plot', type: 'feature', name: 'Scatter plot', aliases: ['visualize/scatterplot'], batch: expect.any(String),
      description: 'Points in two dimensions.', home: 'public/help/visualize/viewers/scatter-plot.md', owner: 'P:askalkin', provenance: 'annotation',
      source_layer: 'public', status: 'active', visibility: 'public',
    });
    expect(features.find((f) => f.id === 'platform/caching')).toMatchObject({name: 'Caching', description: 'The area home; `platform` itself has none, so it is a stub.', home: 'core/docs/CACHING.md', owner: 'P:askalkin', source_layer: 'core', visibility: 'public'});
    expect(features.find((f) => f.id === 'platform')).toMatchObject({name: 'Platform', status: 'proposed', provenance: 'filesystem', source_layer: 'synthetic'});
    expect(features.find((f) => f.id === 'visualize/viewers/old-scatter')).toMatchObject({status: 'removed', owner: 'P:askalkin'});
    expect(rows('nodes/concept').find((c) => c.id === 'C:dataframe')).toMatchObject({name: 'DataFrame', area: 'data', description: expect.stringMatching(/^An in-memory columnar table/), home: 'core/docs/knowledge-graph/concepts/dataframe.yaml', source_layer: 'core', visibility: 'public', provenance: 'annotation'});
    expect(rows('nodes/developer')).toEqual([expect.objectContaining({id: 'P:askalkin', type: 'developer', name: 'askalkin', email: 'a@example.com', company: 'Team:core', bitbucket: 'askalkin', areas: ['visualize/viewers'], tier: 'staff', visibility: 'internal', home: 'core/docs/knowledge-graph/internal/people/askalkin.yaml'})]);
    expect(rows('nodes/team')).toEqual([expect.objectContaining({id: 'Team:core', lead: 'P:askalkin', description: 'The core team.'})]);
  });

  it('reads a migrated scenario: path from the file, manual_only derived, covers to the feature', async () => {
    const {rows} = await build();
    expect(rows('nodes/scenario')).toEqual([expect.objectContaining({
      id: 'TS:scatter-plot-ui', type: 'scenario', name: 'Scatter plot manual checks', path: 'public/packages/UsageAnalysis/files/TestTrack/Viewers/scatter-plot-ui.md',
      manual_only: true, coverage_type: 'smoke', target_layer: 'manual-only', description: 'Migrated scenario: a home although it sits in the Test Track folder.', visibility: 'public',
    })]);
    expect(rows('edges/covers')).toEqual([{type: 'covers', from: 'TS:scatter-plot-ui', to: 'visualize/viewers/scatter-plot', batch: expect.any(String), confidence: 1, derived_by: 'annotation',
      evidence: ['public/packages/UsageAnalysis/files/TestTrack/Viewers/scatter-plot-ui.md'], level: 'exercised'}]);
  });

  it('reads an annotated page: a doc-page stub and documents edges with their properties', async () => {
    const {rows} = await build();
    expect(rows('nodes/doc-page').find((d) => d.id === 'doc:public/help/visualize/viewers/scatter-plot-tips.md')).toMatchObject({
      type: 'doc-page', name: 'Scatter plot tips', path: 'public/help/visualize/viewers/scatter-plot-tips.md', kind: 'help', status: 'proposed', provenance: 'annotation', source_layer: 'public', visibility: 'public',
    });
    expect(rows('edges/documents')).toEqual([expect.objectContaining({from: 'doc:public/help/visualize/viewers/scatter-plot-tips.md', to: 'visualize/viewers/scatter-plot', audience: 'user', derived_by: 'annotation'})]);
  });

  it('claims code: roots at rung 2 after glob expansion into source-file nodes, cited files at rung 3, a path#Anchor as a declaration edge', async () => {
    const repo = makeRepo();
    write(repo, 'core/docs/viewers/README.md', fs.readFileSync(path.join(repo, 'core/docs/viewers/README.md'), 'utf8')
      .replace('tickets:', 'code:\n  - public/js-api/src/viewer.ts#ScatterPlotViewer\n  - {path: core/client/d4/lib, role: definition}\ntickets:'));
    const {rows} = await build(repo);
    const claims = rows('reports/claims.jsonl');
    expect(claims).toEqual([
      {feature: 'visualize/viewers', file: 'core/client/d4/lib/scatter.dart', line: 4, props: {role: 'definition'}, rung: 2, source: 'home'},
      {feature: 'visualize/viewers/scatter-plot', file: 'core/client/d4/lib/scatter.dart', line: 12, props: {role: 'ui'}, rung: 2, source: 'home'},
      {feature: 'visualize/viewers/scatter-plot', file: 'public/js-api/src/viewer.ts', line: 12, props: {}, rung: 2, source: 'home'},
    ]);
    expect(rows('nodes/source-file')).toEqual([
      expect.objectContaining({id: 'file:core/client/d4/lib/scatter.dart', name: 'scatter.dart', path: 'core/client/d4/lib/scatter.dart', language: 'dart', loc: 1, provenance: 'filesystem', source_layer: 'core', visibility: 'dev'}),
      expect.objectContaining({id: 'file:public/js-api/src/viewer.ts', language: 'ts', source_layer: 'public', visibility: 'public'}),
    ]);
    expect(rows('edges/is-implemented-in')).toEqual([expect.objectContaining({from: 'visualize/viewers', to: 'decl:public/js-api/src/viewer.ts#ScatterPlotViewer', derived_by: 'annotation', confidence: 1})]);
    expect(rows('nodes/declaration').find((d) => d.id === 'decl:public/js-api/src/viewer.ts#ScatterPlotViewer')).toMatchObject({name: 'ScatterPlotViewer', language: 'ts', path: 'public/js-api/src/viewer.ts', status: 'proposed'});
  });

  it('claims a cited implementation file at rung 3 unless a code: root already covers it, and takes a cited help page as documentation', async () => {
    const {rows} = await build();
    expect(rows('reports/claims.jsonl')).toEqual([
      {feature: 'visualize/viewers', file: 'core/client/d4/lib', line: 10, props: {}, rung: 3, source: 'home'},
      expect.objectContaining({feature: 'visualize/viewers/scatter-plot', file: 'core/client/d4/lib/scatter.dart', rung: 2}),
      expect.objectContaining({feature: 'visualize/viewers/scatter-plot', file: 'public/js-api/src/viewer.ts', rung: 2}),
    ]);
    // a help page the home's body cites documents the feature (conventions.md §5.4), so it is no mention
    expect(rows('edges/mentions')).toEqual([]);
    expect(rows('edges/documents')).toEqual([expect.objectContaining({from: 'doc:public/help/visualize/viewers/scatter-plot-tips.md', to: 'visualize/viewers/scatter-plot',
      derived_by: 'annotation', evidence: ['public/help/visualize/viewers/scatter-plot-tips.md', 'public/help/visualize/viewers/scatter-plot.md']})]);
    expect(rows('nodes/doc-page').map((d) => d.id)).toEqual(['doc:public/help/visualize/viewers/scatter-plot-tips.md']);
  });

  it('stubs tickets from tickets: with the tracker from the key shape and tracked-in edges with their properties', async () => {
    const {rows} = await build();
    expect(rows('nodes/ticket')).toEqual([
      expect.objectContaining({id: 'GROK-1', type: 'ticket', tracker: 'jira', key: 'GROK-1', status: 'proposed', provenance: 'annotation', source_layer: 'synthetic'}),
      expect.objectContaining({id: 'GROK-20863', tracker: 'jira', key: 'GROK-20863'}),
    ]);
    expect(rows('edges/tracked-in')).toEqual([
      expect.objectContaining({from: 'visualize/viewers', to: 'GROK-1', relation: 'epic'}),
      expect.objectContaining({from: 'visualize/viewers', to: 'GROK-20863', relation: 'other'}),
    ]);
  });

  it('reads defined_by on the concept as defines-concept edges from declaration stubs, and superseded_by through an alias', async () => {
    const {rows} = await build();
    expect(rows('edges/defines-concept')).toEqual([
      expect.objectContaining({from: 'decl:core/client/d4/lib/scatter.dart#DataFrame', to: 'C:dataframe', derived_by: 'annotation'}),
      expect.objectContaining({from: 'decl:public/js-api/src/viewer.ts#DataFrame', to: 'C:dataframe', language: 'ts'}),
    ]);
    expect(rows('nodes/declaration').map((d) => [d.id, d.language, d.path, d.status])).toEqual([
      ['decl:core/client/d4/lib/scatter.dart#DataFrame', 'dart', 'core/client/d4/lib/scatter.dart', 'proposed'],
      ['decl:public/js-api/src/viewer.ts#DataFrame', 'ts', 'public/js-api/src/viewer.ts', 'proposed'],
    ]);
    expect(rows('edges/supersedes')).toEqual([expect.objectContaining({from: 'visualize/viewers/scatter-plot', to: 'visualize/viewers/old-scatter', release: '1.20'})]);
    expect(rows('edges/uses-concept')).toEqual([
      expect.objectContaining({from: 'visualize/viewers/scatter-plot', to: 'C:column', role: 'central'}),
      expect.objectContaining({from: 'visualize/viewers/scatter-plot', to: 'C:dataframe'}),
    ]);
    expect(rows('edges/company')).toEqual([expect.objectContaining({type: 'ref', name: 'company', from: 'P:askalkin', to: 'Team:core'})]);
    expect(rows('edges/areas')).toEqual([expect.objectContaining({from: 'P:askalkin', to: 'visualize/viewers'})]);
    expect(rows('edges/lead')).toEqual([expect.objectContaining({from: 'Team:core', to: 'P:askalkin'})]);
    expect(rows('edges/part-of').map((e) => `${e.from} -> ${e.to}`)).toContain('visualize/viewers/scatter-plot -> visualize/viewers');
    expect(rows('edges/owner').map((e) => e.from).sort()).toEqual(['platform/caching', 'platform/caching/invalidation', 'visualize/viewers', 'visualize/viewers/old-scatter', 'visualize/viewers/scatter-plot']);
  });

  it('marks the homes source partial when check reports errors on a home, and still builds the rest', async () => {
    const repo = makeRepo();
    write(repo, 'core/docs/broken.md', '---\nfeature: govern/permissions\nowner: askalkin\nstatus: bogus\n---\n# P\n');
    const {rows, manifest} = await build(repo);
    // the count comes from the emitter after finalize, so a source cannot call itself ok while its rows were refused
    expect(manifest.sources.homes).toBe('partial(1 rejected)');
    expect(manifest.problems.invalid_rows).toBe(1);
    expect(rows('nodes/feature').some((f) => f.id === 'govern/permissions')).toBe(false);
    expect(rows('reports/invalid.jsonl')).toEqual([expect.objectContaining({id: 'govern/permissions', problems: [expect.stringMatching(/^status: "bogus" is not one of/)]})]);
  });

  it('types a dangling bare id by the expected union, never by sort order, and refuses an ambiguous one', async () => {
    const repo = makeRepo();
    write(repo, 'core/docs/haunted.md', '---\nfeature: govern/haunted\nowner: nobody\nsuperseded_by: [platform/ghost]\n---\n# Haunted\n');
    const {rows, manifest, out} = await build(repo);
    expect(rows('nodes/feature').find((f) => f.id === 'platform/ghost')).toMatchObject({type: 'feature', name: 'Ghost', status: 'proposed'});
    expect(rows('nodes/concept').map((c) => c.id)).toEqual(['C:column', 'C:dataframe']);
    expect(rows('edges/supersedes')).toContainEqual(expect.objectContaining({from: 'platform/ghost', to: 'govern/haunted'}));
    expect(rows('nodes/feature').find((f) => f.id === 'govern/haunted')).not.toHaveProperty('owner');
    expect(rows('edges/owner').some((e) => e.from === 'govern/haunted')).toBe(false);
    expect(manifest.problems.ambiguous_refs).toBe(1);
    expect(JSON.parse(fs.readFileSync(path.join(out, 'reports', 'problems.json'), 'utf8')).ambiguous_refs)
      .toEqual(["core/docs/haunted.md: owner: 'nobody' could be P:nobody or Team:nobody; write the prefix"]);
  });

  it('names every check issue behind a partial homes source in problems.json and counts them', async () => {
    const repo = makeRepo();
    write(repo, 'core/docs/broken.md', '---\nfeature: govern/permissions\nowner: askalkin\nstatus: bogus\n---\n# P\n');
    const {manifest, out} = await build(repo);
    expect(manifest.sources.homes).toBe('partial(1 rejected)');
    expect(manifest.problems.home_issues).toBe(1);
    expect(JSON.parse(fs.readFileSync(path.join(out, 'reports', 'problems.json'), 'utf8')).home_issues)
      .toEqual([expect.stringMatching(/^core\/docs\/broken\.md:4: bad-value: status: "bogus" is not one of/)]);
  });

  it('carries the citation issues to report stale as structured rows, target included', async () => {
    const repo = makeRepo();
    write(repo, 'core/docs/broken.md', '---\nfeature: govern/permissions\nowner: askalkin\n---\n# P\n\nSee `core/nope.dart` and [a](../nope.md).\n');
    const {out} = await build(repo);
    const issues = JSON.parse(fs.readFileSync(path.join(out, 'reports', 'home-issues.json'), 'utf8'));
    expect(issues.map((i: any) => [i.code, i.target])).toEqual([['missing-cited-path', 'core/nope.dart'], ['missing-doc-link', '../nope.md']]);
    const report = makeReport('stale', readGraph(out, repo, system, 'stale'), {system, repoRoot: repo});
    expect(report.sections[0].rows).toContainEqual({kind: 'citation', source: 'core/docs/broken.md:7', target: 'core/nope.dart',
      reason: "missing-cited-path: cited path 'core/nope.dart' does not exist"});
    expect(report.sections[0].rows).toContainEqual({kind: 'citation', source: 'core/docs/broken.md:7', target: '../nope.md',
      reason: "missing-doc-link: linked document '../nope.md' does not exist (resolved to core/nope.md)"});
  });

  it('loads the home documents once per build, whichever extractors ask for them', async () => {
    // the build runs on the transpiled .js, whose callers read the CommonJS exports object; a spy on the ESM namespace is not seen there
    const original = homesModule.loadHomes;
    let calls = 0;
    homesModule.loadHomes = (...args: Parameters<typeof original>) => {
      calls++;
      return original(...args);
    };
    try {
      await build(makeRepo(), {only: 'homes,docs,ts-tests'});
    }
    finally {
      homesModule.loadHomes = original;
    }
    expect(calls).toBe(1);
  });

  it('takes the first prose paragraph after the frontmatter as the description', () => {
    expect(firstParagraph('\n# Title\n\n<!-- note -->\nFirst line\ncontinues here.\n\nSecond paragraph.\n')).toBe('First line continues here.');
    expect(firstParagraph('| a | b |\n|---|---|\n\n```\ncode\n```\n')).toBeUndefined();
  });
});
