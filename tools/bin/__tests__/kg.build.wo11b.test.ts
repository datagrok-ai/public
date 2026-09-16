/// WO-11b, ingestion correctness (build-plan.md, kg-codex-review-3.md #3 to #8): what the build may assert
/// about a row it refused, about two assertions that differ only in a property, about a release record's
/// word, about a test whose name it cannot know, and about the files it never saw. Every case here
/// fails on the pipeline as the third review found it.
import {describe, it, expect} from 'vitest';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {fileURLToPath} from 'url';
import {loadTypeSystem, TypeSystem} from '../utils/kg/types';
import {Emitter} from '../utils/kg/build/emitter';
import {kg} from '../commands/kg';
import {copyFixture, buildFixture, write, Built} from './kg-fixture';

const fixture = path.join(path.dirname(fileURLToPath(import.meta.url)), 'fixtures', 'kg', 'build');
const KG_DIR = path.join('core', 'docs', 'knowledge-graph');
const system: TypeSystem = loadTypeSystem(path.join(fixture, KG_DIR));
const TESTED = 'public/packages/Tested/src';

function copy(prepare?: (repo: string) => void): string {
  return copyFixture('build', prepare);
}

function build(repo: string, only: string): Promise<Built> {
  return buildFixture(repo, only);
}





describe('admission: a refused row asserts nothing (review 3 #3)', () => {
  it('tells an extractor whether the node was taken, and names the reason', () => {
    const e = new Emitter(system, 'b-test');
    const at = (line: number, extra: Record<string, unknown> = {}) => ({type: 'function', id: 'func:Demo:X', name: 'X', language: 'ts',
      path: 'public/packages/Demo/src/package.g.ts', line, package: 'pkg:Demo', provenance: 'annotation', source_layer: 'public', ...extra});
    expect(e.node(at(1))).toEqual({accepted: true, id: 'func:Demo:X'});
    expect(e.node(at(9))).toEqual({accepted: false, id: 'func:Demo:X', reason: 'registration_collision'});
    expect(e.node({type: 'nonsense', id: 'func:Demo:Y'})).toMatchObject({accepted: false, reason: 'invalid'});
    const graph = e.finalize();
    expect(graph.problems.registration_collisions).toBe(1);
    expect(graph.problems.duplicate_ids).toBe(0);
    expect(graph.details.registration_collisions).toEqual(['func:Demo:X: function at public/packages/Demo/src/package.g.ts:9 ignored; function at public/packages/Demo/src/package.g.ts:1 kept']);
  });

  it('counts a cross-extractor merge that cannot be reconciled apart from a name collision', () => {
    const e = new Emitter(system, 'b-test');
    e.node({type: 'feature', id: 'platform/x', name: 'X', provenance: 'annotation', source_layer: 'core'});
    expect(e.node({type: 'ticket', id: 'platform/x', name: 'X', tracker: 'jira', key: 'x', kind: 'bug', state: 'open', provenance: 'external', source_layer: 'process'}))
      .toMatchObject({accepted: false, reason: 'duplicate_id'});
    const graph = e.finalize();
    expect(graph.problems.duplicate_ids).toBe(1);
    expect(graph.problems.registration_collisions).toBe(0);
  });

  it('keeps the semantic types of the registration it retained and none of the one it refused', async () => {
    const {rows, problems} = await build(copy(), 'ts-packages,ts-functions');
    // `Dual` is registered twice: the panel at 153 is kept, the app at 161 refused, and the app's input semType with it
    expect(rows('edges/targets-semtype').filter((e) => e.from === 'func:Demo:Dual')).toEqual([]);
    expect(rows('nodes/panel').find((r) => r.id === 'func:Demo:Dual')).toMatchObject({line: 153});
    expect(problems.registration_collisions).toContain('func:Demo:Dual: app at public/packages/Demo/src/package.g.ts:162 ignored; panel at public/packages/Demo/src/package.g.ts:153 kept');
  });
});

describe('edge identity and ordered lists (review 3 #4)', () => {
  it('keeps one edge per discriminator and merges only within it', () => {
    const e = new Emitter(system, 'b-test');
    e.node({type: 'function', id: 'func:Demo:F', name: 'F', language: 'ts', package: 'pkg:Demo', provenance: 'annotation', source_layer: 'public'});
    e.node({type: 'semantic-type', id: 'semtype:Molecule', name: 'Molecule', language: 'other', provenance: 'annotation', source_layer: 'public'});
    for (const role of ['consumes', 'produces', 'consumes'])
      e.edge({type: 'targets-semtype', from: 'func:Demo:F', to: 'semtype:Molecule', role, derived_by: 'annotation', confidence: 1});
    const edges = e.finalize().edges.filter((x) => x.type === 'targets-semtype');
    expect(edges.map((x) => x.role)).toEqual(['consumes', 'produces']);
  });

  it('keeps every parameter of a signature in order, however often its type repeats, while a set stays a set', () => {
    const e = new Emitter(system, 'b-test');
    const row = {type: 'function', id: 'func:Demo:F', name: 'F', language: 'ts', package: 'pkg:Demo', provenance: 'annotation', source_layer: 'public'};
    // four `list<string>` inputs are four facts about the signature, not one
    e.node({...row, input_types: ['dataframe', 'list<string>', 'list<string>', 'list<string>', 'list<string>'], tags: ['a', 'a', 'b']});
    e.node({...row, input_types: ['dataframe', 'column'], provenance: 'ast'});
    const f = e.finalize().nodes.find((n) => n.id === 'func:Demo:F')!;
    expect(f.input_types).toEqual(['dataframe', 'list<string>', 'list<string>', 'list<string>', 'list<string>']);
    expect(f.tags).toEqual(['a', 'b']);
  });

  it('takes a reference property\'s provenance from the field, not from the node that carries it', () => {
    const e = new Emitter(system, 'b-test');
    e.node({type: 'person', id: 'P:a', name: 'A', handle: 'a', provenance: 'annotation', source_layer: 'core'});
    e.node({type: 'feature', id: 'platform/y', name: 'Y', provenance: 'annotation', source_layer: 'core'});
    e.node({type: 'feature', id: 'platform/y', owner: 'P:a', provenance: 'llm', source_layer: 'core'});
    const ref = e.finalize().edges.find((x) => x.type === 'ref' && x.name === 'owner')!;
    expect(ref).toMatchObject({from: 'platform/y', to: 'P:a', derived_by: 'llm'});
  });
});

describe('a release record speaks for itself (review 3 #5)', () => {
  it('carries the dry run onto every edge the record produces, and never onto the backlog\'s own', async () => {
    const {rows} = await build(copy(), 'homes,process');
    expect(rows('edges/targets-release').filter((e) => e.evidence?.includes('core/docs/release/1.0.1.yaml')).every((e) => e.confidence === 0.7)).toBe(true);
    expect(rows('edges/targets-release').find((e) => e.derived_by === 'external')).toMatchObject({confidence: 1});
  });
});

describe('text extraction says what it cannot know (review 3 #7)', () => {
  it('counts a test whose name is built at run time and marks the row dynamic', async () => {
    const {rows, problems} = await build(copy(), 'homes,ts-tests');
    const dynamic = rows('nodes/test').filter((t) => t.dynamic);
    expect(dynamic.map((t) => t.name)).toEqual(['template…']);
    expect(problems.dynamic_tests).toEqual(['public/packages/Tested/src/tests/demo-tests.ts: Tested: Utils/template… names a registration site, not a runnable test: the title is built at run time']);
  });
});

describe('membership: the inputs and the denominator (review 3 #8)', () => {
  it('inherits a folder that belongs to one feature, however many roots that feature has', async () => {
    // VIEWERS.md has two code: roots; `build/` is outside what its globs expand to, so only rung 4 can reach this file
    const repo = copy((r) => write(r, `${TESTED}/build/helper.ts`, 'export const helper = 1;\n'));
    const {rows} = await build(repo, 'homes,ts-declarations,membership');
    expect(rows('edges/is-implemented-in').filter((e) => e.to === `file:${TESTED}/build/helper.ts`))
      .toEqual([expect.objectContaining({from: 'visualize/viewers', derived_by: 'filesystem', confidence: 0.9})]);
  });

  it('takes a `// ~id` line as participation, and counts one no home declares without inventing it', async () => {
    const repo = copy((r) => {
      write(r, `${TESTED}/marked.ts`, '// ~platform/caching\nexport const a = 1;\n');
      write(r, `${TESTED}/mistyped.ts`, '// ~platform/cashing\nexport const b = 2;\n');
    });
    const {rows, problems} = await build(repo, 'homes,ts-declarations,ts-markers,membership');
    expect(rows('edges/participates-in').filter((e) => e.from === `file:${TESTED}/marked.ts`))
      .toEqual([expect.objectContaining({to: 'platform/caching', derived_by: 'annotation'})]);
    // participation never takes ownership away from the folder
    expect(rows('edges/is-implemented-in').filter((e) => e.to === `file:${TESTED}/marked.ts`)).toEqual([expect.objectContaining({from: 'visualize/viewers'})]);
    expect(problems.unresolved_ids).toContain(`${TESTED}/mistyped.ts:1: ~platform/cashing resolves to no home document`);
    expect(rows('nodes/feature').some((f) => f.id === 'platform/cashing')).toBe(false);
  });

  it('resolves a function\'s //feature: through the home index instead of manufacturing a stub', async () => {
    const repo = copy((r) => {
      const file = path.join(r, ...'public/packages/Plain/src/package.ts'.split('/'));
      fs.appendFileSync(file, '\n//name: Mistyped\n//feature: platform/cashing\nexport function mistyped(): void {}\n');
    });
    const {rows, problems} = await build(repo, 'homes,ts-packages,ts-functions');
    expect(rows('nodes/feature').some((f) => f.id === 'platform/cashing')).toBe(false);
    expect(rows('edges/is-implemented-in').some((e) => e.from === 'platform/cashing')).toBe(false);
    expect(problems.unresolved_ids.some((p) => p.includes('~platform/cashing resolves to no home document'))).toBe(true);
  });

  it('publishes what it observed beside what it owns, so the orphan count has a denominator', async () => {
    const {report, manifest} = await build(copy(), 'homes,ts-packages,ts-functions,ts-declarations,ts-tests,membership');
    const ownership = report('ownership');
    expect(manifest.problems.orphans).toBe(ownership.orphans.length);
    expect(ownership.inventory.observed_files).toBe(ownership.inventory.owned_files + ownership.orphans.length);
    expect(ownership.inventory.observed_loc).toBeGreaterThan(ownership.inventory.owned_loc);
    expect(ownership.inventory.participating_files).toBeGreaterThan(0);
  });
});
