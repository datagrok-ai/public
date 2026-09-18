/// `grok kg build` WO-4 (build-plan.md): membership resolution over the mini monorepo under
/// fixtures/kg/build — which feature owns each file (conventions.md §8), which ones only
/// participate, the tests that follow an owned file, and reports/ownership.json.
import {describe, it, expect} from 'vitest';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {fileURLToPath} from 'url';
import {TypeSystem} from '../utils/kg/types';
import {Emitter, Graph} from '../utils/kg/build/emitter';
import {membershipExtractor} from '../utils/kg/build/extract/membership';
import {kg} from '../commands/kg';
import {copyFixture, buildFixture, kgRoot, fixtureTypes} from './kg-fixture';

const fixture = path.join(path.dirname(fileURLToPath(import.meta.url)), 'fixtures', 'kg', 'build');
const system: TypeSystem = fixtureTypes('build');
const SCATTER = 'core/client/d4/lib/src/viewers/scatterplot/scatter.dart';
const LEGEND = 'core/client/d4/lib/src/legends/legend.dart';
const RENDERER = 'core/client/d4/lib/src/legends/legend_renderer.dart';
const CACHE = 'core/client/d4/lib/src/viewers/legend_cache.dart';
const HISTOGRAM = 'core/client/d4/lib/src/viewers/histogram/histogram.dart';
const VIEWER = 'core/client/d4/lib/src/viewers/viewer.dart';
const RUN = 'public/packages/Demo/scripts/run.js';
const TESTS = 'public/packages/Tested/src/tests/demo-tests.ts';



async function build(): Promise<{rows: (file: string) => any[], ownership: any, manifest: any}> {
  const b = await buildFixture(copyFixture('build'), 'homes,ts-functions,ts-tests,membership');
  return {rows: b.rows, ownership: b.report('ownership'), manifest: b.manifest};
}
const graph = build();
const owners = (rows: any[], file: string) => rows.filter((e) => e.to === `file:${file}`);
const participants = (rows: any[], file: string) => rows.filter((e) => e.from === `file:${file}`).map((e) => e.to);

/** The claims of one file against the fixture's homes, resolved on their own. */
function resolveClaims(claims: any[], file = SCATTER): Graph {
  const repo = copyFixture('build');
  const emitter = new Emitter(system, 'b-test');
  for (const id of ['visualize/viewers', 'visualize/viewers/scatter-plot', 'platform/caching', 'visualize/legends'])
    emitter.node({type: 'feature', id, name: id, provenance: 'annotation', source_layer: 'core'});
  emitter.node({type: 'source-file', id: `file:${file}`, name: path.posix.basename(file), path: file, loc: 3, language: 'dart', provenance: 'filesystem', source_layer: 'core'});
  for (const claim of claims) emitter.claim({file, source: 'marker', props: {}, ...claim});
  membershipExtractor.run({system, kgRoot: kgRoot(repo), repoRoot: repo, mode: 'full'}, emitter);
  return emitter.finalize();
}

describe('membership resolution: the rungs of conventions.md §8 (build-plan.md WO-4)', () => {
  it('lets a marker beat the code: root that also claims the file, which then only participates', async () => {
    const {rows} = await graph;
    expect(owners(rows('edges/is-implemented-in'), RUN)).toEqual([expect.objectContaining({
      from: 'domains/bio', to: `file:${RUN}`, derived_by: 'annotation', confidence: 1, evidence: [RUN],
    })]);
    expect(participants(rows('edges/participates-in'), RUN)).toEqual(['platform/caching']);
  });

  it('resolves nested code: roots by the part-of chain, the leaf owning and the ancestor staying out of participates-in', async () => {
    const {rows, ownership} = await graph;
    expect(owners(rows('edges/is-implemented-in'), SCATTER)).toEqual([expect.objectContaining({
      from: 'visualize/viewers/scatter-plot', to: `file:${SCATTER}`, derived_by: 'annotation', confidence: 1, role: 'definition',
      evidence: ['core/client/d4/lib/src/viewers/scatterplot/CLAUDE.md'],
    })]);
    expect(participants(rows('edges/participates-in'), SCATTER)).toEqual([]);
    expect(ownership.resolved_by_chain).toContainEqual({file: SCATTER, owner: 'visualize/viewers/scatter-plot', over: ['visualize/viewers']});
  });

  // conventions.md §8: between unrelated rung-2 roots the more specific one wins, the exact path of CACHING.md over the
  // viewers glob; the other claimant participates and nothing is ambiguous
  it('gives a file two unrelated roots claim to the more specific root; the other participates', async () => {
    const {rows, ownership, manifest} = await graph;
    expect(owners(rows('edges/is-implemented-in'), CACHE)).toEqual([expect.objectContaining({from: 'platform/caching'})]);
    expect(participants(rows('edges/participates-in'), CACHE)).toEqual(['visualize/viewers']);
    expect(ownership.ambiguous).toEqual([]);
    expect(manifest.problems.ambiguous_owners ?? 0).toBe(0);
  });

  it('lets a citation own only the files whose nearest home it is, and participate in the rest', async () => {
    const {rows} = await graph;
    expect(owners(rows('edges/is-implemented-in'), LEGEND)).toEqual([expect.objectContaining({
      from: 'visualize/legends', derived_by: 'annotation', confidence: 1, evidence: ['core/client/d4/lib/src/legends/README.md'],
    })]);
    expect(owners(rows('edges/is-implemented-in'), 'public/packages/Demo/src/utils.ts')).toEqual([]);
    expect(participants(rows('edges/participates-in'), 'public/packages/Demo/src/utils.ts')).toEqual(['domains/bio']);
  });

  it('takes the folder of a home with no code: root of its own as a rung-2 root, the file it does not cite included', async () => {
    const {rows} = await graph;
    expect(owners(rows('edges/is-implemented-in'), RENDERER)).toEqual([expect.objectContaining({
      from: 'visualize/legends', to: `file:${RENDERER}`, derived_by: 'annotation', confidence: 1,
      evidence: ['core/client/d4/lib/src/legends/README.md'],
    })]);
    expect(participants(rows('edges/participates-in'), RENDERER)).toEqual(['visualize/viewers']);
  });

  it('lets a folder home inside the glob of an ancestor own its own folder, the ancestor keeping the rest', async () => {
    const {rows, ownership} = await graph;
    expect(owners(rows('edges/is-implemented-in'), HISTOGRAM)).toEqual([expect.objectContaining({
      from: 'visualize/viewers/histogram', to: `file:${HISTOGRAM}`, derived_by: 'annotation', confidence: 1,
      evidence: ['core/client/d4/lib/src/viewers/histogram/CLAUDE.md'],
    })]);
    expect(owners(rows('edges/is-implemented-in'), VIEWER).map((e: any) => e.from)).toEqual(['visualize/viewers']);
    expect(ownership.resolved_by_chain).toContainEqual({file: HISTOGRAM, owner: 'visualize/viewers/histogram', over: ['visualize/viewers']});
  });

  it('follows every test of an owned file to its owner', async () => {
    const {rows} = await graph;
    const tests = rows('edges/tests').filter((e) => e.derived_by === 'filesystem');
    expect(tests.every((e) => e.to === 'visualize/viewers' && e.confidence === 0.9 && e.evidence[0] === TESTS)).toBe(true);
    expect(tests).toHaveLength(rows('nodes/test').filter((t) => t.path === TESTS).length);
    expect(rows('edges/tests').some((e) => e.from.startsWith('test:playwright:') && e.derived_by === 'filesystem')).toBe(false);
  });

  it('counts a file no rung reaches as an orphan, sorted by loc, and only under the code roots', async () => {
    const {ownership, manifest} = await graph;
    expect(ownership.orphans[0]).toEqual({file: 'public/packages/Demo/detectors.js', loc: 69});
    expect(ownership.orphans.map((o: any) => o.loc)).toEqual([...ownership.orphans.map((o: any) => o.loc)].sort((a: number, b: number) => b - a));
    expect(ownership.orphans.map((o: any) => o.file)).toContain('public/packages/Demo/src/utils.ts');
    expect(ownership.orphans.every((o: any) => /^(core\/(client|server|shared)|public\/(packages|libraries|js-api))\//.test(o.file))).toBe(true);
    expect(manifest.problems.orphans).toBe(ownership.orphans.length);
  });
});

describe('membership resolution: claims that no extractor makes yet (build-plan.md WO-4, WO-6)', () => {
  it('takes an inline // ~id marker as participation only, leaving the code: root the owner', () => {
    const graph = resolveClaims([
      {feature: 'platform/caching', rung: 1, mode: 'participates'},
      {feature: 'visualize/viewers', rung: 2, source: 'home'},
    ]);
    expect(graph.edges.filter((e) => e.type === 'is-implemented-in')).toEqual([expect.objectContaining({from: 'visualize/viewers', to: `file:${SCATTER}`})]);
    expect(graph.edges.filter((e) => e.type === 'participates-in').map((e) => e.to)).toEqual(['platform/caching']);
  });

  it('inherits the folder of the nearest home at rung 4 when no root claimed the file, with filesystem provenance', () => {
    const graph = resolveClaims([], RENDERER);
    expect(graph.edges.filter((e) => e.type === 'is-implemented-in')).toEqual([expect.objectContaining({
      from: 'visualize/legends', to: `file:${RENDERER}`, derived_by: 'filesystem', confidence: 0.9,
    })]);
  });

  it('refuses to pick between two markers on one file, reporting the ambiguity at rung 1', () => {
    const graph = resolveClaims([{feature: 'platform/caching', rung: 1}, {feature: 'visualize/viewers', rung: 1}]);
    expect(graph.edges.filter((e) => e.type === 'is-implemented-in')).toEqual([]);
    expect(graph.edges.filter((e) => e.type === 'participates-in').map((e) => e.to).sort()).toEqual(['platform/caching', 'visualize/viewers']);
    expect(graph.problems.ambiguous_owners).toBe(1);
    expect(graph.details.ambiguous_owners).toEqual([`${SCATTER}: platform/caching and visualize/viewers claim it at rung 1`]);
  });
});
