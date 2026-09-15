/// The lexical Dart extractor (build-plan.md WO-6): the source files, top-level declarations, tests and
/// `~id` markers a regex pass over `core/**/*.dart` finds, and what the markers make of the ownership
/// the other rungs would have given those files.
import {describe, it, expect, vi} from 'vitest';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {fileURLToPath} from 'url';
import {currentDir} from '../utils/kg/build/write';
import {kg} from '../commands/kg';

const fixture = path.join(path.dirname(fileURLToPath(import.meta.url)), 'fixtures', 'kg', 'build');
const KG_DIR = path.join('core', 'docs', 'knowledge-graph');
const D4 = 'core/client/d4/lib/src';
const LEGEND = `${D4}/legends/legend.dart`;
const RENDERER = `${D4}/legends/legend_renderer.dart`;
const VIEWER = `${D4}/viewers/viewer.dart`;
const GENERATED = `${D4}/viewers/viewer.g.dart`;
const TEST_FILE = `${D4}/legends/test/legend_test.dart`;
const HISTOGRAM = `${D4}/viewers/histogram/histogram.dart`;
const HELP_TABLE = 'core/shared/grok_shared/lib/src/help_url.dart';

interface Built {
  manifest: any;
  rows: (file: string) => any[];
}

function copy(prepare?: (repo: string) => void): string {
  const repo = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-kg-dart-'));
  fs.cpSync(fixture, repo, {recursive: true});
  prepare?.(repo);
  return repo;
}

async function build(repo: string, only = 'homes,dart,membership'): Promise<Built> {
  const log = vi.spyOn(console, 'log').mockImplementation(() => {});
  const error = vi.spyOn(console, 'error').mockImplementation(() => {});
  const before = process.exitCode;
  try {
    await kg({_: ['kg', 'build'], kg: path.join(repo, KG_DIR), only, db: false, output: 'json'});
    expect(error.mock.calls).toEqual([]);
    const out = currentDir(path.join(repo, '.kg'))!;
    const rows = (file: string) => {
      const p = path.join(out, file.startsWith('reports/') ? file : `data/${file}.jsonl`);
      return fs.existsSync(p) ? fs.readFileSync(p, 'utf8').split('\n').filter(Boolean).map((l) => JSON.parse(l)) : [];
    };
    return {manifest: JSON.parse(String(log.mock.calls[0][0])), rows};
  }
  finally {
    process.exitCode = before;
    log.mockRestore();
    error.mockRestore();
  }
}

const graph = build(copy());

describe('the lexical Dart pass (build-plan.md WO-6)', () => {
  it('emits one source file per Dart file, generated ones marked, and counts them per package', async () => {
    const {manifest, rows} = await graph;
    expect(manifest.sources.dart).toBe('ok');
    expect(manifest.dart_depth).toBe('lexical');
    const files = rows('nodes/source-file').filter((f) => f.language === 'dart');
    expect(files.map((f) => f.path).sort()).toEqual([LEGEND, TEST_FILE, RENDERER, `${D4}/viewers/legend_cache.dart`,
      `${D4}/viewers/scatterplot/scatter.dart`, VIEWER, GENERATED, HISTOGRAM, HELP_TABLE].sort());
    expect(files.find((f) => f.path === GENERATED)).toMatchObject({generated: true, provenance: 'filesystem', source_layer: 'core'});
    expect(files.filter((f) => f.generated).map((f) => f.path)).toEqual([GENERATED]);
    expect(manifest.dart_packages).toEqual({d4: 8, grok_shared: 1});
  }, 60_000);

  it('takes the top-level types of a file, with the doc comment and the annotation above them', async () => {
    const {rows} = await graph;
    const decls = rows('nodes/declaration').filter((d) => d.language === 'dart');
    expect(decls.map((d) => d.id).sort()).toEqual([
      `decl:${LEGEND}#Legend`, `decl:${RENDERER}#LegendRenderer`, `decl:${D4}/viewers/legend_cache.dart#LegendCache`,
      `decl:${D4}/viewers/scatterplot/scatter.dart#ScatterPlot`, `decl:${VIEWER}#Viewer`, `decl:${GENERATED}#ViewerProps`,
      `decl:${HISTOGRAM}#Histogram`, `decl:${HELP_TABLE}#HelpUrl`,
    ].sort());
    expect(decls.find((d) => d.id === `decl:${LEGEND}#Legend`)).toMatchObject({name: 'Legend', kind: 'class', exported: true,
      documented: true, deprecated: true, line: 3, path: LEGEND, provenance: 'ast', source_layer: 'core'});
    expect(decls.find((d) => d.id === `decl:${GENERATED}#ViewerProps`)).toMatchObject({generated: true, documented: false});
    expect(rows('edges/declares')).toContainEqual(expect.objectContaining({from: `file:${LEGEND}`, to: `decl:${LEGEND}#Legend`, derived_by: 'ast'}));
  }, 60_000);

  it('takes the tests of a test file under the groups enclosing them, in one suite', async () => {
    const {rows} = await graph;
    expect(rows('nodes/test').filter((t) => t.framework === 'dart')).toMatchObject([
      {id: `test:dart:${TEST_FILE}#a legend measures its labels once`, level: 'unit', suite: `suite:dart:${TEST_FILE}`},
      {id: `test:dart:${TEST_FILE}#placement/a legend takes the slot it is given`, level: 'unit', category: 'placement', suite: `suite:dart:${TEST_FILE}`},
    ]);
    expect(rows('nodes/test').find((t) => t.name === 'a legend measures its labels once').category).toBeUndefined();
    expect(rows('nodes/test-suite').filter((s) => s.framework === 'dart')).toMatchObject([
      {id: `suite:dart:${TEST_FILE}`, name: 'legend_test.dart', path: TEST_FILE},
    ]);
  }, 60_000);

  it('lets a `/// ~id` marker own the file the folder would have, and a `// ~id` line only participate', async () => {
    const {rows} = await graph;
    expect(rows('edges/is-implemented-in').filter((e) => e.to === `file:${RENDERER}`)).toEqual([expect.objectContaining({
      from: 'visualize/viewers/scatter-plot', derived_by: 'annotation', confidence: 1, evidence: [RENDERER],
    })]);
    expect(rows('edges/participates-in').filter((e) => e.from === `file:${VIEWER}`).map((e) => e.to)).toEqual(['platform/caching']);
    expect(rows('edges/is-implemented-in').filter((e) => e.to === `file:${VIEWER}`).map((e) => e.from)).toEqual(['visualize/viewers']);
    expect(rows('reports/claims.jsonl').filter((c) => c.source === 'marker')).toMatchObject([
      {file: RENDERER, feature: 'visualize/viewers/scatter-plot', rung: 1, line: 1},
      {file: VIEWER, feature: 'platform/caching', rung: 1, mode: 'participates', line: 1},
    ]);
  }, 60_000);

  it('follows a Dart test to the feature that owns its file', async () => {
    const {rows} = await graph;
    const tests = rows('edges/tests').filter((e) => e.from.startsWith('test:dart:'));
    expect(tests).toHaveLength(2);
    expect(tests.every((e) => e.to === 'visualize/legends' && e.derived_by === 'filesystem' && e.confidence === 0.9)).toBe(true);
  }, 60_000);

  it('draws documents from the help page a file names to the feature that owns the file, and reports a page that is gone', async () => {
    const {rows, manifest} = await graph;
    // the literal and the doc-comment path of viewer.dart; the page the histogram home cites as well is one edge, annotated
    expect(rows('edges/documents').filter((e) => e.derived_by === 'ast')).toEqual([
      expect.objectContaining({from: 'doc:public/help/datagrok/project.md', to: 'visualize/viewers', confidence: 0.8, evidence: [VIEWER]}),
      expect.objectContaining({from: 'doc:public/help/domains/bio/sequences.md', to: 'visualize/viewers', confidence: 0.8, evidence: [VIEWER]}),
    ]);
    expect(rows('edges/documents')).toContainEqual(expect.objectContaining({from: 'doc:public/help/visualize/viewers/histogram.md',
      to: 'visualize/viewers/histogram', derived_by: 'annotation', confidence: 1,
      evidence: ['core/client/d4/lib/src/viewers/histogram/CLAUDE.md', HISTOGRAM]}));
    // the page HelpUrl.Gone names exists nowhere: stale where it is defined and where it is used
    expect(manifest.problems.unresolved_ids).toBe(2);
  }, 60_000);

  it('counts a marker that names no home and reports the pass as partial', async () => {
    const {manifest, rows} = await build(copy((repo) => fs.writeFileSync(path.join(repo, ...`${D4}/viewers/axes.dart`.split('/')),
      '// ~visualize/viewers/boi\nclass Axes {\n}\n')));
    expect(manifest.sources.dart).toBe('partial');
    // the marker, plus the two files naming the page HelpUrl.Gone points at
    expect(manifest.problems.unresolved_ids).toBe(3);
    expect(rows('edges/participates-in').some((e) => e.to === 'visualize/viewers/boi')).toBe(false);
  }, 60_000);
});
