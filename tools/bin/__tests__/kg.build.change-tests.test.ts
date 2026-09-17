/// What `grok test --recent` walks (change-tests/plan.md, work order A): Dart `imports`, the client's `regTest`
/// cases, `declares` from every test file to its tests, `mentions` from a page to the files it cites, the
/// `entry` flag on package entry points, and headings inside a fence left alone.
import {describe, it, expect} from 'vitest';
import {copyFixture, buildFixture, write} from './kg-fixture';

const D4 = 'core/client/d4/lib/src';
const LEGEND = `${D4}/legends/legend.dart`;
const LEGEND_TEST = `${D4}/legends/test/legend_test.dart`;
const VIEWER = `${D4}/viewers/viewer.dart`;
const HISTOGRAM = `${D4}/viewers/histogram/histogram.dart`;
const REG_TESTS = 'core/client/xamgle/lib/src/tests/browse_test.dart';
const DG_TESTS = 'public/packages/Tested/src/tests/demo-tests.ts';
const PLAYWRIGHT = 'public/packages/Tested/playwright/basic.test.ts';
const FENCED = 'core/docs/FENCED.md';
const DDT = 'core/shared/ddt';
const DDT_LIB = `${DDT}/lib/ddt.dart`;
const DATA_FRAME = `${DDT}/lib/src/data_frame.dart`;
const COLUMN = `${DDT}/lib/src/column.dart`;
const DATA_FRAME_TEST = `${DDT}/test/data_frame_test.dart`;
const SORT_TEST = `${DDT}/test/data_frame_sort_test.dart`;
const ROWS_TEST = `${DDT}/test/rows_test.dart`;
const SETUP = `${DDT}/test/setup.dart`;
const HISTOGRAM_TESTS = 'public/packages/Tested/src/tests/histogram-tests.ts';
const TEST_TRACK = 'public/packages/UsageAnalysis/files/TestTrack/Viewers';

/** A Dart library with two parts, a test named after one part, one named after it with a suffix, one whose name two files share. */
function writeDdt(repo: string): void {
  write(repo, DDT_LIB, "library ddt;\n\npart 'src/data_frame.dart';\npart 'src/column.dart';\n");
  write(repo, DATA_FRAME, 'part of ddt;\n\nclass DataFrame {}\nclass Stats {}\n');
  write(repo, COLUMN, 'part of ddt;\n\nclass Column {}\nclass Stats {}\nclass Row {}\nclass Named {}\n');
  write(repo, `${DDT}/lib/src/util/rows.dart`, 'part of ddt;\n');
  write(repo, `${DDT}/lib/src/grid/rows.dart`, 'part of ddt;\n');
  write(repo, DATA_FRAME_TEST, ["import 'package:ddt/ddt.dart';", "import 'package:test/test.dart';", '', 'void main() {',
    "  test('a DataFrame holds a Column', () { DataFrame df; Column c; Stats s; Row r; Nowhere n; parse('Named([1])'); });",
    '  const raw = r"Named(x)"; const multi = \'\'\'\nNamed\n\'\'\';', '}', ''].join('\n'));
  write(repo, SORT_TEST, ["import 'package:ddt/ddt.dart';", '', 'void main() {', "  test('sorts', () {});", '}', ''].join('\n'));
  write(repo, ROWS_TEST, ['void main() {', "  test('rows', () {});", '}', ''].join('\n'));
  write(repo, SETUP, ["import 'package:ddt/ddt.dart';", '', 'DataFrame frame;', ''].join('\n'));
}

/** A DG test file named after a feature, a Test Track spec whose folder two features spell, and the second of them. */
function writeNamed(repo: string): void {
  write(repo, HISTOGRAM_TESTS, ["import {category, test} from '@datagrok-libraries/utils/src/test';", '', "category('Tested: Histogram', () => {",
    "  test('renders', async () => {});", '});', ''].join('\n'));
  write(repo, `${TEST_TRACK}/Legends/colors-spec.ts`, ["import {test} from '@playwright/test';", '', "test('colors', async () => {});", ''].join('\n'));
  write(repo, 'core/docs/LEGENDS2.md', ['---', 'feature: platform/legends', 'name: Legends', '---', '', '# Legends', ''].join('\n'));
}

const graph = buildFixture(copyFixture('build', (repo) => {
  write(repo, REG_TESTS, [
    'part of xamgle.tests;', '', 'void _initBrowseTests() {',
    "  regTest('Browse | Construction', () async {});",
    "  regTest('Browse | Platform | Users', () async {});",
    "  regTest('Smoke', () async {});",
    "  regTest('${DartLibraryTestCategoryName.ddt} | Ddt benchmark | $test', test.run);",
    "  // regTest('Browse | Gone', () async {});",
    '}', ''].join('\n'));
  write(repo, FENCED, ['Intro line.', '', '```dart', '# not a heading', '```', '', '# Fenced page', '', 'Body.', ''].join('\n'));
  writeDdt(repo);
  writeNamed(repo);
}), 'homes,ts-packages,ts-declarations,ts-tests,docs,dart,membership');
const edges = (rows: any[], from?: string, to?: string) => rows.filter((e) => (from === undefined || e.from === from) && (to === undefined || e.to === to));

describe('the graph rows a change walk needs (change-tests/plan.md, work order A)', () => {
  it('resolves Dart import, export and part directives to the files the pass walked, with the names a show clause lists', async () => {
    const {rows} = await graph;
    const imports = rows('edges/imports').filter((e) => e.from.startsWith('file:core/client/'));
    expect(imports.map((e) => [e.from, e.to, e.symbols]).sort()).toEqual([
      [`file:${HISTOGRAM}`, `file:${D4}/viewers/legend_cache.dart`, ['LegendCache']],
      [`file:${HISTOGRAM}`, `file:${VIEWER}`, undefined],
      [`file:${LEGEND_TEST}`, `file:${LEGEND}`, undefined],
      [`file:${VIEWER}`, `file:${D4}/viewers/viewer.g.dart`, ['part']],
      [`file:${VIEWER}`, `file:${LEGEND}`, undefined],
    ].sort());
    expect(imports[0]).toMatchObject({derived_by: 'ast', confidence: 1, evidence: [imports[0].from.slice('file:'.length)]});
    // dart:html, a pub package and a package path nobody walked are neither edges nor problems
    expect(rows('nodes/source-file').some((f) => f.path.includes('nowhere'))).toBe(false);
    expect((await graph).problems.unresolved_ids ?? []).not.toContainEqual(expect.stringContaining('nowhere.dart'));
  });

  it('turns regTest titles under the client tests folder into xamgle tests, the last segment the name, the rest the category', async () => {
    const {rows} = await graph;
    const tests = rows('nodes/test').filter((t) => t.framework === 'xamgle');
    expect(tests.map((t) => [t.id, t.category, t.name, t.dynamic]).sort()).toEqual([
      [`test:xamgle:${REG_TESTS}#Browse/Construction`, 'Browse', 'Construction', false],
      [`test:xamgle:${REG_TESTS}#Browse | Platform/Users`, 'Browse | Platform', 'Users', false],
      [`test:xamgle:${REG_TESTS}#Smoke`, undefined, 'Smoke', false],
      [`test:xamgle:${REG_TESTS}#ddt | Ddt benchmark/$test`, 'ddt | Ddt benchmark', '$test', true],
    ].sort());
    expect(tests.find((t) => t.name === 'Construction')).toMatchObject({level: 'e2e', path: REG_TESTS, suite: `suite:xamgle:${REG_TESTS}`, provenance: 'ast', source_layer: 'core'});
    expect(rows('nodes/test-suite').filter((s) => s.path === REG_TESTS)).toEqual([expect.objectContaining({id: `suite:xamgle:${REG_TESTS}`, framework: 'xamgle', name: 'browse_test.dart'})]);
    expect((await graph).manifest.dart_packages).toMatchObject({xamgle: 1});
  });

  it('declares every test from the file that holds it, whatever the framework', async () => {
    const {rows} = await graph;
    const declares = rows('edges/declares').filter((e) => e.to.startsWith('test:'));
    expect(declares.every((e) => e.from === `file:${e.to.split(':')[2].split('#')[0]}` && e.derived_by === 'ast' && e.confidence === 1)).toBe(true);
    expect(edges(declares, `file:${LEGEND_TEST}`).map((e) => e.to).sort()).toEqual([
      `test:dart:${LEGEND_TEST}#a legend measures its labels once`, `test:dart:${LEGEND_TEST}#placement/a legend takes the slot it is given`]);
    expect(edges(declares, `file:${REG_TESTS}`)).toHaveLength(4);
    expect(edges(declares, `file:${DG_TESTS}`)).toHaveLength(rows('nodes/test').filter((t) => t.path === DG_TESTS).length);
    expect(edges(declares, `file:${PLAYWRIGHT}`)).toHaveLength(4);
    expect(rows('nodes/source-file').find((f) => f.path === PLAYWRIGHT)).toMatchObject({language: 'ts', package: 'pkg:Tested', loc: 25, provenance: 'filesystem'});
  });

  it('draws mentions from a page to the source files it cites by path, and nothing to a folder or a page', async () => {
    const {rows} = await graph;
    const cited = rows('edges/mentions').filter((e) => e.to.startsWith('file:'));
    expect(cited.map((e) => [e.from, e.to]).sort()).toEqual([
      [`doc:${D4}/legends/README.md`, `file:${LEGEND}`],
      ['doc:core/docs/VIEWERS.md', `file:${D4}/legends/legend_renderer.dart`],
      ['doc:public/help/domains/bio/bio.md', 'file:public/packages/Demo/src/utils.ts'],
    ]);
    expect(cited[0]).toMatchObject({derived_by: 'annotation', confidence: 0.8, evidence: [cited[0].from.slice('doc:'.length)]});
  });

  it('flags the entry points of a package and the Dart library files that declare parts on their file rows', async () => {
    const {rows} = await graph;
    expect(rows('nodes/source-file').filter((f) => f.entry).map((f) => f.path).sort()).toEqual([
      VIEWER, DDT_LIB, 'public/packages/Demo/src/package.ts', 'public/packages/Plain/src/package-test.ts', 'public/packages/Plain/src/package.ts'].sort());
    expect(rows('nodes/source-file').find((f) => f.path === 'public/packages/Demo/src/viewer.ts').entry).toBe(false);
    expect(rows('nodes/source-file').find((f) => f.path === DATA_FRAME).entry).toBe(false);
    // the walk from a part reaches its library and stops there: the library imports the part, as a part
    expect(edges(rows('edges/imports'), `file:${DDT_LIB}`).map((e) => [e.to, e.symbols])).toEqual([[`file:${COLUMN}`, ['part']], [`file:${DATA_FRAME}`, ['part']]]);
  });

  it('draws mirrors from a test file to the one source file of its unit it is named after, exact or before a suffix, never to two', async () => {
    const {rows, problems} = await graph;
    const mirrors = rows('edges/mirrors');
    expect(mirrors.map((e) => [e.from, e.to, e.confidence]).sort()).toEqual([
      [`file:${DATA_FRAME_TEST}`, `file:${DATA_FRAME}`, 0.9],
      [`file:${LEGEND_TEST}`, `file:${LEGEND}`, 0.9],
      [`file:${SORT_TEST}`, `file:${DATA_FRAME}`, 0.8],
    ].sort());
    expect(mirrors[0]).toMatchObject({derived_by: 'filesystem', evidence: [mirrors[0].from.slice('file:'.length)]});
    expect(problems.ambiguous_mirrors).toEqual([`${ROWS_TEST}: ${DDT}/lib/src/grid/rows.dart and ${DDT}/lib/src/util/rows.dart share the name`]);
  });

  it('draws lexical uses from a Dart test file, or a helper of its test folder, to the one type of its package or an imported one that a capitalized word names', async () => {
    const {rows, problems} = await graph;
    const uses = rows('edges/uses').filter((e) => e.derived_by === 'lexical');
    expect(uses.map((e) => [e.from, e.to]).sort()).toEqual([
      [`file:${DATA_FRAME_TEST}`, `decl:${COLUMN}#Column`],
      [`file:${DATA_FRAME_TEST}`, `decl:${DATA_FRAME}#DataFrame`],
      [`file:${SETUP}`, `decl:${DATA_FRAME}#DataFrame`],
    ]);
    expect(rows('nodes/test').some((t) => t.path === SETUP)).toBe(false);
    expect(uses[0]).toMatchObject({kind: 'type', confidence: 0.7, evidence: [DATA_FRAME_TEST]});
    // Stats is declared twice, Row is too short, Nowhere is declared nowhere, Named appears only inside string literals
    expect(problems.ambiguous_uses).toEqual([`${DATA_FRAME_TEST}: Stats is declared in ${COLUMN} and ${DATA_FRAME}`]);
  });

  it('draws tests by name: a file base name, a category segment or a Test Track folder that spells one feature, never one two features share', async () => {
    const {rows, problems} = await graph;
    const named = rows('edges/tests').filter((e) => e.derived_by === 'name');
    expect(named.map((e) => [e.from, e.to]).sort()).toEqual([
      [`test:dg:${HISTOGRAM_TESTS}#Tested: Histogram/renders`, 'visualize/viewers/histogram'],
      [`test:playwright:${TEST_TRACK}/ScatterPlot/scatterplot-legend-spec.ts#Scatter plot legend/color legend`, 'visualize/viewers/scatter-plot'],
      [`test:playwright:${TEST_TRACK}/ScatterPlot/scatterplot-legend-spec.ts#Scatter plot legend/marker legend`, 'visualize/viewers/scatter-plot'],
    ]);
    expect(named[0]).toMatchObject({confidence: 0.7, evidence: [HISTOGRAM_TESTS]});
    // the folder rung already links the Tested file to the viewer area, so the "Viewers" of its category draws no second edge
    expect(rows('edges/tests').filter((e) => e.from === named[0].from && e.to === 'visualize/viewers')).toEqual([expect.objectContaining({derived_by: 'filesystem'})]);
    expect(problems.ambiguous_test_names).toEqual(['legends: platform/legends and visualize/legends']);
    expect(rows('edges/tests').some((e) => e.from.includes('/Legends/colors-spec.ts'))).toBe(false);
  });

  it('names a page by its first heading outside any fence', async () => {
    const {rows} = await graph;
    expect(rows('nodes/doc-page').find((d) => d.path === FENCED)).toMatchObject({name: 'Fenced page', description: 'Intro line.'});
  });
});
