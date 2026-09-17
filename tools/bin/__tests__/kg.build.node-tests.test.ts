/// `grok kg build` over the CLI (change-tests/plan.md, work order D1): public/tools as a source root the ts pass
/// walks, and the ts-node-tests extractor over the vitest and jest suites of the fixture.
import {describe, it, expect} from 'vitest';
import {parseNodeTests} from '../utils/kg/build/extract/ts/node-tests';
import {copyFixture, buildFixture, write} from './kg-fixture';

const TOOLS = 'public/tools/bin';
const IDS = `${TOOLS}/utils/ids.ts`;
const KG = `${TOOLS}/commands/kg.ts`;
const IDS_TEST = `${TOOLS}/__tests__/ids.test.ts`;
const THING_TEST = 'public/libraries/utils/src/thing.test.ts';
const FOLDER_TEST = 'public/libraries/utils/tests/thing.test.js';
const UNIT_TEST = 'public/js-api/scripts/unit/qnum.test.cjs';

/** A node:test suite in a library's tests folder importing the compiled path of a source, one registering through a wrapper,
 * a JS API unit test in CommonJS, and a fixture that is no suite. */
const graph = buildFixture(copyFixture('build', (repo) => {
  write(repo, FOLDER_TEST, ["import {test, describe} from 'node:test';", "import {category} from '../src/test.js';", '',
    'function smoke(name, body) {', '  test(name, async () => { await body(); });', '}', '',
    "describe('folder', () => {", "  smoke('smokes', () => category);", "  test('runs', () => {});", '});', ''].join('\n'));
  write(repo, 'public/libraries/utils/tests/fixtures/sample.test.js', "test('never', () => {});\n");
  write(repo, UNIT_TEST, ["const test = require('node:test');", '', "test('encodes', () => {});", ''].join('\n'));
}), 'ts-packages,ts-declarations,ts-imports,ts-node-tests');
const byId = (rows: any[], id: string) => rows.find((r) => r.id === id);
const edges = (rows: any[], from?: string, to?: string) => rows.filter((e) => (from === undefined || e.from === from) && (to === undefined || e.to === to));

describe('parseNodeTests', () => {
  it('reads test and it under their describe chain, skip and todo, a commented-out test, and the dynamic titles', () => {
    const tests = parseNodeTests([
      "describe('outer', () => {", "  it('one', () => {});", "  describe.skip('parked', () => { test('two', () => {}); });",
      "  test.todo('three');", "  // it('gone', () => {});", "  /* test('gone too', () => {}); */",
      "  test.each([[1, 2]])('adds %i and %i', () => {});", "  it.each`a | b`('$a plus $b', () => {});",
      "  test(`with ${n}`, () => {});", "  test('head ' + name, () => {});", '});', "it('alone', () => {});",
    ].join('\n'));
    expect(tests.map((t) => [t.describes.join(' > '), t.title, t.skipped, t.dynamic])).toEqual([
      ['outer', 'one', false, false], ['outer > parked', 'two', true, false], ['outer', 'three', true, false],
      ['outer', 'adds…', false, true], ['outer', '…', false, true], ['outer', 'with…', false, true], ['outer', 'head…', false, true],
      ['', 'alone', false, false],
    ]);
  });
});

describe('public/tools as a source root', () => {
  it('walks bin/**/*.ts of the CLI into source files declared by lib:tools, leaving out Babel output and fixtures', async () => {
    const {rows} = await graph;
    expect(byId(rows('nodes/library'), 'lib:tools')).toMatchObject({name: 'tools', npm: 'datagrok-tools', path: 'public/tools'});
    const files = rows('nodes/source-file').filter((f) => f.path.startsWith('public/tools/')).map((f) => f.path);
    expect(files).toEqual([IDS_TEST, KG, `${TOOLS}/grok.ts`, IDS]);
    expect(byId(rows('nodes/source-file'), `file:${IDS}`)).toMatchObject({language: 'ts', loc: 14, source_layer: 'public', visibility: 'public'});
    expect(byId(rows('nodes/source-file'), `file:${IDS}`).package).toBeUndefined();
    expect(edges(rows('edges/declares'), 'lib:tools').map((e) => e.to).sort()).toEqual(files.map((f) => `file:${f}`).sort());
    expect(edges(rows('edges/declares'), `file:${IDS}`).map((e) => e.to).sort()).toEqual([`decl:${IDS}#Ids`, `decl:${IDS}#ParsedId`]);
  });

  it('resolves the imports of the CLI relative to the file and to the JS API', async () => {
    const {rows} = await graph;
    const imports = rows('edges/imports').filter((e) => e.from.startsWith('file:public/tools/'));
    expect(imports.map((e) => [e.from, e.to, e.symbols])).toEqual([
      [`file:${IDS_TEST}`, `file:${IDS}`, ['Ids', 'fileId']],
      [`file:${KG}`, `file:${IDS}`, ['Ids', 'fileId']],
      [`file:${KG}`, 'lib:js-api', ['*']],
      [`file:${TOOLS}/grok.ts`, `file:${KG}`, ['kg']],
    ]);
  });
});

describe('ts-node-tests extractor', () => {
  it('emits one suite per vitest or jest file and a test per case, the describe chain as the category', async () => {
    const {rows, manifest} = await graph;
    expect(manifest.sources['ts-node-tests']).toBe('ok');
    expect(rows('nodes/test-suite').map((s) => [s.id, s.name, s.framework, s.path]).sort()).toEqual([
      [`suite:node:${UNIT_TEST}`, 'qnum.test.cjs', 'node', UNIT_TEST],
      [`suite:node:${THING_TEST}`, 'thing.test.ts', 'node', THING_TEST],
      [`suite:node:${FOLDER_TEST}`, 'thing.test.js', 'node', FOLDER_TEST],
      [`suite:node:${IDS_TEST}`, 'ids.test.ts', 'node', IDS_TEST],
    ].sort());
    const tests = rows('nodes/test');
    expect(tests.map((t) => [t.id, t.category, t.skipped, t.dynamic]).sort()).toEqual([
      [`test:node:${UNIT_TEST}#encodes`, undefined, false, false],
      [`test:node:${THING_TEST}#Thing/goes`, 'Thing', false, false],
      [`test:node:${FOLDER_TEST}#folder/runs`, 'folder', false, false],
      [`test:node:${FOLDER_TEST}#folder/smokes`, 'folder', false, false],
      [`test:node:${IDS_TEST}#ids > parse/reads a prefix`, 'ids > parse', true, false],
      [`test:node:${IDS_TEST}#ids > parse/reads the scheme`, 'ids > parse', false, false],
      [`test:node:${IDS_TEST}#ids/fileId(…`, 'ids', false, true],
      [`test:node:${IDS_TEST}#ids/prefixes a file`, 'ids', false, false],
      [`test:node:${IDS_TEST}#stands alone…`, undefined, false, true],
    ].sort());
    expect(byId(tests, `test:node:${IDS_TEST}#ids/prefixes a file`)).toMatchObject({name: 'prefixes a file', path: IDS_TEST, framework: 'node', level: 'unit', suite: `suite:node:${IDS_TEST}`, provenance: 'ast', source_layer: 'public'});
    expect(rows('nodes/source-file').some((f) => f.path.includes('/fixtures/'))).toBe(false);
  });

  it('observes the test folders of a library and the JS API as source files, a compiled-path import resolving to the source', async () => {
    const {rows} = await graph;
    expect(byId(rows('nodes/source-file'), `file:${FOLDER_TEST}`)).toMatchObject({language: 'js', loc: 11, source_layer: 'public'});
    expect(byId(rows('nodes/source-file'), `file:${UNIT_TEST}`)).toMatchObject({language: 'js'});
    expect(edges(rows('edges/declares'), 'lib:utils', `file:${FOLDER_TEST}`)).toHaveLength(1);
    expect(edges(rows('edges/imports'), `file:${FOLDER_TEST}`).map((e) => [e.to, e.symbols])).toEqual([['file:public/libraries/utils/src/test.ts', ['category']]]);
    expect(edges(rows('edges/declares'), `file:${FOLDER_TEST}`).map((e) => e.to).sort()).toEqual([`test:node:${FOLDER_TEST}#folder/runs`, `test:node:${FOLDER_TEST}#folder/smokes`]);
  });

  it('declares every case from its file and counts a run-time title as a dynamic test', async () => {
    const {rows, problems} = await graph;
    const declares = rows('edges/declares').filter((e) => e.to.startsWith('test:node:'));
    expect(declares.every((e) => e.from === `file:${e.to.split(':')[2].split('#')[0]}` && e.derived_by === 'ast' && e.confidence === 1)).toBe(true);
    expect(edges(declares, `file:${IDS_TEST}`)).toHaveLength(5);
    expect(edges(declares, `file:${THING_TEST}`)).toHaveLength(1);
    expect(problems.dynamic_tests).toEqual([
      `${IDS_TEST}: ids/fileId(… names a registration site, not a runnable test: the title is built at run time`,
      `${IDS_TEST}: stands alone… names a registration site, not a runnable test: the title is built at run time`,
    ]);
  });
});
