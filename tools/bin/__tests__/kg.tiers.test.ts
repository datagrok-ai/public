/// The tiers of `tests-for` and `impact`, the change set and the runner rows (change-tests/plan.md WO-B).
import {describe, it, expect, beforeAll, afterAll} from 'vitest';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {fileURLToPath} from 'url';
import {load, loadKuzu, open, run} from '../utils/kg/kuzu';
import {testsFor, impact, resolveTarget, resolveTargets, runRows, parseTiers} from '../utils/kg/ops';
import {testUnitsOf} from '../utils/kg/ids';
import {changeSet} from '../utils/kg/changes';
import {replay, Commit} from '../utils/kg/report';
import {loadQuestions, ask} from '../utils/kg/questions';
import {copyFixture, buildFixture, fixtureTypes, git, write} from './kg-fixture';

const LIMIT = {limit: 50};
const ONLY = 'homes,ts-packages,ts-declarations,ts-imports,ts-uses,ts-tests,ts-node-tests,membership';
const UTILS = 'file:public/packages/Demo/src/utils.ts';
const TESTS = '@datagrok-libraries/utils/src/test';
const TEST_LIB = 'file:public/libraries/utils/src/test.ts';
const questionsRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..', '..', '..', '..', 'core', 'docs', 'knowledge-graph');
const withKuzu = loadKuzu() ? it : it.skip;

/** The fixture with a test file importing a source file, a mirror test file, an import chain of three hops, one
 * that only passes through an entry file, a test using a JS API declaration, and a playwright config. */
function prepare(repo: string): void {
  write(repo, 'public/packages/Demo/src/tests/utils-tests.ts', `import {category, test} from '${TESTS}';\nimport * as DG from 'datagrok-api/dg';\n` +
    `import {callThings} from '../utils';\n\ncategory('Demo: Utils', () => {\n  test('calls things', async () => { await callThings(); });\n` +
    `  test('calls again', async () => { DG.DataFrame.fromCsv('a'); });\n});\n`);
  write(repo, 'public/packages/Demo/src/tests/base-tests.ts', `import {category, test} from '${TESTS}';\n\ncategory('Demo: Base', () => {\n` +
    `  test('mirrors base', async () => {});\n});\n`);
  write(repo, 'public/packages/Demo/src/deep.ts', `import {DemoViewer} from './viewer';\nexport const deep = DemoViewer;\n`);
  write(repo, 'public/packages/Demo/src/tests/deep-tests.ts', `import {category, test} from '${TESTS}';\nimport {deep} from '../deep';\n\n` +
    `category('Demo: Deep', () => {\n  test('reaches utils', async () => { console.log(deep); });\n});\n`);
  write(repo, 'public/packages/Demo/src/package-api.ts', `import {callThings} from './utils';\nexport const api = callThings;\n`);
  write(repo, 'public/packages/Demo/src/tests/api-tests.ts', `import {category, test} from '${TESTS}';\nimport {api} from '../package-api';\n\n` +
    `category('Demo: Api', () => {\n  test('through the entry', async () => { console.log(api); });\n});\n`);
  write(repo, 'public/packages/Tested/playwright/playwright.config.ts', 'export default {};\n');
  // a helper of the tests folder using a JS API declaration, the test importing it, and a barrel over a plain file
  write(repo, 'public/packages/Tested/src/tests/setup.ts', `import * as DG from 'datagrok-api/dg';\nexport const df = DG.DataFrame.fromCsv('a');\n`);
  write(repo, 'public/packages/Tested/src/lib.ts', 'export const lib = 1;\n');
  write(repo, 'public/packages/Tested/src/all.ts', `export * from './lib';\n`);
  write(repo, 'public/packages/Tested/src/everything.ts', `export * from './all';\n`);
  write(repo, 'public/packages/Tested/src/tests/setup-tests.ts', `import {category, test} from '${TESTS}';\nimport {df} from './setup';\nimport {lib} from '../lib';\n\n` +
    `category('Tested: Setup', () => {\n  test('boots', async () => { console.log(df, lib); });\n});\n`);
}

describe('the tiers of tests-for and impact (change-tests/plan.md)', () => {
  let index: {repo: string, opened: Awaited<ReturnType<typeof open>>};

  beforeAll(async () => {
    if (!loadKuzu()) return;
    const built = await buildFixture(copyFixture('build', prepare), ONLY);
    await load(built.out, fixtureTypes('build'));
    index = {repo: built.repo, opened: await open(built.out, true)};
  });

  afterAll(async () => {
    if (!index?.opened) return;
    await index.opened.conn.close();
    await index.opened.db.close();
  });

  const conn = () => index.opened!.conn;
  const titles = (sections: {title: string}[]) => sections.map((s) => s.title);
  const rows = (answer: {sections: {title: string, rows: Record<string, unknown>[]}[]}, title: string) => answer.sections.find((s) => s.title === title)!.rows;

  withKuzu('puts the tests importing a file in the immediate tier and those reaching it through imports in the reachable one, never through an entry file', async () => {
    const answer = await testsFor(conn(), (await resolveTarget(conn(), UTILS))!, {...LIMIT, repoRoot: index.repo});
    expect(titles(answer.sections)).toEqual(['features', 'immediate', 'reachable', 'feature', 'scenarios', 'automations', 'run']);
    expect(rows(answer, 'immediate').map((r) => [r.name, r.via, r.tier])).toEqual([
      ['calls again', `${UTILS} → imports ← file:public/packages/Demo/src/tests/utils-tests.ts`, 'immediate'],
      ['calls things', `${UTILS} → imports ← file:public/packages/Demo/src/tests/utils-tests.ts`, 'immediate'],
    ]);
    expect(rows(answer, 'immediate')[0]).toMatchObject({framework: 'dg', level: 'unit', path: 'public/packages/Demo/src/tests/utils-tests.ts', category: 'Demo: Utils', feature: null});
    expect(rows(answer, 'reachable').map((r) => [r.name, r.via])).toEqual([
      ['reaches utils', `${UTILS} → imports ← file:public/packages/Demo/src/viewer.ts → imports ← file:public/packages/Demo/src/deep.ts → imports ← file:public/packages/Demo/src/tests/deep-tests.ts`],
    ]);
    // api-tests.ts reaches utils.ts only through package-api.ts, an entry file the walk does not expand
    expect(JSON.stringify(answer.sections)).not.toContain('through the entry');
    // the feature tier is the rest: what the owning feature carries beyond the two tiers above, each once
    const feature = rows(answer, 'feature');
    expect(feature.length).toBeGreaterThan(0);
    expect(feature.every((r) => r.tier === 'feature' && String(r.via).includes(' → tests ← '))).toBe(true);
    expect(feature.map((r) => r.name)).not.toContain('calls things');
    // two of Demo's four DG test files are selected across the tiers: half, so the DG row is the whole suite rather than one call per category
    expect(rows(answer, 'run')[0]).toEqual({framework: 'dg', cwd: 'public/packages/Demo', command: 'grok test', tests: 3,
      names: ['calls again', 'calls things', 'reaches utils'], note: 'whole suite: 2 of 4 test files selected'});
    const only = await testsFor(conn(), (await resolveTarget(conn(), UTILS))!, {...LIMIT, tiers: ['immediate']});
    expect(titles(only.sections)).toEqual(['features', 'immediate', 'run']);
    expect(rows(only, 'run')).toEqual([{framework: 'dg', cwd: 'public/packages/Demo', command: 'grok test --category "Demo: Utils"', tests: 2, names: ['calls again', 'calls things']}]);
  });

  withKuzu('finds the mirror test file through the mirrors edge, and the tests using what a file declares', async () => {
    const base = await testsFor(conn(), (await resolveTarget(conn(), 'public/packages/Demo/src/base.ts'))!, LIMIT);
    expect(rows(base, 'immediate').map((r) => [r.name, r.via])).toEqual([
      ['mirrors base', 'file:public/packages/Demo/src/base.ts → mirrors ← file:public/packages/Demo/src/tests/base-tests.ts'],
    ]);
    // the JS API is another unit than the package using it: the link holds, the tier is reachable
    const dataframe = await testsFor(conn(), (await resolveTarget(conn(), 'public/js-api/src/dataframe.ts'))!, LIMIT);
    expect(rows(dataframe, 'immediate')).toEqual([]);
    const users = rows(dataframe, 'reachable').filter((r) => r.path === 'public/packages/Demo/src/tests/utils-tests.ts');
    expect(users.map((r) => [r.name, r.tier])).toEqual([['calls again', 'reachable'], ['calls things', 'reachable']]);
    expect(users[0].via).toMatch(/^file:public\/js-api\/src\/dataframe\.ts → declares → decl:[^ ]+ → uses ← file:public\/packages\/Demo\/src\/tests\/utils-tests\.ts$/);
    // a test reaches the declaration through the helper of its tests folder that uses it: one more imports hop
    const helped = rows(dataframe, 'reachable').filter((r) => r.path === 'public/packages/Tested/src/tests/setup-tests.ts');
    expect(helped.map((r) => r.name)).toEqual(['boots']);
    expect(helped[0].via).toMatch(/^file:public\/js-api\/src\/dataframe\.ts → declares → decl:[^ ]+ → uses ← file:public\/packages\/Tested\/src\/tests\/setup\.ts → imports ← file:public\/packages\/Tested\/src\/tests\/setup-tests\.ts$/);
  });

  withKuzu('lets a barrel stand for the files it re-exports, in tests-for and in impact, and says so', async () => {
    const barrel = 'file:public/packages/Tested/src/all.ts';
    const reexport = await run(conn(), `MATCH (b)-[e:IMPORTS]->(m) WHERE b.id = '${barrel}' RETURN m.id AS target, e.reexport AS reexport`);
    expect(reexport.rows).toEqual([{target: 'file:public/packages/Tested/src/lib.ts', reexport: true}]);
    const answer = await testsFor(conn(), (await resolveTarget(conn(), barrel))!, LIMIT);
    expect(rows(answer, 'immediate').map((r) => [r.name, r.via])).toEqual([
      ['boots', `${barrel} → imports → file:public/packages/Tested/src/lib.ts → imports ← file:public/packages/Tested/src/tests/setup-tests.ts`]]);
    expect(answer.notes).toEqual(['barrel: public/packages/Tested/src/all.ts expanded to 1 re-exported file']);
    const reach = await impact(conn(), (await resolveTarget(conn(), barrel))!, LIMIT);
    expect(rows(reach, 'importers').map((r) => r.caller)).toEqual(['file:public/packages/Tested/src/everything.ts', 'file:public/packages/Tested/src/tests/setup-tests.ts']);
    expect(reach.notes).toEqual(answer.notes);
    // a barrel of barrels is followed down to the files, with one note for the target
    const nested = await testsFor(conn(), (await resolveTarget(conn(), 'public/packages/Tested/src/everything.ts'))!, LIMIT);
    expect(rows(nested, 'immediate').map((r) => r.via)).toEqual([
      `file:public/packages/Tested/src/everything.ts → imports → ${barrel} → imports → file:public/packages/Tested/src/lib.ts → imports ← file:public/packages/Tested/src/tests/setup-tests.ts`]);
    expect(nested.notes).toEqual(['barrel: public/packages/Tested/src/everything.ts expanded to 2 re-exported files']);
    // a file that declares something is no barrel, whatever it re-exports
    const plain = await impact(conn(), (await resolveTarget(conn(), 'public/packages/Tested/src/lib.ts'))!, LIMIT);
    expect(plain.notes).toBeUndefined();
  });

  withKuzu('puts a test file of another unit in the reachable tier even when it imports the changed file directly', async () => {
    const answer = await testsFor(conn(), (await resolveTarget(conn(), TEST_LIB))!, LIMIT);
    expect(rows(answer, 'immediate').map((r) => [r.name, r.via, r.framework])).toEqual([
      ['goes', `${TEST_LIB} → imports ← file:public/libraries/utils/src/thing.test.ts`, 'node']]);
    const dg = rows(answer, 'reachable').filter((r) => r.path === 'public/packages/Demo/src/tests/utils-tests.ts');
    expect(dg.map((r) => [r.name, r.via, r.tier])).toEqual([
      ['calls again', `${TEST_LIB} → imports ← file:public/packages/Demo/src/tests/utils-tests.ts`, 'reachable'],
      ['calls things', `${TEST_LIB} → imports ← file:public/packages/Demo/src/tests/utils-tests.ts`, 'reachable']]);
    expect(rows(answer, 'reachable').every((r) => String(r.path).startsWith('public/packages/'))).toBe(true);
  });

  withKuzu('answers the tests-for-change question with the four links as one table', async () => {
    const question = loadQuestions(questionsRoot).questions.find((q) => q.id === 'tests-for-change')!;
    expect(question.params.file.default).toBe('core/shared/ddt/lib/src/data_frame/data_frame.dart');
    const base = await ask(conn(), question, {file: 'public/packages/Demo/src/base.ts'});
    expect(base.rows).toEqual([{id: 'test:dg:public/packages/Demo/src/tests/base-tests.ts#Demo: Base/mirrors base', framework: 'dg', name: 'mirrors base',
      path: 'public/packages/Demo/src/tests/base-tests.ts', file: 'file:public/packages/Demo/src/tests/base-tests.ts', via: 'mirrors'}]);
    // utils-tests.ts both imports and mirrors utils.ts: one row per link
    const utils = await ask(conn(), question, {file: 'public/packages/Demo/src/utils.ts'});
    expect(utils.rows.map((r) => [r.via, r.name]).sort()).toEqual([['imports', 'calls again'], ['imports', 'calls things'], ['mirrors', 'calls again'], ['mirrors', 'calls things']]);
    const declared = await ask(conn(), question, {file: 'public/packages/Demo/src/tests/utils-tests.ts'});
    expect(declared.rows.map((r) => [r.via, r.name, r.file]).sort()).toEqual([
      ['declares', 'calls again', 'file:public/packages/Demo/src/tests/utils-tests.ts'], ['declares', 'calls things', 'file:public/packages/Demo/src/tests/utils-tests.ts']]);
    const api = await ask(conn(), question, {file: 'public/js-api/src/dataframe.ts'});
    expect(api.rows.filter((r) => r.via === 'uses').map((r) => r.name)).toContain('calls again');
  });

  withKuzu('answers several targets at once, leads with the changes, and restricts the tiers and the run rows to --tier', async () => {
    const resolved = await resolveTargets(conn(), ['public/packages/Demo/src/utils.ts', 'public/packages/Demo/src/base.ts', 'public/nowhere.ts']);
    expect(resolved.get('public/nowhere.ts')).toBeNull();
    const targets = [resolved.get('public/packages/Demo/src/utils.ts')!, resolved.get('public/packages/Demo/src/base.ts')!];
    const changes = [{path: 'public/packages/Demo/src/utils.ts', known: true, repo: 'public'}, {path: 'public/nowhere.ts', known: false, repo: 'public'}];
    const answer = await testsFor(conn(), targets, {...LIMIT, tiers: ['immediate']}, changes);
    expect(answer.target).toBeUndefined();
    expect(answer.targets).toHaveLength(2);
    expect(titles(answer.sections)).toEqual(['changes', 'features', 'immediate', 'run']);
    expect(rows(answer, 'changes')).toEqual(changes);
    expect(rows(answer, 'immediate').map((r) => r.name)).toEqual(['calls again', 'calls things', 'mirrors base']);
    // two of Demo's four DG test files: half, so the whole suite in one row
    expect(rows(answer, 'run')).toEqual([{framework: 'dg', cwd: 'public/packages/Demo', command: 'grok test', tests: 3, names: ['calls again', 'calls things', 'mirrors base'],
      note: 'whole suite: 2 of 4 test files selected'}]);
    const feature = await testsFor(conn(), (await resolveTarget(conn(), 'domains/bio'))!, {...LIMIT, tiers: ['feature']});
    expect(titles(feature.sections)).toEqual(['features', 'feature', 'scenarios', 'automations', 'run']);
    const whole = await testsFor(conn(), (await resolveTarget(conn(), 'domains/bio'))!, LIMIT);
    expect(whole.sections.find((s) => s.title === 'immediate')!.empty).toBe('immediate (0): none of the targets is a source file or a declaration');
  });

  withKuzu('gives impact the direct importers, the reachable files and the citing documents of a file', async () => {
    const answer = await impact(conn(), (await resolveTarget(conn(), UTILS))!, LIMIT);
    expect(titles(answer.sections)).toEqual(['features', 'owners', 'evidence', 'cites', 'documents', 'work', 'importers', 'reachable']);
    expect(rows(answer, 'importers').map((r) => r.caller)).toEqual([
      'file:public/packages/Demo/src/package-api.ts', 'file:public/packages/Demo/src/tests/utils-tests.ts', 'file:public/packages/Demo/src/viewer.ts']);
    expect(rows(answer, 'reachable').map((r) => r.importer)).toEqual(['file:public/packages/Demo/src/deep.ts', 'file:public/packages/Demo/src/tests/deep-tests.ts']);
    expect(rows(answer, 'reachable')[0].via).toBe(`${UTILS} → imports ← file:public/packages/Demo/src/viewer.ts → imports ← file:public/packages/Demo/src/deep.ts`);
    expect(rows(answer, 'cites')).toEqual([]);
    const decl = (await run(conn(), `MATCH (n)-[:DECLARES]->(d) WHERE n.id STARTS WITH 'file:public/packages/Demo/src/' AND d.id STARTS WITH 'decl:' RETURN d.id AS id ORDER BY id LIMIT 1`)).rows[0];
    const called = await impact(conn(), (await resolveTarget(conn(), String(decl.id)))!, LIMIT);
    expect(titles(called.sections)).toContain('cites');
    expect(titles(called.sections).slice(-1)).toEqual(['callers']);
  });

  withKuzu('groups the playwright run rows by the folder holding playwright.config.ts', async () => {
    const tests = (await run(conn(), `MATCH (t:Artifact) WHERE t.type = 'test' AND t.framework = 'playwright' AND t.skipped = false AND t.path STARTS WITH 'public/packages/Tested/' ` +
      'RETURN t.framework AS framework, t.name AS name, t.path AS path, t.category AS category, t.dynamic AS dynamic ORDER BY name')).rows;
    expect(tests.length).toBeGreaterThan(1);
    expect(runRows(tests, index.repo).rows).toEqual([
      {framework: 'playwright', cwd: 'public/packages/Tested/playwright', command: 'npx playwright test basic.test.ts', tests: tests.length, names: tests.map((t) => t.name)},
    ]);
    expect(runRows(tests.slice(0, 1), index.repo).rows[0].command).toBe(`npx playwright test basic.test.ts --grep "${tests[0].name}"`);
    expect(runRows(tests.slice(0, 1)).rows[0].cwd).toBe('public/packages/Tested');
  });

  withKuzu('scores the history replay: a changed test file in the immediate tier is a hit, one only reachable a linked hit, the rest a miss', async () => {
    const sha = (c: string) => c.repeat(40);
    const commits: Commit[] = [
      {sha: sha('a'), subject: 'utils and its test', files: ['public/packages/Demo/src/utils.ts', 'public/packages/Demo/src/tests/utils-tests.ts']},
      {sha: sha('b'), subject: 'utils and a test three imports away', files: ['public/packages/Demo/src/utils.ts', 'public/packages/Demo/src/tests/deep-tests.ts']},
      {sha: sha('c'), subject: 'base and a test behind the entry file', files: ['public/packages/Demo/src/base.ts', 'public/packages/Demo/src/tests/api-tests.ts', 'README.md']},
      {sha: sha('d'), subject: 'a source the index does not know', files: ['public/packages/Demo/src/new.ts', 'public/packages/Demo/src/tests/utils-tests.ts']},
      {sha: sha('e'), subject: 'docs and a test', files: ['public/help/a.md', 'public/packages/Demo/src/tests/utils-tests.ts']},
      {sha: sha('f'), subject: 'a source alone', files: ['public/packages/Demo/src/utils.ts']},
    ];
    const report = await replay(conn(), commits, {dart: 'ok'});
    expect(report.summary).toMatch(/^3 of 6 commits changed a source file the index knows and a test file; every changed test file is in the immediate tier for 33% \(1\/3\), in immediate or reachable for 67% \(2\/3\); \d+\.\d s\.$/);
    expect(report.notes).toEqual(['1 commit with a source and a test file skipped: none of the source files is in the index (new or deleted since the build); 0 of the scored commits have some source files the index does not know']);
    expect(rows(report, 'frameworks')).toEqual([{framework: 'dg', test_files: 3, immediate: 1, linked: 2, immediate_rate: '33%', linked_rate: '67%'}]);
    expect(rows(report, 'misses')).toEqual([
      {commit: 'bbbbbbbbbb', subject: 'utils and a test three imports away', sources: 'public/packages/Demo/src/utils.ts', missed: 'public/packages/Demo/src/tests/deep-tests.ts'},
      {commit: 'cccccccccc', subject: 'base and a test behind the entry file', sources: 'public/packages/Demo/src/base.ts', missed: 'public/packages/Demo/src/tests/api-tests.ts'},
    ]);
  });
});

describe('the test-unit table (conventions.md §8.1)', () => {
  it('names the units a test file stands for: ApiTests for the JS API, a client test for the library its category leads with, datlas tests for grok_shared too', () => {
    expect(testUnitsOf('public/packages/ApiTests/src/dataframe/dataframe.ts', 'Dataframe')).toEqual(['public/js-api', 'public/packages/ApiTests']);
    expect(testUnitsOf('core/client/xamgle/lib/src/tests/viewers_test.dart', 'd4 | Viewers')).toEqual(['core/client/d4']);
    expect(testUnitsOf('core/client/xamgle/lib/src/tests/frames_test.dart', 'ddt | Frames')).toEqual(['core/shared/ddt']);
    expect(testUnitsOf('core/client/xamgle/lib/src/tests/browse_test.dart', 'Browse | Tree')).toEqual(['core/client/xamgle']);
    expect(testUnitsOf('core/client/xamgle/lib/src/tests/smoke_test.dart')).toEqual(['core/client/xamgle']);
    expect(testUnitsOf('core/server/datlas/test/services/domain_model_test.dart')).toEqual(['core/server/datlas', 'core/shared/grok_shared']);
    expect(testUnitsOf('public/packages/Chem/src/tests/a-tests.ts', 'Chem: A')).toEqual(['public/packages/Chem']);
    expect(testUnitsOf('public/playwright-public/Viewers/grid.test.ts')).toEqual([]);
  });
});

describe('the tier list', () => {
  it('reads a tier, linked, all or a comma list, in tier order, and refuses anything else', () => {
    expect(parseTiers('immediate')).toEqual(['immediate']);
    expect(parseTiers('linked')).toEqual(['immediate', 'reachable']);
    expect(parseTiers('all')).toEqual(['immediate', 'reachable', 'feature']);
    expect(parseTiers('feature, immediate')).toEqual(['immediate', 'feature']);
    expect(parseTiers('immediate,linked')).toBeUndefined();
    expect(parseTiers('nope')).toBeUndefined();
    expect(parseTiers('')).toBeUndefined();
  });
});

describe('the run rows', () => {
  const dg = (name: string, category = 'Chem: A', dynamic = false) => ({framework: 'dg', name, path: 'public/packages/Chem/src/tests/a-tests.ts', category, dynamic});

  it('runs one DG test by name and a category whole when it has several or a dynamic one', () => {
    expect(runRows([dg('one')]).rows).toEqual([{framework: 'dg', cwd: 'public/packages/Chem', command: 'grok test --category "Chem: A" --test "one"', tests: 1, names: ['one']}]);
    expect(runRows([dg('one'), dg('two'), dg('three', 'Chem: B')]).rows.map((r) => r.command)).toEqual(['grok test --category "Chem: A"', 'grok test --category "Chem: B" --test "three"']);
    expect(runRows([dg('template…', 'Chem: A', true)]).rows[0].command).toBe('grok test --category "Chem: A"');
    expect(runRows([dg('say "hi"')]).rows[0].command).toBe('grok test --category "Chem: A" --test "say \\"hi\\""');
  });

  it('sends client tests to DevTools, Dart tests to their package, and vitest cases to the nearest package.json', () => {
    const xamgle = (category: string | null, name = 'renders', dynamic = false) => ({framework: 'xamgle', name, path: 'core/client/xamgle/lib/src/tests/grid_tests.dart', category, dynamic});
    expect(runRows([xamgle('Browse | Tree')]).rows[0]).toEqual({framework: 'xamgle', cwd: 'public/packages/DevTools', command: 'grok test --category "Core: Browse: Tree" --test "renders"', tests: 1, names: ['renders']});
    expect(runRows([xamgle('d4 | Grid | Rendering')]).rows[0]).toMatchObject({command: 'grok test --category "Core: d4: Grid: Rendering" --test "renders"'});
    expect(runRows([xamgle('ddt | Frames'), xamgle('ddt | Frames', 'again')]).rows[0]).toMatchObject({command: 'grok test --category "Core: ddt: Frames"', tests: 2});
    // one row per DevTools category; a test in no category sits right under Core
    expect(runRows([xamgle('d4 | Grid'), xamgle('d4 | Legend'), xamgle(null, 'smoke')]).rows.map((r) => r.command)).toEqual([
      'grok test --category "Core: d4: Grid" --test "renders"', 'grok test --category "Core: d4: Legend" --test "renders"', 'grok test --category "Core" --test "smoke"']);
    // a client test registered in a loop has no name to run by: skipped, and said so
    const mixed = runRows([xamgle('d4 | Grid'), xamgle('d4 | Grid', 'loop…', true)]);
    expect(mixed.rows[0]).toMatchObject({command: 'grok test --category "Core: d4: Grid" --test "renders"', tests: 1, note: '1 dynamic test skipped: registered in a loop, no runnable name'});
    expect(mixed.notes).toEqual([]);
    const only = runRows([xamgle('d4 | Grid', 'loop…', true)]);
    expect(only.rows).toEqual([]);
    expect(only.notes).toEqual(['run: Core: d4: Grid: 1 dynamic test skipped: registered in a loop, no runnable name']);
    const dart = (p: string, name: string) => ({framework: 'dart', name, path: p, category: 'legend'});
    expect(runRows([dart('core/client/d4/test/legend/a_test.dart', 'fits')]).rows[0]).toEqual({framework: 'dart', cwd: 'core/client/d4', command: 'pub run test test/legend/a_test.dart -n "fits"', tests: 1, names: ['fits']});
    expect(runRows([dart('core/client/libs/dock_spawn/test/a_test.dart', 'x'), dart('core/client/libs/dock_spawn/test/b_test.dart', 'y')]).rows[0])
      .toMatchObject({cwd: 'core/client/libs/dock_spawn', command: 'pub run test test/a_test.dart test/b_test.dart', tests: 2});
    const node = {framework: 'node', name: 'parses', path: 'public/tools/bin/__tests__/kg.test.ts', category: null};
    const repo = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-kg-run-'));
    write(repo, 'public/tools/package.json', '{}');
    expect(runRows([node], repo).rows[0]).toMatchObject({cwd: 'public/tools', command: 'npx vitest run bin/__tests__/kg.test.ts'});
    const java = runRows([{framework: 'java', name: 'Postgres', path: 'public/connectors/src/test/PostgresTest.java', category: null}]);
    expect(java.rows[0]).toMatchObject({framework: 'java', cwd: '', command: '?', tests: 1});
    expect(java.notes).toEqual(["run: no runner known for framework 'java'"]);
  });

  it('runs the whole suite of a unit when at least half of its test files of a framework are selected, and says so', () => {
    const dart = (file: string, name: string) => ({framework: 'dart', name, path: `core/shared/ddt/test/${file}`, category: null});
    const suites = new Map([['dart core/shared/ddt', 4], ['dg public/packages/Chem', 10], ['node public/tools', 2], ['playwright public/packages/Chem', 1]]);
    const selected = [dart('a_test.dart', 'a'), dart('b_test.dart', 'b'), dart('b_test.dart', 'b2'), dart('c_test.dart', 'c'), dg('one'), dg('two', 'Chem: B'),
      {framework: 'node', name: 'parses', path: 'public/tools/bin/__tests__/kg.test.ts', category: null},
      {framework: 'playwright', name: 'opens', path: 'public/packages/Chem/playwright/a.test.ts', category: null}];
    const {rows, notes} = runRows(selected, undefined, suites);
    expect(rows).toEqual([
      {framework: 'dart', cwd: 'core/shared/ddt', command: 'pub run test', tests: 4, names: ['a', 'b', 'b2'], note: 'whole suite: 3 of 4 test files selected'},
      {framework: 'dg', cwd: 'public/packages/Chem', command: 'grok test --category "Chem: A" --test "one"', tests: 1, names: ['one']},
      {framework: 'dg', cwd: 'public/packages/Chem', command: 'grok test --category "Chem: B" --test "two"', tests: 1, names: ['two']},
      {framework: 'node', cwd: 'public/tools', command: 'npx vitest run', tests: 1, names: ['parses'], note: 'whole suite: 1 of 2 test files selected'},
      {framework: 'playwright', cwd: 'public/packages/Chem', command: 'npx playwright test', tests: 1, names: ['opens'], note: 'whole suite: 1 of 1 test files selected'},
    ]);
    expect(notes).toEqual([]);
    // the same selection without the sizes, or with a unit the index does not know, runs by file as before
    expect(runRows(selected.slice(0, 1)).rows[0].command).toBe('pub run test test/a_test.dart -n "a"');
    expect(runRows(selected.slice(0, 1), undefined, new Map()).rows[0].command).toBe('pub run test test/a_test.dart -n "a"');
  });
});

describe('the change set', () => {
  /** A monorepo with `public/` as a repository of its own, recorded as a gitlink, both on `master` with one commit. */
  function makeRepo(): string {
    const repo = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-kg-changes-'));
    write(repo, 'core/a.dart', 'a\n');
    write(repo, 'public/packages/P/src/p.ts', 'p\n');
    for (const dir of [path.join(repo, 'public'), repo]) {
      git(dir, 'init', '-q');
      git(dir, 'symbolic-ref', 'HEAD', 'refs/heads/master');
      git(dir, 'add', '-A');
      git(dir, 'commit', '-q', '-m', 'base');
    }
    return repo;
  }

  it('diffs both repositories against the merge base with master, working tree and untracked files included', () => {
    const repo = makeRepo();
    const base = git(repo, 'rev-parse', 'HEAD');
    const publicBase = git(path.join(repo, 'public'), 'rev-parse', 'HEAD');
    git(repo, 'checkout', '-q', '-b', 'work');
    fs.appendFileSync(path.join(repo, 'core', 'a.dart'), 'more\n');
    git(repo, 'commit', '-q', '-am', 'committed');
    fs.appendFileSync(path.join(repo, 'public', 'packages', 'P', 'src', 'p.ts'), 'more\n');
    write(repo, 'core/new.md', 'new\n');
    write(repo, 'public/packages/P/src/q.ts', 'q\n');
    const set = changeSet(repo);
    expect(set.files).toEqual([
      {path: 'core/a.dart', repo: 'core'}, {path: 'core/new.md', repo: 'core'},
      {path: 'public/packages/P/src/p.ts', repo: 'public'}, {path: 'public/packages/P/src/q.ts', repo: 'public'}]);
    expect(set.base).toEqual({core: base, public: publicBase});
    expect(set.notes).toEqual([]);
    // against a ref: the public base is the gitlink that ref recorded, so the committed core change drops out and the public one stays
    const at = changeSet(repo, 'HEAD');
    expect(at.files.map((f) => f.path)).toEqual(['core/new.md', 'public/packages/P/src/p.ts', 'public/packages/P/src/q.ts']);
    expect(at.base).toEqual({core: 'HEAD', public: publicBase});
  });

  it('falls back to the parent commit, and says so, when no master exists', () => {
    const repo = makeRepo();
    git(repo, 'branch', '-m', 'master', 'trunk');
    fs.appendFileSync(path.join(repo, 'core', 'a.dart'), 'more\n');
    git(repo, 'commit', '-q', '-am', 'second');
    const set = changeSet(repo);
    expect(set.files.map((f) => f.path)).toEqual(['core/a.dart']);
    expect(set.base.core).toBe('HEAD~1');
    expect(set.notes).toEqual(['core: neither master nor origin/master exists; the base is HEAD~1']);
  });
});
