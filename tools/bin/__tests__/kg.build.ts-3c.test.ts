/// `grok kg build` WO-3c (build-plan.md): the ts-tests, ts-samples, ts-changelog and docs extractors against the mini
/// monorepo under fixtures/kg/build, and the text parsers they share.
import {describe, it, expect} from 'vitest';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {fileURLToPath} from 'url';
import {parseDgTests, parsePlaywrightTests, blankComments} from '../utils/kg/build/extract/ts/tests';
import {parseChangelog} from '../utils/kg/build/extract/ts/changelog';
import {parseSampleHeader} from '../utils/kg/build/extract/ts/samples';
import {idTokens, ticketTokens, leadingId} from '../utils/kg/build/extract/markers';
import {kebab} from '../utils/kg/ids';
import {kg} from '../commands/kg';
import {copyFixture, buildFixture} from './kg-fixture';

const fixture = path.join(path.dirname(fileURLToPath(import.meta.url)), 'fixtures', 'kg', 'build');
const TESTS = 'public/packages/Tested/src/tests/demo-tests.ts';
const PLAYWRIGHT = 'public/packages/Tested/playwright/basic.test.ts';
const SPEC = 'public/packages/UsageAnalysis/files/TestTrack/Viewers/ScatterPlot/scatterplot-legend-spec.ts';
const TRACK_TEST = 'public/packages/UsageAnalysis/files/TestTrack/Connections/initial runs/basic.test.ts';
const LEGACY = 'public/packages/UsageAnalysis/files/TestTrack/Viewers/ScatterPlot/scatterplot-legend.md';
const BIO = 'public/help/domains/bio/bio.md';
const PROJECT = 'public/help/datagrok/project.md';
const SEQUENCES = 'public/help/domains/bio/sequences.md';


const graph = buildFixture(copyFixture('build'), 'homes,ts-packages,ts-declarations,ts-tests,ts-samples,ts-changelog,docs');
const byId = (rows: any[], id: string) => rows.find((r) => r.id === id);
const edges = (rows: any[], from?: string, to?: string) => rows.filter((e) => (from === undefined || e.from === from) && (to === undefined || e.to === to));

describe('text parsers (build-plan.md WO-3c)', () => {
  it('reads DG tests in order under the last category, the trailing options object, a run-time title and a conditional skip', () => {
    const tests = parseDgTests([
      'test(\'unregistered\', async () => {});', '// test(\'out\', async () => {});', '/* test(\'out too\', async () => {}); */',
      'category(\'A\', () => {', '  test(\'one\', async () => { if (re.test(\'x\')) call(1, {timeout: 5}); }, {timeout: 60000, skipReason: \'flaky\'});',
      '  test("two", async () => {}, {benchmark: true, tags: ["~x/y", \'z\']});', '  test(`three ${n}`, async () => {});', '});',
      'category(\'B\', () => { test(\'four\', fn, {skipReason: reason}); });',
      'category(\'C \', () => { test(\'five: \' + name, async () => {}); });',
    ].join('\n'));
    const plain = {dynamic: false, skipConditional: false, skipReason: undefined, benchmark: false, tags: []};
    expect(tests).toEqual([
      {...plain, category: 'A', name: 'one', skipReason: 'flaky'},
      {...plain, category: 'A', name: 'two', benchmark: true, tags: ['~x/y', 'z']},
      // a template and a concatenation both name a registration site: the literal head with an ellipsis, marked dynamic
      {...plain, category: 'A', name: 'three…', dynamic: true},
      // skipReason is an expression, so the row may not say the test is skipped
      {...plain, category: 'B', name: 'four', skipConditional: true},
      {...plain, category: 'C', name: 'five:…', dynamic: true},
    ]);
  });

  it('nests Playwright tests under their describes and marks skip, fixme and describe.skip', () => {
    const tests = parsePlaywrightTests([
      'test.describe(\'Outer\', () => {', '  test.beforeEach(async ({page}) => {});', '  test(\'a\', async ({page}) => { const s = "}"; });',
      '  test.describe.serial(\'Inner\', () => {', '    test.fixme(\'b\', async () => {});', '  });', '});',
      'test.skip(\'c\', async () => {});', 'test.describe.skip(\'Parked\', () => { test(\'d\', async () => {}); });', 'test.describe.configure({mode: \'serial\'});',
      'test(\'e\', {tag: \'@smoke\'}, async ({page}) => {});',
    ].join('\n'));
    expect(tests).toEqual([
      {describes: ['Outer'], title: 'a', skipped: false},
      {describes: ['Outer', 'Inner'], title: 'b', skipped: true},
      {describes: [], title: 'c', skipped: true},
      {describes: ['Parked'], title: 'd', skipped: true},
      {describes: [], title: 'e', skipped: false},
    ]);
    expect(blankComments('a // b\nc /* d\ne */ f "//x" \'/*y\'')).toBe('a     \nc     \n     f "//x" \'/*y\'');
  });

  it('reads changelog bullets under version headings, v.next as version 0, continuation lines joined', () => {
    expect(parseChangelog('# X changelog\n\n## v.next\n\n* one\n  continued\n- two\n\n## 1.2.3 (2026-01-31)\n\n* three\n\n## 1.2.3-rc (2026-01-30)\n\n* four\n\n## Unreleased\n\n* ignored\n\n## 1.0.0\n\n* five\n\n## 0.9.0 (2025-27-01)\n\n* six\n')).toEqual([
      {version: '0', date: undefined, n: 1, text: 'one continued'}, {version: '0', date: undefined, n: 2, text: 'two'},
      {version: '1.2.3', date: '2026-01-31', n: 1, text: 'three'}, {version: '1.2.3', date: '2026-01-30', n: 2, text: 'four'},
      {version: '1.0.0', date: undefined, n: 1, text: 'five'}, {version: '0.9.0', date: undefined, n: 1, text: 'six'},
    ]);
  });

  it('reads a sample header: keys, the first prose line, the block text; python uses #', () => {
    expect(parseSampleHeader('\n//api: DG.A, ui.b\n//help-url: https://datagrok.ai/help/x\n// https://example.com is not a key\n// Shows ~a/b.\nlet x = 1;\n', '//'))
      .toEqual({keys: {api: 'DG.A, ui.b', 'help-url': 'https://datagrok.ai/help/x'}, description: 'https://example.com is not a key',
        text: 'api: DG.A, ui.b\nhelp-url: https://datagrok.ai/help/x\nhttps://example.com is not a key\nShows ~a/b.', body: 'let x = 1;\n'});
    const python = parseSampleHeader('#api: grok.shell.info\nprint(1)\n', '#');
    expect(python.keys).toEqual({api: 'grok.shell.info'});
    expect(python.description).toBeUndefined();
  });

  it('finds ~id and ticket tokens and leaves paths, strikethrough and numbers alone', () => {
    expect([...idTokens('see ~visualize/viewers/scatter-plot#legend and `~C:dataframe`, ~/.claude, ~~gone~~, ~1.5s, foo~bar/baz, ~visualize/viewers/scatter-plot').keys()])
      .toEqual(['visualize/viewers/scatter-plot', 'C:dataframe']);
    expect([...ticketTokens('GROK-1 and GROK-1, [#12](https://github.com/datagrok-ai/public/issues/12), #13, GROK-123')]).toEqual([['GROK-1', 2], ['GROK-123', 1], ['gh:public#12', 1]]);
    expect([leadingId('~domains/bio renders'), leadingId(' ~C:x'), leadingId('plain ~a/b')]).toEqual(['domains/bio', 'C:x', undefined]);
    expect(['ScatterPlot', 'MLMethods', 'initial runs', 'EDA', 'scatterplot-legend', 'User groups'].map(kebab)).toEqual(['scatter-plot', 'ml-methods', 'initial-runs', 'eda', 'scatterplot-legend', 'user-groups']);
  });
});

describe('ts-tests extractor (build-plan.md WO-3c)', () => {
  it('emits DG tests with ids, levels, options and suites; a test before any category or in a comment is not one', async () => {
    const {rows} = await graph;
    const tests = rows('nodes/test').filter((t) => t.framework === 'dg');
    expect(tests.map((t) => t.id)).toEqual([
      'test:dg:public/packages/ApiTests/src/tests/shell.ts#Shell/windows',
      `test:dg:${TESTS}#Tested: Utils/conditional skip`, `test:dg:${TESTS}#Tested: Utils/tagged`, `test:dg:${TESTS}#Tested: Utils/template…`,
      `test:dg:${TESTS}#~domains/bio/detects sequences`, `test:dg:${TESTS}#~domains/bio/renders slowly`, `test:dg:${TESTS}#~domains/bio/skipped one`,
    ]);
    expect(byId(tests, 'test:dg:public/packages/ApiTests/src/tests/shell.ts#Shell/windows')).toMatchObject({level: 'api', category: 'Shell', skipped: false, benchmark: false, provenance: 'ast'});
    expect(byId(tests, `test:dg:${TESTS}#~domains/bio/detects sequences`)).toEqual({
      id: `test:dg:${TESTS}#~domains/bio/detects sequences`, type: 'test', name: 'detects sequences', batch: expect.any(String), benchmark: false, category: '~domains/bio',
      framework: 'dg', level: 'unit', path: TESTS, provenance: 'ast', dynamic: false, skipped: false, skip_conditional: false, suite: 'suite:dg:Tested:~domains/bio',
      source_layer: 'public', status: 'active', visibility: 'public',
    });
    expect(byId(tests, `test:dg:${TESTS}#~domains/bio/renders slowly`)).toMatchObject({benchmark: true, tags: ['render', 'slow'], skipped: false});
    expect(byId(tests, `test:dg:${TESTS}#~domains/bio/skipped one`)).toMatchObject({skipped: true, skip_reason: 'GROK-100: flaky on CI'});
    const conditional = byId(tests, `test:dg:${TESTS}#Tested: Utils/conditional skip`);
    expect(conditional).toMatchObject({skipped: false, skip_conditional: true});
    expect(conditional).not.toHaveProperty('skip_reason');
    expect(byId(tests, `test:dg:${TESTS}#Tested: Utils/template…`)).toMatchObject({dynamic: true, name: 'template…'});
    expect(rows('nodes/test-suite').filter((s) => s.framework === 'dg')).toEqual([
      expect.objectContaining({id: 'suite:dg:ApiTests:Shell', name: 'Shell', package: 'pkg:ApiTests', provenance: 'ast'}),
      expect.objectContaining({id: 'suite:dg:Tested:Tested: Utils', name: 'Tested: Utils', package: 'pkg:Tested'}),
      expect.objectContaining({id: 'suite:dg:Tested:~domains/bio', name: '~domains/bio', package: 'pkg:Tested'}),
    ]);
    expect(byId(tests, `test:dg:${TESTS}#Tested: Utils/tagged`).suite).toBe('suite:dg:Tested:Tested: Utils');
    expect(rows('edges/suite')).toContainEqual(expect.objectContaining({type: 'ref', name: 'suite', from: `test:dg:${TESTS}#Tested: Utils/tagged`, to: 'suite:dg:Tested:Tested: Utils', derived_by: 'ast', confidence: 1}));
    expect(rows('edges/package').filter((e) => e.from.startsWith('suite:')).map((e) => [e.from, e.to])).toContainEqual(['suite:dg:Tested:~domains/bio', 'pkg:Tested']);
  });

  it('turns a ~id category or first tag into a tests edge', async () => {
    const {rows} = await graph;
    expect(rows('edges/tests').map((e) => [e.from, e.to, e.derived_by, e.confidence])).toEqual([
      [`test:dg:${TESTS}#Tested: Utils/tagged`, 'domains/bio', 'annotation', 1],
      [`test:dg:${TESTS}#~domains/bio/detects sequences`, 'domains/bio', 'annotation', 1],
      [`test:dg:${TESTS}#~domains/bio/renders slowly`, 'domains/bio', 'annotation', 1],
      [`test:dg:${TESTS}#~domains/bio/skipped one`, 'domains/bio', 'annotation', 1],
      [`test:playwright:${PLAYWRIGHT}#Basic/~domains/bio sequence view opens`, 'domains/bio', 'annotation', 1],
    ]);
    expect(rows('edges/tests').every((e) => e.kind === undefined)).toBe(true);
  });

  it('never makes a feature out of a ~id marker no home declares: the token is counted and the edge dropped', async () => {
    const {rows, problems} = await graph;
    expect(problems.unresolved_ids).toContain(`${TESTS}: ~nowhere/thing resolves to no home document`);
    expect(rows('nodes/feature').map((f) => f.id)).not.toContain('nowhere/thing');
    expect(edges(rows('edges/tests'), `test:dg:${TESTS}#Tested: Utils/template \${name}`)).toEqual([]);
  });

  it('emits Playwright tests named describe > title, one suite per file, skipped from skip and describe.skip', async () => {
    const {rows} = await graph;
    const tests = rows('nodes/test').filter((t) => t.framework === 'playwright');
    expect(tests.map((t) => [t.id, t.name, t.skipped])).toEqual([
      [`test:playwright:${PLAYWRIGHT}#Basic > Inner/nested title`, 'Basic > Inner > nested title', false],
      [`test:playwright:${PLAYWRIGHT}#Basic/~domains/bio sequence view opens`, 'Basic > ~domains/bio sequence view opens', false],
      [`test:playwright:${PLAYWRIGHT}#Parked/parked test`, 'Parked > parked test', true],
      [`test:playwright:${PLAYWRIGHT}#basic.test/later`, 'later', true],
      [`test:playwright:${TRACK_TEST}#Initial runs/opens a connection`, 'Initial runs > opens a connection', false],
      [`test:playwright:${SPEC}#Scatter plot legend/color legend`, 'Scatter plot legend > color legend', false],
      [`test:playwright:${SPEC}#Scatter plot legend/marker legend`, 'Scatter plot legend > marker legend', false],
    ]);
    expect(tests[0]).toMatchObject({level: 'e2e', category: 'Basic > Inner', path: PLAYWRIGHT, provenance: 'ast'});
    expect(tests[3].category).toBeUndefined();
    expect(rows('nodes/test-suite').filter((s) => s.framework === 'playwright')).toEqual([
      expect.objectContaining({id: `suite:playwright:${PLAYWRIGHT}`, name: 'basic.test.ts', path: PLAYWRIGHT, package: 'pkg:Tested', provenance: 'filesystem'}),
      expect.objectContaining({id: `suite:playwright:${TRACK_TEST}`, name: 'basic.test.ts', package: 'pkg:UsageAnalysis'}),
      expect.objectContaining({id: `suite:playwright:${SPEC}`, name: 'scatterplot-legend-spec.ts', package: 'pkg:UsageAnalysis'}),
    ]);
    expect(edges(rows('edges/suite'), undefined, `suite:playwright:${PLAYWRIGHT}`)).toHaveLength(4);
    expect(rows('edges/suite')).toHaveLength(14);
  });
});

describe('ts-samples extractor (build-plan.md WO-3c)', () => {
  it('reads api_members from an //api: header, else from DG./ui./grok. usages, as a set', async () => {
    const {rows} = await graph;
    expect(byId(rows('nodes/sample'), 'sample:dapi/projects-list')).toEqual({
      id: 'sample:dapi/projects-list', type: 'sample', name: 'projects-list', api_members: ['DG.HttpDataSource.filter', 'DG.HttpDataSource.list'], batch: expect.any(String),
      description: 'Lists demo projects; the canonical sample for ~domains/bio project browsing.', folder: 'dapi', language: 'js', path: 'public/packages/ApiSamples/scripts/dapi/projects-list.js',
      provenance: 'annotation', source_layer: 'public', status: 'active', visibility: 'public',
    });
    expect(byId(rows('nodes/sample'), 'sample:ui/hello')).toMatchObject({api_members: ['DG.Viewer.fromType', 'grok.data.demo.demog', 'grok.shell.newView', 'ui.divText'], folder: 'ui', language: 'js', provenance: 'ast'});
    expect(byId(rows('nodes/sample'), 'sample:ui/info')).toMatchObject({api_members: ['grok.shell.info'], language: 'python', description: 'A Python sample.'});
    expect(rows('nodes/sample')).toHaveLength(4);
  });

  it('mentions the help page of //help-url: (absolute or /help/ form, anchor dropped), demonstrates a ~id in the header, counts a missing page', async () => {
    const {rows, problems, manifest} = await graph;
    expect(edges(rows('edges/mentions'), 'sample:dapi/projects-list')).toEqual([expect.objectContaining({to: `doc:${PROJECT}`, derived_by: 'annotation', evidence: ['public/packages/ApiSamples/scripts/dapi/projects-list.js']})]);
    expect(edges(rows('edges/mentions'), 'sample:ui/info').map((e) => e.to)).toEqual([`doc:${PROJECT}`]);
    expect(edges(rows('edges/mentions'), 'sample:misc/missing')).toEqual([]);
    expect(rows('edges/demonstrates')).toEqual([expect.objectContaining({from: 'sample:dapi/projects-list', to: 'domains/bio', derived_by: 'annotation', confidence: 1})]);
    expect(problems.unresolved_ids).toContain('public/packages/ApiSamples/scripts/misc/missing.js: help-url https://datagrok.ai/help/nowhere/at-all names no page under public/help');
    expect(manifest.sources['ts-samples']).toBe('partial');
  });

  it('never makes a feature out of a ~id no home declares: the token is counted and the edge dropped', async () => {
    const {rows, problems} = await graph;
    expect(problems.unresolved_ids).toContain('public/packages/ApiSamples/scripts/misc/missing.js: ~item resolves to no home document');
    expect(rows('nodes/feature').map((f) => f.id)).not.toContain('item');
    expect(edges(rows('edges/demonstrates'), 'sample:misc/missing')).toEqual([]);
  });

  it('uses the JS API declarations the code calls, the header left out of the count', async () => {
    const {rows} = await graph;
    expect(rows('edges/uses').filter((e) => e.from.startsWith('sample:')).map((e) => [e.from, e.to, e.kind, e.count])).toEqual([
      ['sample:dapi/projects-list', 'decl:public/js-api/src/shell.ts#Shell.info', 'function', 1],
      ['sample:misc/missing', 'decl:public/js-api/src/shell.ts#Shell.info', 'function', 1],
      ['sample:ui/hello', 'decl:public/js-api/src/viewer.ts#Viewer', 'class', 1],
      ['sample:ui/info', 'decl:public/js-api/src/shell.ts#Shell.info', 'function', 1],
    ]);
  });
});

describe('ts-changelog extractor (build-plan.md WO-3c)', () => {
  it('emits one entry per bullet with package, version, date and text; v.next is version 0', async () => {
    const {rows} = await graph;
    expect(rows('nodes/changelog-entry').map((e) => [e.id, e.version, e.date])).toEqual([
      ['chg:Tested:0:1', '0', undefined], ['chg:Tested:0:2', '0', undefined], ['chg:Tested:1.0.0:1', '1.0.0', '2026-01-01'], ['chg:Tested:1.1.0:1', '1.1.0', '2026-02-01'], ['chg:Tested:1.1.0:2', '1.1.0', '2026-02-01'],
    ]);
    expect(byId(rows('nodes/changelog-entry'), 'chg:Tested:0:1')).toMatchObject({
      name: 'GROK-20753: Added `~domains/bio` sequence rendering for [#12](https://github.com/datagrok-ai/public…',
      text: 'GROK-20753: Added `~domains/bio` sequence rendering for [#12](https://github.com/datagrok-ai/public/issues/12)',
      package: 'pkg:Tested', path: 'public/packages/Tested/CHANGELOG.md', provenance: 'annotation', source_layer: 'public',
    });
    expect(byId(rows('nodes/changelog-entry'), 'chg:Tested:0:2').text).toBe('Fixed a crash in ~C:dataframe handling when the table is empty');
  });

  it('mentions GROK-n and linked GitHub issues as ticket stubs, and a ~id that is a feature becomes changes with the verb as kind', async () => {
    const {rows, problems, manifest} = await graph;
    expect(rows('edges/mentions').filter((e) => e.from.startsWith('chg:')).map((e) => [e.from, e.to, e.count])).toEqual([
      ['chg:Tested:0:1', 'GROK-20753', 1], ['chg:Tested:0:1', 'gh:public#12', 1], ['chg:Tested:1.1.0:1', 'GROK-1', 1],
    ]);
    expect(rows('edges/changes').map((e) => [e.from, e.to, e.kind, e.derived_by, e.confidence])).toEqual([
      ['chg:Tested:0:1', 'domains/bio', 'added', 'annotation', 1], ['chg:Tested:1.1.0:2', 'domains/bio', 'removed', 'annotation', 1],
    ]);
    expect(byId(rows('nodes/ticket'), 'gh:public#12')).toMatchObject({tracker: 'github', key: '#12', kind: 'unknown', state: 'open', status: 'proposed', provenance: 'annotation'});
    expect(byId(rows('nodes/ticket'), 'GROK-20753')).toMatchObject({tracker: 'jira', key: 'GROK-20753'});
    expect(problems.unresolved_ids).toContain('public/packages/Tested/CHANGELOG.md (0 #2): ~C:dataframe resolves to no home document');
    expect(manifest.sources['ts-changelog']).toBe('partial');
  });
});

describe('docs extractor (build-plan.md WO-3c)', () => {
  it('emits a doc-page per markdown file with kind by folder, title, keywords, mdx and unlisted', async () => {
    const {rows} = await graph;
    expect(rows('nodes/doc-page').map((d) => [d.id, d.kind, d.provenance])).toEqual([
      ['doc:core/client/d4/lib/src/legends/README.md', 'readme', 'annotation'], ['doc:core/client/d4/lib/src/viewers/histogram/CLAUDE.md', 'agent', 'annotation'],
      ['doc:core/client/d4/lib/src/viewers/scatterplot/CLAUDE.md', 'agent', 'annotation'],
      ['doc:core/docs/CACHING.md', 'core-doc', 'annotation'], ['doc:core/docs/NOTES.md', 'core-doc', 'filesystem'], ['doc:core/docs/VIEWERS.md', 'core-doc', 'annotation'],
      [`doc:${PROJECT}`, 'help', 'annotation'], [`doc:${BIO}`, 'help', 'annotation'], [`doc:${SEQUENCES}`, 'help', 'annotation'],
      ['doc:public/help/visualize/viewers/histogram.md', 'help', 'annotation'],
      ['doc:public/packages/Tested/README.md', 'readme', 'filesystem'],
      ['doc:public/packages/UsageAnalysis/files/TestTrack/Connections/initial runs/basic.md', 'other', 'annotation'],
      ['doc:public/packages/UsageAnalysis/files/TestTrack/Viewers/ScatterPlot/scatter-plot-ui.md', 'other', 'annotation'], [`doc:${LEGACY}`, 'other', 'annotation'],
    ]);
    expect(byId(rows('nodes/doc-page'), 'doc:core/docs/NOTES.md')).toMatchObject({name: 'Notes', description: 'Loose notes on the fixture; see ~domains/bio and GROK-42.', mdx: false, unlisted: false, visibility: 'dev', source_layer: 'core'});
    expect(byId(rows('nodes/doc-page'), `doc:${PROJECT}`)).toMatchObject({name: 'Projects', keywords: ['project', 'sharing'], description: 'A project is a collection of entities saved and shared together.', status: 'active'});
    expect(byId(rows('nodes/doc-page'), `doc:${SEQUENCES}`)).toMatchObject({name: 'Sequences', unlisted: true, mdx: true});
  });

  it('merges the page of a home into the node the homes layer stubbed: one node, the home\'s description, status active', async () => {
    const {rows} = await graph;
    const pages = rows('nodes/doc-page').filter((d) => d.id === `doc:${BIO}`);
    expect(pages).toHaveLength(1);
    const feature = byId(rows('nodes/feature'), 'domains/bio');
    expect(feature).toMatchObject({home: BIO, description: 'Sequence analysis for biologics: notation conversion, MSA and activity cliffs on macromolecules.', status: 'active'});
    expect(pages[0]).toMatchObject({name: 'Bioinformatics', kind: 'help', description: feature.description, keywords: ['macromolecules', 'sequences'], status: 'active', provenance: 'annotation'});
    // the help page the home's body cites documents the feature instead of being mentioned by it
    expect(edges(rows('edges/mentions'), `doc:${BIO}`).map((e) => e.to)).toEqual([]);
  });

  it('emits a doc-anchor per #..#### heading with GitHub slugs, -1 for a duplicate and the explicit {#id}; deeper, fenced and unsluggable lines are not anchors', async () => {
    const {rows} = await graph;
    expect(rows('nodes/doc-anchor').filter((a) => a.page === `doc:${PROJECT}`).map((a) => [a.slug, a.depth, a.name])).toEqual([
      ['code-data-d42', 2, 'Code & data (`d42`)'], ['links', 3, 'Links'], ['sharing', 2, 'Sharing'], ['sharing-1', 2, 'Sharing'],
    ]);
    expect(rows('nodes/doc-anchor').some((a) => a.name === 'Соглашения' || a.name === 'Deep')).toBe(false);
    expect(byId(rows('nodes/doc-anchor'), `doc:${PROJECT}#sharing-1`)).toEqual({
      id: `doc:${PROJECT}#sharing-1`, type: 'doc-anchor', name: 'Sharing', batch: expect.any(String), depth: 2, page: `doc:${PROJECT}`, path: PROJECT,
      provenance: 'annotation', slug: 'sharing-1', source_layer: 'public', status: 'active', visibility: 'public',
    });
    expect(rows('nodes/doc-anchor').filter((a) => a.page === `doc:${BIO}`).map((a) => a.slug)).toEqual(['bioinformatics', 'notation-conversion', 'overview']);
    expect(rows('edges/page').filter((e) => e.from === `doc:${PROJECT}#links`)).toEqual([expect.objectContaining({type: 'ref', name: 'page', to: `doc:${PROJECT}`})]);
  });

  it('mentions resolved ~id tokens and tickets outside fences with a count, counts unresolved ids, and leaves documents: to the homes layer', async () => {
    const {rows, problems, manifest} = await graph;
    expect(edges(rows('edges/mentions'), `doc:${SEQUENCES}`).map((e) => [e.to, e.count])).toEqual([['GROK-777', 2], ['domains/bio', 1]]);
    expect(edges(rows('edges/mentions'), 'doc:core/docs/NOTES.md').map((e) => [e.to, e.count])).toEqual([['GROK-42', 1], ['domains/bio', 1]]);
    expect(edges(rows('edges/mentions'), 'doc:public/packages/Tested/README.md').map((e) => e.to)).toEqual(['GROK-42', 'GROK-43']);
    expect(problems.unresolved_ids).toContain(`${SEQUENCES}: ~C:nothing resolves to no home document`);
    expect(rows('edges/documents')).toEqual([
      expect.objectContaining({from: `doc:${PROJECT}`, to: 'domains/bio', derived_by: 'annotation', evidence: [BIO]}),
      expect.objectContaining({from: `doc:${SEQUENCES}`, to: 'domains/bio', audience: 'user', derived_by: 'annotation'}),
      expect.objectContaining({from: 'doc:public/help/visualize/viewers/histogram.md', to: 'visualize/viewers/histogram', derived_by: 'annotation'}),
    ]);
    expect(manifest.sources).toEqual({docs: 'partial', homes: 'ok', 'ts-changelog': 'partial', 'ts-declarations': 'ok', 'ts-packages': 'ok', 'ts-samples': 'partial', 'ts-tests': 'ok'});
  });

  it('maps a legacy Test Track scenario to TS:<kebab folders>/<stem> with its frontmatter, steps, automates from the sibling spec and mentions from related_bugs', async () => {
    const {rows, problems} = await graph;
    const id = 'TS:viewers/scatter-plot/scatterplot-legend';
    expect(byId(rows('nodes/scenario'), id)).toEqual({
      id, type: 'scenario', name: 'Scatter plot legend', batch: expect.any(String), coverage_type: 'smoke', description: 'Legend entries follow the Color and Marker columns.', manual_only: false,
      path: LEGACY, priority: 'p0', provenance: 'annotation', source_layer: 'public', status: 'active', steps: 3, target_layer: 'playwright', visibility: 'public',
    });
    expect(byId(rows('nodes/scenario'), 'TS:connections/initial-runs/basic')).toMatchObject({target_layer: 'apitest', coverage_type: 'regression', steps: 1, manual_only: false});
    expect(byId(rows('nodes/scenario'), 'TS:viewers/scatter-plot/ui')).toMatchObject({manual_only: true, home: 'public/packages/UsageAnalysis/files/TestTrack/Viewers/ScatterPlot/scatter-plot-ui.md'});
    expect(rows('nodes/scenario').map((s) => [s.id, s.status])).toEqual([
      ['TS:connections', 'proposed'], ['TS:connections/initial-runs', 'proposed'], ['TS:connections/initial-runs/basic', 'active'], ['TS:viewers', 'proposed'],
      ['TS:viewers/scatter-plot', 'proposed'], [id, 'active'], ['TS:viewers/scatter-plot/ui', 'active'],
    ]);
    expect(rows('edges/automates').map((e) => [e.from, e.to, e.derived_by])).toEqual([
      [`test:playwright:${TRACK_TEST}#Initial runs/opens a connection`, 'TS:connections/initial-runs/basic', 'annotation'],
      [`test:playwright:${SPEC}#Scatter plot legend/color legend`, id, 'annotation'], [`test:playwright:${SPEC}#Scatter plot legend/marker legend`, id, 'annotation'],
    ]);
    expect(edges(rows('edges/mentions'), id).map((e) => e.to)).toEqual(['GROK-17227', 'GROK-19083']);
    expect(rows('edges/covers').map((e) => [e.from, e.to])).toEqual([['TS:viewers/scatter-plot/ui', 'domains/bio']]);
    expect(problems.unresolved_ids).toContain(`${LEGACY}: realized_as missing-spec.ts is not beside the scenario`);
  });

  it('emits tutorials from the Tutorials tracks with the track name', async () => {
    const {rows} = await graph;
    expect(rows('nodes/tutorial')).toEqual([expect.objectContaining({
      id: 'tutorial:chem/activity-cliffs', name: 'Activity Cliffs', track: 'Cheminformatics', package: 'pkg:Tutorials', path: 'public/packages/Tutorials/src/tracks/chem/tutorials/activity-cliffs.ts',
      description: 'Detects pairs of molecules with similar structures but different activity.', provenance: 'ast',
    })]);
  });
});
