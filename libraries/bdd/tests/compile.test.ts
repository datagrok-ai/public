import assert from 'node:assert/strict';
import {beforeEach, test} from 'node:test';
import {compileFeature} from '../src/compile.js';
import type {Bindings} from '../src/discover.js';
import {parseFeature} from '../src/gherkin.js';
import {StepMatcher} from '../src/match.js';
import {context, dataset, defineParameterType, element, Given, kind, resetRegistry, StepFn, Then, When} from '../src/registry.js';

const ROOT = 'C:/pkg/bdd';
const STEPS = '@datagrok-libraries/bdd/bindings/common/steps';

function bindings(fns: Record<string, StepFn>, registryModules: string[] = []): Bindings {
  const module = {file: 'C:/lib/bindings/common/steps.ts', specifier: STEPS, exports: fns, stepDefs: [], registers: false};
  const modules = [module, ...registryModules.map((s) => ({file: s, specifier: s, exports: {}, stepDefs: [], registers: true}))];
  const exportOf = new Map<StepFn, {module: typeof module; name: string}>();
  for (const [name, fn] of Object.entries(fns))
    exportOf.set(fn, {module, name});
  return {modules, exportOf};
}

const FEATURE = `@demo @realizes:u2.dialog
Feature: Toolbox
  Background:
    Given user opens spgi dataset

  Scenario: Add a viewer
    When user clicks on scatter plot icon on toolbox
    Then scatter plot viewer should be visible

  Scenario Outline: Several viewers
    When user clicks on <viewer> icon on toolbox
    Then <viewer> viewer should be visible
    Examples:
      | viewer    |
      | histogram |
      | bar chart |
`;

let fns: Record<string, StepFn>;

beforeEach(() => {
  resetRegistry();
  defineParameterType({name: 'element', regexp: /.+?/});
  defineParameterType({name: 'dataset', regexp: /[\w./:-]+?/});
  defineParameterType({name: 'state', regexp: /visible|hidden/});
  kind('icon', {selector: '[name^="icon-"]', match: ['dart'], dartNames: ['icon-{q}']});
  kind('viewer', {selector: '[name^="viewer-"]', match: ['dart'], dartNames: ['viewer-{q}']});
  element('toolbox', {selector: '.d4-toolbox'});
  dataset('spgi', {path: 'System:AppData/Chem/tests/spgi-100.csv'});
  fns = {
    openDataset: Given('user opens {dataset} dataset', async () => undefined),
    clickOn: When('user clicks (on ){element}', async () => undefined),
    shouldBe: Then('{element} should be {state}', async () => undefined),
  };
});

function compile(source = FEATURE, registryModules: string[] = []) {
  const feature = parseFeature(`${ROOT}/features/viewers/toolbox.feature`, source);
  return compileFeature(feature, {root: ROOT, matcher: new StepMatcher(), bindings: bindings(fns, registryModules)});
}

test('one test per Gherkin scenario and outline row on the feature page, background inlined', () => {
  const {code, diagnostics, outFile} = compile();
  assert.equal(diagnostics.filter((d) => d.level === 'error').length, 0);
  assert.ok(outFile.replace(/\\/g, '/').endsWith('bdd/generated/viewers/toolbox.test.ts'));
  assert.match(code, /test\.describe\("Toolbox", \(\) => \{\n  const session = feature\(test\);/);
  assert.equal((code.match(/^  test\(/gm) ?? []).length, 3);
  assert.match(code, /test\("Several viewers \[viewer=bar chart\]", \{tag: \["@demo", "@realizes:u2.dialog"\]\}, async \(\{browser\}\) => \{\n    const page = await session\.page\(browser\);/);
  assert.equal((code.match(/openDataset\(page, ds\("spgi"\)\)/g) ?? []).length, 3);
  assert.match(code, /clickOn\(page, el\("scatter plot icon on toolbox"\)\)/);
  assert.match(code, /shouldBe\(page, el\("bar chart viewer"\), "visible"\)/);
  assert.match(code, /sub_features_covered: \[u2.dialog\]/);
  assert.match(code, /import \{clickOn, openDataset, shouldBe\} from '@datagrok-libraries\/bdd\/bindings\/common\/steps';/);
  assert.match(code, /import \{ds, el, feature\} from '@datagrok-libraries\/bdd\/runtime';/);
});

test('registry modules are side-effect imports: library ones by package subpath, project ones relative', () => {
  const {code} = compile(FEATURE, ['@datagrok-libraries/bdd/bindings/common/kinds', `${ROOT}/bindings/elements.ts`]);
  const lines = code.split('\n').filter((l) => l.startsWith("import '"));
  assert.deepEqual(lines, ["import '../../bindings/elements.js';", "import '@datagrok-libraries/bdd/bindings/common/kinds';"]);
});

test('the output is deterministic', () => {
  assert.equal(compile().code, compile().code);
});

test('unbound steps, unknown elements and unknown datasets are errors with pointers', () => {
  const {code, diagnostics} = compile(`Feature: Broken
  Scenario: Oops
    Given user opens nope dataset
    When user frobs the grid
    Then the frobnicator should be visible
`);
  const errors = diagnostics.filter((d) => d.level === 'error');
  assert.equal(errors.length, 3);
  assert.match(errors[0].message, /dataset "nope" is not registered/);
  assert.match(errors[1].message, /no step definition matches "user frobs the grid" — nearest: .*"user clicks \(on \)\{element\}"/);
  assert.match(errors[2].message, /element "the frobnicator"/);
  assert.equal(errors[1].line, 4);
  assert.match(code, /throw new Error\('no step definition matches this step'\)/);
});

test('ambiguous definitions are reported, not picked', () => {
  fns.clickOn2 = When('user clicks on {element}', async () => undefined);
  const {diagnostics} = compile(`Feature: A
  Scenario: B
    When user clicks on toolbox
`);
  assert.match(diagnostics.find((d) => d.level === 'error')?.message ?? '', /ambiguous step/);
});

test('a literal definition beats a parameterized one', () => {
  fns.openToolbox = When('user clicks on toolbox', async () => undefined);
  const {code, diagnostics} = compile(`Feature: A
  Scenario: B
    When user clicks on toolbox
`);
  assert.equal(diagnostics.filter((d) => d.level === 'error').length, 0);
  assert.match(code, /openToolbox\(page\)/);
});

test('data tables and doc strings are passed as trailing arguments', () => {
  fns.fillIn = When('user fills in:', async () => undefined);
  const {code} = compile(`Feature: A
  Scenario: B
    When user fills in:
      | name input | Batch 7 |
      | method     | clustal |
`);
  assert.match(code, /fillIn\(page, \[\["name input","Batch 7"\],\["method","clustal"\]\]\)/);
});

test('a step that enters a context switches the vocabulary and emits the switch', () => {
  const workbench = context('MSA workbench', {selector: '[data-u2-name="msaWorkbench"]'});
  workbench.element('results', {selector: '[data-u2-name="results"]'});
  fns.openWorkbench = Given('user opens the MSA workbench', async () => undefined, {enters: 'MSA workbench'});
  const {code, diagnostics} = compile(`Feature: A
  Background:
    Given user opens the MSA workbench
  Scenario: B
    Then results should be visible
  Scenario: C
    Then toolbox should be visible
`);
  assert.equal(diagnostics.filter((d) => d.level === 'error').length, 0);
  assert.match(code, /openWorkbench\(page\)\);\n    enter\(page, "MSA workbench"\);/);
  assert.match(code, /import \{el, enter, feature\} from/);
  assert.equal((code.match(/enter\(page/g) ?? []).length, 2);
});

test('phrases before a context is entered do not see its names', () => {
  const workbench = context('MSA workbench', {selector: '.w'});
  workbench.element('results', {selector: '.r'});
  fns.openWorkbench = Given('user opens the MSA workbench', async () => undefined, {enters: 'MSA workbench'});
  const {diagnostics} = compile(`Feature: A
  Scenario: B
    Then results should be visible
    Given user opens the MSA workbench
    Then results should be visible
`);
  const errors = diagnostics.filter((d) => d.level === 'error');
  assert.equal(errors.length, 1);
  assert.equal(errors[0].line, 3);
});

test('entering an unregistered context is an error', () => {
  fns.openThing = Given('user opens the thing', async () => undefined, {enters: 'the thing'});
  const {diagnostics} = compile(`Feature: A
  Scenario: B
    Given user opens the thing
`);
  assert.match(diagnostics.find((d) => d.level === 'error')?.message ?? '', /not a registered context/);
});

test('"of" names a part of a generic kind', () => {
  kind('input', {selector: '[data-u2$="-input"]', match: ['name'], parts: {label: '[data-u2-part="label"]'}});
  const {diagnostics} = compile(`Feature: A
  Scenario: B
    Then label of name input should be visible
`);
  assert.equal(diagnostics.filter((d) => d.level === 'error').length, 0);
  assert.match(diagnostics[0]?.message ?? '', /part "label" \[\[data-u2-part="label"\]\] within kind "input" qualified "name"/);
});
