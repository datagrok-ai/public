/* The vocabulary the Queries, Scripts and Connections round promoted out of its package bindings:
   every phrase resolves to a library step, and the ones a package could have kept (the visual query
   builder, the Results table) are deliberately not here. A phrase that stops resolving here is a
   feature of three suites that stops compiling. */
import assert from 'node:assert/strict';
import {join, resolve} from 'node:path';
import {before, test} from 'node:test';
import {compileFeature} from '../src/compile.js';
import {type Bindings, loadBindings} from '../src/discover.js';
import {parseFeature} from '../src/gherkin.js';
import {StepMatcher} from '../src/match.js';
import {librarySource} from '../src/project.js';

const ROOT = resolve('/pkg/bdd').split('\\').join('/');
let bindings: Bindings;

before(async () => {
  bindings = await loadBindings([librarySource('common'), librarySource('platform'), librarySource(join('tiers', 'viewers'))]);
});

function errorsOf(steps: string[]): string[] {
  const source = `Feature: Promoted\n  Scenario: Phrases\n${steps.map((s) => `    ${s}`).join('\n')}\n`;
  const feature = parseFeature(`${ROOT}/features/promoted.feature`, source);
  return compileFeature(feature, {root: ROOT, matcher: new StepMatcher(), bindings})
    .diagnostics.filter((d) => d.level === 'error').map((d) => d.message);
}

test('the entities a Queries, Scripts or Connections feature saves are named on the server', () => {
  assert.deepEqual(errorsOf([
    'Given no query named "BDD-Q-{run}" is on the server',
    'And no script named "BddRun{time}" is on the server',
    'And no connection named "BDD-Conn-{run}" is on the server',
    'And a script "BddRun{time}" is on the server:',
    '  """',
    '  #language: javascript',
    '  """',
    'And a "Postgres" connection named "BDD-Conn-{run}" is on the server',
    'Then 1 query named "BDD-Q-{run}" should be on the server',
    'And 0 scripts named "BddRun{time}" should be on the server',
    'And 1 connection named "BDD-Conn-{run}" should be on the server',
    'And the script "BddRun{time}" on the server should contain "count"',
    'And the script "BddRun{time}" on the server should have an output "count" of type "int"',
    'And the query "BDD-Q-{run}" on the server should have a post-process containing "PP77"',
    'And the query "BDD-Q-{run}" on the server should have a layout',
  ]), []);
});

test('the editor, the console, the panes and the shell phrases of the round resolve', () => {
  assert.deepEqual(errorsOf([
    'Given browser alerts are recorded',
    'And the toolbox pane is hidden',
    'When user replaces the code of code editor with "select 1"',
    'And user appends "r = 1" to code editor',
    'And user notes the console output',
    'And user calls the script "BddRun{time}" from the console with \'"cars"\'',
    'And user opens the "BDD-P-{run}" project and waits for its table',
    'Then code editor should hold the code "select 1"',
    'And the text area of Query pane should hold "select 1"',
    'And the console should show "count: 510"',
    'And the console should show "count: 510" 1 time',
    'And the "Activity" pane of the context panel should count at least 1',
    'And the current view should be a DataQueryView view',
    'And the browser should have shown the alert "Hello World!"',
  ]), []);
});

test('the Scripts view and the gallery kinds the round added resolve', () => {
  assert.deepEqual(errorsOf([
    'Given user opens the Scripts view',
    'When user clicks on "BddRun{time}" gallery card',
    'And user clicks on "Run query..." action in toolbox',
    'Then "BddRun{time}" gallery card should be visible',
  ]), []);
});
