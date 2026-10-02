/* The legend phrases of the viewers tier against the real bindings: each resolves to its own step,
   none of them takes over a phrase an existing feature already uses, and the runtime refuses a
   splitter direction it cannot drag. */
import assert from 'node:assert/strict';
import {join, resolve} from 'node:path';
import {before, test} from 'node:test';
import {compileFeature} from '../src/compile.js';
import {type Bindings, loadBindings} from '../src/discover.js';
import {parseFeature} from '../src/gherkin.js';
import {StepMatcher} from '../src/match.js';
import {librarySource} from '../src/project.js';
import {dragLegendSplitter} from '../src/runtime/viewer-legend.js';

const ROOT = resolve('/pkg/bdd').split('\\').join('/');
let bindings: Bindings;

before(async () => {
  bindings = await loadBindings([librarySource('common'), librarySource('platform'), librarySource(join('tiers', 'viewers'))]);
});

function compile(steps: string[]) {
  const source = `Feature: Legend\n  Scenario: Phrases\n${steps.map((s) => `    ${s}`).join('\n')}\n`;
  const feature = parseFeature(`${ROOT}/features/legend.feature`, source);
  return compileFeature(feature, {root: ROOT, matcher: new StepMatcher(), bindings});
}

test('the legend mode, size and item phrases resolve to their own steps', () => {
  const {code, diagnostics} = compile([
    'Then the legend of scatter plot viewer should be in a corner',
    'And the legend of scatter plot viewer should be collapsed to the mini icon',
    'And the legend of scatter plot viewer should be placed nowhere',
    'When user drags the legend splitter of bar chart viewer by 60 pixels to the left',
    'Then the legend of bar chart viewer should be wider than before',
    'And the legend of bar chart viewer should be narrower than before',
    'And the legend of bar chart viewer should be taller than before',
    'And the legend of bar chart viewer should be shorter than before',
    'And every item in the legend of pie chart viewer should be drawn as a structure',
    'And every item in the legend of pie chart viewer should be drawn as text',
  ]);
  assert.deepEqual(diagnostics.filter((d) => d.level === 'error'), []);
  for (const call of ['legendInCorner(page, el("scatter plot viewer"))', 'legendMiniIcon(page, el("scatter plot viewer"))',
    'legendPlacedNowhere(page, el("scatter plot viewer"))', 'dragLegendSplitter(page, el("bar chart viewer"), 60, "left")',
    'legendWider(page, el("bar chart viewer"))', 'legendNarrower(page, el("bar chart viewer"))', 'legendTaller(page, el("bar chart viewer"))',
    'legendShorter(page, el("bar chart viewer"))', 'legendItemsAsStructures(page, el("pie chart viewer"))', 'legendItemsAsText(page, el("pie chart viewer"))'])
    assert.ok(code.includes(call), `the generated spec calls ${call}`);
});

test('the phrases existing features use keep their steps', () => {
  const {code, diagnostics} = compile([
    'Then legend of scatter plot viewer should be hidden',
    'And the legend of line chart viewer should be hidden',
    'And the legend of scatter plot viewer should be docked',
    'And the legend of trellis plot viewer should be in the "left" slot',
    'When user drags the "splitter" area of grid by 40 pixels to the right',
  ]);
  assert.deepEqual(diagnostics.filter((d) => d.level === 'error'), []);
  assert.ok(code.includes('shouldBe(page, el("legend of scatter plot viewer"), "hidden")'));
  assert.ok(code.includes('shouldBe(page, el("the legend of line chart viewer"), "hidden")'));
  assert.ok(code.includes('legendDocked(page, el("scatter plot viewer"))'));
  assert.ok(code.includes('legendSlot(page, el("trellis plot viewer"), "left")'));
  assert.ok(code.includes('dragAreaBy(page, "splitter", el("grid"), 40, "right")'));
  assert.ok(!code.includes('dragLegendSplitter'));
});

test('a legend splitter drags left, right, up or down, and nowhere else', async () => {
  await assert.rejects(dragLegendSplitter({} as never, {phrase: 'bar chart viewer'} as never, 40, 'sideways'),
    /left, right, up or down, not "sideways"/);
});
