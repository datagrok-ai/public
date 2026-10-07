/* Steps against a static page with a stand-in `grok`/`DG`: what each claim reads and when it
   fails. Runs in the library's Chromium; no stand. */
import assert from 'node:assert/strict';
import {after, before, test} from 'node:test';
import {type Browser, chromium, type Page} from '@playwright/test';

import '../bindings/common/kinds.js';
import '../bindings/platform/elements.js';
import {fillsParent, shouldOffer, visibleCount} from '../bindings/common/steps.js';
import {mouseOverRowIs} from '../bindings/platform/columns.js';
import {newestMatchingDistinct, newestMatchingFilled} from '../bindings/platform/commands.js';
import {setTableTag, tableTagIsFile} from '../bindings/platform/data.js';
import {taskBarFinished, taskBarShown, watchTaskBar} from '../bindings/platform/events.js';
import {openTableOf, SKETCHER_CONTROLS_TAG, sketcherPin} from '../bindings/platform/steps.js';
import {fixtureFamilies, isStaleFixture, serviceGap} from '../src/runtime/server.js';
import {el} from '../src/runtime/args.js';
import {select} from '../src/runtime/gestures.js';
import {hitArea, readValue} from '../src/runtime/viewers.js';
import {locate} from '../src/runtime/locate.js';
import {knownFailure} from '../src/runtime/harness.js';
import {whileExpectedToFail} from '../src/runtime/patience.js';

let missing = '';
let browser: Browser | undefined;
let page: Page | undefined;

const scenario = (name: string, fn: () => Promise<void>): void => void test(name, async (t) => {
  if (page === undefined) {
    t.skip(missing || 'no page');
    return;
  }
  await fn();
});

/** A step that must fail, failing fast: the claims poll, and a stand-in page never changes. */
const fails = (body: () => Promise<unknown>): Promise<void> => assert.rejects(() => whileExpectedToFail(async () => { await body(); }));

before(async () => {
  try {
    browser = await chromium.launch();
  }
  catch (e) {
    missing = `Chromium does not start here: ${String(e).slice(0, 80)}`;
    return;
  }
  page = await browser.newPage();
  // tsx names the functions it compiles with a helper the page does not have
  await page.addInitScript(() => { (window as any).__name = (fn: unknown) => fn; });
  await page.goto('about:blank');
});

after(async () => {
  await browser?.close();
});

/** A current table of named columns (each a list of values, null for missing) with tags. */
async function standIn(p: Page, columns: Record<string, (string | null)[]>): Promise<void> {
  await p.evaluate((cols) => {
    const names = Object.keys(cols);
    const tags: Record<string, string> = {};
    const col = (n: string) => ({name: n, get: (i: number) => cols[n][i], isNone: (i: number) => cols[n][i] === null});
    const t = {
      rowCount: cols[names[0]].length,
      columns: {names: () => names},
      col: (n: string) => names.includes(n) ? col(n) : null,
      getTag: (k: string) => tags[k] ?? null,
      setTag: (k: string, v: string) => { tags[k] = v; },
    };
    // the event streams the viewer runtime subscribes to when a step changes the table through it
    const stream = {subscribe: () => ({unsubscribe: () => undefined})};
    (window as any).grok = {shell: {t, tables: [t], tableViews: []},
      events: {onViewerAdded: stream, onViewerClosed: stream, onEvent: () => stream, onCustomEvent: () => stream},
      functions: {onBeforeRunAction: stream, onAfterRunAction: stream}};
  }, columns);
}

scenario('the choices a dropdown offers are read in order, exactly', async () => {
  await page!.setContent(`<div class="ui-input-root" name="input-host-Linkage"><select name="input-Linkage"><option>single</option><option>ward</option></select></div>`);
  await shouldOffer(page!, el('Linkage input'), 'single, ward');
  await fails(() => shouldOffer(page!, el('Linkage input'), 'ward, single'));
  await fails(() => shouldOffer(page!, el('Linkage input'), 'single'));
});

scenario('visible elements are counted, hidden ones are not', async () => {
  await page!.setContent(`<i class="grok-icon" aria-label="Assign Clusters">A</i><i class="grok-icon" aria-label="Assign Clusters" style="display:none">A</i>`);
  await visibleCount(page!, 1, el('"Assign Clusters" icon'));
  await fails(() => visibleCount(page!, 2, el('"Assign Clusters" icon')));
});

scenario('the newest column matching a pattern: distinct values and missing values', async () => {
  await page!.setContent('<div></div>');
  await standIn(page!, {'Cluster (1.00)': ['1', '1', '1'], 'Cluster (2.00)': ['1', '2', null]});
  await newestMatchingDistinct(page!, '^Cluster \\(', 2);
  await fails(() => newestMatchingDistinct(page!, '^Cluster \\(', 3));
  await fails(() => newestMatchingFilled(page!, '^Cluster \\('));
  await newestMatchingFilled(page!, '1\\.00');
  await assert.rejects(() => newestMatchingDistinct(page!, '^Missing', 1), /no column matching/);
});

scenario('a table tag is set on the current table', async () => {
  await page!.setContent('<div></div>');
  await standIn(page!, {node: ['a']});
  await setTableTag(page!, '.newick-alt', '(a,b);');
  assert.equal(await page!.evaluate(() => (window as any).grok.shell.t.getTag('.newick-alt')), '(a,b);');
});

scenario('a table written in the feature is parsed as a CSV into a view of its own', async () => {
  await page!.setContent('<div></div>');
  await page!.evaluate(() => {
    const w = window as any;
    w.DG = {DataFrame: {fromCsv: (csv: string) => ({csv, name: ''})}};
    w.grok = {shell: {addTableView: (df: any) => { w.grok.shell.tv = {dataFrame: df, grid: {}}; }}};
  });
  await openTableOf(page!, 'leaves', [['leaf', 'note'], ['A', 'x, y'], ['B', 'say "hi"']]);
  assert.equal(await page!.evaluate(() => (window as any).grok.shell.tv.dataFrame.csv), 'leaf,note\nA,"x, y"\nB,"say ""hi"""');
  await assert.rejects(() => openTableOf(page!, 'empty', [['leaf']]), /at least one row/);
});

scenario('a Dart property row opens its editor on the value before a choice is made, and names its label and value', async () => {
  await page!.setContent(`<table><tr class="property-grid-item" name="prop-newick-tag"><td class="property-grid-item-name"><div class="property-grid-item-name-text"><span>Newick Tag</span></div></td>
    <td class="property-grid-item-value"><div class="property-grid-item-view-label" name="prop-view-newick-tag"></div></td></tr></table>`);
  await page!.evaluate(() => {
    const cell = document.querySelector('.property-grid-item-value')!;
    cell.addEventListener('click', () => {
      if (!cell.querySelector('select'))
        cell.insertAdjacentHTML('beforeend', '<select><option value=""></option><option value=".newick-alt">.newick-alt</option></select>');
    });
  });
  await select(page!, el('"Newick Tag" property'), '.newick-alt');
  assert.equal(await page!.locator('select').inputValue(), '.newick-alt');
  assert.equal(await (await locate(page!, el('value of "Newick Tag" property'))).count(), 1);
  assert.equal((await (await locate(page!, el('label of "Newick Tag" property'))).textContent())?.trim(), 'Newick Tag');
});

scenario('the task bar record keeps an entry the bar has already removed', async () => {
  await page!.setContent('<div class="d4-task-bar"></div>');
  await watchTaskBar(page!);
  await page!.evaluate(async () => {
    const bar = document.querySelector('.d4-task-bar')!;
    bar.insertAdjacentHTML('beforeend', '<div class="d4-task-bar-entry">Creating dendrogram ...</div>');
    await new Promise((r) => setTimeout(r, 0));
    bar.innerHTML = '';
  });
  await taskBarShown(page!, 'Creating dendrogram');
  await fails(() => taskBarShown(page!, 'Loading'));
  await taskBarFinished(page!, 'Creating dendrogram');
  await fails(() => taskBarFinished(page!, 'Loading'));
  await page!.evaluate(() => { document.querySelector('.d4-task-bar')!.insertAdjacentHTML('beforeend', '<div class="d4-task-bar-entry">Loading ...</div>'); });
  await fails(() => taskBarFinished(page!, 'Loading'));
});

scenario('the mouse-over row counts from 1, and a filled element is told from one inside its padding', async () => {
  await page!.setContent(`<div id="host" style="width:200px;height:100px;padding:10px"><div name="viewer-Full" style="width:100%;height:100%">f</div></div>
    <div id="host2" style="width:200px;height:100px"><div name="viewer-Half" style="width:50%;height:100%">h</div></div>`);
  await page!.evaluate(() => { (window as any).grok = {shell: {t: {mouseOverRowIdx: 4}}}; });
  await mouseOverRowIs(page!, 5);
  await fails(() => mouseOverRowIs(page!, 4));
  await fillsParent(page!, el('Full viewer'));
  await fails(() => fillsParent(page!, el('Half viewer')));
});

scenario('a widget outside any viewer that reports a status is found up from its element: its hit areas and readings', async () => {
  await page!.setContent(`<div style="padding:30px"><div data-widget="true" name="Probe" style="position:relative;left:5px;width:200px;height:100px">
    <div id="inner" style="width:50px;height:50px"></div></div><div data-widget="true" name="Plain" style="width:20px;height:20px"></div></div>`);
  await page!.evaluate(() => {
    const w = window as any;
    const root = document.querySelector('[name="Probe"]')!;
    const probe = {root, type: 'Probe', isRenderPending: false,
      getWidgetStatus: () => ({parts: {}, hitAreas: {'atom 0': {x: 10, y: 20, width: 8, height: 8}}, values: {smiles: 'CC'}})};
    w.grok = {shell: {tableViews: []}};
    w.DG = {Widget: {find: (e: Element) => e === root ? probe : null}};
  });
  const r = await page!.locator('[name="Probe"]').boundingBox();
  assert.deepEqual(await hitArea(page!, el('Probe widget'), 'atom 0'), {x: r!.x + 10, y: r!.y + 20, width: 8, height: 8});
  assert.equal(await readValue(page!, el('Probe widget'), 'smiles'), 'CC');
  await assert.rejects(() => readValue(page!, el('Plain widget'), 'smiles'), /not a viewer/);
});

scenario('a table tag is compared with a file byte for byte', async () => {
  await page!.setContent('<div></div>');
  await standIn(page!, {node: ['a']});
  await page!.evaluate(() => {
    const w = window as any;
    w.grok.dapi = {files: {readAsText: async (p: string) => p === 'System:A/tree.nwk' ? '(a,b);\r\n' : ''}};
    w.grok.shell.t.setTag('.newick', '(a,b);\r\n');
  });
  await tableTagIsFile(page!, '.newick', 'System:A/tree.nwk');
  await page!.evaluate(() => (window as any).grok.shell.t.setTag('.newick', '(a,b);'));
  await assert.rejects(() => tableTagIsFile(page!, '.newick', 'System:A/tree.nwk'));
});

test('a known failure passes when its steps fail and fails when they pass', async () => {
  await knownFailure(async () => { throw new Error('the defect'); });
  await assert.rejects(() => knownFailure(async () => undefined), /tagged @known-failure and passed/);
});

test('a fixture of a dead run is stale after an hour; a live, a foreign or an unsuffixed one is not', () => {
  const families = fixtureFamilies(['BDD-GL-Group-1789927107045', 'BDD-Share-Model-2f1c9c3e-0b1a-4c2d-8e9f-0123456789ab', 'BDD-CP-Root']);
  assert.deepEqual(families, ['BDD-GL-Group', 'BDD-Share-Model']);
  const now = 1789930000000;
  const old = now - 2 * 3600 * 1000;
  const at = (name: string, createdOn: number) => ({name, friendlyName: name, createdOn});
  assert.equal(isStaleFixture(at('BDD-GL-Group-1789900000000', old), families, now), true);
  assert.equal(isStaleFixture(at('BDD-Share-Model-9f8e7d6c-5b4a-4321-8765-ba9876543210', old), families, now), true);
  assert.equal(isStaleFixture(at('BDD-GL-Group-1789900000000', now - 5 * 60 * 1000), families, now), false);
  assert.equal(isStaleFixture(at('BDD-GL-Renamed-1789900000000', old), families, now), false);
  assert.equal(isStaleFixture(at('BDD-GL-Group', old), families, now), false);
  assert.equal(isStaleFixture(at('BDD-GL-Group-1789900000000', 0), families, now), false);
});

test('a service gate skips on a service the stand reports missing or down, and lets a stand that reports none go on', () => {
  const services = [{name: 'Jupyter', enabled: true, status: 'Running'}, {name: 'Grok Spawner', enabled: false, status: 'Running'},
    {name: 'Grok Connect', enabled: true, status: 'Failed'}];
  assert.equal(serviceGap(services, 'Jupyter'), '');
  assert.equal(serviceGap(services, 'Grok Spawner'), 'disabled, Running');
  assert.equal(serviceGap(services, 'Grok Connect'), 'Failed');
  assert.equal(serviceGap(services, 'No Such Service'), 'absent');
  assert.equal(serviceGap([], 'Jupyter'), '');
  assert.equal(serviceGap([], 'Jupyter', true), 'the stand reports no service health');
});

test('the molecule sketcher pinned is the one the feature names, or the run override; a feature about the controls of one sketcher skips under another', () => {
  assert.deepEqual(sketcherPin('OpenChemLib', undefined, []), {pin: 'OpenChemLib'});
  assert.deepEqual(sketcherPin('OpenChemLib', '', []), {pin: 'OpenChemLib'});
  assert.deepEqual(sketcherPin('OpenChemLib', ' ', [SKETCHER_CONTROLS_TAG]), {pin: 'OpenChemLib'});
  assert.deepEqual(sketcherPin('OpenChemLib', 'Crux', ['@journey']), {pin: 'Crux'});
  assert.deepEqual(sketcherPin('Crux', 'Crux', [SKETCHER_CONTROLS_TAG]), {pin: 'Crux'});
  const skipped = sketcherPin('Ketcher', ' Crux ', ['@journey', SKETCHER_CONTROLS_TAG]);
  assert.ok('skip' in skipped);
  assert.match((skipped as {skip: string}).skip, /Ketcher's own controls .*"Crux" \(BDD_MOLECULE_SKETCHER\)/);
});
