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
import {tableTagIsFile} from '../bindings/platform/data.js';
import {taskBarFinished, taskBarShown, watchTaskBar} from '../bindings/platform/events.js';
import {el} from '../src/runtime/args.js';
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
    (window as any).grok = {shell: {t}};
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
