/* What guide mode records of the pointer, against a static page in the library's Chromium: every
   press where the browser delivered it — a locator's own click as much as the page's mouse — with
   its button and the stop of the step it belongs to. No stand. */
import assert from 'node:assert/strict';
import {mkdtempSync, readFileSync, rmSync} from 'node:fs';
import {tmpdir} from 'node:os';
import {join} from 'node:path';
import {after, before, test} from 'node:test';
import {type Browser, chromium, type Page, type TestInfo} from '@playwright/test';
import {attach, begin, end, hop, type GuideManifest} from '../src/runtime/guide.js';

const PAGE = `<body style="margin:0">
  <button name="button-A" style="position:absolute;left:40px;top:40px;width:120px;height:30px">A</button>
  <button name="button-B" style="position:absolute;left:40px;top:120px;width:120px;height:30px">B</button>
  <div name="div-C" style="position:absolute;left:300px;top:300px;width:200px;height:40px;background:#ddd">C</div>
</body>`;

let browser: Browser | undefined;
let page: Page | undefined;
let missing = '';
const dir = mkdtempSync(join(tmpdir(), 'bdd-guide-capture-'));

before(async () => {
  process.env.BDD_GUIDE = dir;
  process.env.BDD_GUIDE_SETTLE = '0';
  try {
    browser = await chromium.launch();
    page = await browser.newPage({viewport: {width: 800, height: 600}});
    await attach(page);
    await page.goto(`data:text/html,${encodeURIComponent(PAGE)}`);
  }
  catch (e) {
    missing = `no Chromium here: ${(e as Error).message.split('\n')[0]}`;
  }
});

after(async () => {
  await browser?.close();
  delete process.env.BDD_GUIDE;
  delete process.env.BDD_GUIDE_SETTLE;
  rmSync(dir, {recursive: true, force: true});
});

const info = {titlePath: ['capture.test.ts', 'Capture', 'Every press'], tags: []} as unknown as TestInfo;

function manifest(): GuideManifest {
  return JSON.parse(readFileSync(join(dir, 'capture', 'every-press', 'steps.json'), 'utf8'));
}

test('every press is recorded where it landed, with its button, on the stop it was made at', async (t) => {
  if (!page) {
    t.skip(missing);
    return;
  }
  await begin(page, info, 1, 'When user clicks on everything');
  await page.locator('[name="button-A"]').click();
  // a click of its own a moment later, a few pixels off, is no double-click
  await page.locator('[name="button-A"]').click({position: {x: 63, y: 15}});
  await page.mouse.click(600, 200, {button: 'right'});
  await page.locator('[name="button-B"]').dblclick();
  await hop(page, page.locator('[name="div-C"]'));
  await page.locator('[name="div-C"]').click({position: {x: 5, y: 5}});
  await end(page);
  const [step] = manifest().steps;
  assert.equal(step.legs.length, 2, 'the page the step began with, then the stop');
  const presses = step.legs.map((leg) => leg.pointer.filter((p) => p.type === 'down')
    .map((p) => `${p.button} ${p.count} ${p.x},${p.y} ${p.on}`));
  assert.deepEqual(presses[0], ['left 1 100,55 button-A', 'left 1 105,57 button-A', 'right 1 600,200 body',
    'left 1 100,135 button-B', 'left 2 100,135 button-B']);
  assert.deepEqual(presses[1], ['left 1 305,305 div-C']);
  for (const leg of step.legs) {
    for (let i = 1; i < leg.pointer.length; i++)
      assert.ok(!(leg.pointer[i].type === 'move' && leg.pointer[i - 1].type === 'move'), 'a pointer at rest is one move');
  }
  assert.deepEqual(step.target, {x: 300, y: 300, width: 200, height: 40});
});

test('a drag keeps the points the button was held over, and pictures the page along it', async (t) => {
  if (!page) {
    t.skip(missing);
    return;
  }
  await begin(page, info, 2, 'When user drags across the page');
  await page.mouse.move(400, 450);
  await page.mouse.down();
  await page.mouse.move(640, 540);
  await page.mouse.up();
  await end(page);
  const step = manifest().steps[1];
  const pointer = step.legs[0].pointer;
  assert.equal(pointer.find((p) => p.type === 'down')?.x, 400);
  const held = pointer.filter((p) => p.type === 'move' && p.held);
  assert.ok(held.length >= 2, 'the way the button was held over');
  assert.deepEqual([held[held.length - 1].x, held[held.length - 1].y], [640, 540]);
  assert.ok(held.some((p) => p.shot?.endsWith('-drag1.png')), 'a picture of what the drag drew');
  assert.equal(pointer[pointer.length - 1].type, 'up');
});
