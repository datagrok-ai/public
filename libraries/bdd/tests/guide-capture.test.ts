/* What guide mode records of the pointer, against a static page in the library's Chromium: every
   press where the browser delivered it — a locator's own click as much as the page's mouse — with
   its button and the stop of the step it belongs to. No stand. */
import assert from 'node:assert/strict';
import {mkdtempSync, readFileSync, rmSync} from 'node:fs';
import {tmpdir} from 'node:os';
import {join} from 'node:path';
import {after, before, test} from 'node:test';
import {type Browser, chromium, type Page, type TestInfo} from '@playwright/test';
import {attach, begin, end, hop, keyPressed, moveUnfilmed, pressShown, say, waypoint, type GuideManifest}
  from '../src/runtime/guide.js';

const PAGE = `<body style="margin:0">
  <button name="button-A" style="position:absolute;left:40px;top:40px;width:120px;height:30px">A</button>
  <button name="button-B" style="position:absolute;left:40px;top:120px;width:120px;height:30px">B</button>
  <div name="div-C" style="position:absolute;left:300px;top:300px;width:200px;height:40px;background:#ddd">C</div>
  <input name="input-D" style="position:absolute;left:40px;top:200px;width:200px">
  <div class="d4-accordion-pane-header" name="div-section--Explore"
    style="position:absolute;left:40px;top:500px;width:600px;height:24px;font:14px sans-serif">Explore</div>
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

const info = {titlePath: ['capture.test.ts', 'Capture', 'Every press'], tags: [], project: {use: {}}} as unknown as TestInfo;

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

test('typing through a locator is pictured as the text goes in', async (t) => {
  if (!page) {
    t.skip(missing);
    return;
  }
  await begin(page, info, 3, 'When user types into D');
  const input = page.locator('[name="input-D"]');
  await input.fill('twelve chars');
  await input.pressSequentially(' more');
  await end(page);
  const stills = manifest().steps[2].legs.flatMap((l) => l.pointer).filter((p) => p.type === 'still');
  assert.ok(stills.length >= 7, `a picture after each piece of the text, ${stills.length} taken`);
  assert.ok(stills.every((p) => p.shot?.startsWith('03-beat')));
  assert.equal(await input.inputValue(), 'twelve chars more');
});

test('a corner of the way is kept, a press only the video shows is on its box, held and pressed keys are named',
  async (t) => {
  if (!page) {
    t.skip(missing);
    return;
  }
  await begin(page, info, 4, 'When user walks into a menu');
  await hop(page, page.locator('[name="div-C"]'), 'menu');
  await page.mouse.move(100, 300);
  await waypoint(page);
  await page.mouse.move(150, 300);
  await page.mouse.move(150, 400);
  await pressShown(page, {x: 300, y: 300, width: 200, height: 40}, 'option');
  await page.keyboard.down('Shift');
  await page.locator('[name="button-B"]').click();
  await page.keyboard.up('Shift');
  keyPressed(page, 'Escape');
  await end(page);
  const step = manifest().steps[3];
  assert.equal(step.legs[0].kind, 'menu', 'the stop is paced as a menu\'s');
  const moves = step.legs[0].pointer.filter((p) => p.type === 'move')
    .map((p) => `${p.x},${p.y}${p.kept ? ' kept' : ''}`);
  assert.deepEqual(moves.slice(0, 3), ['100,300 kept', '150,400', '400,320'],
    'the corner, where the pointer came to rest, the shown press');
  const presses = step.legs[0].pointer.filter((p) => p.type === 'down');
  assert.deepEqual(presses.map((p) => `${p.on} ${p.shown ? 'shown' : 'real'} ${(p.mods ?? []).join('+')}`),
    ['option shown ', 'button-B real Shift']);
  assert.deepEqual(presses[0].box, {x: 300, y: 300, width: 200, height: 40});
  assert.deepEqual(step.pressed, ['Escape']);
});

test('a dialog the step opened is counted by its inputs; one open before the step is not the step\'s', async (t) => {
  if (!page) {
    t.skip(missing);
    return;
  }
  const dialog = (inputs: number): Promise<void> => page!.evaluate((n) => {
    const d = document.createElement('div');
    d.className = 'd4-dialog';
    d.innerHTML = Array.from({length: n}, (_, i) => `<div class="ui-input-root"><input name="i${i}"></div>`).join('');
    document.body.appendChild(d);
  }, inputs);
  await begin(page, info, 5, 'When user opens a dialog');
  await dialog(3);
  await end(page);
  await begin(page, info, 6, 'When user does nothing with the dialog open');
  await end(page);
  const [opened, kept] = manifest().steps.slice(4);
  assert.deepEqual(opened.dialog, {inputs: 3});
  assert.equal(kept.dialog, undefined);
  await page.evaluate(() => document.querySelectorAll('.d4-dialog').forEach((d) => d.remove()));
});

test('a row as wide as its panel is pressed on its caption, a move left unfilmed is not in the video, a note says where it is pinned',
  async (t) => {
    if (!page) {
      t.skip(missing);
      return;
    }
    await begin(page, info, 7, 'When user expands "Explore" section');
    await page.mouse.move(700, 300);
    await moveUnfilmed(page, 600, 100);
    await page.locator('[name="div-section--Explore"]').click();
    await end(page);
    await begin(page, info, 8, 'When the video says "Here" at the top');
    say(page, 'Here', 'top');
    await end(page);
    const [step, note] = manifest().steps.slice(6);
    const pointer = step.legs.flatMap((l) => l.pointer);
    assert.ok(!pointer.some((p) => p.x === 600 && p.y === 100), 'the unfilmed move is not in the video');
    const press = pointer.find((p) => p.type === 'down')!;
    assert.ok(press.x < 120, `the header pressed on its caption, ${press.x} px from the page's left (its middle is 340)`);
    assert.equal(press.y, 512);
    assert.equal(note.kind, 'note');
    assert.equal(note.noteAt, 'top');
  });

test('a demo takes a hover\'s tooltip down before set-up that follows a note, not before set-up right after the hover',
  async (t) => {
    if (!page) {
      t.skip(missing);
      return;
    }
    process.env.BDD_GUIDE_DEMO = '1';
    try {
      // a string, not a function: tsx would wrap the stub's arrow in a `__name` call the page lacks
      await page.evaluate('window.__bddHides = 0; window.ui = {tooltip: {hide() { window.__bddHides++; }}};');
      const hides = (): Promise<number> => page!.evaluate(() => (window as any).__bddHides as number);
      await begin(page, info, 9, 'When user hovers over C');
      await page.mouse.move(400, 320);
      await end(page);
      await begin(page, info, 10, 'Given the rows of the grid are 90 px tall');
      await end(page);
      assert.equal(await hides(), 0, 'set-up right after a hover is part of what the hover led to');
      await begin(page, info, 11, 'When the video says "The rows"');
      say(page, 'The rows');
      await end(page);
      await begin(page, info, 12, 'Given the rows of the grid are 60 px tall');
      await end(page);
      assert.equal(await hides(), 1, 'set-up after a note starts with the tooltip down');
    }
    finally {
      delete process.env.BDD_GUIDE_DEMO;
    }
  });
