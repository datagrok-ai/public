/* The layout audit and the geometry claims against static pages: each audit rule produces exactly its
   finding on a broken page and nothing on its twin; every geometry phrase passes, and fails with its
   sentence. Runs in the library's Chromium; no stand. */
import assert from 'node:assert/strict';
import {after, before, test} from 'node:test';
import {type Browser, chromium, type Page} from '@playwright/test';

import '../bindings/common/kinds.js';
import {above, asTallAs, asWideAs, atLeastWidePx, atMostWidePx, below, centeredWith, clearOf, inside, insideWindow, leftAligned, leftOf,
  noLayoutProblems, notClipped, rightOf, tallPx, topAligned, widePx} from '../bindings/common/geometry.js';
import {el} from '../src/runtime/args.js';
import {type AllowEntry, auditLayout, type LayoutResult, rectOf} from '../src/runtime/layout.js';
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

const STYLE = `<style>
  body { margin: 0; font: 14px sans-serif; }
  [data-u2$="-input"] { display: flex; align-items: center; height: 28px; }
  [data-u2-part="label"] { display: block; width: 100px; white-space: nowrap; }
  [data-u2-part="editor"] { display: block; }
  [data-u2-part="editor"] input { display: block; width: 150px; height: 24px; box-sizing: border-box; margin: 0; }
</style>`;

const attr = (name: string, value: string) => value ? ` ${name}="${value}"` : '';
const INPUT = (label: string, opts: {root?: string; label?: string; editor?: string; attrs?: string} = {}) =>
  `<div data-u2="text-input"${attr('style', opts.root ?? '')}${opts.attrs ?? ''}><label data-u2-part="label"${attr('style', opts.label ?? '')}>${label}</label>` +
  `<span data-u2-part="editor"${attr('style', opts.editor ?? '')}><input></span></div>`;
const OVERLAPPING = INPUT('Name') + INPUT('Notes', {root: 'margin-top: -20px', attrs: ' data-u2-name="notes"'});

async function auditOf(html: string, allow: AllowEntry[] = []): Promise<LayoutResult> {
  await page!.setContent(STYLE + html);
  return auditLayout(page!, {allow});
}

const rules = (r: LayoutResult) => r.findings.map((f) => f.rule);

scenario('a page with no u2 control is not audited', async () => {
  const r = await auditOf('<div>plain text</div>');
  assert.deepEqual(r, {u2: false, controls: 0, findings: [], ms: 0, capped: false, busy: 0});
});

scenario('overlap: two controls lying on one another; a badge positioned over its icon is not one', async () => {
  const broken = await auditOf(OVERLAPPING);
  assert.deepEqual(rules(broken), ['overlap']);
  assert.equal(broken.findings[0].detail, '"Name" text input overlaps "Notes" text input by 1280×20 px at (0, 8)');
  assert.equal(broken.findings[0].what, '"Name" text input and "Notes" text input');
  const twin = await auditOf(`<div style="padding: 20px"><div data-u2="icon-button" aria-label="Inbox" style="position: relative; width: 24px; height: 24px">✉` +
    `<span data-u2="badge" style="position: absolute; top: -6px; right: -6px; width: 12px; height: 12px; font-size: 9px">3</span></div></div>`);
  assert.deepEqual(twin.findings, []);
  assert.equal(twin.controls, 2);
});

scenario('a fixed dialog over the page is another layer: nothing under it overlaps it', async () => {
  const r = await auditOf(INPUT('Name') + INPUT('Notes') +
    `<div data-u2="dialog" style="position: fixed; left: 0; top: 0; width: 400px; height: 200px; background: white">${INPUT('Title')}</div>`);
  assert.deepEqual(r.findings, []);
});

scenario('label-misaligned: a label beside its editor 5 px off its middle; aligned, nothing', async () => {
  const broken = await auditOf(INPUT('Dose, mg', {label: 'position: relative; top: 5px'}));
  assert.deepEqual(rules(broken), ['label-misaligned']);
  assert.equal(broken.findings[0].detail, 'label of "Dose, mg" text input sits 5 px below the middle of its editor');
  assert.equal(broken.findings[0].what, 'label of "Dose, mg" text input');
  assert.deepEqual((await auditOf(INPUT('Dose, mg'))).findings, []);
});

scenario('label-misaligned: beside a two-row editor (48 px) a 28 px label on its first line or its middle is aligned; off both, reported', async () => {
  const list = (top: number) => INPUT('Sites', {root: 'height: auto; align-items: flex-start', editor: 'height: 48px',
    label: `height: 28px; position: relative; top: ${top}px`});
  assert.deepEqual((await auditOf(list(0))).findings, []);
  assert.deepEqual((await auditOf(list(10))).findings, []);
  const nearFirst = await auditOf(list(3));
  assert.deepEqual(rules(nearFirst), ['label-misaligned']);
  assert.equal(nearFirst.findings[0].detail, 'label of "Sites" text input sits 3 px below the first line of its editor');
  assert.equal((await auditOf(list(16))).findings[0]?.detail, 'label of "Sites" text input sits 6 px below the middle of its editor');
});

scenario('label-misaligned: a label above its editor starting 6 px off its left edge; aligned, nothing', async () => {
  const top = {root: 'display: block; height: auto'};
  const broken = await auditOf(INPUT('Name', {...top, editor: 'margin-left: 6px'}));
  assert.deepEqual(rules(broken), ['label-misaligned']);
  assert.equal(broken.findings[0].detail, 'label of "Name" text input starts 6 px left of its editor');
  assert.deepEqual((await auditOf(INPUT('Name', top))).findings, []);
});

scenario('label-column: a form whose third editor starts 6 px right of the others', async () => {
  const form = (third: string) => `<div data-u2="form" aria-label="Patient">${INPUT('Name')}${INPUT('Age')}${INPUT('Dose', {editor: third})}</div>`;
  const broken = await auditOf(form('margin-left: 6px'));
  assert.deepEqual(rules(broken), ['label-column']);
  assert.equal(broken.findings[0].detail,
    'the editors of "Patient" form do not start in one column: "Dose" text input at 106 px; the other 2 at 100 px');
  assert.deepEqual((await auditOf(form(''))).findings, []);
});

scenario('label-column: inputs inside another input\'s editor (a columns map\'s key rows) are not rows of the form', async () => {
  const key = {label: 'width: 20px'};
  const map = `<div data-u2="columns-map-input" style="height: auto; align-items: flex-start"><label data-u2-part="label">Mapping</label>` +
    `<span data-u2-part="editor">${INPUT('x', key)}${INPUT('y', key)}</span></div>`;
  const r = await auditOf(`<div data-u2="form" aria-label="Plot">${INPUT('Table')}${map}${INPUT('Color')}</div>`);
  assert.deepEqual(r.findings, []);
});

scenario('clipped-text and truncated-text: a label cut, a label ellipsized; an ellipsized item label is a list\'s design', async () => {
  const cut = 'width: 40px; overflow: hidden';
  const clipped = await auditOf(INPUT('Dose, mg per kilogram', {label: cut}));
  assert.deepEqual(rules(clipped), ['clipped-text']);
  assert.match(clipped.findings[0].detail, /^label of "Dose, mg per kilogram" text input is cut to "Dose[^"…]*" \(needs \d+ px, has 40\)$/);
  const truncated = await auditOf(INPUT('Dose, mg per kilogram', {label: `${cut}; text-overflow: ellipsis`}));
  assert.deepEqual(rules(truncated), ['truncated-text']);
  assert.match(truncated.findings[0].detail, /^label of "Dose, mg per kilogram" text input is cut to "D[^"]*…" \(needs \d+ px, has 40\)$/);
  const item = await auditOf(`<div data-u2="list" aria-label="Fruits"><div data-u2-item="a" role="option">` +
    `<span data-u2-part="item-label" style="display: block; ${cut}; white-space: nowrap; text-overflow: ellipsis">Pineapple and mango</span></div></div>`);
  assert.deepEqual(item.findings, []);
});

scenario('clipped-text: a button cut by the control that clips it; one wholly hidden in a collapsed pane is put away, not cut', async () => {
  const cut = await auditOf('<div data-u2="card" aria-label="Actions" style="width: 60px; overflow: hidden"><button style="width: 100px">Save</button></div>');
  assert.deepEqual(rules(cut), ['clipped-text']);
  assert.match(cut.findings[0].detail, /^"Save" button of "Actions" card is cut off by "Actions" card: 60×\d+ of 100×\d+ px show$/);
  const away = await auditOf('<div data-u2="section" aria-label="More"><div>More</div>' +
    '<div style="height: 0; overflow: hidden"><button>Hidden away</button></div></div>');
  assert.deepEqual(away.findings, []);
});

scenario('clipped-text: a label partly scrolled out of its scroller is not cut; one in the strip a scroller overrunning its panel hides is', async () => {
  const panel = (sizing: string) => `<div style="height: 100px; overflow: hidden"><div style="height: 100%; padding: 8px; box-sizing: ${sizing}; ` +
    `overflow: auto"><div style="height: 80px"></div>${INPUT('Cyl')}<div style="height: 200px"></div></div></div>`;
  assert.deepEqual((await auditOf(panel('border-box'))).findings, []);
  const cut = await auditOf(panel('content-box'));
  assert.deepEqual(rules(cut), ['clipped-text']);
  assert.match(cut.findings[0].detail,
    /^label of "Cyl" text input is cut off by its container: 100×\d+ of 100×\d+ px show; the scroller inside it runs 16 px past its edge$/);
});

scenario('overflow: content spilling out of its control; content that fits, nothing', async () => {
  const card = (width: number) => `<div data-u2="card" aria-label="Summary" style="width: 100px; height: 40px"><div style="width: ${width}px; height: 20px"></div></div>`;
  const broken = await auditOf(card(200));
  assert.deepEqual(rules(broken), ['overflow']);
  assert.equal(broken.findings[0].detail, 'content spills out of "Summary" card: it needs 200×40 px, has 100×40');
  assert.deepEqual((await auditOf(card(80))).findings, []);
});

scenario('zero-size: a control with content and no width; an empty container is not one', async () => {
  const broken = await auditOf('<div data-u2="button" aria-label="Go" style="width: 0; overflow: hidden">Go</div>');
  assert.deepEqual(rules(broken), ['zero-size']);
  assert.match(broken.findings[0].detail, /^"Go" button has no size on the page \(0×\d+ px\)$/);
  assert.deepEqual((await auditOf('<div data-u2="panel"></div>')).findings, []);
});

scenario('zero-size: a form section whose header and fields are all hidden shows nothing; one squeezed to no height over a shown field does', async () => {
  const section = (body: string, fields: string) => `<div data-u2="section" aria-label="Battery" style="width: 406px${body}">` +
    `<div data-u2-part="header" hidden><span data-u2-part="title">Battery</span></div><div data-u2-part="body">${fields}</div></div>`;
  assert.deepEqual((await auditOf(section('', INPUT('Battery Capacity', {root: 'display: none', attrs: ' hidden'})))).findings, []);
  const squeezed = await auditOf(section('; height: 0; overflow: hidden', INPUT('Battery Capacity')));
  assert.deepEqual(rules(squeezed), ['zero-size']);
  assert.equal(squeezed.findings[0].detail, '"Battery" section has no size on the page (406×0 px)');
});

scenario('offscreen: a control outside the window that nothing scrolls to; one inside a scroller is reachable', async () => {
  const broken = await auditOf('<div data-u2="button" aria-label="Away" style="position: relative; left: -500px; width: 100px">Away</div>');
  assert.deepEqual(rules(broken), ['offscreen']);
  assert.equal(broken.findings[0].detail, '"Away" button lies outside the window at (-500, 0) and nothing scrolls to it');
  const twin = await auditOf('<div style="height: 100px; overflow: auto"><div style="height: 2000px"></div>' +
    '<div data-u2="button" aria-label="Far">Far</div></div>');
  assert.deepEqual(twin.findings, []);
  assert.equal(twin.controls, 1);
});

scenario('an allow entry drops the findings of its rule on the elements under its selector, and only those', async () => {
  assert.deepEqual((await auditOf(OVERLAPPING, [{rule: 'overlap', selector: '[data-u2-name="notes"]', reason: 'test'}])).findings, []);
  assert.deepEqual(rules(await auditOf(OVERLAPPING, [{rule: 'overflow', selector: '[data-u2-name="notes"]', reason: 'test'}])), ['overlap']);
});

scenario('a busy control is waited for, then left out with everything in it and counted', async () => {
  const busy = await auditOf(`<div data-u2="form" aria-label="Loading" aria-busy="true">${OVERLAPPING}</div>`);
  assert.deepEqual(busy.findings, []);
  assert.equal(busy.busy, 1);
  assert.equal(busy.controls, 0);
  await page!.setContent(STYLE + `<div data-u2="form" aria-label="Loading" aria-busy="true">${OVERLAPPING}</div>`);
  await page!.evaluate(() => setTimeout(() => document.querySelector('[aria-busy]')!.removeAttribute('aria-busy'), 300));
  const done = await auditLayout(page!);
  assert.deepEqual(rules(done), ['overlap']);
  assert.equal(done.busy, 0);
});

scenario('a page of more than 400 controls is capped at 400 and audited in under 100 ms', async () => {
  const r = await auditOf(`<div data-u2="form" aria-label="Many">${Array.from({length: 420}, (_, i) => INPUT(`Field ${i}`)).join('')}</div>`);
  assert.equal(r.controls, 400);
  assert.equal(r.capped, true);
  assert.deepEqual(r.findings, []);
  assert.ok(r.ms < 100, `audited in ${r.ms} ms`);
});

const BOXES = `<div data-u2-name="alpha" style="position: absolute; left: 10px; top: 10px; width: 100px; height: 30px"></div>
  <div data-u2-name="beta" style="position: absolute; left: 120px; top: 10px; width: 100px; height: 30px"></div>
  <div data-u2-name="gamma" style="position: absolute; left: 10px; top: 50px; width: 210px; height: 100px">
    <div data-u2-name="delta" style="position: absolute; left: 10px; top: 10px; width: 50px; height: 30px"></div></div>
  <div data-u2-name="epsilon" style="position: absolute; left: 60px; top: 20px; width: 100px; height: 30px"></div>
  <div data-u2-name="frame" style="position: absolute; left: 300px; top: 10px; width: 60px; height: 30px; overflow: hidden">
    <div data-u2-name="wide" style="width: 100px; height: 20px"></div></div>
  <div data-u2-name="caption" style="position: absolute; left: 300px; top: 60px; width: 40px; overflow: hidden; white-space: nowrap;
    text-overflow: ellipsis">A long caption text</div>
  <div data-u2-name="twin" style="position: absolute; left: 400px; top: 10px; width: 20px; height: 20px"></div>
  <div data-u2-name="twin" style="position: absolute; left: 430px; top: 10px; width: 20px; height: 20px"></div>`;

const e = (name: string) => el(`${name} element`);

/** Claims that must fail, run side by side under the narrowed budget, each with its sentence. */
async function failWith(checks: [() => Promise<unknown>, RegExp][]): Promise<void> {
  const settled = await whileExpectedToFail(() => Promise.allSettled(checks.map(([check]) => check())));
  for (const [i, s] of settled.entries()) {
    assert.equal(s.status, 'rejected', `check ${i} should have failed`);
    assert.match(String((s as PromiseRejectedResult).reason?.message), checks[i][1]);
  }
}

scenario('the geometry claims hold where the boxes stand so', async () => {
  await page!.setContent(BOXES);
  const p = page!;
  await leftOf(p, e('alpha'), e('beta'));
  await rightOf(p, e('beta'), e('alpha'));
  await above(p, e('alpha'), e('gamma'));
  await below(p, e('gamma'), e('alpha'));
  await inside(p, e('delta'), e('gamma'));
  await clearOf(p, e('alpha'), e('beta'));
  await leftAligned(p, e('alpha'), e('gamma'));
  await topAligned(p, e('alpha'), e('beta'));
  await centeredWith(p, e('alpha'), e('beta'));
  await asTallAs(p, e('alpha'), e('beta'));
  await tallPx(p, e('alpha'), 30);
  await asWideAs(p, e('alpha'), e('beta'));
  await widePx(p, e('alpha'), 100);
  await atLeastWidePx(p, e('alpha'), 80);
  await atMostWidePx(p, e('alpha'), 120);
  await notClipped(p, e('alpha'));
  assert.deepEqual(await rectOf(p, e('alpha')), {x: 10, y: 10, width: 100, height: 30, what: '"alpha"'});
});

scenario('the geometry claims fail with both boxes in one sentence', async () => {
  await page!.setContent(BOXES);
  const p = page!;
  await failWith([
    [() => leftOf(p, e('beta'), e('alpha')),
      /beta element is not left of alpha element — its right edge is 210 px past the left edge of the other: beta element spans 120\.\.220 × 10\.\.40, alpha element 10\.\.110 × 10\.\.40/],
    [() => rightOf(p, e('alpha'), e('beta')), /alpha element is not right of beta element — its left edge is 210 px short of the right edge of the other/],
    [() => above(p, e('gamma'), e('alpha')), /gamma element is not above alpha element — its bottom is 140 px below the top of the other/],
    [() => below(p, e('alpha'), e('gamma')), /alpha element is not below gamma element — its top is 140 px above the bottom of the other/],
    [() => inside(p, e('alpha'), e('gamma')), /alpha element is not inside gamma element — it sticks out 40 px at the top/],
    [() => clearOf(p, e('alpha'), e('epsilon')), /alpha element is not clear of epsilon element — they overlap by 50×20 px/],
    [() => leftAligned(p, e('beta'), e('alpha')), /beta element is not left-aligned with alpha element — its left edge is 110 px right of the other's/],
    [() => topAligned(p, e('alpha'), e('gamma')), /alpha element is not top-aligned with gamma element — its top is 40 px above the other's/],
    [() => centeredWith(p, e('gamma'), e('alpha')), /gamma element is not vertically centered with alpha element — its middle is 75 px below the other's/],
    [() => asTallAs(p, e('gamma'), e('alpha')), /gamma element is not as tall as alpha element — it is 100 px tall, the other 30 px/],
    [() => tallPx(p, e('alpha'), 40), /alpha element is 30 px tall, not 40/],
    [() => asWideAs(p, e('gamma'), e('alpha')), /gamma element is not as wide as alpha element — it is 210 px wide, the other 100 px/],
    [() => widePx(p, e('alpha'), 120), /alpha element is 100 px wide, not 120/],
    [() => atLeastWidePx(p, e('alpha'), 150), /alpha element is 100 px wide, not at least 150/],
    [() => atMostWidePx(p, e('alpha'), 60), /alpha element is 100 px wide, not at most 60/],
    [() => notClipped(p, e('wide')), /wide element is clipped: it is cut off by a div it sits in: 40 px on the right/],
    [() => notClipped(p, e('caption')), /caption element is clipped: its text needs \d+ px and has 40, ending in an ellipsis/],
    [() => leftOf(p, e('omega'), e('alpha')), /omega element is not shown/],
    [() => leftOf(p, e('twin'), e('alpha')), /^twin element names 2 shown elements: "twin" ×2$/],
  ]);
});

scenario('inside the window: a fixed box within the window holds; one past its bottom fails with the window\'s size', async () => {
  const size = page!.viewportSize()!;
  await page!.setViewportSize({width: 800, height: 600});
  try {
    await page!.setContent(`<div data-u2-name="shown" style="position: fixed; left: 10px; top: 10px; width: 100px; height: 30px"></div>
      <div data-u2-name="low" style="position: fixed; left: 10px; top: 580px; width: 100px; height: 40px"></div>`);
    await insideWindow(page!, e('shown'));
    await failWith([[() => insideWindow(page!, e('low')),
      /low element is not inside the window — it sticks out 20 px at the bottom: it spans 10\.\.110 × 580\.\.620, the window 800×600/]]);
  }
  finally {
    await page!.setViewportSize(size);
  }
});

scenario('no layout problems: the audit over the element\'s subtree, failing with its findings or when it holds no u2 control', async () => {
  await page!.setContent(STYLE + `<div data-u2="form" aria-label="Clean">${INPUT('Name')}${INPUT('Age')}</div>` +
    `<div data-u2="form" aria-label="Broken">${OVERLAPPING}</div><div data-u2-name="plain">text</div>`);
  await noLayoutProblems(page!, el('Clean form'));
  await failWith([
    [() => noLayoutProblems(page!, el('Broken form')),
      /Broken form has 1 layout problem\(s\): overlap: "Name" text input overlaps "Notes" text input by \d+×20 px/],
    [() => noLayoutProblems(page!, el('plain element')), /plain element holds no u2 control the layout audit reads/],
  ]);
});
