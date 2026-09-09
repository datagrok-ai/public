/* The filter builder and the query input in a real browser: row editing through the native
   selects, the horizontal wrap, completion at the caret, advanced mode (nesting, `not`, the mode
   guard), the drag gesture and template locks — against the gallery's literal schema — then the
   U2Demo page, where the tree filters a DataFrame and rides the app URL as `?q=`. Gallery-hosted
   (`e2e/gallery-server.mjs`), so it ends by rebooting the client for the suites after it. */
import {BOOT_TIMEOUT, LOCAL_URL, consoleErrors, note, ok, pageErrors, shot} from '../local.mjs';
import {APP, reopenApp} from '../lib.mjs';
import {startGalleryServer} from '../gallery-server.mjs';

const APP_PATH = '/apps/U2demo/U2Demo';
const ROW = '[data-u2="filter-row"]';
const GROUP = '[data-u2="filter-group"]';
const POPUP = '.u2-fq-popup';

let server;

const builder = (page, name) => page.locator(`[data-u2="filter-builder"][data-u2-name="${name}"]`).first();
const rows = (page, name) => builder(page, name).locator(ROW);
const queryOf = (page, name) => builder(page, name).locator('[data-u2-part="query"]').first().textContent();
const modeLink = (page, name) => builder(page, name).locator('[data-u2-part="mode"]').first();
const hidden = (locator) => locator.evaluate((el) => el.hidden || el.offsetParent == null);

/** The rows' ids in document order — a group precedes its children, so the order is the shape. */
const nodeIds = (page, name) => builder(page, name).evaluate((b) =>
  [...b.querySelectorAll('[data-u2-node]')].map((el) => `${el.dataset.u2 === 'filter-group' ? 'G' : 'r'}:${el.dataset.u2Node}`));
const errors = () => consoleErrors.length + pageErrors.length;

/** `page.mouse` speaks viewport coordinates: a section under the fold has to be brought up first. */
async function bringUp(locator) {
  await locator.evaluate((el) => el.scrollIntoView({block: 'center'}));
  await locator.page().waitForTimeout(150);
}

async function openPage(page, hash, name) {
  await page.goto(`${server.url}/gallery/#${hash}`);
  await page.waitForSelector(`[data-u2="filter-builder"][data-u2-name="${name}"]`, {timeout: 30000});
}

async function setRow(row, property, operator, value) {
  await row.locator('[data-u2-part="prop"] select').selectOption(property);
  await row.locator('[data-u2-part="op"] select').selectOption(operator);
  const input = row.locator('[data-u2-part="value"] input').first();
  await input.fill(String(value));
  await input.press('Enter');
  await input.press('Tab');
}

/** A press-move-release on a row's grip: the four-pixel threshold, then the target. */
async function drag(page, handle, to, {release = true} = {}) {
  const box = await handle.boundingBox();
  await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
  await page.mouse.down();
  await page.mouse.move(box.x + box.width / 2 + 12, box.y + box.height / 2 + 12);
  await page.mouse.move(to.x, to.y, {steps: 6});
  await page.waitForTimeout(100);
  if (release)
    await page.mouse.up();
  await page.waitForTimeout(150);
}

const center = async (locator) => {
  const box = await locator.boundingBox();
  return {x: box.x + box.width / 2, y: box.y + box.height / 2};
};

export async function fixture(page) {
  server = await startGalleryServer();
  await openPage(page, 'filters', 'empty');
}

async function checkAddRow(page) {
  await builder(page, 'empty').locator('.u2-fb-header [data-u2-part="add"]').click();
  await page.waitForTimeout(150);
  const added = await rows(page, 'empty').count();
  await setRow(rows(page, 'empty').first(), 'age', '>', 30);
  await page.waitForTimeout(150);
  const query = await queryOf(page, 'empty');
  await shot(page, 'filters-1-row-added');
  ok('filters/1a/plus-adds-one-row-and-property-operator-value-write-the-query',
    added === 1 && query === 'age > 30', `rows=${added} query="${query}"`);
}

/** Two rows and the chip between them: at 720px the chip shares the first row's line and the
 * second row wraps; at 400px the chip drops onto a line of its own. A row is never split. */
async function checkHorizontalWrap(page) {
  const measure = () => builder(page, 'criteriaH').evaluate((b) => {
    const host = b.parentElement.getBoundingClientRect();
    const parts = [...b.querySelectorAll('[data-u2="filter-row"], .u2-fb-join')].map((el) => el.getBoundingClientRect())
      .sort((a, c) => a.top - c.top);
    // the chip is centred on its row, so a shared line is an overlap in y, not an equal top
    let lines = 0;
    let bottom = -Infinity;
    for (const r of parts) {
      if (r.top >= bottom - 2) {
        lines++;
        bottom = r.bottom;
      }
    }
    return {host: Math.round(host.width), lines, widest: Math.round(Math.max(...parts.map((r) => r.width)))};
  });
  await page.getByRole('button', {name: '400px', exact: true}).first().click();
  await page.waitForTimeout(200);
  const narrow = await measure();
  await shot(page, 'filters-2-horizontal-wrap');
  await page.getByRole('button', {name: '720px', exact: true}).first().click();
  await page.waitForTimeout(200);
  const wide = await measure();
  ok('filters/2a/horizontal-rows-and-their-chip-wrap-onto-more-lines-as-the-host-narrows',
    wide.host === 720 && wide.lines === 2 && narrow.host === 400 && narrow.lines === 3,
    `720px=${JSON.stringify(wide)} 400px=${JSON.stringify(narrow)}`);
  if (narrow.widest > narrow.host)
    note('filters/2b/a-row-wider-than-the-host-overflows-it', `row ${narrow.widest}px in a ${narrow.host}px host`);
}

async function checkQueryInputCompletion(page) {
  await page.goto(`${server.url}/gallery/#filter-query`);
  const input = page.locator('[data-u2="filter-query-input"][data-u2-name="query"] input').first();
  await input.waitFor({timeout: 30000});
  await input.click();
  await page.keyboard.press('Control+a');
  await page.keyboard.type('na', {delay: 20});
  await page.waitForSelector(`${POPUP} .u2-fq-option`, {timeout: 5000});
  const offered = await page.locator(`${POPUP} .u2-fq-option`).allTextContents();
  await page.keyboard.press('Enter');
  await page.waitForTimeout(150);
  const afterProperty = await input.inputValue();
  await page.keyboard.type('=', {delay: 20});
  await page.waitForTimeout(150);
  await page.keyboard.press('Enter');
  await page.waitForTimeout(150);
  const afterOperator = await input.inputValue();
  await page.keyboard.type('"x"', {delay: 20});
  await page.waitForTimeout(150);
  await shot(page, 'filters-3-completion-before-commit');
  await page.keyboard.press('Enter');
  await page.waitForTimeout(200);
  const committed = await input.inputValue();
  const count = await rows(page, 'criteria').count();
  const prop = await rows(page, 'criteria').first().locator('[data-u2-part="prop"] select').inputValue();
  await shot(page, 'filters-3-completion-committed');
  ok('filters/3a/na-offers-name-then-the-operator-then-enter-commits-one-builder-row',
    offered.some((t) => t.includes('Name')) && afterProperty === 'name ' && afterOperator === 'name = ' &&
    committed === 'name = "x"' && count === 1 && prop === 'name',
    `offered=${JSON.stringify(offered)} "${afterProperty}" → "${afterOperator}" → "${committed}" rows=${count} prop=${prop}`);
}

/** Advanced mode on the "Empty, then add" builder, which check 1 left holding `age > 30`. */
async function checkNesting(page) {
  await openPage(page, 'filters', 'empty');
  const b = builder(page, 'empty');
  await b.locator('.u2-fb-header [data-u2-part="add"]').click();
  await setRow(rows(page, 'empty').first(), 'age', '>', 30);
  await modeLink(page, 'empty').click();
  await page.waitForTimeout(150);
  await b.locator('.u2-fb-header [data-u2-part="add-group"]').click();
  await page.waitForTimeout(150);
  const g1 = b.locator(GROUP).first();
  await g1.locator('.u2-fb-group-header [data-u2-part="add-group"]').first().click();
  await page.waitForTimeout(150);
  const g2 = g1.locator(GROUP).first();
  await g2.locator('.u2-fb-group-header [data-u2-part="add"]').first().click();
  await page.waitForTimeout(150);
  await setRow(g2.locator(ROW).first(), 'age', '>', 5);
  await page.waitForTimeout(150);
  const ids = await nodeIds(page, 'empty');
  const nested = await g1.evaluate((el) => {
    const inner = el.querySelector('[data-u2="filter-group"]');
    return inner != null && inner.querySelector('[data-u2="filter-row"]') != null;
  });
  await shot(page, 'filters-4-nested-groups');
  ok('filters/4a/advanced-nests-a-group-in-a-group-and-a-group-precedes-its-children',
    ids.length === 4 && ids[0].startsWith('r:') && ids[1].startsWith('G:') && ids[2].startsWith('G:') &&
    ids[3].startsWith('r:') && nested, `nodes=${JSON.stringify(ids)} nested=${nested}`);
}

async function checkNot(page) {
  const g1 = builder(page, 'empty').locator(GROUP).first();
  await g1.locator('.u2-fb-group-header [data-u2-part="not"]').first().click();
  await page.waitForTimeout(150);
  const negated = await g1.evaluate((el) => el.classList.contains('u2-fb-negated'));
  const query = await queryOf(page, 'empty');
  await shot(page, 'filters-5-not-group');
  ok('filters/5a/the-not-chip-negates-the-group-and-the-query-reads-not',
    negated && query.includes('not ('), `negated=${negated} query="${query}"`);
  await g1.locator('.u2-fb-group-header [data-u2-part="not"]').first().click();
  await page.waitForTimeout(150);
}

/** The guard: `or [G2, cond]` under an `and` root cannot flatten; aligning the connector can. */
async function checkModeGuard(page) {
  const b = builder(page, 'empty');
  const g1 = b.locator(GROUP).first();
  const hint = b.locator('[data-u2-part="hint"]').first();
  const flatten = hint.locator('a').first();
  await g1.locator('.u2-fb-group-header [data-u2-part="add"]').first().click();
  await page.waitForTimeout(150);
  await modeLink(page, 'empty').click();
  await page.waitForTimeout(150);
  const refused = !(await hidden(hint)) && (await modeLink(page, 'empty').textContent()) === 'Simple';
  const noFlatten = await hidden(flatten);
  await shot(page, 'filters-6-simple-refused');
  await g1.locator('[data-u2-part="connector"] button').filter({hasText: /^and$/}).first().click();
  await page.waitForTimeout(150);
  await modeLink(page, 'empty').click();
  await page.waitForTimeout(150);
  const offered = !(await hidden(hint)) && !(await hidden(flatten));
  await flatten.click();
  await page.waitForTimeout(200);
  const groups = await b.locator(GROUP).count();
  const mode = await modeLink(page, 'empty').textContent();
  const count = await rows(page, 'empty').count();
  await shot(page, 'filters-6-flattened');
  ok('filters/6a/simple-is-refused-on-mixed-connectors-and-flatten-works-once-they-align',
    refused && noFlatten && offered && groups === 0 && mode === 'Advanced' && count === 3,
    `refused=${refused} noFlatten=${noFlatten} offered=${offered} groups=${groups} mode=${mode} rows=${count}`);
}

/** The "Advanced" section: `age > 30`, `or [status in …, and [sex = F, mw between]]`, `created > -1w`. */
async function checkDrag(page) {
  const b = builder(page, 'advanced');
  await bringUp(b);
  const before = await nodeIds(page, 'advanced');
  const last = b.locator(`.u2-fb-rows > ${ROW}`).last();
  const lastId = await last.getAttribute('data-u2-node');
  const orGroup = b.locator(GROUP).first();
  await drag(page, last.locator('[data-u2-part="handle"]'), await center(orGroup.locator('.u2-fb-group-header').first()));
  const into = await orGroup.evaluate((el, id) => {
    const children = [...el.querySelector('[data-u2-part="rows"]').children];
    return children.length > 0 && children[children.length - 1].dataset.u2Node === id;
  }, lastId);
  await shot(page, 'filters-7-dragged-into-group');

  const first = b.locator(`.u2-fb-rows > ${ROW}`).first();
  const box = await first.boundingBox();
  const moved = b.locator(`[data-u2-node="${lastId}"]`).first();
  await drag(page, moved.locator('[data-u2-part="handle"]'), {x: box.x + box.width / 2, y: box.y + 3});
  const line = await b.locator('.u2-fb-rows > [data-u2-node]').first().getAttribute('data-u2-node') === lastId;
  await shot(page, 'filters-7-dragged-before-first');

  const order = await nodeIds(page, 'advanced');
  const indicator = b.locator('.u2-fb-drop').first();
  await drag(page, first.locator('[data-u2-part="handle"]'),
    await center(orGroup.locator('.u2-fb-group-header').first()), {release: false});
  const shown = !(await hidden(indicator));
  await page.keyboard.press('Escape');
  await page.waitForTimeout(100);
  const gone = await hidden(indicator);
  await page.mouse.up();
  await page.waitForTimeout(150);
  const unchanged = JSON.stringify(await nodeIds(page, 'advanced')) === JSON.stringify(order);
  ok('filters/7a/a-row-drags-into-a-group-header-and-between-rows-and-escape-cancels',
    into && line && shown && gone && unchanged,
    `before=${JSON.stringify(before)} into=${into} line=${line} indicator shown=${shown} gone=${gone} unchanged=${unchanged}`);
}

async function checkTemplateLocks(page) {
  const b = builder(page, 'templateAdvanced');
  await bringUp(b);
  const locked = b.locator(`${GROUP}.u2-fb-locked`).first();
  // a locked property/operator is plain text (`.u2-fb-prop-text` / `.u2-fb-op-text`), the pickers hidden
  const lockedRows = await locked.locator(ROW).evaluateAll((els) => els.map((el) => ({
    locked: el.classList.contains('u2-fb-locked'),
    text: ['prop', 'op'].every((part) => {
      const text = el.querySelector(`.u2-fb-${part}-text[data-u2-part="${part}-text"]`);
      return text != null && !text.hidden && text.textContent !== '' &&
        el.querySelector(`[data-u2-part="${part}"]`)?.hidden === true;
    }),
    grip: !(el.querySelector('[data-u2-part="handle"]')?.hidden ?? true),
  })));
  const free = b.locator(`${GROUP}:not(.u2-fb-locked)`).first();
  const freeRow = free.locator(ROW).first();
  const freeId = await freeRow.getAttribute('data-u2-node');
  const indicator = b.locator('.u2-fb-drop').first();
  await drag(page, freeRow.locator('[data-u2-part="handle"]'),
    await center(locked.locator('.u2-fb-group-header').first()), {release: false});
  const refused = await hidden(indicator);
  await page.mouse.up();
  await page.waitForTimeout(150);
  const stayed = await free.locator(`[data-u2-node="${freeId}"]`).count() === 1;
  const simpleOnlyLink = await hidden(modeLink(page, 'simpleOnly'));
  await shot(page, 'filters-8-template-locks');
  ok('filters/8a/a-locked-group-renders-its-rows-as-text-hides-grips-refuses-drops-and-simple-only-hides-the-mode-link',
    lockedRows.length === 2 && lockedRows.every((r) => r.locked && r.text && !r.grip) && refused && stayed &&
    simpleOnlyLink,
    `rows=${JSON.stringify(lockedRows)} dropRefused=${refused} stayed=${stayed} simpleOnlyLinkHidden=${simpleOnlyLink}`);
}

/** The gallery lives on its own origin, so the client has to be re-booted, not just reopened. */
async function restoreClient(page) {
  await server.close();
  await page.goto(`${LOCAL_URL}/login.html?mode=local`, {waitUntil: 'load', timeout: BOOT_TIMEOUT});
  await page.waitForFunction(() => {
    try {
      return !!(window.DG && DG.Func && grok.shell.user);
    } catch (e) {
      return false;
    }
  }, null, {timeout: BOOT_TIMEOUT});
}

/** Opens U2 Demo with an explicit `path` through the staged bundle's export (u2demo-routing). */
const openDemoAt = (page, path) => page.evaluate(async ({pkg, path}) => {
  const sibling = DG.Func.find({package: pkg})[0];
  if (!sibling)
    throw new Error(`package ${pkg} is not in the local fixture — stage it and refresh the fixture`);
  await sibling.package.load();
  const module = sibling.package.getModule('package.js');
  grok.shell.addView(await module.u2DemoApp(path));
}, {pkg: APP.package, path});

const rowsReadout = (page) => page.evaluate(() => [...document.querySelectorAll('.u2demo-content .u2demo-status')]
  .find((el) => el.firstElementChild?.textContent === 'rows = ')?.lastElementChild?.textContent ?? '');

async function checkDemoPage(page) {
  await restoreClient(page);
  const q = 'city = "Basel"';
  const errorsAt = [errors()];
  await openDemoAt(page, `/platform/filters?q=${encodeURIComponent(q)}`);
  await page.waitForSelector('.u2demo-content [data-u2="filter-builder"][data-u2-name="orders"]', {timeout: 30000});
  await page.waitForTimeout(300);
  errorsAt.push(errors());
  const inbound = await page.evaluate(() => grok.shell.v.path);
  const seeded = await rows(page, 'orders').count();
  // the inbound query is applied as soon as the pane is up: no Apply needed
  const basel = await rowsReadout(page);
  errorsAt.push(errors());
  await shot(page, 'filters-9-demo-applied');

  const input = page.locator('.u2demo-content [data-u2-name="ordersQuery"] input').first();
  await input.click();
  await page.keyboard.press('Control+a');
  await page.keyboard.type('total > 1000', {delay: 20});
  await page.keyboard.press('Enter');
  await page.waitForTimeout(300);
  const path = await page.evaluate(() => grok.shell.v.path);
  await page.locator('.u2demo-content').getByRole('button', {name: 'Apply', exact: true}).first().click();
  await page.waitForTimeout(400);
  const big = await rowsReadout(page);
  errorsAt.push(errors());
  await shot(page, 'filters-9-demo-requeried');
  // another tab: the path names it, and carries that tab's (empty) query — the DataFrame query stays
  await page.locator('.u2demo-content .u2-tabs-tab', {hasText: /^Entity type$/}).first().click();
  await page.waitForTimeout(300);
  const entityPath = await page.evaluate(() => grok.shell.v.path);
  await page.locator('.u2demo-content .u2-tabs-tab', {hasText: /^DataFrame$/}).first().click();
  await page.waitForTimeout(300);
  const backPath = await page.evaluate(() => grok.shell.v.path);
  await shot(page, 'filters-9-demo-tab-path');
  await reopenApp(page);
  errorsAt.push(errors());
  ok('filters/9a/the-inbound-q-seeds-the-tree-and-applies-the-path-follows-the-query-and-the-tab',
    inbound === `${APP_PATH}/platform/filters?q=${encodeURIComponent(q)}` && seeded === 1 &&
    basel === '2 of 6' && big === '3 of 6' &&
    path === `${APP_PATH}/platform/filters?q=${encodeURIComponent('total > 1000')}` &&
    entityPath === `${APP_PATH}/platform/filters?tab=entity` && backPath === path &&
    errorsAt.every((n) => n === errorsAt[0]),
    `inbound="${inbound}" seeded=${seeded} rows ${basel} (auto-applied) → ${big} path="${path}" ` +
    `entity="${entityPath}" back="${backPath}" errors at boot/open/requery/close=${JSON.stringify(errorsAt)}`);
}

export const checks = [
  {id: 'filters/1 + adds a row; property, operator and value write the query', run: checkAddRow},
  {id: 'filters/2 horizontal rows wrap', run: checkHorizontalWrap},
  {id: 'filters/3 query input completion commits into the builder', run: checkQueryInputCompletion},
  {id: 'filters/4 advanced mode nests groups two levels deep', run: checkNesting},
  {id: 'filters/5 the not chip', run: checkNot},
  {id: 'filters/6 the simple/advanced guard and Flatten', run: checkModeGuard},
  {id: 'filters/7 drag into a group, between rows, and Escape', run: checkDrag},
  {id: 'filters/8 template locks', run: checkTemplateLocks},
  {id: 'filters/9 U2Demo: ?q= in, df.filter out, path follows', run: checkDemoPage},
];
