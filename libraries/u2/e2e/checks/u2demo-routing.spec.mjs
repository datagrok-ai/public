/* Routing, Lane 1 — what local mode can prove: path → leaf resolution through the app function's
   `path` argument, and leaf → `grok.shell.v.path` on in-app navigation. The URL → app direction
   needs a server route and lives in e2e/u2demo.stand.spec.mjs. */
import {ok, shot} from '../local.mjs';
import {APP, demoPath} from '../lib.mjs';

const APP_PATH = '/apps/U2demo/U2Demo';
const NAV_ROW = '.u2demo-nav .u2-tree-row';

/** Opens U2 Demo with an explicit `path` argument through the staged bundle's export — the same
 * module fallback `openApp` takes, but with the argument the fixture-registered function may not
 * know yet. */
const openDemoAt = (page, path) => page.evaluate(async ({pkg, path}) => {
  const sibling = DG.Func.find({package: pkg})[0];
  if (!sibling)
    throw new Error(`package ${pkg} is not in the local fixture — stage it and refresh the fixture`);
  await sibling.package.load();
  const module = sibling.package.getModule('package.js');
  grok.shell.addView(await module.u2DemoApp(path));
}, {pkg: APP.package, path});

/** Clicks a leaf row in the visible nav tree (stacked hidden views keep their rows in the DOM). */
const clickLeaf = (page, label) => page.evaluate(({sel, label}) => {
  const row = [...document.querySelectorAll(sel)]
    .find((r) => (r.querySelector('.u2-tree-label')?.textContent ?? '').trim() === label &&
      r.offsetParent != null);
  if (row == null)
    return false;
  row.click();
  return true;
}, {sel: NAV_ROW, label});

const waitNav = (page) => page.waitForFunction((sel) =>
  [...document.querySelectorAll(sel)].some((r) => r.offsetParent != null), NAV_ROW, {timeout: 30000});

export async function fixture(page) {
  await page.evaluate(() => localStorage.removeItem('u2demo.tab'));
}

async function checkTwoSegmentPath(page) {
  await openDemoAt(page, '/collections/trees');
  await waitNav(page);
  const path = await demoPath(page);
  const treesShown = await page.evaluate(() =>
    [...document.querySelectorAll('.u2demo-content [data-u2="tree"]')].some((e) => e.offsetParent != null));
  await shot(page, 'u2demo-routing-1-trees');
  ok('u2demo-routing/1 open-at-two-segment-path', path === `${APP_PATH}/collections/trees` && treesShown,
    `path="${path}" treesShown=${treesShown}`);
}

async function checkOneSegmentPath(page) {
  await openDemoAt(page, '/funcs');
  await waitNav(page);
  await page.waitForSelector('.u2demo-funcs', {timeout: 30000});
  const path = await demoPath(page);
  ok('u2demo-routing/2 one-segment-path-resolves-to-funcs', path === `${APP_PATH}/forms/funcs`,
    `path="${path}"`);
}

async function checkNavigationWritesPath(page) {
  const clicked = await clickLeaf(page, 'Lists');
  await page.waitForTimeout(500);
  const path = await demoPath(page);
  ok('u2demo-routing/3 in-app-navigation-updates-view-path',
    clicked && path === `${APP_PATH}/collections/lists`, `clicked=${clicked} path="${path}"`);
}

async function checkReplaceStateNoHistoryGrowth(page) {
  const before = await page.evaluate(() => history.length);
  for (const label of ['Trees', 'Files', 'Form']) {
    await clickLeaf(page, label);
    await page.waitForTimeout(400);
  }
  const after = await page.evaluate(() => history.length);
  const path = await demoPath(page);
  ok('u2demo-routing/4 three-navigations-cost-no-history-entries',
    after === before && path === `${APP_PATH}/forms/form`,
    `history ${before} -> ${after}, path="${path}"`);
}

export const checks = [
  {id: 'u2demo-routing/1 u2DemoApp(\'/collections/trees\') opens on Trees', run: checkTwoSegmentPath},
  {id: 'u2demo-routing/2 a one-segment path resolves to Function forms', run: checkOneSegmentPath},
  {id: 'u2demo-routing/3 in-app navigation writes /apps/U2demo/U2Demo/<group>/<leaf>', run: checkNavigationWritesPath},
  {id: 'u2demo-routing/4 navigation replaces rather than pushes history', run: checkReplaceStateNoHistoryGrowth},
];
