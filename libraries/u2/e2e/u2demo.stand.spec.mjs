/* Lane 2 — what only a stand can answer for U2 Demo: the platform URL route reaching a sub-demo,
   the Browse sub-tree under Apps ▸ Dev ▸ U2 Demo (appTreeBrowser resolves from the server's func
   registry), and a Browse leaf click navigating the open view instead of stacking a new one.
   Everything about in-app behaviour is Lane 1's job (checks/u2demo-*.spec.mjs). Run with
   `npm run e2e:stand-u2demo` on a stand nobody else is using; publish first
   (`npm run debug-u2demo-<host>` in packages/U2Demo). */
import {ARTIFACTS, note, ok, report, results, shot} from './local.mjs';
import {STAND_URL, openStand} from './stand.mjs';

const APP_PATH = '/apps/U2demo/U2Demo';
const GROUPS = ['Start', 'Inputs', 'Containers', 'Collections', 'Display', 'Forms', 'Platform'];

const demoState = (page) => page.evaluate(() => ({
  view: grok.shell.v?.name,
  path: grok.shell.v?.path,
  views: [...grok.shell.views].length,
  nav: !!document.querySelector('.u2demo-nav'),
}));

/** S1 — the deep link a shared URL is: the platform route fills the app function's `path` param. */
async function checkUrlRoute(page) {
  await page.goto(`${STAND_URL}${APP_PATH}/platform/dataframes`, {waitUntil: 'load', timeout: 120000});
  await page.waitForSelector('.u2demo-nav', {timeout: 120000});
  await page.waitForTimeout(2000);
  const state = await demoState(page);
  await shot(page, 'u2demo-stand-1-url');
  ok('S1/url-deep-link-opens-dataframes',
    state.nav && (state.path ?? '').endsWith('/platform/dataframes'), JSON.stringify(state));
}

const browseLabel = (page, text) => page.locator(
  '.d4-tree-view-group-label:visible, .d4-tree-view-item-label:visible',
  {hasText: new RegExp(`^${text}$`)}).first();

const DEMO_PATH = ['Apps', 'Dev', 'U2 Demo'];

/* A group expands on twistie click or caption double-click only — a label click selects and
   previews (tree_view.dart:605-611). `annotate` names the twistie after the pipe-joined path with
   every separator character turned into '-' (html_utils.dart:461). */
const expanderName = (path) => `tree-expander-${path.join(' | ')}`.replace(/[:_; *\\[\]{}|]/g, '-');

/* Expansion is persisted per user (tree_view.dart:321), so a blind toggle would collapse a node
   a previous run left open. */
async function expandTri(tri) {
  await tri.waitFor({timeout: 60000});
  const cls = (await tri.getAttribute('class')) ?? '';
  if (!cls.includes('d4-tree-view-tri-expanded'))
    await tri.click();
}

const expand = (page, path) => expandTri(page.locator(`[name="${expanderName(path)}"]`).first());

/** S2 — `Apps ▸ Dev ▸ U2 Demo` is a group (meta.role: appTreeBrowser + meta.app match) and its
 * lazy expansion yields the demo's group nodes. */
async function checkBrowseTree(page) {
  for (let i = 0; i < 3; i++) {
    if (await page.locator('.d4-tree-view-group-label', {hasText: /^Apps$/}).count() > 0)
      break;
    await page.locator('[name="icon-compass"]').first().click().catch(() => {});
    await page.waitForTimeout(2000);
  }
  await expand(page, ['Apps']);
  await page.waitForTimeout(2000);
  await expand(page, ['Apps', 'Dev']);
  await page.waitForTimeout(2000);
  // expand, don't open: the tri toggles the group; a label dblclick would launch the app instead
  const appGroup = page.locator('.d4-tree-view-group',
    {has: page.locator('.d4-tree-view-group-label', {hasText: /^U2 Demo$/})}).last();
  await expandTri(appGroup.locator('.d4-tree-view-tri').first());
  await page.waitForTimeout(3000);
  const labels = await page.evaluate(() =>
    [...document.querySelectorAll('.d4-tree-view-group-label')].map((e) => e.textContent.trim()));
  await shot(page, 'u2demo-stand-2-browse');
  ok('S2/browse-node-expands-into-the-demo-groups', GROUPS.every((g) => labels.includes(g)),
    `expected ${JSON.stringify(GROUPS)} among ${JSON.stringify(labels)}`);
}

/** S3 — a Browse leaf click opens the demo on that sub-demo and the URL says so. */
async function checkBrowseLeafOpens(page) {
  await expand(page, [...DEMO_PATH, 'Inputs']);
  await page.waitForTimeout(1500);
  await browseLabel(page, 'All inputs').click();
  await page.waitForTimeout(4000);
  const state = await demoState(page);
  await shot(page, 'u2demo-stand-3-browse-leaf');
  ok('S3/browse-leaf-click-lands-on-all-inputs',
    (state.path ?? '').endsWith('/inputs/all-inputs') &&
    page.url().toLowerCase().endsWith('/inputs/all-inputs'), `${page.url()} ${JSON.stringify(state)}`);
}

/** S4 — a second Browse leaf reuses the live view (openDemo) instead of stacking another one. */
async function checkBrowseReuse(page) {
  const before = (await demoState(page)).views;
  await expand(page, [...DEMO_PATH, 'Collections']);
  await page.waitForTimeout(1500);
  await browseLabel(page, 'Trees').click();
  await page.waitForTimeout(4000);
  const state = await demoState(page);
  await shot(page, 'u2demo-stand-4-reuse');
  ok('S4/second-browse-leaf-reuses-the-open-view',
    state.views === before && (state.path ?? '').endsWith('/collections/trees'),
    `views ${before} -> ${state.views}, ${JSON.stringify(state)}`);
}

const started = Date.now();
const {browser, page} = await openStand();
try {
  for (const check of [checkUrlRoute, checkBrowseTree, checkBrowseLeafOpens, checkBrowseReuse]) {
    try {
      await check(page);
    } catch (e) {
      ok(check.name, false, `check threw: ${String(e).slice(0, 300)}`);
      await shot(page, `u2demo-stand-${check.name}-failure`);
    }
  }
} finally {
  await browser.close();
}

note('stand/target', STAND_URL);
const code = report('u2demo.stand');
console.log(`${((Date.now() - started) / 1000).toFixed(1)}s · ${results.length} checks · artifacts in ${ARTIFACTS}`);
process.exit(code);
