/* The U2Demo nav shell (packages/U2Demo/src/demo.ts + nav.ts): the left VirtualTree over the
   2-level registry, per-page scope disposal on navigation, the shell status bar text and the
   demo ribbon commands. The headline count is 7 groups / 26 leaves — the ask-A table, the Display
   group the UX audit's M8 asked for, the Start group holding the Overview front door, and the two
   value editors (Range slider, Multi-select) the acceptance pass moved out of Display. */
import {consoleErrors, ok, pageErrors, shot} from '../local.mjs';
import {openDemoPage} from '../lib.mjs';

const NAV = '.u2demo-nav';
const CONTENT = '.u2demo-content';

/** The tree in order — every leaf label is its id's noun, so the tree, the URL and the Browse
 * sub-tree speak one vocabulary and nothing has to truncate. What a page covers lives in the
 * leaf's description (nav.ts), which the source panel shows. */
const LABELS = [
  'Start', 'Overview',
  'Inputs', 'All inputs', 'Basic inputs', 'Range slider', 'Multi-select', 'Async',
  'Containers', 'Layout', 'Popups',
  'Collections', 'Lists', 'Trees',
  'Display', 'Cards', 'Feedback', 'Tables', 'Sections & wizard', 'Message input',
  'Forms', 'Form', 'Property grid', 'Object form', 'Functions', 'FuncCalls',
  'Platform', 'Dataframes', 'Files', 'Entities', 'Spaces', 'Molecules', 'Bridge',
];

/** The Display group's leaves — the home the orphaned display controls got. */
const DISPLAY = ['Cards', 'Feedback', 'Tables', 'Sections & wizard', 'Message input'];

/** The nav pane the UX audit measured, where 8 of 18 labels truncated: 178 px, minus the twistie
 * and the row's indent. Viewport-independent, so the check means the same in both lanes. */
const PANE = 178;
const TWISTIE = 14;

const navRow = (page, label) => page.locator(`${NAV} .u2-tree-row`)
  .filter({has: page.locator('.u2-tree-label', {hasText: new RegExp(`^${label}$`)})}).first();

async function goTo(page, label) {
  await navRow(page, label).click();
  await page.waitForTimeout(400);
}

export async function fixture(page) {
  await openDemoPage(page, 'lists');
}

async function checkTree(page) {
  // the tree is virtualized: the row elements arrive over a layout pass or two, and `aria-setsize`
  // is what the list itself says the row count is — the DOM row count is only ever a lower bound
  await page.waitForFunction(({nav, total}) => {
    const rows = document.querySelectorAll(`${nav} .u2-list-row`);
    return rows.length === total &&
      [...rows].every((r) => r.getAttribute('aria-setsize') === String(total));
  }, {nav: NAV, total: LABELS.length}, {timeout: 10000}).catch(() => {});
  const counts = await page.evaluate((nav) => {
    // in list order, not DOM order: a recycled row is appended wherever the pool had it
    const rows = [...document.querySelectorAll(`${nav} .u2-list-row`)]
      .sort((a, b) => Number(a.dataset.index) - Number(b.dataset.index))
      .map((r) => r.querySelector('.u2-tree-row'))
      .filter((r) => r != null);
    const groups = rows.filter((r) =>
      r.querySelector('.u2-tree-twistie:not(.u2-tree-twistie-hidden)'));
    return {
      total: rows.length, groups: groups.length,
      claimed: Number(document.querySelector(`${nav} .u2-list-row`)?.getAttribute('aria-setsize') ?? 0),
      labels: rows.map((r) => r.querySelector('.u2-tree-label')?.textContent ?? ''),
      // the natural text width of each row: label + twistie + the row's own indent
      widths: rows.map((r) => Math.ceil(parseFloat(getComputedStyle(r).paddingLeft) +
        (r.querySelector('.u2-tree-label')?.scrollWidth ?? 0))),
    };
  }, NAV);
  const wide = counts.labels.filter((_, i) => counts.widths[i] + TWISTIE > PANE);
  await shot(page, 'u2demo-nav-1-tree');
  ok('u2demo-nav/1a the tree shows the 7 expanded groups and their 26 leaves',
    counts.groups === 7 && counts.claimed === LABELS.length && counts.total === LABELS.length,
    `groups=${counts.groups} rows=${counts.total} claimed=${counts.claimed} ` +
    `first=${JSON.stringify(counts.labels.slice(0, 6))}`);
  ok('u2demo-nav/1b every label is its id\'s noun and fits the nav pane without truncating',
    counts.labels.join('|') === LABELS.join('|') && wide.length === 0,
    `labels=${JSON.stringify(counts.labels)} tooWide=${JSON.stringify(wide)}`);
}

async function checkTreesLeaf(page) {
  await goTo(page, 'Trees');
  const shape = await page.evaluate((content) => {
    const host = document.querySelector(content);
    const strayLists = host ? [...host.querySelectorAll('.u2-list')]
      .filter((el) => !el.closest('[data-u2="tree"]')) : [];
    return {
      hasPage: host?.querySelector('.u2demo-page') != null,
      trees: host?.querySelectorAll('[data-u2="tree"]').length ?? 0,
      strayLists: strayLists.length,
    };
  }, CONTENT);
  await shot(page, 'u2demo-nav-2-trees');
  ok('u2demo-nav/2 clicking Trees renders the tree page, not the list page',
    shape.hasPage && shape.trees >= 1 && shape.strayLists === 0,
    `page=${shape.hasPage} trees=${shape.trees} strayLists=${shape.strayLists}`);
}

async function checkRebuild(page) {
  await page.evaluate((content) =>
    window.__u2NavPrev = document.querySelector(content).firstElementChild, CONTENT);
  await goTo(page, 'Lists');
  await goTo(page, 'Trees');
  const detached = await page.evaluate(() => window.__u2NavPrev?.isConnected === false);
  ok('u2demo-nav/3 navigating away and back rebuilds the page (the old root is detached)',
    detached, `detached=${detached}`);
}

async function checkStatus(page) {
  const text = await page.evaluate(() => [...document.querySelectorAll('div')]
    .find((x) => x.children.length === 0 &&
      (x.textContent ?? '').startsWith('Collections / Trees'))?.textContent ?? '');
  ok('u2demo-nav/4 the shell status bar text starts with the group label',
    text.startsWith('Collections / Trees'), `text="${text}"`);
}

async function checkRibbon(page) {
  const counts = await page.evaluate(() => ({
    icons: document.querySelectorAll('.d4-ribbon .u2-icon-btn').length,
    dropdowns: document.querySelectorAll('.d4-ribbon .u2-dropdown-btn').length,
  }));
  await shot(page, 'u2demo-nav-5-ribbon');
  ok('u2demo-nav/5 the ribbon holds the 3 icon commands and one drop-down',
    counts.icons === 3 && counts.dropdowns === 1,
    `icons=${counts.icons} dropdowns=${counts.dropdowns}`);
}

async function checkDisplayLeaves(page) {
  const before = consoleErrors.length + pageErrors.length;
  const rendered = [];
  for (const label of DISPLAY) {
    await goTo(page, label);
    rendered.push(await page.evaluate((content) =>
      document.querySelector(`${content} .u2demo-page`)?.childElementCount ?? 0, CONTENT));
  }
  const errors = consoleErrors.length + pageErrors.length - before;
  await shot(page, 'u2demo-nav-6-display');
  ok('u2demo-nav/6 every Display leaf renders, none of them with a console error',
    rendered.every((n) => n > 3) && errors === 0,
    `children=${JSON.stringify(rendered)} newErrors=${errors}`);
}

/** The front door: Overview is the first leaf of the first group and the app's default, it names
 * the six content areas, and its "start here" pointers navigate. */
async function checkOverview(page) {
  await goTo(page, 'Overview');
  const state = await page.evaluate((content) => {
    const page_ = document.querySelector(`${content} .u2demo-page`);
    return {
      areas: [...(page_?.querySelectorAll('.u2demo-overview-area') ?? [])].map((e) => e.textContent),
      links: [...(page_?.querySelectorAll('a.u2-link') ?? [])].map((e) => e.textContent),
    };
  }, CONTENT);
  const clicked = await page.locator(`${CONTENT} a.u2-link`).first().click()
    .then(() => true).catch(() => false);
  await page.waitForTimeout(600);
  const landed = await page.evaluate(() => grok.shell.v?.path ?? '');
  await shot(page, 'u2demo-nav-9-overview');
  ok('u2demo-nav/9 Overview names every area and its start-here pointers navigate',
    state.areas.length === 6 && !state.areas.includes('Start') && state.links.length === 3 &&
    clicked && landed.endsWith('/inputs/all-inputs'),
    `areas=${JSON.stringify(state.areas)} links=${JSON.stringify(state.links)} path="${landed}"`);
}

/** The acceptance pass at 900 px: the demo column measured ~85 px, words broke mid-word, the
 * range slider collapsed onto its two handles and the nav vanished. The shell reads its OWN width
 * through a container query (css/u2demo.css), so the check narrows the window, measures, and puts
 * it back — every check after this one runs at the standard viewport again. */
async function checkNarrowLayout(page) {
  await goTo(page, 'Range slider');
  const measure = () => {
    const panels = [...document.querySelectorAll('.u2demo-shell > .u2-splitter-panel')];
    const track = document.querySelector('.u2demo-content .u2-range-track');
    return {
      shell: Math.round(document.querySelector('.u2demo-root')?.clientWidth ?? 0),
      nav: panels[0]?.offsetParent == null ? 0 : Math.round(panels[0].getBoundingClientRect().width),
      host: Math.round(document.querySelector('.u2demo-content')?.clientWidth ?? 0),
      // what the reader actually reads: the host scrolls under it rather than squeezing it
      page: Math.round(document.querySelector('.u2demo-content .u2demo-page')?.clientWidth ?? 0),
      track: Math.round(track?.getBoundingClientRect().width ?? 0),
      panel: document.querySelector('.grok-context-panel, .d4-accordion') != null,
    };
  };
  const wide = await page.evaluate(measure);
  await page.setViewportSize({width: 900, height: 1000});
  await page.waitForTimeout(800);
  const narrow = await page.evaluate(measure);
  await page.setViewportSize({width: 1600, height: 1000});
  await page.waitForTimeout(800);
  await shot(page, 'u2demo-nav-8-narrow');
  // folded away (0) or a rail wide enough to read — never a squeezed sliver
  const navOk = narrow.nav === 0 || narrow.nav >= 140;
  ok('u2demo-nav/8 at 900px the nav folds deliberately and the demo keeps a readable column',
    navOk && narrow.page >= 260 && narrow.track >= 60 && wide.host > narrow.host,
    `wide=${JSON.stringify(wide)} narrow=${JSON.stringify(narrow)}`);
}

/** The Tour is a document-level dim layer and the compose box's mentions are an overlay popup:
 * neither is inside the page, so only the page's own disposal can take them down. A stuck one is
 * both a dead app (every click lands on the dim layer) and the known OOM trap for this runner. */
async function checkOverlaysReleased(page) {
  await goTo(page, 'Feedback');
  await page.locator(`${CONTENT} button`).filter({hasText: 'Start tour'}).first().click();
  const tourUp = await page.waitForSelector('.u2-tour-overlay', {timeout: 8000})
    .then(() => true).catch(() => false);
  await goTo(page, 'Message input');
  const tourGone = await page.evaluate(() => document.querySelector('.u2-tour-overlay') == null &&
    document.querySelector('.u2-tour-popup') == null);

  await page.locator(`${CONTENT} .u2-msg-editor`).first().click();
  await page.keyboard.type('@a', {delay: 30});
  await page.waitForTimeout(800);
  const mentionsUp = await page.evaluate(() => document.querySelector('.u2-msg-popup') != null);
  await goTo(page, 'Cards');
  const mentionsGone = await page.evaluate(() => document.querySelector('.u2-msg-popup') == null);
  await shot(page, 'u2demo-nav-7-overlays');
  ok('u2demo-nav/7 leaving the Feedback and Message input leaves closes their overlays',
    tourUp && tourGone && mentionsGone,
    `tourUp=${tourUp} tourGone=${tourGone} mentionsOpened=${mentionsUp} mentionsGone=${mentionsGone}`);
  // the verdict is recorded: whatever is still up must not poison every check after this one
  await page.evaluate(() => {
    for (const el of document.querySelectorAll('.u2-tour-overlay, .u2-tour-popup, .u2-msg-popup'))
      el.remove();
  });
}

export const checks = [
  {id: 'u2demo-nav/1 the tree shows the groups expanded, labelled by id noun', run: checkTree},
  {id: 'u2demo-nav/2 clicking Trees renders the tree page', run: checkTreesLeaf},
  {id: 'u2demo-nav/3 navigation disposes and rebuilds pages', run: checkRebuild},
  {id: 'u2demo-nav/4 the shell status bar reports group / leaf', run: checkStatus},
  {id: 'u2demo-nav/5 the ribbon holds the demo commands', run: checkRibbon},
  {id: 'u2demo-nav/6 every Display leaf renders cleanly', run: checkDisplayLeaves},
  {id: 'u2demo-nav/7 the tour and the mention popup die with their page', run: checkOverlaysReleased},
  {id: 'u2demo-nav/8 the shell degrades sensibly at 900px', run: checkNarrowLayout},
  {id: 'u2demo-nav/9 Overview is the front door and points at three sub-demos', run: checkOverview},
];
