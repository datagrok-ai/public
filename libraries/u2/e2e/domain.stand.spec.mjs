/* Lane 2 — the real-browser interaction gate for the domain apps (consolidation 2-6). What no
   headless test can answer: geometry at the width every L3 failure was found at (1335 × 900, the
   context panel OPEN), the browser's own history, a physical hit target, an autocomplete popup's
   rectangles, a virtualizer that was measured at zero height, and the staged restore end to end
   against the real server. Six checks over the real Stockroom tables through the platform's
   `/domains/<schema>/<table>` route.

   Every check asserts STATE and GEOMETRY, never only that a selector exists — a shim's rectangles
   are supplied data, which is why these checks are here at all.

   Run with `npm run e2e:stand-domain` against a stand nobody else is using and that is running the
   batch's code (the lead publishes U2Demo/Stockroom/Grit first); `U2_STAND_URL` points at the
   browser-facing composite, e.g. `http://localhost:56690` for the ems-u2 worktree. */
import {execSync} from 'child_process';
import {ARTIFACTS, launch, note, ok, report, results, shot} from './local.mjs';
import {STAND_URL, openStand} from './stand.mjs';

/** The width every L3 geometry failure was found at. */
const VIEWPORT = {width: 1335, height: 900};
const SUBSTANCE = '/domains/stockroom/substance';
const CONTAINER = '/domains/stockroom/container';
const TAG = `${Date.now()}`.slice(-9);
const PROBE = `E2E-${TAG}`;

/** The probe row every write-bearing check acts on, so nothing here touches a row a person made. */
let probeId = null;

/* ------------------------------------------------------------------ the page-side primitives */

/** The element's rect, its host's rect and the viewport — one round trip, JSON-safe. */
const rectsOf = (page, selector, host) => page.evaluate(({selector, host}) => {
  const box = (el) => {
    const r = el.getBoundingClientRect();
    return {x: r.x, y: r.y, width: r.width, height: r.height, right: r.right, bottom: r.bottom};
  };
  const el = document.querySelector(selector);
  if (el === null)
    return null;
  const owner = host === null ? el.parentElement : document.querySelector(host);
  return {rect: box(el), host: owner === null ? null : box(owner),
    view: {width: window.innerWidth, height: window.innerHeight}};
}, {selector, host});

const visible = (r) => r !== null && r.width > 1 && r.height > 1;
const inViewport = (r, view) => r.x >= 0 && r.y >= 0 && r.right <= view.width + 1 && r.bottom <= view.height + 1;
/** Fully inside its host, give or take a rounding pixel — a control clipped by its own pane is
 * the defect class check 1 exists for. */
const within = (r, host) => host !== null && r.x >= host.x - 1 && r.y >= host.y - 1 &&
  r.right <= host.right + 1 && r.bottom <= host.bottom + 1;

/** Whether the first list row is reachable at its own coordinates, and what is over it when not. */
const topRowHit = (page) => page.evaluate(() => {
  const row = document.querySelector('[data-u2="domain-list"] [data-u2-row]');
  if (row === null)
    return 'no rows';
  const r = row.getBoundingClientRect();
  const hit = document.elementFromPoint(r.x + Math.min(40, r.width / 2), r.y + r.height / 2);
  if (hit === null)
    return 'nothing at the row point';
  if (row.contains(hit))
    return '';
  return hit.closest('.u2-suggest-popup') === null ? `covered by ${hit.className}` : 'covered by the suggestion popup';
});

/** What the app says about itself, read from the DOM alone — no handle on the app object. */
const appState = (page) => page.evaluate(() => ({
  app: !!document.querySelector('[data-u2="domain-app"]'),
  crumbs: document.querySelector('[data-u2="breadcrumbs"]')?.textContent?.trim() ?? '',
  rows: document.querySelectorAll('[data-u2="domain-list"] [data-u2-row]').length,
  deleted: document.querySelectorAll('[data-u2="domain-list"] .u2-domain-list-deleted').length,
  restored: document.querySelectorAll('[data-u2="domain-list"] .u2-domain-list-restored').length,
  saveDisabled: document.querySelector('[data-u2="save-button"]')?.disabled ?? null,
  url: location.pathname + location.search,
}));

/* --------------------------------------------------------------------------- the way around */

const address = (path, query) => `${STAND_URL}${path}${query === undefined ? '' :
  `?q=${encodeURIComponent(query)}`}`;

/** A domain app open at `path`, its list settled. */
async function openDomain(page, path, query, extra = '') {
  await page.goto(`${address(path, query)}${extra}`, {waitUntil: 'load', timeout: 180000});
  await page.waitForSelector('[data-u2="domain-app"]', {timeout: 180000});
  await page.waitForTimeout(3000);
}

/** The ⋯ menu's item by label. */
async function menuPick(page, label) {
  await page.locator('[data-u2="actions-menu"]').first().click();
  await page.waitForSelector('.u2-menu', {timeout: 15000});
  await page.locator('.u2-menu .u2-menu-item', {hasText: new RegExp(`^${label}$`)}).first().click();
}

/** A row's context menu item by label — the menu is the only place every row action is offered. */
async function rowPick(page, label) {
  await page.locator('[data-u2="domain-list"] [data-u2-row]').first().click({button: 'right'});
  await page.waitForSelector('.u2-menu', {timeout: 15000});
  await page.locator('.u2-menu .u2-menu-item', {hasText: new RegExp(`^${label}$`)}).first().click();
}

/* --------------------------------------------------------------------------------- the checks */

/** 1 — the list at 1335 px with the context panel open: the filter row and NEW are both whole,
 * inside the viewport and inside their own host. The geometry class of defect no headless test
 * can see (U29's real question, a check instead of a screenshot). */
async function checkGeometry(page) {
  await openDomain(page, SUBSTANCE);
  await page.evaluate(() => grok.shell.windows.showContextPanel = true);
  await page.waitForTimeout(2500);
  const filters = await rectsOf(page, '[data-u2="domain-filters"]', '[data-u2-part="list-page"]');
  const add = await rectsOf(page, '[data-u2="new-button"]', null);
  await shot(page, 'domain-stand-1-geometry');
  const pass = filters !== null && add !== null &&
    visible(filters.rect) && within(filters.rect, filters.host) && inViewport(filters.rect, filters.view) &&
    visible(add.rect) && within(add.rect, add.host) && inViewport(add.rect, add.view);
  ok('domain-stand/1 the filter row and NEW are whole at 1335 px with the context panel open',
    pass, JSON.stringify({filters, add}));
}

/** 2 — ⋯ → Trash, Back, Forward: one history entry per move, and the list's mode following the
 * address in both directions. */
async function checkTrashHistory(page) {
  await openDomain(page, SUBSTANCE);
  const before = await appState(page);
  await menuPick(page, 'Trash');
  await page.waitForTimeout(3000);
  const trash = await appState(page);
  await shot(page, 'domain-stand-2-trash');
  await page.goBack();
  await page.waitForTimeout(3000);
  const back = await appState(page);
  await page.goForward();
  await page.waitForTimeout(3000);
  const forward = await appState(page);
  await shot(page, 'domain-stand-2-forward');
  const inTrash = (s) => s.url.includes('trash=1') && /Trash/.test(s.crumbs);
  ok('domain-stand/2 ⋯ → Trash, Back and Forward move the list with the address',
    before.app && !inTrash(before) && inTrash(trash) && !inTrash(back) && back.app && inTrash(forward),
    JSON.stringify({before, trash, back, forward}));
}

/** 3 — a PHYSICAL click on a bulk-edit include checkbox, then OK. The hit target U23 is about:
 * the box is reachable at its own coordinates rather than covered by the label cell, the click
 * enables its editor, and the edit lands on the row the list is narrowed to. */
async function checkBulkCheckbox(page) {
  const value = `E2E${TAG}`;
  await openDomain(page, SUBSTANCE, `cas = "${PROBE}"`);
  const state = await appState(page);
  await menuPick(page, 'Bulk edit…');
  await page.waitForSelector('.u2-domain-bulk', {timeout: 20000});
  const box = page.locator('[data-u2-include="molecular_formula"]').first();
  const editor = page.locator('.u2-domain-bulk-row:has([data-u2-include="molecular_formula"]) input.u2-input-editor')
    .first();
  const before = await editor.isEnabled();
  // a real click at the element's own coordinates — dispatchEvent would prove nothing about U23
  await box.click();
  await page.waitForTimeout(500);
  const after = await editor.isEnabled();
  await editor.fill(value);
  await page.waitForTimeout(300);
  await shot(page, 'domain-stand-3-bulk');
  await page.locator('.u2-dialog button', {hasText: /^OK$/}).first().click();
  await page.waitForTimeout(4000);
  const landed = await page.evaluate(async (id) => {
    const row = await grok.dapi.domains.table('stockroom.substance').get(id);
    return String(row?.molecular_formula ?? '');
  }, probeId);
  ok('domain-stand/3 a real click on a bulk include checkbox enables its editor, and OK writes it',
    state.rows === 1 && before === false && after === true && landed === value,
    JSON.stringify({rows: state.rows, before, after, landed, expected: value}));
}

/** 4 — a multi-word entity value typed into the filter query box (U33/U34): the suggestion list
 * opens over the two words, ArrowDown + Enter applies the pick, and the popup that applied it is
 * not left covering the first rows of the fresh result. */
async function checkAutocomplete(page) {
  await openDomain(page, CONTAINER);
  const box = page.locator('[data-u2="filter-query-input"] input').first();
  await box.click();
  await box.fill('');
  await page.keyboard.type('location_id under "Main ca', {delay: 50});
  await page.waitForTimeout(2500);
  const options = await page.locator('.u2-suggest-popup .u2-suggest-option').count();
  const popup = await rectsOf(page, '.u2-suggest-popup', null);
  await shot(page, 'domain-stand-4-suggest');
  await page.keyboard.press('ArrowDown');
  await page.keyboard.press('Enter');
  await page.waitForTimeout(3000);
  const applied = await box.inputValue();
  const covering = await topRowHit(page);
  await shot(page, 'domain-stand-4-applied');
  ok('domain-stand/4 a two-word entity value completes, the pick applies, and the popup leaves the rows alone',
    options > 0 && popup !== null && inViewport(popup.rect, popup.view) &&
    applied !== 'location_id under "Main ca' && applied.includes('under') && covering === '',
    JSON.stringify({options, popup: popup && popup.rect, applied, covering}));
}

/** 5 — the zero-height virtualizer recovery (`VirtualRows._retry`): the view is hidden behind
 * another one and shown again, and the rows come back with a height. */
async function checkHiddenPane(page) {
  await openDomain(page, SUBSTANCE);
  const before = await appState(page);
  await page.evaluate(() => {
    const view = DG.View.create();
    view.name = 'e2e blank';
    window._u2e2e = {domain: grok.shell.v, blank: view};
    grok.shell.addView(view);
  });
  await page.waitForTimeout(2500);
  await page.evaluate(() => {
    grok.shell.v = window._u2e2e.domain;
  });
  await page.waitForTimeout(3000);
  const after = await appState(page);
  const row = await rectsOf(page, '[data-u2="domain-list"] [data-u2-row]', null);
  await shot(page, 'domain-stand-5-reshown');
  await page.evaluate(() => {
    window._u2e2e.blank.close();
  });
  ok('domain-stand/5 the list survives being hidden behind another view and shown again',
    before.rows > 0 && after.rows === before.rows && row !== null && visible(row.rect),
    JSON.stringify({before: before.rows, after: after.rows, row: row && row.rect}));
}

/** 6 — R-c end to end: in the trash, Restore STAGES (the row is marked, Save comes alive, nothing
 * is written), Ctrl+S lands it, and the row leaves the trash. */
async function checkStagedRestore(page) {
  await page.evaluate(async (id) => {
    await grok.dapi.domains.table('stockroom.substance').delete(id);
  }, probeId);
  await openDomain(page, SUBSTANCE, `cas = "${PROBE}"`, '&trash=1');
  const inTrash = await appState(page);
  await rowPick(page, 'Restore');
  await page.waitForTimeout(1500);
  const staged = await appState(page);
  // a staged restore writes NOTHING: the row is still the trash's until Save
  const untouched = await page.evaluate(async (id) => {
    const row = await grok.dapi.domains.table('stockroom.substance').get(id, {deleted: 'only'});
    return row !== null && row['~is_deleted'] === true;
  }, probeId);
  await shot(page, 'domain-stand-6-staged');
  await page.keyboard.press('Control+s');
  await page.waitForTimeout(6000);
  const saved = await appState(page);
  const live = await page.evaluate(async (id) => {
    const rows = await grok.dapi.domains.table('stockroom.substance').query({filter: `id = "${id}"`});
    return rows.length;
  }, probeId);
  await shot(page, 'domain-stand-6-saved');
  ok('domain-stand/6 Restore in the trash stages, Ctrl+S lands it, and the row leaves the trash',
    inTrash.rows === 1 && inTrash.deleted === 1 && staged.restored === 1 && staged.saveDisabled === false &&
    untouched && saved.rows === 0 && live === 1,
    JSON.stringify({inTrash, staged, untouched, saved, live}));
}

/* ------------------------------------------------------------------------------------ the run */

/** The way in. A worktree stand answers the dev-key login and not the login form (brief.md), so
 * that is tried first; a stand without dev keys falls through to the shared form login. */
async function openDomainStand() {
  const {browser, page} = await launch({viewport: VIEWPORT});
  await page.goto(`${STAND_URL}/login.html`, {waitUntil: 'load', timeout: 180000});
  await page.evaluate(() => fetch('/api/users/login/dev/admin', {method: 'POST'}).catch(() => {}));
  await page.goto(`${STAND_URL}/`, {waitUntil: 'load', timeout: 180000});
  try {
    await page.waitForFunction(() => {
      try {
        return !!(window.DG && DG.Func && grok.shell.user && grok.shell.user.login);
      } catch (e) {
        return false;
      }
    }, null, {timeout: 90000});
  } catch (e) {
    note('stand/login', 'the dev key did not answer; falling back to the login form');
    await browser.close();
    return openStand({viewport: VIEWPORT});
  }
  await page.waitForTimeout(4000);
  return {browser, page};
}

/** The `/domains` route is a Beta flag, and a fresh browser profile has it off. */
async function enableDomains(page) {
  await page.evaluate(() => {
    try {
      grok.shell.settings.enableDomainDatabases = true;
      grok.shell.settings.domainsDartUi = false;
    } catch (e) {
      // the settings proxy answered nothing useful; the localStorage copy below is the fallback
    }
    const stored = JSON.parse(window.localStorage.getItem('grok-settings') ?? '{}');
    stored.enableDomainDatabases = true;
    stored.domainsDartUi = false;
    window.localStorage.setItem('grok-settings', JSON.stringify(stored));
  });
  await page.reload({waitUntil: 'load', timeout: 180000});
  await page.waitForFunction(() => {
    try {
      return !!(window.DG && DG.Func && grok.shell.user);
    } catch (e) {
      return false;
    }
  }, null, {timeout: 180000});
  await page.waitForTimeout(3000);
}

const started = Date.now();
const {browser, page} = await openDomainStand();
try {
  await enableDomains(page);
  probeId = await page.evaluate(async ({cas, tag}) => {
    const [row] = await grok.dapi.domains.table('stockroom.substance')
      .insert([{name: `E2E probe ${tag}`, cas, molecular_formula: 'probe'}]);
    return String(row.id);
  }, {cas: PROBE, tag: TAG});
  note('stand/probe', `${PROBE} → ${probeId}`);
  for (const check of [checkGeometry, checkTrashHistory, checkBulkCheckbox, checkAutocomplete,
    checkHiddenPane, checkStagedRestore]) {
    try {
      await check(page);
    } catch (e) {
      ok(check.name, false, `check threw: ${String(e).slice(0, 300)}`);
      await shot(page, `domain-stand-${check.name}-failure`);
    }
  }
} finally {
  if (probeId !== null) {
    await page.evaluate(async (id) => {
      const table = grok.dapi.domains.table('stockroom.substance');
      try {
        await table.restore(id);
      } catch (e) {
        // already live — the delete below is what the probe is cleaned up by
      }
      await table.delete(id);
    }, probeId).catch(() => {});
  }
  await browser.close();
}

note('stand/target', `${STAND_URL} at ${VIEWPORT.width}×${VIEWPORT.height}`);
try {
  note('stand/code', execSync('git rev-parse HEAD', {encoding: 'utf8'}).trim());
} catch (e) {
  note('stand/code', 'unknown');
}
const code = report('domain.stand');
console.log(`${((Date.now() - started) / 1000).toFixed(1)}s · ${results.length} checks · artifacts in ${ARTIFACTS}`);
process.exit(code);
