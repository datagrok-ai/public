// Ordering context menu: the "Order by" section on the Spaces node, on a space, and in a
// DataSourceCardView. Checks the field list (model props + schema props + meta params), the
// direction items, and the EditSettings-gated "for all users" icon.
const { chromium } = require('playwright');
const BASE = process.env.DG_BASE || 'http://localhost:61006';
const API = process.env.DG_API || 'http://localhost:61002';
let failures = 0;
const errors = [];
function check(label, ok, detail) {
  if (!ok) failures++;
  console.log(`${ok ? 'PASS' : 'FAIL'}  ${label}${detail ? `\n        ${detail}` : ''}`);
}
const NOISE = /Failed to load resource|packages\/published|detectors\.js|package\.js|MIME type|401 \(Unauthorized\)|favicon|Provider \w+ not found|ConnectorsService|^Stack trace \w+|Translating stack trace/i;

async function login(page) {
  await page.goto(`${BASE}/login.html`, { waitUntil: 'domcontentloaded', timeout: 180000 });
  const u = page.locator('input[placeholder="Login or Email"]:visible').first();
  await u.waitFor({ state: 'visible', timeout: 300000 });
  await u.fill('admin');
  await page.locator('input[placeholder="Password"]:visible').first().fill('admin');
  await page.locator('button:has-text("LOGIN"), .ui-btn:has-text("LOGIN")').first().click();
  await page.locator('.d4-tree-view-group-label', { hasText: /^Spaces$/ }).first().waitFor({ timeout: 300000 });
  await page.waitForTimeout(8000);
}

// Menu items are rendered into the document body; read their labels.
const menuLabels = (page) => page.evaluate(() =>
  [...document.querySelectorAll('.d4-menu-item-label')].map(e => e.innerText.trim()).filter(Boolean));

// Escape does not dismiss a d4 context menu, and clicking a "neutral" corner lands on the left
// sidebar, whose menu groups expand on hover and then swallow every later click. Each right-click
// builds a fresh Menu, so dropping the old popup is safe and is the only reliable close.
async function closeMenus(page) {
  await page.keyboard.press('Escape');
  await page.evaluate(() =>
    document.querySelectorAll('.d4-menu-popup').forEach(e => e.remove()));
  await page.mouse.move(800, 600);
  await page.waitForTimeout(800);
}

// Hovers "Order by" so its submenu renders, then returns every label on screen.
async function openOrderBy(page) {
  const item = page.locator('.d4-menu-item-label', { hasText: /^Order by$/ }).first();
  if (await item.count() === 0) return null;
  await item.hover();
  await page.waitForTimeout(2500);
  return await menuLabels(page);
}

const APPLY_ALL = 'Apply for all users';

(async () => {
  const browser = await chromium.launch({ headless: true });
  const ctx = await browser.newContext({ viewport: { width: 1600, height: 1000 } });
  const page = await ctx.newPage();
  page.on('console', m => { const t = m.text().slice(0, 300);
    if (m.type() === 'error' && !NOISE.test(t)) errors.push('CONSOLE ' + t); });
  page.on('pageerror', e => errors.push('PAGEERROR ' + String(e).slice(0, 300)));

  await login(page);
  check('client compiled and app booted', true);

  // ---------- 1. Spaces node
  const spacesNode = page.locator('.d4-tree-view-group-label', { hasText: /^Spaces$/ }).first();
  await spacesNode.click({ button: 'right' });
  await page.waitForTimeout(2500);
  let labels = await menuLabels(page);
  console.log('    spaces menu:', JSON.stringify(labels.slice(0, 12)));
  check('Spaces node context menu has "Order by"', labels.includes('Order by'), JSON.stringify(labels));

  const spacesOrder = await openOrderBy(page);
  console.log('    Order by (Spaces):', JSON.stringify((spacesOrder || []).slice(0, 25)));
  check('field list offers model properties',
    spacesOrder && spacesOrder.includes('friendlyName'), JSON.stringify(spacesOrder));
  check('direction items present',
    spacesOrder && spacesOrder.includes('Ascending') && spacesOrder.includes('Descending'),
    JSON.stringify(spacesOrder));
  check('"' + APPLY_ALL + '" shown for an EditSettings holder',
    spacesOrder && spacesOrder.includes(APPLY_ALL), JSON.stringify(spacesOrder));
  const metaAt = (spacesOrder || []).indexOf('Meta parameters');
  const plainAt = (spacesOrder || []).indexOf('friendlyName');
  console.log('    Meta group at', metaAt, '/ first model field at', plainAt);
  check('"Meta parameters" group comes before the model fields',
    metaAt >= 0 && plainAt >= 0 && metaAt < plainAt,
    JSON.stringify((spacesOrder || []).slice(0, 14)));
  await page.screenshot({ path: 'om-01-spaces.png' });
  await closeMenus(page);

  // ---------- 2. a space
  await page.evaluate(() => {
    const tri = document.querySelector('[name="tree-expander-Spaces"]');
    if (tri) tri.click();
  });
  await page.waitForTimeout(6000);
  const firstSpace = page.locator('.d4-tree-view-group-label, .d4-tree-view-item-label')
    .filter({ hasText: /^(alpha|middle|ordered|plain|zebra)$/ }).first();
  const haveSpace = await firstSpace.count() > 0;
  check('a space node is present to right-click', haveSpace);
  if (haveSpace) {
    await firstSpace.click({ button: 'right' });
    await page.waitForTimeout(2500);
    labels = await menuLabels(page);
    check('space context menu has "Order by"', labels.includes('Order by'), JSON.stringify(labels.slice(0, 20)));
    check('the old "Order Children By..." command is gone',
      !labels.some(l => /Order Children By/i.test(l)), JSON.stringify(labels.slice(0, 20)));
    const spaceOrder = await openOrderBy(page);
    console.log('    Order by (space):', JSON.stringify((spaceOrder || []).slice(0, 25)));
    check('space ordering offers fields',
      spaceOrder && spaceOrder.includes('friendlyName'), JSON.stringify(spaceOrder));
    // the tree can show every bucket under a space, so the menu covers all their fields
    check('space menu in Browse offers other types’ fields too',
      spaceOrder && spaceOrder.length > 12, JSON.stringify((spaceOrder || []).slice(0, 20)));
    await page.screenshot({ path: 'om-02-space.png' });
    await closeMenus(page);
  }

  // ---------- 3. a DataSourceCardView
  await page.goto(`${BASE}/connections`, { waitUntil: 'domcontentloaded', timeout: 180000 });
  await page.waitForTimeout(12000);
  // The handler is on the gallery grid itself, not the view root.
  const grid = page.locator('.grok-gallery-grid').first();
  const haveGrid = await grid.count() > 0;
  check('DataSourceCardView rendered', haveGrid,
    haveGrid ? '' : await page.evaluate(() => document.body.innerText.slice(0, 200)));
  if (haveGrid) {
    // A real right-click lands on a card, whose own handler swallows it; dispatch at the grid.
    await page.evaluate(() => {
      const g = document.querySelector('.grok-gallery-grid');
      const r = g.getBoundingClientRect();
      g.dispatchEvent(new MouseEvent('contextmenu', { bubbles: true, cancelable: true,
        view: window, clientX: r.left + 30, clientY: r.bottom - 30, button: 2 }));
    });
    await page.waitForTimeout(3500);
    labels = await menuLabels(page);
    check('card view context menu has "Order by"', labels.includes('Order by'), JSON.stringify(labels.slice(0, 20)));
    const viewOrder = await openOrderBy(page);
    console.log('    Order by (card view):', JSON.stringify((viewOrder || []).slice(0, 25)));
    // the view's own sortableBy labels win over raw field names, so this reads "Name", not
    // "friendlyName"; a field the view has no label for (updatedOn) still shows raw
    check('card view ordering offers fields under the view own labels',
      viewOrder && viewOrder.includes('Name') && viewOrder.includes('Author') &&
      viewOrder.includes('Created'), JSON.stringify(viewOrder));
    check('card view ordering also offers fields the view does not label',
      viewOrder && viewOrder.includes('updatedOn'), JSON.stringify(viewOrder));
    check('"' + APPLY_ALL + '" shown in the card view too',
      viewOrder && viewOrder.includes(APPLY_ALL), JSON.stringify(viewOrder));
    await page.screenshot({ path: 'om-03-cardview.png' });
    await closeMenus(page);
  }

  // ---------- 4. picking a field actually persists a rule and reorders the tree
  await page.goto(`${BASE}/login.html`, { waitUntil: 'domcontentloaded', timeout: 180000 });
  await page.locator('.d4-tree-view-group-label', { hasText: /^Spaces$/ }).first().waitFor({ timeout: 300000 });
  await page.waitForTimeout(8000);
  await page.evaluate(() => { window.grok.settings.orderRules = '[]'; });

  const treeOrder = () => page.evaluate(() =>
    [...document.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-item-label')]
      .map(e => e.innerText.trim())
      .filter(l => ['alpha', 'middle', 'ordered', 'plain', 'zebra'].includes(l)));

  await page.evaluate(() => {
    const tri = document.querySelector('[name="tree-expander-Spaces"]');
    if (tri) tri.click();
  });
  await page.waitForTimeout(8000);
  const before = await treeOrder();
  console.log('    tree before:', JSON.stringify(before));

  await closeMenus(page);
  await page.locator('.d4-tree-view-group-label', { hasText: /^Spaces$/ }).first().click({ button: 'right' });
  await page.waitForTimeout(2500);
  await openOrderBy(page);
  await page.locator('.d4-menu-item-label', { hasText: /^createdOn$/ }).first().click();
  await page.waitForTimeout(6000);

  const rules = await page.evaluate(() => window.grok.settings.orderRules);
  console.log('    orderRules after picking createdOn:', JSON.stringify(rules));
  check('picking a field writes the personal rule',
    rules && rules.indexOf('createdOn') >= 0 && rules.indexOf('@root') >= 0, JSON.stringify(rules));

  await page.evaluate(() => {
    const tri = document.querySelector('[name="tree-expander-Spaces"]');
    if (tri) tri.click();
  });
  await page.waitForTimeout(8000);
  const after = await treeOrder();
  console.log('    tree after:', JSON.stringify(after));
  check('the tree reorders by the picked field',
    after.length === before.length && after.length > 1 && after.join(',') !== before.join(','),
    `before=${JSON.stringify(before)} after=${JSON.stringify(after)}`);

  // ---------- 5. the users toggle writes the All-Users override
  const overrides = () => page.evaluate(async (api) => {
    const r = await fetch(`${api}/admin/settings/override`, { credentials: 'include' });
    return r.ok ? await r.json() : { __status: r.status };
  }, API);

  await closeMenus(page);
  await closeMenus(page);
  await page.locator('.d4-tree-view-group-label', { hasText: /^Spaces$/ }).first().click({ button: 'right' });
  await page.waitForTimeout(2500);
  await openOrderBy(page);
  await page.locator('.d4-menu-item-label', { hasText: /^updatedOn$/ }).first().click();
  await page.waitForTimeout(3000);

  await closeMenus(page);
  await page.locator('.d4-tree-view-group-label', { hasText: /^Spaces$/ }).first().click({ button: 'right' });
  await page.waitForTimeout(2500);
  await openOrderBy(page);
  const applyItem = page.locator('.d4-menu-item-label', { hasText: new RegExp('^' + APPLY_ALL + '$') }).first();
  check('the "' + APPLY_ALL + '" item is present', await applyItem.count() > 0);
  await applyItem.click();
  await page.waitForTimeout(8000);

  const ov = await overrides();
  const stored = ov && ov.orderRules;
  console.log('    All-Users override:', JSON.stringify(stored).slice(0, 250));
  check('"' + APPLY_ALL + '" writes the All-Users rule',
    stored && JSON.stringify(stored).indexOf('updatedOn') >= 0, JSON.stringify(ov).slice(0, 250));

  // cleanup so a later run starts clean
  await page.evaluate(async (api) => {
    await fetch(`${api}/admin/settings/override`, { method: 'DELETE', credentials: 'include' });
  }, API);
  await page.evaluate(() => { window.grok.settings.orderRules = '[]'; });

  // ---------- 6. a space's own override drives the space's view, not just the tree
  await closeMenus(page);
  const spaceName = 'plain';
  const spaceId = await page.evaluate(async ([api, name]) => {
    const r = await fetch(`${api}/projects/namespaces?include=metaParams`, { credentials: 'include' });
    const l = r.ok ? await r.json() : [];
    const p = (Array.isArray(l) ? l : []).find(x => (x.friendlyName || x.name) === name);
    return p ? p.id : null;
  }, [API, spaceName]);
  check('found the space to scope the view to', !!spaceId, String(spaceId));

  // Order its children by id — a field no view exposes as a sortableBy label.
  const setOverride = (field) => page.evaluate(async ([api, id, f]) => {
    const r = await fetch(`${api}/projects/${id}`, { credentials: 'include' });
    const p = await r.json();
    p.metaParams = Object.assign({}, p.metaParams, { orderBy: f, orderDescending: 'false' });
    await fetch(`${api}/projects`, { method: 'POST', credentials: 'include',
      headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(p) });
  }, [API, spaceId, field]);

  const viewOrderOf = async () => {
    await page.goto(`${BASE}/s/${spaceName}`, { waitUntil: 'domcontentloaded', timeout: 180000 });
    await page.waitForTimeout(15000);
    return await page.evaluate(() =>
      [...document.querySelectorAll('.grok-gallery-grid .d4-link-label label')]
        .map(e => e.innerText.trim()).filter(Boolean).slice(0, 10));
  };

  await setOverride('friendlyName');
  const viewByName = await viewOrderOf();
  console.log('    space view ordered by friendlyName:', JSON.stringify(viewByName));
  await setOverride('id');
  const viewById = await viewOrderOf();
  console.log('    space view ordered by id:', JSON.stringify(viewById));
  check('the space view lists its contents', viewByName.length > 1, JSON.stringify(viewByName));
  check('the space\'s own ordering override changes the view',
    viewByName.length > 1 && viewById.join(',') !== viewByName.join(','),
    `byName=${JSON.stringify(viewByName)} byId=${JSON.stringify(viewById)}`);

  // ---------- 7. a space's override orders everything inside it.
  // Driven through the menu: writing metaParams over REST silently no-ops, because `metaParams`
  // is a read projection and the write shape is `entityMetaParams`.
  const ALL = 'plain';
  await closeMenus(page);
  await page.goto(BASE + '/login.html', { waitUntil: 'domcontentloaded', timeout: 180000 });
  await page.locator('.d4-tree-view-group-label', { hasText: /^Spaces$/ }).first().waitFor({ timeout: 300000 });
  await page.waitForTimeout(8000);
  await page.evaluate(() => { window.grok.settings.orderRules = '[]'; });

  const expand = (path) => page.evaluate((v) => {
    const t = document.querySelector('[name="tree-expander-' + v + '"]');
    if (t) { t.click(); return true; }
    return false;
  }, path);
  const childrenOf = (nm) => page.evaluate((s2) => {
    const all = [...document.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-item-label')]
      .map(e => e.innerText.trim()).filter(Boolean);
    const i = all.indexOf(s2);
    return i < 0 ? [] : all.slice(i + 1, i + 5);
  }, nm);
  const spaceMeta = () => page.evaluate(async (args) => {
    const r = await fetch(args[0] + '/projects/namespaces?include=metaParams', { credentials: 'include' });
    const l = r.ok ? await r.json() : [];
    const a2 = (Array.isArray(l) ? l : []).find(x => (x.friendlyName || x.name) === args[1]);
    return a2 ? a2.metaParams : null;
  }, [API, ALL]);
  // dispatch, not click: after a reload the expanded subtree covers the node
  const pickOn = async (nm, label) => {
    await closeMenus(page);
    const found = await page.evaluate((n) => {
      const el = [...document.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-item-label')]
        .find(e => e.innerText.trim() === n);
      if (!el) return false;
      const r = el.getBoundingClientRect();
      el.dispatchEvent(new MouseEvent('contextmenu', { bubbles: true, cancelable: true,
        view: window, clientX: r.left + 5, clientY: r.top + 5, button: 2 }));
      return true;
    }, nm);
    if (!found) return false;
    await page.waitForTimeout(2500);
    await page.locator('.d4-menu-item-label', { hasText: /^Order by$/ }).first().hover();
    await page.waitForTimeout(2500);
    const it = page.locator('.d4-menu-item-label', { hasText: new RegExp('^' + label + '$') }).first();
    if (await it.count() === 0) return false;
    await it.click({ timeout: 20000 });
    await page.waitForTimeout(7000);
    return true;
  };
  const snapshot = async () => {
    await page.reload({ waitUntil: 'domcontentloaded', timeout: 180000 });
    await page.locator('.d4-tree-view-group-label', { hasText: /^Spaces$/ }).first().waitFor({ timeout: 300000 });
    await page.waitForTimeout(9000);
    await expand('Spaces'); await page.waitForTimeout(6000);
    await expand('Spaces---' + ALL); await page.waitForTimeout(8000);
    return (await childrenOf(ALL)).join(',');
  };

  await expand('Spaces'); await page.waitForTimeout(6000);
  await expand('Spaces---' + ALL); await page.waitForTimeout(8000);
  const mixed = await page.evaluate(async (args) => {
    const r = await fetch(args[0] + '/projects/namespaces', { credentials: 'include' });
    const l = r.ok ? await r.json() : [];
    const sp = (Array.isArray(l) ? l : []).find(x => (x.friendlyName || x.name) === args[1]);
    if (!sp) return [];
    const c = await fetch(args[0] + '/spaces/' + sp.id + '/children?includeLinked=true', { credentials: 'include' });
    const cs = c.ok ? await c.json() : [];
    return [...new Set((Array.isArray(cs) ? cs : []).map(x => x['#type'] || '?'))];
  }, [API, ALL]);
  console.log('    types inside ' + ALL + ':', JSON.stringify(mixed));
  check('the space holds more than one entity type', mixed.length > 1, JSON.stringify(mixed));

  check('picked an order field on the space', await pickOn(ALL, 'friendlyName'));
  const metaAsc = await spaceMeta();
  check('the space override was stored', metaAsc && metaAsc.orderBy === 'friendlyName', JSON.stringify(metaAsc));
  const ascOrder = await snapshot();
  check('picked Descending on the space', await pickOn(ALL, 'Descending'));
  const metaDesc = await spaceMeta();
  check('the direction was stored', metaDesc && metaDesc.orderDescending === 'true', JSON.stringify(metaDesc));
  const descOrder = await snapshot();
  console.log('    asc :', ascOrder);
  console.log('    desc:', descOrder);
  check('the space override reorders its children',
    ascOrder.length > 0 && ascOrder !== descOrder, 'asc=' + ascOrder + ' desc=' + descOrder);

  // "Default" is the only way back out of an override
  check('picked Default on the space', await pickOn(ALL, 'Default'));
  const metaCleared = await spaceMeta();
  console.log('    meta after Default:', JSON.stringify(metaCleared));
  check('Default clears the space override',
    metaCleared && metaCleared.orderBy == null && metaCleared.orderDescending == null,
    JSON.stringify(metaCleared));
  const clearedOrder = await snapshot();
  console.log('    cleared:', clearedOrder);
  check('cleared space falls back to the name order',
    clearedOrder === ascOrder, 'cleared=' + clearedOrder + ' expected=' + ascOrder);

  // ---------- 8. the space view has its own ordering button carrying the same section
  await closeMenus(page);
  await page.goto(BASE + '/s/plain', { waitUntil: 'domcontentloaded', timeout: 180000 });
  await page.waitForTimeout(16000);
  const btnCount = await page.locator('i[name*="icon-sort"], i[class*="fa-sort"]').count();
  check('SpaceView shows an ordering button', btnCount > 0, 'found=' + btnCount);

  // the listener sits on the icon element itself (sortButton.root), not on any wrapper
  const opened = await page.evaluate(() => {
    const i = document.querySelector('i[name*="icon-sort-alt"], i[class*="fa-sort-alt"]');
    if (!i) return false;
    i.click();
    return true;
  });
  await page.waitForTimeout(5000);
  const btnMenu = await menuLabels(page);
  console.log('    SpaceView order button menu:', JSON.stringify(btnMenu.slice(0, 22)));
  check('the ordering button opens a menu', opened && btnMenu.length > 0,
    JSON.stringify(btnMenu.slice(0, 10)));
  check('button menu uses Ascending/Descending, not high/low',
    btnMenu.includes('Ascending') && btnMenu.includes('Descending') &&
    !btnMenu.some(l => /high to low|low to high/i.test(l)), JSON.stringify(btnMenu.slice(0, 22)));
  check('button menu carries "' + APPLY_ALL + '"',
    btnMenu.includes(APPLY_ALL), JSON.stringify(btnMenu.slice(0, 22)));
  await closeMenus(page);

  console.log('\n--- errors (' + errors.length + '):');
  errors.slice(0, 10).forEach(e => console.log('      ' + e));
  check('no unexpected console errors', errors.length === 0, errors.slice(0, 3).join('\n        '));

  await browser.close();
  console.log(failures === 0 ? '\nALL ORDER-MENU CHECKS PASSED' : `\n${failures} CHECK(S) FAILED`);
  process.exit(failures === 0 ? 0 : 1);
})().catch(e => { console.error('HARNESS ERROR', e); process.exit(2); });
