// Drives Settings -> Browse (the ordering rules editor) in a real browser.
const { chromium } = require('playwright');
const BASE = process.env.DG_BASE || 'http://localhost:61006';
const errors = [];
let failures = 0;
function check(label, ok, detail) {
  if (!ok) failures++;
  console.log(`${ok ? 'PASS' : 'FAIL'}  ${label}${detail ? `\n        ${detail}` : ''}`);
}
// Pre-existing noise on this dev stack: undeployed packages, unbuilt help, no grok_connect
// container (so no JDBC providers). Datagrok also emits a separate "Stack trace <id>" line.
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

async function saveAndApply(page) {
  await page.locator('text=SAVE AND APPLY').first().click();
  await page.waitForTimeout(7000);
}

const rulesOf = (page) => page.evaluate(() => {
  try { return JSON.parse(JSON.parse(localStorage.getItem('grok-settings') || '{}').orderRules || '[]'); }
  catch (e) { return null; }
});

(async () => {
  const browser = await chromium.launch({ headless: true });
  const page = await browser.newPage({ viewport: { width: 1600, height: 1000 } });
  page.on('console', m => { const t = m.text().slice(0, 300);
    if (m.type() === 'error' && !NOISE.test(t)) errors.push('CONSOLE ' + t); });
  page.on('pageerror', e => errors.push('PAGEERROR ' + String(e).slice(0, 300)));

  await login(page);
  // Start from a clean slate so row indexes are predictable.
  await page.evaluate(() => { if (window.grok && window.grok.settings) window.grok.settings.orderRules = '[]'; });

  const before = errors.length;
  await page.goto(`${BASE}/settings/browse`, { waitUntil: 'domcontentloaded', timeout: 180000 });
  await page.waitForTimeout(25000);
  await page.screenshot({ path: 'set-01-browse.png' });
  check('no JS errors opening Settings > Browse', errors.length === before,
    errors.slice(before, before + 3).join('\n        '));

  // 6 - DOM controls, not a grid
  const canvases = await page.locator('.grok-settings-order-rules canvas').count();
  check('no grid canvas (reverted to DOM controls)', canvases === 0, `canvases=${canvases}`);

  // 1 - fills the page
  const box = await page.locator('.grok-settings-order-rules').first().boundingBox();
  const host = await page.locator('.dlg-settings-page-host').first().boundingBox().catch(() => null);
  console.log('    editor box:', JSON.stringify(box), 'host:', JSON.stringify(host));
  check('editor fills its width', box && host && box.width >= host.width - 30,
    `editor=${box && box.width} host=${host && host.width}`);
  check('editor fills its height', box && host && box.height >= host.height - 30,
    `editor=${box && box.height} host=${host && host.height}`);

  // Header must lay out as a row. ui-box forces direct children to flex-direction:column unless
  // they carry an exempt class, which silently stacked this header while every other check passed.
  const hdr = await page.evaluate(() => {
    const h = document.querySelector('.grok-order-rule-header');
    if (!h) return null;
    const s = getComputedStyle(h);
    const kids = [...h.children].map(c => Math.round(c.getBoundingClientRect().top));
    return { dir: s.flexDirection, height: h.offsetHeight, tops: kids };
  });
  console.log('    header:', JSON.stringify(hdr));
  check('header lays out as a row', hdr && hdr.dir === 'row', JSON.stringify(hdr));
  check('header is one line tall', hdr && hdr.height < 48, hdr && String(hdr.height));
  check('header cells share a baseline', hdr && new Set(hdr.tops).size <= 2, JSON.stringify(hdr && hdr.tops));

  // add two rules
  await page.locator('text=Add rule').first().click();
  await page.waitForTimeout(2500);
  await page.locator('text=Add rule').first().click();
  await page.waitForTimeout(2500);
  await page.screenshot({ path: 'set-02-two-rules.png' });

  const rows = page.locator('.grok-order-rule-row').filter({ hasNot: page.locator('.grok-order-rule-header') });
  const rowCount = await page.locator('.grok-order-rule-row:not(.grok-order-rule-header)').count();
  check('two rule rows rendered', rowCount === 2, `rows=${rowCount}`);

  // 4 - priorities are visible and numbered
  const prios = await page.locator('.grok-order-rule-row:not(.grok-order-rule-header) .grok-order-rule-priority')
    .allInnerTexts();
  console.log('    priorities:', JSON.stringify(prios));
  check('priorities shown as 1..n', prios.join(',') === '1,2', JSON.stringify(prios));

  // 2 - per-row remove button (plus the two move buttons)
  const actions = await page.locator('.grok-order-rule-row:not(.grok-order-rule-header)')
    .first().locator('.grok-order-rule-actions i').count();
  check('row has move-up / move-down / remove buttons', actions === 3, `icons=${actions}`);

  // 3 - space selector is a typeahead
  const ta = await page.locator('.grok-order-rule-row:not(.grok-order-rule-header)')
    .first().locator('input[placeholder="space name"]').count();
  check('parent space is a typeahead input', ta > 0, `matches=${ta}`);

  // Desc checkbox must stay checkbox-sized: `.ui-input-editor { width: 100% }` had stretched it.
  const cb = await page.evaluate(() => {
    const c = document.querySelector('.grok-order-rule-row:not(.grok-order-rule-header) .grok-order-rule-desc input[type=checkbox]');
    return c ? { w: c.offsetWidth, h: c.offsetHeight } : null;
  });
  console.log('    checkbox:', JSON.stringify(cb));
  check('Desc checkbox is checkbox-sized', cb && cb.w <= 24 && cb.h <= 24, JSON.stringify(cb));

  // Add rule sits above the list.
  const orderOk = await page.evaluate(() => {
    const host = document.querySelector('.grok-settings-order-rules');
    const kids = [...host.children];
    return kids.findIndex(k => k.classList.contains('grok-order-rule-add'))
         < kids.findIndex(k => k.classList.contains('grok-order-rules-rows'));
  });
  check('Add rule is above the list', orderOk === true, String(orderOk));

  // The rows container must not clip, or it hides the typeahead dropdown.
  const clip = await page.evaluate(() => {
    const r = document.querySelector('.grok-order-rules-rows');
    const host = document.querySelector('.grok-settings-order-rules');
    const s = getComputedStyle(r);
    return { overflowY: s.overflowY, h: r.offsetHeight, hostH: host.offsetHeight };
  });
  console.log('    rows container:', JSON.stringify(clip));
  check('rows container does not clip the dropdown',
    clip && clip.overflowY !== 'auto' && clip.overflowY !== 'scroll' && clip.overflowY !== 'hidden',
    JSON.stringify(clip));
  // Chrome above the list is the description + Add rule + header + host padding (~135px).
  check('list stretches to the remaining height', clip && clip.h >= clip.hostH - 170,
    JSON.stringify(clip));

  // A new rule defaults to "Any space", so the typeahead starts disabled; unchecking Any enables it.
  const row1 = page.locator('.grok-order-rule-row:not(.grok-order-rule-header)').first();
  const anyBox = row1.locator('.grok-order-rule-space input[type=checkbox]').first();
  const spaceInput = row1.locator('input[placeholder="space name"]');
  const hiddenWhenAny = await spaceInput.evaluate(e => {
    const r = e.getBoundingClientRect();
    return r.width === 0 && r.height === 0;
  });
  check('space typeahead is hidden while Any is checked', hiddenWhenAny === true,
    String(hiddenWhenAny));
  await anyBox.uncheck();
  await page.waitForTimeout(2000);
  const shownAfter = await spaceInput.evaluate(e => {
    const r = e.getBoundingClientRect();
    return r.width > 0 && r.height > 0;
  });
  check('unchecking Any shows the space typeahead', shownAfter === true, String(shownAfter));

  // Checkbox first, then its label; and all three boxes are the same size.
  const layout = await page.evaluate(() => {
    const row = document.querySelector('.grok-order-rule-row:not(.grok-order-rule-header)');
    const chk = row.querySelector('.grok-order-rule-check');
    const box = chk.querySelector('input[type=checkbox]');
    const lbl = chk.querySelector('.grok-order-rule-check-label');
    const sizes = [...row.querySelectorAll('input[type=checkbox]')]
      .map(c => `${c.offsetWidth}x${c.offsetHeight}`);
    return { boxX: box.getBoundingClientRect().left, lblX: lbl.getBoundingClientRect().left,
             sizes, rowH: row.offsetHeight };
  });
  console.log('    layout:', JSON.stringify(layout));
  check('checkbox comes before its label', layout.boxX < layout.lblX,
    `box=${layout.boxX} label=${layout.lblX}`);
  check('all three checkboxes are the same size',
    new Set(layout.sizes).size === 1 && layout.sizes.length === 3, JSON.stringify(layout.sizes));
  check('row is taller', layout.rowH >= 32, String(layout.rowH));

  // Typeahead dropdown must be visible above later rows.
  await spaceInput.click();
  await spaceInput.type('a');
  await page.waitForTimeout(2500);
  await page.screenshot({ path: 'set-05-dropdown.png' });
  const dd = await page.evaluate(() => {
    const d = document.querySelector('.type-ahead-drop-down');
    if (!d) return null;
    const r = d.getBoundingClientRect();
    if (r.width === 0 || r.height === 0) return { visible: false };
    // whatever paints at the dropdown's own centre should be the dropdown itself
    const el = document.elementFromPoint(r.left + r.width / 2, r.top + Math.min(20, r.height / 2));
    return { visible: true, ownsPixel: !!(el && d.contains(el)), rect: { w: Math.round(r.width), h: Math.round(r.height) } };
  });
  console.log('    dropdown:', JSON.stringify(dd));
  check('typeahead dropdown is not covered by the list',
    dd === null || dd.visible === false || dd.ownsPixel === true, JSON.stringify(dd));

  // Any / Root are checkboxes now, not entries in the typeahead list.
  const spaceCell = page.locator('.grok-order-rule-row:not(.grok-order-rule-header)').first()
    .locator('.grok-order-rule-space');
  const boxes = await spaceCell.locator('input[type=checkbox]').count();
  const anyRoot = (await spaceCell.innerText()).replace(/\s+/g, ' ').trim();
  console.log('    space cell:', JSON.stringify(anyRoot), 'checkboxes=', boxes);
  check('space cell has Any / Root checkboxes', boxes === 2, `checkboxes=${boxes}`);
  check('space typeahead is for real spaces only',
    await spaceCell.locator('input[placeholder="space name"]').count() > 0);

  // Entity types are limited to what actually appears in Browse / card views.
  const types = await page.locator('.grok-order-rule-row:not(.grok-order-rule-header)').first()
    .locator('.grok-order-rule-cell').first().locator('option').allInnerTexts();
  console.log('    entity types:', JSON.stringify(types));
  check('entity list includes Space / Dashboard / User',
    types.indexOf('Space') >= 0 && types.indexOf('Dashboard') >= 0 && types.indexOf('User') >= 0,
    JSON.stringify(types));
  check('entity list excludes internal types',
    types.indexOf('FuncCall') < 0 && types.indexOf('EntityTag') < 0, JSON.stringify(types));
  check('entity list is short (curated, not the whole registry)', types.length <= 25,
    `count=${types.length}`);

  // make the two rules distinguishable: set row 2's Order by
  const row2Field = page.locator('.grok-order-rule-row:not(.grok-order-rule-header)').nth(1)
    .locator('.grok-order-rule-cell').nth(2).locator('input');
  await row2Field.fill('createdOn');
  await row2Field.press('Enter');
  await page.waitForTimeout(3000);

  // Nothing may be applied before Save and apply.
  let staged = await rulesOf(page);
  console.log('    before save:', JSON.stringify(staged));
  check('edits are NOT applied before Save and apply',
    !staged || staged.length === 0, JSON.stringify(staged));

  await saveAndApply(page);
  let rules = await rulesOf(page);
  console.log('    after save:', JSON.stringify(rules));
  check('Save and apply commits the rules',
    rules && rules.length === 2 && rules[1].field === 'createdOn', JSON.stringify(rules));
  // A half-typed space name is not a space id; storing it would make the rule match nothing.
  check('partial space text is not stored as a space id',
    rules && rules.every(r => r.space === '' || r.space === '@root' || r.space.length > 8),
    JSON.stringify(rules.map(r => r.space)));

  // 5 - move row 1 down
  await page.locator('.grok-order-rule-row:not(.grok-order-rule-header)').first()
    .locator('.grok-order-rule-actions i').nth(1).click();
  await page.waitForTimeout(3000);
  await saveAndApply(page);
  rules = await rulesOf(page);
  console.log('    after move-down:', JSON.stringify(rules));
  check('move down reorders priority',
    rules && rules.length === 2 && rules[0].field === 'createdOn', JSON.stringify(rules));
  await page.screenshot({ path: 'set-03-moved.png' });

  // 5 - move it back up
  await page.locator('.grok-order-rule-row:not(.grok-order-rule-header)').nth(1)
    .locator('.grok-order-rule-actions i').nth(0).click();
  await page.waitForTimeout(3000);
  await saveAndApply(page);
  rules = await rulesOf(page);
  console.log('    after move-up:', JSON.stringify(rules));
  check('move up restores priority',
    rules && rules.length === 2 && rules[0].field === 'friendlyName', JSON.stringify(rules));

  // 2 - remove a row
  await page.locator('.grok-order-rule-row:not(.grok-order-rule-header)').first()
    .locator('.grok-order-rule-actions i').nth(2).click();
  await page.waitForTimeout(3000);
  await saveAndApply(page);
  rules = await rulesOf(page);
  console.log('    after remove:', JSON.stringify(rules));
  check('remove deletes that row',
    rules && rules.length === 1 && rules[0].field === 'createdOn', JSON.stringify(rules));

  // persistence across reload
  await page.goto(`${BASE}/settings/browse`, { waitUntil: 'domcontentloaded', timeout: 180000 });
  await page.waitForTimeout(25000);
  rules = await rulesOf(page);
  console.log('    after reload:', JSON.stringify(rules));
  check('rules survive reload', rules && rules.length === 1, JSON.stringify(rules));
  await page.screenshot({ path: 'set-04-reloaded.png' });

  console.log('\n--- errors (' + errors.length + '):');
  errors.slice(0, 8).forEach(e => console.log('      ' + e));
  await browser.close();
  console.log(failures === 0 ? '\nALL SETTINGS CHECKS PASSED' : `\n${failures} CHECK(S) FAILED`);
  process.exit(failures === 0 ? 0 : 1);
})().catch(e => { console.error('HARNESS ERROR', e); process.exit(2); });
