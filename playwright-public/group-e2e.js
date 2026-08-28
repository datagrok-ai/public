// Group / All-Users defaults for the ordering rules: set them in the scope selector as admin,
// confirm the server stored the override, then confirm a fresh session inherits them and the
// Browse tree actually orders by them.
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

// Reads the stored All-Users override straight from the server, through the page's session.
const overrides = (page) => page.evaluate(async (api) => {
  const r = await fetch(`${api}/admin/settings/override`, { credentials: 'include' });
  return r.ok ? await r.json() : { __status: r.status };
}, API);

(async () => {
  const browser = await chromium.launch({ headless: true });

  // ---------- 1. admin sets the All-Users default through the scope selector
  const ctx1 = await browser.newContext({ viewport: { width: 1600, height: 1000 } });
  const page = await ctx1.newPage();
  page.on('console', m => { const t = m.text().slice(0, 300);
    if (m.type() === 'error' && !NOISE.test(t)) errors.push('CONSOLE ' + t); });
  page.on('pageerror', e => errors.push('PAGEERROR ' + String(e).slice(0, 300)));

  await login(page);
  await page.evaluate(() => { if (window.grok) window.grok.settings.orderRules = '[]'; });
  await page.goto(`${BASE}/settings/browse`, { waitUntil: 'domcontentloaded', timeout: 180000 });
  await page.waitForTimeout(25000);

  const scopeBtns = page.locator('.grok-settings-scope-btn');
  const nBtns = await scopeBtns.count();
  check('scope selector present (personal / group)', nBtns >= 2, `buttons=${nBtns}`);
  if (nBtns < 2) { await browser.close(); process.exit(1); }

  const errBefore = errors.length;
  await scopeBtns.nth(1).click();           // switch to group / shared defaults
  await page.waitForTimeout(9000);
  await page.screenshot({ path: 'grp-01-scope.png' });
  const scopeName = await page.locator('.grok-settings-scope-name').first().innerText().catch(() => '');
  console.log('    scope:', JSON.stringify(scopeName));
  check('switched to a shared scope', scopeName && scopeName.trim().length > 0, scopeName);
  check('no errors switching scope', errors.length === errBefore,
    errors.slice(errBefore, errBefore + 3).join('\n        '));

  // add a rule in this scope and set it to something visibly different from the name default
  await page.locator('text=Add rule').first().click();
  await page.waitForTimeout(3000);
  const field = page.locator('.grok-order-rule-row:not(.grok-order-rule-header)').first()
    .locator('.grok-order-rule-cell').nth(2).locator('input');
  await field.fill('createdOn');
  await field.press('Enter');
  await page.waitForTimeout(3000);
  await page.screenshot({ path: 'grp-02-rule.png' });

  // Save and apply -> writes the scope override server-side
  const save = page.locator('text=SAVE AND APPLY').first();
  check('Save and apply button present', await save.count() > 0);
  await save.click();
  await page.waitForTimeout(12000);
  await page.screenshot({ path: 'grp-03-saved.png' });

  const ov = await overrides(page);
  console.log('    server overrides:', JSON.stringify(ov).slice(0, 400));
  const stored = ov && ov.orderRules;
  check('server stored the orderRules override for the scope',
    !!stored, JSON.stringify(ov).slice(0, 200));
  check('stored override carries the rule',
    stored && JSON.stringify(stored).indexOf('createdOn') >= 0, JSON.stringify(stored));

  // ---------- 2. a fresh session (no personal settings) inherits it
  const ctx2 = await browser.newContext({ viewport: { width: 1600, height: 1000 } });
  const page2 = await ctx2.newPage();
  await login(page2);
  const inherited = await page2.evaluate(() => window.grok && window.grok.settings
    ? window.grok.settings.orderRules : null);
  console.log('    fresh session orderRules:', JSON.stringify(inherited));
  check('fresh session inherits the group default',
    !!inherited && inherited.indexOf('createdOn') >= 0, JSON.stringify(inherited));

  // ---------- 3. and the tree actually orders by it
  await page2.evaluate(() => {
    const tri = document.querySelector('[name="tree-expander-Spaces"]');
    if (tri) tri.click();
  });
  await page2.waitForTimeout(10000);
  await page2.screenshot({ path: 'grp-04-tree.png' });
  const order = (await page2.evaluate(() =>
    [...document.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-item-label')]
      .map(e => e.innerText.trim()).join(String.fromCharCode(10))))
    .split(String.fromCharCode(10)).map(x => x.trim())
    .filter(l => ['alpha', 'middle', 'ordered', 'plain', 'zebra'].includes(l));
  console.log('    root spaces under the group default (createdOn):', JSON.stringify(order));
  check('group default reorders the Browse tree',
    order.length === 5 && order.join(',') !== 'alpha,middle,ordered,plain,zebra',
    JSON.stringify(order));

  // ---------- cleanup: drop the override so later runs start clean
  await page.evaluate(async (api) => {
    await fetch(`${api}/admin/settings/override`, { method: 'DELETE', credentials: 'include' });
  }, API);
  const after = await overrides(page);
  console.log('    overrides after cleanup:', JSON.stringify(after).slice(0, 200));

  console.log('\n--- errors (' + errors.length + '):');
  errors.slice(0, 8).forEach(e => console.log('      ' + e));
  await browser.close();
  console.log(failures === 0 ? '\nALL GROUP-DEFAULT CHECKS PASSED' : `\n${failures} CHECK(S) FAILED`);
  process.exit(failures === 0 ? 0 : 1);
})().catch(e => { console.error('HARNESS ERROR', e); process.exit(2); });
