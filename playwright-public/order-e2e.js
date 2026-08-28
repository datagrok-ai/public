// End-to-end check of the Browse ordering feature against the worktree dev stack.
// Usage: node order-e2e.js
const { chromium } = require('playwright');

const BASE = process.env.DG_BASE || 'http://localhost:61006';
const errors = [];
let failures = 0;

function check(label, ok, detail) {
  if (!ok) failures++;
  console.log(`${ok ? 'PASS' : 'FAIL'}  ${label}${detail ? `\n        ${detail}` : ''}`);
}

(async () => {
  const browser = await chromium.launch({ headless: true });
  const page = await browser.newPage({ viewport: { width: 1600, height: 1000 } });
  const NOISE = /packages\/published|detectors\.js|package\.js|MIME type|401 \(Unauthorized\)|favicon/i;
  page.on('console', m => {
    if (m.type() !== 'error') return;
    const t = m.text().slice(0, 300);
    if (!NOISE.test(t)) errors.push('CONSOLE ' + t);
  });
  page.on('pageerror', e => errors.push('PAGEERROR ' + String(e).slice(0, 300)));

  console.log(`--- opening ${BASE}/login.html`);
  await page.goto(`${BASE}/login.html`, { waitUntil: 'domcontentloaded', timeout: 180000 });

  // Datagrok keeps hidden Sign-Up inputs with the same placeholders, so pin to the visible ones.
  const user = page.locator('input[placeholder="Login or Email"]:visible').first();
  await user.waitFor({ state: 'visible', timeout: 300000 });
  await user.fill('admin');
  await page.locator('input[placeholder="Password"]:visible').first().fill('admin');
  await page.locator('button:has-text("LOGIN"), .ui-btn:has-text("LOGIN")').first().click();

  // Datagrok is a SPA served from login.html - the URL never changes, so wait for the shell.
  const tree = page.locator('.d4-tree-view-group-label:visible', { hasText: /^Spaces$/ }).first();
  await tree.waitFor({ state: 'visible', timeout: 300000 }).catch(() => {});
  await page.waitForTimeout(10000);
  await page.screenshot({ path: 'shot-01-after-login.png' });

  const loggedIn = await tree.count() > 0;
  check('logged in, Browse tree rendered', loggedIn);
  if (!loggedIn) { await browser.close(); process.exit(1); }

  // ---- expand the Spaces node
  // Each group carries name="tree-expander-<caption>"; a real click misses the 12px triangle,
  // and the collapsed host is display:none so :visible filters never match.
  async function expand(name, waitMs) {
    const ok = await page.evaluate((n) => {
      const tri = document.querySelector(`[name="tree-expander-${n}"]`);
      if (!tri) return false;
      tri.click();
      return true;
    }, name);
    await page.waitForTimeout(waitMs || 7000);
    return ok;
  }
  const spaces = page.locator('.d4-tree-view-group-label', { hasText: /^Spaces$/ }).first();
  check('Spaces node present', await spaces.count() > 0);
  const before = errors.length;
  await expand('Spaces', 9000);
  await page.screenshot({ path: 'shot-02-spaces.png' });
  check('no new console errors expanding Spaces', errors.length === before,
    errors.slice(before, before + 3).join('\n        '));

  // ---- read the order of the root spaces
  const treeText = await page.evaluate(() =>
    [...document.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-item-label')]
      .map(e => e.innerText.trim()).join(String.fromCharCode(10)));
  const lines = treeText.split('\n').map(s => s.trim()).filter(Boolean);
  console.log('--- browse tree:\n' + lines.map(l => '      ' + l).join('\n'));

  const rootOrder = lines.filter(l => ['alpha', 'middle', 'ordered', 'plain', 'zebra'].includes(l));
  console.log('    root spaces in tree order:', JSON.stringify(rootOrder));
  check('all five root spaces listed', rootOrder.length === 5, JSON.stringify(rootOrder));

  // ---- expand 'ordered' (override: rank desc -> aaa, ccc, bbb) and 'plain' (name -> aaa, bbb, ccc)
  // Expander names are path-based: Spaces---ordered
  for (const name of ['ordered', 'plain'])
    check(`${name} expanded`, await expand(`Spaces---${name}`, 10000));
  await page.screenshot({ path: 'shot-03-children.png' });
  const afterText = await page.evaluate(() =>
    [...document.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-item-label')]
      .map(e => e.innerText.trim()).join(String.fromCharCode(10)));
  const afterLines = afterText.split('\n').map(s => s.trim()).filter(Boolean);
  console.log('--- tree after expanding:\n' + afterLines.map(l => '      ' + l).join('\n'));

  function childrenAfter(parent) {
    const i = afterLines.indexOf(parent);
    if (i < 0) return [];
    const out = [];
    for (let j = i + 1; j < afterLines.length; j++) {
      const l = afterLines[j];
      if (['alpha', 'middle', 'ordered', 'plain', 'zebra'].includes(l)) break;
      if (['aaa', 'bbb', 'ccc'].includes(l)) out.push(l);
    }
    return out;
  }
  const oc = childrenAfter('ordered');
  const pc = childrenAfter('plain');
  console.log('    ordered children:', JSON.stringify(oc));
  console.log('    plain   children:', JSON.stringify(pc));
  check("'ordered' children use its rank-desc override (aaa, ccc, bbb)",
    oc.join(',') === 'aaa,ccc,bbb', JSON.stringify(oc));
  check("'plain' children fall back to name (aaa, bbb, ccc)",
    pc.join(',') === 'aaa,bbb,ccc', JSON.stringify(pc));

  // The reported failure: a Project rule on `id` did nothing at any level. Set it through the
  // public JS API - the same path the settings page takes (property set -> change event) - rather
  // than poking localStorage, which the boot sequence rewrites.
  const applied = await page.evaluate(() => {
    if (!window.grok || !window.grok.settings) return 'no grok.settings';
    window.grok.settings.orderRules =
      JSON.stringify([{ type: 'Project', space: '', field: 'id', desc: false }]);
    return window.grok.settings.orderRules;
  });
  console.log('    orderRules now:', JSON.stringify(applied));
  check('rule set through grok.settings', typeof applied === 'string' && applied.indexOf('"id"') >= 0,
    JSON.stringify(applied));
  await page.waitForTimeout(12000);
  await expand('Spaces', 8000);
  await page.waitForTimeout(6000);
  await page.screenshot({ path: 'shot-04-by-id.png' });
  const idOrder = (await page.evaluate(() =>
    [...document.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-item-label')]
      .map(e => e.innerText.trim()).join(String.fromCharCode(10))))
    .split(String.fromCharCode(10)).map(x => x.trim())
    .filter(l => ['alpha', 'middle', 'ordered', 'plain', 'zebra'].includes(l));
  console.log('    root spaces ordered by id:', JSON.stringify(idOrder));
  check("ordering by 'id' differs from name order",
    idOrder.length === 5 && idOrder.join(',') !== 'alpha,middle,ordered,plain,zebra',
    JSON.stringify(idOrder));

  // Reported: the Desc switch did nothing for a root-level rule.
  const descApplied = await page.evaluate(() => {
    window.grok.settings.orderRules = JSON.stringify(
      [{ type: '', space: '@root', field: 'friendlyName', desc: true }]);
    return window.grok.settings.orderRules;
  });
  console.log('    root desc rule:', JSON.stringify(descApplied));
  await page.waitForTimeout(12000);
  await expand('Spaces', 8000);
  await page.waitForTimeout(6000);
  await page.screenshot({ path: 'shot-05-root-desc.png' });
  const descOrder = (await page.evaluate(() =>
    [...document.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-item-label')]
      .map(e => e.innerText.trim()).join(String.fromCharCode(10))))
    .split(String.fromCharCode(10)).map(x => x.trim())
    .filter(l => ['alpha', 'middle', 'ordered', 'plain', 'zebra'].includes(l));
  console.log('    root spaces, root-level rule desc=true:', JSON.stringify(descOrder));
  check('Desc reverses a root-level rule',
    descOrder.join(',') === 'zebra,plain,ordered,middle,alpha', JSON.stringify(descOrder));

  console.log('\n--- all console errors (' + errors.length + '):');
  errors.slice(0, 10).forEach(e => console.log('      ' + e));

  await browser.close();
  console.log(failures === 0 ? '\nALL BROWSE CHECKS PASSED' : `\n${failures} CHECK(S) FAILED`);
  process.exit(failures === 0 ? 0 : 1);
})().catch(e => { console.error('HARNESS ERROR', e); process.exit(2); });
