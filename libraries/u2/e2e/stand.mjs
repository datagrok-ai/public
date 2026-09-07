/* Lane 2 — a real stand: login, a published package, real routing. Everything else (capture,
   reporting, screenshots) is shared with local.mjs; only the way in differs. Requires a stand that
   is free: the login lands as `admin`, and a run publishes nothing but does drive somebody's UI. */
import {launch, note} from './local.mjs';

/** The browser-facing composite (client + API); the API port alone will not serve login.html. */
export const STAND_URL = process.env.U2_STAND_URL ?? 'http://localhost:8888';
export const STAND_LOGIN = process.env.U2_STAND_LOGIN ?? 'admin';
export const STAND_PASSWORD = process.env.U2_STAND_PASSWORD ?? 'admin';
export const LOGIN_TIMEOUT = Number(process.env.U2_LOGIN_TIMEOUT ?? 300000);

export async function openStand(options = {}) {
  const {browser, page} = await launch(options);
  await page.goto(`${STAND_URL}/login.html`, {waitUntil: 'load', timeout: LOGIN_TIMEOUT});
  await page.waitForSelector('input[placeholder="Login or Email"]',
    {timeout: LOGIN_TIMEOUT, state: 'attached'});
  await page.waitForTimeout(2500);
  await page.locator('input[placeholder="Login or Email"]:visible').first().fill(STAND_LOGIN);
  await page.locator('input[type="password"]:visible').first().fill(STAND_PASSWORD);
  await page.locator('button:visible', {hasText: /^Login$/i}).first().click();
  await page.waitForFunction(() => {
    try {
      return !!(window.DG && DG.Func && grok.shell.user && grok.shell.user.login);
    } catch (e) {
      return false;
    }
  }, null, {timeout: LOGIN_TIMEOUT});
  await page.waitForTimeout(4000);
  return {browser, page};
}

/** Opens an app the way a user does — the platform route, which only resolves for a published
 * package. Answers whether the view came up. */
export async function openAppByRoute(page, pkg, appName, selector, timeout = 120000) {
  await page.goto(`${STAND_URL}/apps/${pkg}/${appName}`, {waitUntil: 'load', timeout});
  try {
    await page.waitForSelector(selector, {timeout});
    return true;
  } catch (e) {
    note('stand/route', `${pkg}/${appName} did not render ${selector}: ${String(e).slice(0, 200)}`);
    return false;
  }
}

/** Browse → Apps → <path> → <app>, the sidebar route a `meta.browsePath` app is filed under. */
export async function openAppFromBrowse(page, browsePath, appName) {
  for (let i = 0; i < 3; i++) {
    if (await page.locator('.d4-tree-view-group-label', {hasText: /^Apps$/}).count() > 0)
      break;
    await page.locator('[name="icon-compass"]').first().click().catch(() => {});
    await page.waitForTimeout(2000);
  }
  await page.locator('.d4-tree-view-group-label', {hasText: /^Apps$/}).first().click();
  await page.waitForTimeout(2000);
  await page.locator('.d4-tree-view-group-label:visible, .d4-tree-view-item-label:visible',
    {hasText: new RegExp(`^${browsePath}$`)}).first().click();
  await page.waitForTimeout(2000);
  // :visible matters: collapsed tree groups keep hidden item labels in the DOM, and .first()
  // would pick one of those over the visible gallery card, then time out waiting to click it.
  // Double-click: a single click only previews the app in the context panel.
  await page.locator(
    '.d4-tree-view-group-label:visible, .d4-tree-view-item-label:visible, .grok-gallery-grid *:visible',
    {hasText: new RegExp(`^${appName}$`)}).first().dblclick();
  await page.waitForTimeout(4000);
}
