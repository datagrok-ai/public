import { Page } from '@playwright/test';
import { RIBBON, BROWSE_HEADER_HOME } from '../browse/selectors';
import { BASE, ensureBrowsePanelOpen } from '../browse/helpers';
import {
  WELCOME_VIEW,
  WIDGETS_HOST,
  WIDGET_HOST,
  customizeLink,
  customizeToggle,
  WIDGET_FUNC_KEY,
  widgetByTitle,
  widgetCloseIcon,
} from './selectors';

// Server endpoint that the welcome view's settings sync hits. User widget settings are
// kept in memory and flushed to the server on a ~10s timer (see UserSettingsStorage docs),
// so persistence checks must wait for this PUT before reloading — a blind delay is flaky.
const WIDGETS_SYNC_URL = '/api/user_settings_storage/widgets';

/** Re-export so test files can use a single import surface. */
export { watchErrors, expectNoErrors } from '../browse/helpers';

// dev can be slow — and occasionally hang — booting the app shell. Give the first load a
// long window (a reload can interrupt a slow-but-progressing boot); if it still hangs, a
// stuck SPA usually recovers after a reload.
const RIBBON_TIMEOUTS = [70_000, 40_000, 40_000];
const RIBBON_BACKOFF = 3_000;
// Cap for the navigation call itself (goto/reload); the real wait is on the ribbon.
const NAV_TIMEOUT = 70_000;

/** Navigate (goto/reload) and wait for the app ribbon, retrying the navigation if it hangs. */
async function loadAppShell(page: Page, navigate: () => Promise<unknown>): Promise<void> {
  let lastErr: unknown;
  for (let attempt = 0; attempt < RIBBON_TIMEOUTS.length; attempt++) {
    try {
      await navigate();
      await page.waitForSelector(RIBBON, { timeout: RIBBON_TIMEOUTS[attempt] });
      return;
    } catch (e) {
      lastErr = e;
      // Let a transient dev blip clear before re-navigating.
      await page.waitForTimeout(RIBBON_BACKOFF);
    }
  }
  throw lastErr;
}

/** Open the app and make sure the Home page (welcome view) with widgets is shown. */
export async function openHomePage(page: Page): Promise<void> {
  await loadAppShell(page, () => page.goto(BASE, { waitUntil: 'domcontentloaded', timeout: NAV_TIMEOUT }));
  await ensureHomePage(page);
}

/** Reload the current page and wait for the Home page to be ready again. */
export async function reloadHome(page: Page): Promise<void> {
  await loadAppShell(page, () => page.reload({ waitUntil: 'domcontentloaded', timeout: NAV_TIMEOUT }));
  await ensureHomePage(page);
}

/**
 * Reset to the Home page between tests on a shared (worker-scoped) page.
 * First test in the worker pays for a full app-shell boot; afterwards we just close all
 * open views in-app (`grok.shell.closeAll()` returns to Home) — far cheaper and avoids the
 * flaky SPA reboot.
 */
export async function resetHome(page: Page): Promise<void> {
  if (!(await page.locator(RIBBON).isVisible().catch(() => false))) {
    await openHomePage(page);
    return;
  }
  await page.evaluate(() => (window as any).grok?.shell?.closeAll());
  await ensureHomePage(page);
}

/** Ensure the welcome view is the current view (click Home if some other view is open). */
export async function ensureHomePage(page: Page): Promise<void> {
  const welcome = page.locator(WELCOME_VIEW);
  if (!(await welcome.isVisible().catch(() => false))) {
    await ensureBrowsePanelOpen(page);
    await page.locator(BROWSE_HEADER_HOME).click();
    await welcome.waitFor({ state: 'visible', timeout: 20_000 });
  }
  await page.locator(WIDGETS_HOST).waitFor({ state: 'visible', timeout: 20_000 });
}

/** Returns the widget titles in their DOM order (top to bottom). */
export async function widgetTitlesInOrder(page: Page): Promise<string[]> {
  return page.locator(`${WIDGET_HOST}[widget-title]`).evaluateAll((hosts) =>
    hosts.map((h) => h.getAttribute('widget-title') || ''));
}

/** Reads the persisted `ignored` flag for a widget from user settings (verification read). */
export async function readWidgetIgnored(page: Page, friendlyName: string): Promise<boolean | undefined> {
  const key = WIDGET_FUNC_KEY[friendlyName];
  return page.evaluate((fn) => {
    const grok = (window as any).grok;
    const raw = grok?.userSettings?.get('widgets') || {};
    if (raw[fn] == null) return undefined;
    try { return JSON.parse(raw[fn]).ignored; }
    catch { return undefined; }
  }, key);
}

/**
 * Arm a wait for the next settings-sync PUT, run `action`, then await the sync.
 * Use whenever a test changes widget visibility and then reloads — guarantees the
 * change reached the server first.
 */
export async function withSettingsSync(page: Page, action: () => Promise<void>): Promise<void> {
  const sync = page.waitForResponse(
    (r) => r.url().includes(WIDGETS_SYNC_URL) && r.request().method() === 'PUT',
    { timeout: 20_000 },
  );
  await action();
  await sync;
}

/** Hover the widget to reveal its close icon, then click it to remove the widget. */
export async function closeWidgetViaIcon(page: Page, title: string): Promise<void> {
  await widgetByTitle(page, title).hover();
  const icon = widgetCloseIcon(page, title);
  // The icon is hover-revealed; force-click as a fallback if hover visibility lags.
  await icon.click({ timeout: 5_000 }).catch(() => icon.click({ force: true }));
  await widgetByTitle(page, title).waitFor({ state: 'detached', timeout: 5_000 });
}

/** Open the "Customize widgets..." form in the Context Panel. */
export async function openCustomizeForm(page: Page): Promise<void> {
  await customizeLink(page).click();
  // Wait until at least the Community toggle is rendered in the panel.
  await customizeToggle(page, 'Community').waitFor({ state: 'visible', timeout: 10_000 });
}

/** Set a widget's visibility through the Customize form and wait for it to take effect. */
export async function setWidgetVisible(page: Page, title: string, visible: boolean): Promise<void> {
  const toggle = customizeToggle(page, title);
  const checked = await toggle.isChecked();
  if (checked === visible) return;
  await toggle.click();
  if (visible)
    await widgetByTitle(page, title).waitFor({ state: 'visible', timeout: 5_000 });
  else
    await widgetByTitle(page, title).waitFor({ state: 'detached', timeout: 5_000 });
}

/**
 * Restore a widget to visible state in-app (no wait for the server sync). Cheap self-heal for
 * beforeEach — the shared worker page keeps the in-memory setting, so subsequent non-reload
 * tests see the widget immediately. Server-side cleanliness is handled once by cleanWidgetSettings().
 */
export async function restoreWidgetVisible(page: Page, title: string): Promise<void> {
  if (await widgetByTitle(page, title).isVisible().catch(() => false)) return;
  await openCustomizeForm(page);
  await setWidgetVisible(page, title, true);
}

/**
 * Reset all widget `ignored` flags to false and wait for the sync to land. Called once in
 * afterAll so the next run boots with every widget shown — the mutating tests only persist
 * their changes in-memory, and Nav-01 reloads (reading server state), so the server must be
 * left clean across runs.
 */
export async function cleanWidgetSettings(page: Page): Promise<void> {
  const synced = page.waitForResponse(
    (r) => r.url().includes(WIDGETS_SYNC_URL) && r.request().method() === 'PUT',
    { timeout: 15_000 },
  ).catch(() => undefined);
  await page.evaluate(() => {
    const grok = (window as any).grok;
    const raw = grok?.userSettings?.get('widgets') || {};
    const clean: Record<string, string> = {};
    for (const k of Object.keys(raw)) clean[k] = JSON.stringify({ ignored: false });
    if (Object.keys(clean).length) grok.userSettings.addAll('widgets', clean);
  });
  await synced;
}
