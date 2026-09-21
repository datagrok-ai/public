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

const WIDGETS_SYNC_URL = '/api/user_settings_storage/widgets';

export { watchErrors, expectNoErrors } from '../browse/helpers';

const RIBBON_TIMEOUTS = [70_000, 40_000, 40_000];
const RIBBON_BACKOFF = 3_000;

const NAV_TIMEOUT = 70_000;

async function loadAppShell(page: Page, navigate: () => Promise<unknown>): Promise<void> {
  let lastErr: unknown;
  for (let attempt = 0; attempt < RIBBON_TIMEOUTS.length; attempt++) {
    try {
      await navigate();
      await page.waitForSelector(RIBBON, { timeout: RIBBON_TIMEOUTS[attempt] });
      return;
    } catch (e) {
      lastErr = e;

      await page.waitForTimeout(RIBBON_BACKOFF);
    }
  }
  throw lastErr;
}

export async function openHomePage(page: Page): Promise<void> {
  await loadAppShell(page, () => page.goto(BASE, { waitUntil: 'domcontentloaded', timeout: NAV_TIMEOUT }));
  await ensureHomePage(page);
}

export async function reloadHome(page: Page): Promise<void> {
  await loadAppShell(page, () => page.reload({ waitUntil: 'domcontentloaded', timeout: NAV_TIMEOUT }));
  await ensureHomePage(page);
}

export async function resetHome(page: Page): Promise<void> {
  if (!(await page.locator(RIBBON).isVisible().catch(() => false))) {
    await openHomePage(page);
    return;
  }
  await page.evaluate(() => (window as any).grok?.shell?.closeAll());
  await ensureHomePage(page);
}

export async function ensureHomePage(page: Page): Promise<void> {
  const welcome = page.locator(WELCOME_VIEW);
  if (!(await welcome.isVisible().catch(() => false))) {
    await ensureBrowsePanelOpen(page);
    await page.locator(BROWSE_HEADER_HOME).click();
    await welcome.waitFor({ state: 'visible', timeout: 20_000 });
  }
  await page.locator(WIDGETS_HOST).waitFor({ state: 'visible', timeout: 20_000 });
}

export async function widgetTitlesInOrder(page: Page): Promise<string[]> {
  return page.locator(`${WIDGET_HOST}[widget-title]`).evaluateAll((hosts) =>
    hosts.map((h) => h.getAttribute('widget-title') || ''));
}

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

export async function withSettingsSync(page: Page, action: () => Promise<void>): Promise<void> {
  const sync = page.waitForResponse(
    (r) => r.url().includes(WIDGETS_SYNC_URL) && r.request().method() === 'PUT',
    { timeout: 20_000 },
  );
  await action();
  await sync;
}

export async function closeWidgetViaIcon(page: Page, title: string): Promise<void> {
  await widgetByTitle(page, title).hover();
  const icon = widgetCloseIcon(page, title);

  await icon.click({ timeout: 5_000 }).catch(() => icon.click({ force: true }));
  await widgetByTitle(page, title).waitFor({ state: 'detached', timeout: 5_000 });
}

export async function openCustomizeForm(page: Page): Promise<void> {
  await customizeLink(page).click();

  await customizeToggle(page, 'Community').waitFor({ state: 'visible', timeout: 10_000 });
}

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

export async function restoreWidgetVisible(page: Page, title: string): Promise<void> {
  if (await widgetByTitle(page, title).isVisible().catch(() => false)) return;
  await openCustomizeForm(page);
  await setWidgetVisible(page, title, true);
}

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
