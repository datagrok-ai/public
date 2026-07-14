import { test as base, expect, type Page, type BrowserContext } from '@playwright/test';

export const BASE = process.env.DATAGROK_URL!;

/** Selector for a host wrapping any DG input — `[name="input-host-<safeName>"]` is set in input_base.dart. */
export const inputHost = (safeName: string): string => `[name="input-host-${safeName}"]`;

/** Editor inside an input host. */
export const inputEditor = (safeName: string): string => `${inputHost(safeName)} input.ui-input-editor`;

/**
 * Shared, already-authenticated context + page reused across every test in the worker.
 * `global-setup` logs in and writes `e2e/.auth.json`; we boot the Datagrok SPA exactly
 * once here, then each test re-enters Diff Studio through a single model-load navigation
 * (`openModelFromLibrary`) instead of the previous two full reloads (`openDiffStudio` +
 * model). Reusing one context keeps the HTTP cache warm, so those navigations are fast.
 */
export const test = base.extend<{}, { sharedContext: BrowserContext; sharedPage: Page }>({
  sharedContext: [async ({ browser }, use) => {
    const context = await browser.newContext({ storageState: 'e2e/.auth.json' });
    await use(context);
    await context.close();
  }, { scope: 'worker' }],
  sharedPage: [async ({ sharedContext }, use) => {
    const page = await sharedContext.newPage();
    await page.goto(BASE, { waitUntil: 'domcontentloaded', timeout: 120_000 });
    await page.waitForFunction(() => document.querySelector('.grok-preloader') == null,
      undefined, { timeout: 120_000 });
    await page.locator('.d4-ribbon').first().waitFor({ timeout: 60_000 });
    await use(page);
  }, { scope: 'worker' }],
  // Override the built-in test-scoped `context`/`page` so existing `async ({ page }) => …`
  // and `async ({ page, context }) => …` bodies transparently receive the shared worker pair.
  context: async ({ sharedContext }, use) => {
    await use(sharedContext);
  },
  page: async ({ sharedPage }, use) => {
    await use(sharedPage);
  },
});

export { expect };

/** Navigate to Datagrok home and wait for the main ribbon to mount. */
export async function openHome(page: Page): Promise<void> {
  await page.goto(BASE);
  await page.waitForSelector('.d4-ribbon', { timeout: 30_000 });
}

/**
 * Ensure the Datagrok SPA is up (ribbon mounted). The worker `sharedPage` fixture already
 * booted the platform once, so this is a fast readiness check rather than a full
 * `/apps/DiffStudio` reload — the single per-test page load happens in `openModelFromLibrary`,
 * which navigates straight to the model. Kept as a named step so test bodies read unchanged.
 */
export async function openDiffStudio(page: Page): Promise<void> {
  await page.locator('.d4-ribbon').first().waitFor({ timeout: 30_000 });
}

/**
 * Map of TITLE → EDITOR_STATE values (lifted from DiffStudio/src/app.ts).
 * Used for direct URL navigation to a loaded model.
 */
export const STATE_BY_TITLE: Record<string, string> = {
  'Basic': 'basic',
  'Advanced': 'advanced',
  'Extended': 'extended',
  'Chem reactions': 'chem-react',
  "Robertson's model": 'robertson',
  'Fermentation': 'fermentation',
  'PK': 'pk',
  'PK-PD': 'pk-pd',
  'Acid production': 'ga-production',
  'Nimotuzumab': 'nimotuzumab',
  'Bioreactor': 'bioreactor',
  'Pollution': 'pollution',
};

/**
 * Load a built-in model directly via URL navigation (the address-bar equivalent of
 * clicking through Library → Model in the UI). DiffStudio's `runSolverApp` parses
 * `startingPath` and only enters the model-load branch when the path contains `?params:`
 * (see `DiffStudio/src/app.ts` runSolverApp:287). Without it the parser treats the last
 * segment as a folder name and falls back to "last called model".
 *
 * `modelTitle` matches TITLE values (e.g. 'PK-PD', 'Bioreactor', 'Acid production').
 */
export async function openModelFromLibrary(page: Page, modelTitle: string): Promise<void> {
  const state = STATE_BY_TITLE[modelTitle];
  if (!state) throw new Error(`Unknown DiffStudio model title: "${modelTitle}"`);
  await page.goto(`${BASE}/apps/DiffStudio/Library/${state}?params:`);
  await page.waitForSelector('.d4-ribbon', { timeout: 30_000 });
  // `.diff-studio-ribbon-widget` is the folder/Open icon mounted on the model-view ribbon.
  await page.waitForSelector('.diff-studio-ribbon-widget', { timeout: 60_000 });
  // Wait for the inputs panel to populate (model-specific input hosts appear).
  await page.waitForFunction(() => document.querySelectorAll('[name^="input-host-"]').length >= 3,
    null, { timeout: 60_000 });
  // Solving runs asynchronously after state set — give it a moment to render charts.
  await page.waitForTimeout(2000);
}

/**
 * Reset the shell between scenarios — Escape to dismiss popups, then `grok.shell.closeAll()`
 * to drop views/tables in-place. No navigation: the shared worker page stays booted, and the
 * next test's `openModelFromLibrary` re-enters Diff Studio with a single model-load reload.
 */
export async function resetShell(page: Page): Promise<void> {
  await page.keyboard.press('Escape').catch(() => {});
  await page.evaluate(() => {
    const g = (window as any).grok;
    g?.shell?.closeAll?.();
  }).catch(() => {});
  await page.waitForTimeout(200);
}

/** Click a ribbon item by its visible text label. */
export async function clickRibbonText(page: Page, text: string): Promise<void> {
  const item = page.locator('.d4-ribbon-item', { hasText: text }).first();
  await item.waitFor({ timeout: 10_000 });
  await item.click();
}

/** Toggle a boolean ribbon switch (e.g. the "Edit" toggle). Returns the new on/off state. */
export async function toggleRibbonSwitch(page: Page, label: string): Promise<boolean> {
  const item = page.locator('.d4-ribbon-item', { hasText: label }).first();
  await item.waitFor({ timeout: 10_000 });
  await item.locator('.ui-input-bool-switch .ui-input-editor').click();
  await page.waitForTimeout(800);
  return await item.locator('.ui-input-switch').evaluate(el => el.classList.contains('ui-input-switch-on'));
}

/** Read the current on/off state of a boolean ribbon switch. */
export async function ribbonSwitchOn(page: Page, label: string): Promise<boolean> {
  const sw = page.locator('.d4-ribbon-item', { hasText: label }).first().locator('.ui-input-switch');
  if (await sw.count() === 0) return false;
  return await sw.evaluate(el => el.classList.contains('ui-input-switch-on'));
}

/** Hover an input's label and read the resulting tooltip text. */
export async function readInputTooltip(page: Page, safeName: string): Promise<string> {
  const label = page.locator(`${inputHost(safeName)} label`).first();
  await label.waitFor({ timeout: 5_000 });
  await label.hover();
  await page.waitForTimeout(900);
  const tooltip = page.locator('.d4-tooltip').first();
  const text = (await tooltip.count() > 0) ? (await tooltip.textContent() ?? '').trim() : '';
  // move pointer away so the tooltip dismisses before the next hover
  await page.mouse.move(0, 0);
  await page.waitForTimeout(200);
  return text;
}

/** Set the value of a text/number input by replacing its content and confirming with Tab. */
export async function setInputValue(page: Page, safeName: string, value: string): Promise<void> {
  const ed = page.locator(inputEditor(safeName)).first();
  // `fill()` clears the existing value and types in one atomic operation — more reliable than
  // click+Ctrl+A+type which sometimes loses Ctrl+A to a global "Select all rows" binding when
  // a TableView grid steals focus.
  await ed.fill(value);
  await page.keyboard.press('Tab');
  await page.waitForTimeout(1500);
}

/**
 * Resolve the real `input-host-<safeName>` whose safeName equals `caption` ignoring case, and
 * return that safeName (or '' if absent). The native library view names hosts by the lowercase
 * variable (e.g. `dose`), while a model reopened from Model Hub names them by caption (e.g.
 * `Dose`). Callers that touch Model-Hub-opened models use this so the casing difference doesn't
 * break selectors.
 */
export async function resolveInputHostName(page: Page, caption: string): Promise<string> {
  return await page.evaluate((cap) => {
    const want = `input-host-${cap}`.toLowerCase();
    for (const h of Array.from(document.querySelectorAll('[name^="input-host-"]'))) {
      const name = h.getAttribute('name') ?? '';
      if (name.toLowerCase() === want) return name.replace(/^input-host-/, '');
    }
    return '';
  }, caption);
}

/** Click the +/- clicker for a numeric input. The icons live inside the input host. */
export async function clickerIncrement(page: Page, safeName: string, times = 1): Promise<void> {
  for (let i = 0; i < times; i++) {
    await page.locator(`${inputHost(safeName)} [name="icon-plus"]`).click({ force: true });
    await page.waitForTimeout(300);
  }
}

/** Click the -- clicker for a numeric input `times` times. */
export async function clickerDecrement(page: Page, safeName: string, times = 1): Promise<void> {
  for (let i = 0; i < times; i++) {
    await page.locator(`${inputHost(safeName)} [name="icon-minus"]`).click({ force: true });
    await page.waitForTimeout(300);
  }
}

/** Choice input that renders as <select>. */
export async function selectChoice(page: Page, safeName: string, value: string): Promise<void> {
  const select = page.locator(`${inputHost(safeName)} select`);
  await select.waitFor({ timeout: 10_000 });
  await select.selectOption({ label: value });
  await page.waitForTimeout(1500);
}

/** Click a dock tab header by its title text. */
export async function clickTab(page: Page, text: string): Promise<void> {
  const tab = page.locator('.tab-handle', {
    has: page.locator('.tab-handle-text', { hasText: new RegExp(`^${escapeRegex(text)}$`) }),
  }).first();
  await tab.waitFor({ timeout: 10_000 });
  await tab.click();
}

/** Returns titles of all currently rendered dock tabs (non-empty). */
export async function listTabs(page: Page): Promise<string[]> {
  const titles = await page.locator('.tab-handle-text').allInnerTexts();
  return titles.map(t => t.trim()).filter(t => t.length > 0);
}

/** True iff a dock tab with the given title is currently rendered. */
export async function hasTab(page: Page, text: string): Promise<boolean> {
  const tabs = await listTabs(page);
  return tabs.includes(text);
}

function escapeRegex(s: string): string {
  return s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}
