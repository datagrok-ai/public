import { Page } from '@playwright/test';

export const BASE = process.env.DATAGROK_URL!;

/** Selector for a host wrapping any DG input — `[name="input-host-<safeName>"]` is set in input_base.dart. */
export const inputHost = (safeName: string): string => `[name="input-host-${safeName}"]`;

/** Editor inside an input host. */
export const inputEditor = (safeName: string): string => `${inputHost(safeName)} input.ui-input-editor`;

/** Navigate to Datagrok home and wait for the main ribbon to mount. */
export async function openHome(page: Page): Promise<void> {
  await page.goto(BASE);
  await page.waitForSelector('.d4-ribbon', { timeout: 30_000 });
}

/**
 * Open the Diff Studio app via URL navigation; waits for the ribbon to mount.
 * Whether the landing page or a previously-loaded model is shown depends on platform
 * routing state — both are acceptable signals that the app is up.
 */
export async function openDiffStudio(page: Page): Promise<void> {
  await page.goto(`${BASE}/apps/DiffStudio`);
  await page.waitForSelector('.d4-ribbon', { timeout: 30_000 });
  await page.waitForTimeout(2000);
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
 * Reset the shell between scenarios — close any dialogs, balloons, then go home.
 * Pure UI: presses Escape to dismiss popups, then navigates away.
 */
export async function resetShell(page: Page): Promise<void> {
  await page.keyboard.press('Escape').catch(() => {});
  await page.goto(BASE);
  await page.waitForSelector('.d4-ribbon', { timeout: 30_000 });
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
