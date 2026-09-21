import { test as base, expect, type Page, type BrowserContext } from '@playwright/test';

export const BASE = process.env.DATAGROK_URL!;

export const inputHost = (safeName: string): string => `[name="input-host-${safeName}"]`;

export const inputEditor = (safeName: string): string => `${inputHost(safeName)} input.ui-input-editor`;

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

  context: async ({ sharedContext }, use) => {
    await use(sharedContext);
  },
  page: async ({ sharedPage }, use) => {
    await use(sharedPage);
  },
});

export { expect };

export async function openHome(page: Page): Promise<void> {
  await page.goto(BASE);
  await page.waitForSelector('.d4-ribbon', { timeout: 30_000 });
}

export async function openDiffStudio(page: Page): Promise<void> {
  await page.locator('.d4-ribbon').first().waitFor({ timeout: 30_000 });
}

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

export async function openModelFromLibrary(page: Page, modelTitle: string): Promise<void> {
  const state = STATE_BY_TITLE[modelTitle];
  if (!state) throw new Error(`Unknown DiffStudio model title: "${modelTitle}"`);
  await page.goto(`${BASE}/apps/DiffStudio/Library/${state}?params:`);
  await page.waitForSelector('.d4-ribbon', { timeout: 30_000 });

  await page.waitForSelector('.diff-studio-ribbon-widget', { timeout: 60_000 });

  await page.waitForFunction(() => document.querySelectorAll('[name^="input-host-"]').length >= 3,
    null, { timeout: 60_000 });

  await page.waitForTimeout(2000);
}

export async function resetShell(page: Page): Promise<void> {
  await page.keyboard.press('Escape').catch(() => {});
  await page.evaluate(() => {
    const g = (window as any).grok;
    g?.shell?.closeAll?.();
  }).catch(() => {});
  await page.waitForTimeout(200);
}

export async function clickRibbonText(page: Page, text: string): Promise<void> {
  const item = page.locator('.d4-ribbon-item', { hasText: text }).first();
  await item.waitFor({ timeout: 10_000 });
  await item.click();
}

export async function toggleRibbonSwitch(page: Page, label: string): Promise<boolean> {
  const item = page.locator('.d4-ribbon-item', { hasText: label }).first();
  await item.waitFor({ timeout: 10_000 });
  await item.locator('.ui-input-bool-switch .ui-input-editor').click();
  await page.waitForTimeout(800);
  return await item.locator('.ui-input-switch').evaluate(el => el.classList.contains('ui-input-switch-on'));
}

export async function ribbonSwitchOn(page: Page, label: string): Promise<boolean> {
  const sw = page.locator('.d4-ribbon-item', { hasText: label }).first().locator('.ui-input-switch');
  if (await sw.count() === 0) return false;
  return await sw.evaluate(el => el.classList.contains('ui-input-switch-on'));
}

export async function readInputTooltip(page: Page, safeName: string): Promise<string> {
  const label = page.locator(`${inputHost(safeName)} label`).first();
  await label.waitFor({ timeout: 5_000 });
  await label.hover();
  await page.waitForTimeout(900);
  const tooltip = page.locator('.d4-tooltip').first();
  const text = (await tooltip.count() > 0) ? (await tooltip.textContent() ?? '').trim() : '';

  await page.mouse.move(0, 0);
  await page.waitForTimeout(200);
  return text;
}

export async function setInputValue(page: Page, safeName: string, value: string): Promise<void> {
  const ed = page.locator(inputEditor(safeName)).first();

  await ed.fill(value);
  await page.keyboard.press('Tab');
  await page.waitForTimeout(1500);
}

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

export async function clickerIncrement(page: Page, safeName: string, times = 1): Promise<void> {
  for (let i = 0; i < times; i++) {
    await page.locator(`${inputHost(safeName)} [name="icon-plus"]`).click({ force: true });
    await page.waitForTimeout(300);
  }
}

export async function clickerDecrement(page: Page, safeName: string, times = 1): Promise<void> {
  for (let i = 0; i < times; i++) {
    await page.locator(`${inputHost(safeName)} [name="icon-minus"]`).click({ force: true });
    await page.waitForTimeout(300);
  }
}

export async function selectChoice(page: Page, safeName: string, value: string): Promise<void> {
  const select = page.locator(`${inputHost(safeName)} select`);
  await select.waitFor({ timeout: 10_000 });
  await select.selectOption({ label: value });
  await page.waitForTimeout(1500);
}

export async function clickTab(page: Page, text: string): Promise<void> {
  const tab = page.locator('.tab-handle', {
    has: page.locator('.tab-handle-text', { hasText: new RegExp(`^${escapeRegex(text)}$`) }),
  }).first();
  await tab.waitFor({ timeout: 10_000 });
  await tab.click();
}

export async function listTabs(page: Page): Promise<string[]> {
  const titles = await page.locator('.tab-handle-text').allInnerTexts();
  return titles.map(t => t.trim()).filter(t => t.length > 0);
}

export async function hasTab(page: Page, text: string): Promise<boolean> {
  const tabs = await listTabs(page);
  return tabs.includes(text);
}

function escapeRegex(s: string): string {
  return s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}
