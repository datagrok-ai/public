// Shared helpers for Projects regression specs.
// All specs run against https://dev.datagrok.ai using the storageState from
// ../auth-dev.json — no interactive login. Helpers transcribed from
// .claude/skills/grok-browser/references/projects.md and
// .claude/skills/grok-browser/references/widgets/dialog.md.
import {Page, expect} from '@playwright/test';
import * as path from 'path';

export const BASE_URL = 'https://dev.datagrok.ai';
export const AUTH_FILE = path.resolve(__dirname, '..', 'auth-dev.json');

export const projectsTestOptions = {
  storageState: AUTH_FILE,
  baseURL: BASE_URL,
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
  actionTimeout: 15_000,
  navigationTimeout: 120_000,
};

export async function evalJs<T = any>(page: Page, script: string): Promise<T> {
  return page.evaluate(script) as Promise<T>;
}

export async function gotoApp(page: Page) {
  await page.goto(BASE_URL);
  await page.locator('[name="Browse"]').waitFor({timeout: 60_000});
}

export async function closeAll(page: Page) {
  await evalJs(page, 'grok.shell.closeAll()');
}

export async function setupSession(page: Page) {
  await evalJs(page, `(() => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
  })()`);
}

export async function openCsv(page: Page, path: string) {
  await evalJs(page, `(async () => {
    const df = await grok.dapi.files.readCsv('${path}');
    grok.shell.addTableView(df);
    return df.rowCount;
  })()`);
  await page.waitForTimeout(1500);
}

// Drive the Save Project dialog. Per the platform's own UITests pattern
// (scripts-layout.test.ts:415) the toolbar Save button is reachable as
// `button[name="button-Save"]:visible` — explicit `button` element +
// `:visible` filter + `.first()`. Bare `[name="button-Save"]` returns
// undefined in Playwright (DOM has hidden duplicates that confuse the
// resolver).
//
// Match the dialog in two stages:
//   * `[name="dialog-Save-project"]` — name attribute (preferred when present)
//   * `.d4-dialog:has-text("Save project")` — text-based fallback
export async function saveProjectViaDialog(page: Page, name: string) {
  // Wait for the toolbar to settle — under Playwright the toolbar render
  // can lag the dataframe render by several seconds. Also guard against
  // overlay banners that intercept pointer events on the toolbar.
  const saveBtn = page.locator('button[name="button-Save"]:visible').first();
  await saveBtn.waitFor({timeout: 30_000, state: 'visible'});
  await page.waitForTimeout(500);  // settle render
  // Try regular click first; fall back to JS DOM click if intercepted.
  try {
    await saveBtn.click({timeout: 5_000});
  } catch (_) {
    await page.evaluate(() => {
      const candidates = Array.from(document.querySelectorAll('button[name="button-Save"]'));
      const visible = candidates.find(b => (b as HTMLElement).offsetParent !== null);
      if (visible) (visible as HTMLElement).click();
    });
  }
  const dialog = page.locator(
    '.d4-dialog[name="dialog-Save-project"], .d4-dialog:has-text("Save project")',
  ).first();
  await dialog.waitFor({timeout: 15_000});
  const nameInput = dialog.locator(
    'input[name="input-Name"], input[type="text"].ui-input-editor',
  ).first();
  await nameInput.fill(name);
  await dialog.locator('button.ui-btn-ok, [name="button-OK"]').first().click();
  // Save dialog auto-opens a Share dialog after OK. Match by name-attribute
  // prefix `dialog-Share-*` (preferred) or by title-text prefix `^Share /`
  // — the user-typed `name` is server-normalized (PascalCase, no dashes) in
  // the share dialog title, so exact `Share ${name}` text match doesn't fire.
  const shareDialog = page.locator(
    '.d4-dialog[name^="dialog-Share-"], .d4-dialog:has-text("Share ")',
  ).first();
  if (await shareDialog.isVisible({timeout: 10_000}).catch(() => false)) {
    const cancel = shareDialog.locator(
      '[name="button-CANCEL"], button.ui-btn-cancel, button:has-text("Cancel")',
    ).first();
    if (await cancel.isVisible({timeout: 2_000}).catch(() => false))
      await cancel.click();
    else
      await page.keyboard.press('Escape');
    await expect(shareDialog).toBeHidden({timeout: 10_000});
  }
}

// Save via JS API (no UI). Use for setup state where the scenario doesn't
// own the Save dialog.
export async function saveProjectViaApi(page: Page, name: string): Promise<string> {
  return await evalJs<string>(page, `(async () => {
    const tv = grok.shell.tv;
    const layout = tv.saveLayout();
    const proj = grok.shell.project;
    proj.name = '${name}';
    const saved = await grok.dapi.projects.save(proj);
    return saved.id;
  })()`);
}

export async function projectExists(page: Page, name: string): Promise<boolean> {
  return await evalJs<boolean>(page,
    `(async () => (await grok.dapi.projects.filter('name = "${name}"').first()) != null)()`);
}

export async function deleteProjectByName(page: Page, name: string) {
  await evalJs(page, `(async () => {
    const p = await grok.dapi.projects.filter('name = "${name}"').first();
    if (p) await grok.dapi.projects.delete(p);
  })()`).catch(() => {});
}

export async function pollUntilProjectExists(page: Page, name: string, timeout = 60_000) {
  await expect.poll(async () => projectExists(page, name),
    {timeout, intervals: [500, 1000, 2000, 5000]}).toBe(true);
}

export async function navigateToDashboards(page: Page) {
  await page.goto(BASE_URL + '/projects');
  await page.waitForSelector('.grok-gallery-grid', {timeout: 30_000});
}
