import { test, expect, Page } from '@playwright/test';
import * as H from './helpers';
import { expandDbProvider, expandDbConnection, expandDbGroupWrapper, expandChildByLabel } from '../connections/helpers';

// Database meta (sticky meta on DB entities) for a Postgres table and column: fill Comment /
// LLM Comment, Save, verify persistence across a full page reload, clear, and verify removal.
// Connection-level Database meta is not surfaced in the UI, so it is out of scope.
// The table and column tests are independent (not serial) so one failing never skips the other.

// The Schemas group renders under a wrapper keyed by the connection's server-side name
// (`NorthwindTest` friendly name -> `PostgresTest` server name).
const CONN_FRIENDLY = 'NorthwindTest';
const CONN_SERVER = 'PostgresTest';
const SCHEMAS_WRAPPER_NAME = `div-Postgres-${CONN_SERVER}-Schemas`;

/** Click a tree node's label (selects the entity) found by exact text inside a parent subtree. */
async function selectChildByLabel(page: Page, parentName: string, label: string): Promise<void> {
  await page.waitForSelector(`[name="${parentName}"]`, { state: 'attached', timeout: 15_000 });
  await page.evaluate(({ p, l }) => {
    const parent = document.querySelector(`[name="${p}"]`)!;
    const el = Array.from(parent.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-item-label'))
      .find((e) => e.textContent?.trim() === l) as HTMLElement | undefined;
    if (!el) throw new Error(`tree node "${l}" not found under ${p}`);
    el.scrollIntoView({ block: 'center' });
    el.click();
  }, { p: parentName, l: label });
  await page.waitForTimeout(2500);
}

/** Expand Databases > Postgres > NorthwindTest > Schemas > public (idempotent). */
async function openPublicSchema(page: Page): Promise<void> {
  await H.gotoHome(page);
  await page.evaluate(() => document.body.classList.add('selenium'));
  await expandDbProvider(page, 'Postgres');
  await expandDbConnection(page, 'Postgres', CONN_FRIENDLY);
  await expandDbGroupWrapper(page, 'Postgres', CONN_SERVER, 'Schemas');
  await expandChildByLabel(page, SCHEMAS_WRAPPER_NAME, 'public');
}

/**
 * Select the categories table so its Database meta entity panel is shown.
 *
 * SCOPE NOTE: clicking a table node opens its grid in the main area, and the context panel then
 * follows the grid (replacing the DbTable entity panel). There is no pure-UI gesture that keeps the
 * DbTable entity selected, so the entity captured on click is restored as the current object via the
 * API. The Database meta fill / save / verify that this test asserts is performed entirely via the UI.
 */
async function selectTable(page: Page): Promise<void> {
  await openPublicSchema(page);
  await page.evaluate((parentName) => {
    const parent = document.querySelector(`[name="${parentName}"]`)!;
    const el = Array.from(parent.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-item-label'))
      .find((e) => e.textContent?.trim() === 'categories') as HTMLElement | undefined;
    if (!el) throw new Error('categories table node not found');
    el.scrollIntoView({ block: 'center' });
    el.click();
  }, SCHEMAS_WRAPPER_NAME);
  await page.waitForTimeout(700);
  await page.evaluate(() => { (window as any)._dbTable = (window as any).grok.shell.o; });
  await page.waitForTimeout(2500); // let the grid finish opening
  await page.evaluate(() => { (window as any).grok.shell.o = (window as any)._dbTable; });
  await page.waitForTimeout(2000);
}

/** Expand categories and select the categoryid column. */
async function selectColumn(page: Page): Promise<void> {
  await openPublicSchema(page);
  await expandChildByLabel(page, SCHEMAS_WRAPPER_NAME, 'categories');
  await selectChildByLabel(page, SCHEMAS_WRAPPER_NAME, 'categoryid');
}

// Scope to the entity property panel: when a table is selected its grid opens too, and the grid's
// own property panel would otherwise shadow the entity's Database meta section.
const PANEL = '.grok-entity-prop-panel';
const COMMENT_SEL = `${PANEL} [name="input-host-Comment"] input, ${PANEL} [name="input-host-Comment"] textarea`;
const LLM_SEL = `${PANEL} [name="input-host-LLM-Comment"] input, ${PANEL} [name="input-host-LLM-Comment"] textarea`;

/** Whether the Context Panel currently shows a Database meta section. */
async function dbMetaPresent(page: Page): Promise<boolean> {
  return (await page.locator(`${PANEL} .d4-accordion-pane-header`).filter({ hasText: /^Database\s*meta/ }).count()) > 0;
}

/** Expand the Database meta accordion pane if it is collapsed. */
async function ensureDbMetaExpanded(page: Page): Promise<void> {
  const header = page.locator(`${PANEL} .d4-accordion-pane-header`).filter({ hasText: /^Database\s*meta/ }).first();
  await header.waitFor({ state: 'visible', timeout: 10_000 });
  await header.scrollIntoViewIfNeeded();
  if (!(await header.evaluate((el) => el.classList.contains('expanded')))) {
    await header.click();
    await page.waitForTimeout(500);
  }
}

/** Read the Comment / LLM Comment values from the Database meta pane. */
async function readDbMeta(page: Page): Promise<{ comment: string; llm: string }> {
  await ensureDbMetaExpanded(page);
  const comment = page.locator(COMMENT_SEL).first();
  await comment.waitFor({ state: 'visible', timeout: 10_000 });
  return { comment: await comment.inputValue(), llm: await page.locator(LLM_SEL).first().inputValue() };
}

/** Type a value into a field (selecting any existing text first), then commit with Tab. */
async function typeField(page: Page, selector: string, value: string): Promise<void> {
  const input = page.locator(selector).first();
  await input.waitFor({ state: 'visible', timeout: 10_000 });
  await input.click({ clickCount: 3 });
  await page.keyboard.press('Control+a');
  await page.keyboard.press('Delete');
  if (value) await page.keyboard.type(value);
  await page.keyboard.press('Tab');
  await page.waitForTimeout(300);
}

/** Set Comment + LLM Comment via real typing, then click the pane's Save button. */
async function fillDbMetaAndSave(page: Page, comment: string, llm: string): Promise<void> {
  await ensureDbMetaExpanded(page);
  await typeField(page, COMMENT_SEL, comment);
  await typeField(page, LLM_SEL, llm);
  await page.evaluate((panelSel) => {
    const panel = document.querySelector(panelSel)!;
    const head = Array.from(panel.querySelectorAll('.d4-accordion-pane-header'))
      .find((h) => /^Database\s*meta/.test((h.textContent || '').trim())) as HTMLElement | undefined;
    (head?.parentElement?.querySelector('[name="button-Save"]') as HTMLElement | null)?.click();
  }, PANEL);
  await page.waitForTimeout(1800);
}

test('Sticky Meta: database meta on a table', async ({ page }) => {
  test.setTimeout(240_000);
  const val = `pwsm_t_${H.uniqueSuffix()}`;
  try {
    await selectTable(page);
    expect(await dbMetaPresent(page)).toBe(true);

    await fillDbMetaAndSave(page, val, val);
    // Persist across a full reload.
    await selectTable(page);
    await expect.poll(async () => (await readDbMeta(page)).comment, { timeout: 15_000, intervals: [800] }).toBe(val);
    expect((await readDbMeta(page)).llm).toBe(val);

    // Clear and verify removal across reload.
    await fillDbMetaAndSave(page, '', '');
    await selectTable(page);
    await expect.poll(async () => (await readDbMeta(page)).comment, { timeout: 15_000, intervals: [800] }).toBe('');
    expect((await readDbMeta(page)).llm).toBe('');
  } finally {
    // Best-effort: ensure no value is left behind.
    await fillDbMetaAndSave(page, '', '').catch(() => {});
  }
});

test('Sticky Meta: database meta on a column', async ({ page }) => {
  test.setTimeout(240_000);
  const val = `pwsm_c_${H.uniqueSuffix()}`;
  try {
    await selectColumn(page);
    expect(await dbMetaPresent(page)).toBe(true);

    await fillDbMetaAndSave(page, val, val);
    await selectColumn(page);
    await expect.poll(async () => (await readDbMeta(page)).comment, { timeout: 15_000, intervals: [800] }).toBe(val);
    expect((await readDbMeta(page)).llm).toBe(val);

    await fillDbMetaAndSave(page, '', '');
    await selectColumn(page);
    await expect.poll(async () => (await readDbMeta(page)).comment, { timeout: 15_000, intervals: [800] }).toBe('');
    expect((await readDbMeta(page)).llm).toBe('');
  } finally {
    await fillDbMetaAndSave(page, '', '').catch(() => {});
  }
});
