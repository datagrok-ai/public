import { test, expect, Page } from '@playwright/test';
import {
  CONTEXT_MENU,
  TREE_EXPAND_ARROW,
  TREE_EXPAND_ARROW_EXPANDED,
  contextMenuItem,
  treeNodeByPath,
} from './selectors';
import {
  goHome,
  ensureBrowsePanelOpen,
  ensureContextPanelOpen,
  watchErrors,
  expectNoErrors,
  expandTreeGroup,
} from './helpers';

const COMPUTE_PATH = ['Apps', 'Compute'];
const MH_PATH = [...COMPUTE_PATH, 'Model-Hub'];
const UNCAT_PATH = [...MH_PATH, 'Uncategorized'];
// A fresh stand has no models, and Model Hub only shows Uncategorized when some model has
// no department, so the specs seed their own.
const SEED_MODEL = 'PwModelHubSeed';
const SEED_PATH = [...UNCAT_PATH, SEED_MODEL];
// dev lists ~3,300 models: Model Hub's groups and the seed can each take longer than 10 s to appear.
const SEED_LIST_TIMEOUT = 30_000;

async function apiEnsureSeedModel(page: Page): Promise<void> {
  await page.evaluate(async (name) => {
    const grok = (window as any).grok;
    if ((await grok.dapi.scripts.filter(`name = "${name}"`).list()).length > 0)
      return;
    await grok.dapi.scripts.save((window as any).DG.Script.create(
      `//name: ${name}\n//language: javascript\n//meta.role: model\n//output: int result\nresult = 1;\n`));
  }, SEED_MODEL);
}

async function apiDeleteSeedModel(page: Page): Promise<void> {
  await page.evaluate(async (name) => {
    const grok = (window as any).grok;
    for (const s of await grok.dapi.scripts.filter(`name = "${name}"`).list())
      await grok.dapi.scripts.delete(s);
  }, SEED_MODEL);
}

test.describe('Browse Model Hub (Browse-ModelHub-*)', () => {
  test.beforeEach(async ({ page }) => {
    await goHome(page);
    await ensureBrowsePanelOpen(page);
    await ensureContextPanelOpen(page);
    // Model Hub ships in the Compute plugin, which isn't built on the minimal CI stack
    // (and Compute currently fails to build there). Skip when Apps > Compute is absent.
    await expandTreeGroup(page, 'Apps').catch(() => undefined);
    const hasCompute = await treeNodeByPath(page, COMPUTE_PATH)
      .waitFor({ state: 'visible', timeout: 5_000 }).then(() => true, () => false);
    test.skip(!hasCompute, 'Apps > Compute (Model Hub) not deployed on this stack');
    await apiEnsureSeedModel(page);
  });

  test.afterEach(async ({ page }) => {
    await apiDeleteSeedModel(page).catch(() => undefined);
  });

  test('Browse-ModelHub-01 — Apps > Compute > Model Hub group is reachable without errors', async ({ page }) => {
    const sink = watchErrors(page);

    await expandTreeGroup(page, 'Apps');
    await expandTreeGroup(page, COMPUTE_PATH);
    const mh = treeNodeByPath(page, MH_PATH);
    await expect(mh, 'Model Hub node must be present').toBeVisible({ timeout: 10_000 });
    // Click expands it.
    await mh.click();
    await page.waitForTimeout(1500);

    await expectNoErrors(page, sink);
  });

  test('Browse-ModelHub-02 — single click on a model in the tree updates Context Panel without errors', async ({ page }) => {
    // Was blocked by platform regression GROK-19740 (model click threw
    // `TypeError: p.append is not a function` in Dart code). That fix has landed — the
    // test now passes, so the former `test.fail` annotation is removed.
    const sink = watchErrors(page);

    await expandTreeGroup(page, 'Apps');
    await expandTreeGroup(page, COMPUTE_PATH);
    await expandTreeGroup(page, MH_PATH);
    await treeNodeByPath(page, UNCAT_PATH).waitFor({ state: 'visible', timeout: SEED_LIST_TIMEOUT });
    await expandTreeGroup(page, UNCAT_PATH);

    const model = treeNodeByPath(page, SEED_PATH);
    await expect(model, 'The seeded model must be listed under Uncategorized').toBeVisible({ timeout: SEED_LIST_TIMEOUT });
    await model.click();
    await page.waitForTimeout(1500);
    await model.hover();
    await page.waitForTimeout(800);

    await expectNoErrors(page, sink);
  });

  test('Browse-ModelHub-03 — double click on a model does not crash, right-click Run is reachable', async ({ page }) => {
    // Was blocked by platform regression GROK-19965 (double-click on a model threw a Dart
    // `TypeError: p.append is not a function`, same root cause as ModelHub-02). That fix has
    // landed — the test now passes, so the former `test.fail` annotation is removed.
    const sink = watchErrors(page);

    await expandTreeGroup(page, 'Apps');
    await expandTreeGroup(page, COMPUTE_PATH);
    await expandTreeGroup(page, MH_PATH);
    await treeNodeByPath(page, UNCAT_PATH).waitFor({ state: 'visible', timeout: SEED_LIST_TIMEOUT });
    await expandTreeGroup(page, UNCAT_PATH);

    const model = treeNodeByPath(page, SEED_PATH);
    await expect(model, 'The seeded model must be listed under Uncategorized').toBeVisible({ timeout: SEED_LIST_TIMEOUT });
    await model.dblclick();
    await page.waitForTimeout(2000);
    await expectNoErrors(page, sink);

    // Opening the model swaps the left pane from Browse to Toolbox.
    await ensureBrowsePanelOpen(page);
    const label = model.locator('.d4-tree-view-item-label, .d4-tree-view-group-label').first();
    await label.click({ button: 'right' });
    await expect(page.locator(CONTEXT_MENU)).toBeVisible({ timeout: 5_000 });
    const runItem = contextMenuItem(page, 'Run...');
    expect(await runItem.count(), 'A "Run..." item should be present in the seeded model context menu')
      .toBeGreaterThanOrEqual(1);

    await page.keyboard.press('Escape');
    await expectNoErrors(page, sink);
  });

  test('Browse-ModelHub-04 — Uncategorized models group expands without errors', async ({ page }) => {
    const sink = watchErrors(page);

    await expandTreeGroup(page, 'Apps');
    await expandTreeGroup(page, COMPUTE_PATH);
    await expandTreeGroup(page, MH_PATH);

    const uncat = treeNodeByPath(page, UNCAT_PATH);
    await uncat.waitFor({ state: 'visible', timeout: SEED_LIST_TIMEOUT });
    const tri = uncat.locator(TREE_EXPAND_ARROW).first();
    if (!(await tri.evaluate((el) => el.classList.contains('d4-tree-view-tri-expanded')).catch(() => false))) {
      await tri.click();
      await page.waitForTimeout(1500);
    }
    await expect(tri, 'Uncategorized arrow should be expanded (ref: GROK-19628)')
      .toHaveClass(/d4-tree-view-tri-expanded/);

    await expectNoErrors(page, sink);
  });

});
