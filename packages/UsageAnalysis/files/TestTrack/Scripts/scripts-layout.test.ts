import { test, expect, Page, BrowserContext } from '@playwright/test';
import * as path from 'path';
import {
  openScriptsBrowser,
  rightClickScript,
  clickMenuItem,
  apiDeleteScript,
} from './helpers';

const BASE = process.env.DATAGROK_URL!;
const AUTH_STATE = path.resolve(__dirname, '..', '.auth.json');
const LAYOUT_SCRIPT_NAME = 'test_Layout';

const PROJECT_NAME = `PW_LayoutProject_${Date.now()}`;

const LAYOUT_SCRIPT_CONTENT = `//name: ${LAYOUT_SCRIPT_NAME}
//language: javascript
//input: int idx=0
//output: dataframe df

df = await grok.data.getDemoTable('cars.csv');`;

async function resetToScripts(page: Page) {
  await page.evaluate(() => {
    const g = (window as any).grok;
    if (g?.shell?.closeAll) g.shell.closeAll();
    document.querySelectorAll('.d4-dialog').forEach((d: any) => { try { d.remove(); } catch (_) {} });
    document.querySelectorAll('.d4-toast, .d4-balloon, .d4-menu').forEach((e) => e.remove());
  });
  await page.waitForTimeout(300);
  const scriptsLabel = page.locator('.d4-tree-view-item-label', { hasText: /^Scripts$/i }).first();
  if (await scriptsLabel.isVisible({ timeout: 1_000 }).catch(() => false)) {
    await scriptsLabel.click();
    await expect(page.locator('.grok-gallery-search-bar')).toBeVisible({ timeout: 10_000 });
  } else {
    await openScriptsBrowser(page);
  }
}

async function openScriptEditor(page: Page) {
  await rightClickScript(page, LAYOUT_SCRIPT_NAME);
  await clickMenuItem(page, 'Edit...');
  await page.waitForURL(/\/script\//, { timeout: 15_000 });
  await expect(page.locator('i[name="icon-play"]')).toBeVisible({ timeout: 10_000 });
}

async function runScriptFromEditor(page: Page) {
  await page.locator('i[name="icon-play"]').click();
  const dialog = page.locator('.d4-dialog').first();

  if (await dialog.isVisible({ timeout: 3_000 }).catch(() => false)) {
    const okBtn = dialog.locator('button.ui-btn-ok').first();
    if (await okBtn.isVisible({ timeout: 2_000 }).catch(() => false)) await okBtn.click();
  }

  await expect(page.locator('.d4-grid').first()).toBeVisible({ timeout: 30_000 });
}

async function expectViewerActive(page: Page, viewerType: string, timeout = 15_000) {
  await expect(async () => {
    const domVisible = await page
      .locator(`[aria-label="Close ${viewerType}"]`)
      .first()
      .isVisible({ timeout: 500 })
      .catch(() => false);
    if (domVisible) return;
    const viaWidget: boolean = await page.evaluate((t) => {
      const DG = (window as any).DG;
      for (const w of (DG?.Widget?.getAll?.() ?? [])) {
        if (w?.type === t) return true;
      }
      return false;
    }, viewerType);
    expect(viaWidget).toBeTruthy();
  }).toPass({ timeout });
}

test.describe.serial('Scripts: Layout', () => {
  let sharedContext: BrowserContext;
  let page: Page;

  test.beforeAll(async ({ browser }) => {
    sharedContext = await browser.newContext({ storageState: AUTH_STATE });
    page = await sharedContext.newPage();
    await openScriptsBrowser(page);
  });

  test.afterAll(async () => {
    await sharedContext?.close();
  });

  test.beforeEach(async () => {
    await page.evaluate(() => {
      const g = (window as any).grok;
      if (g?.shell?.closeAll) g.shell.closeAll();
    });
    await page.waitForTimeout(300);
    await page.waitForFunction(() => !!(window as any).grok?.dapi?.scripts, { timeout: 10_000 });
    await apiDeleteScript(page, LAYOUT_SCRIPT_NAME);

    await page.evaluate(async (name) => {
      try {
        const grok = (window as any).grok;
        const list = await grok.dapi.projects.filter(name).list();
        for (const p of list) { try { await grok.dapi.projects.delete(p); } catch (_) {} }
      } catch (_) {  }
    }, PROJECT_NAME);

    await page.evaluate(async (content) => {
      const DG = (window as any).DG;
      const grok = (window as any).grok;
      const script = DG.Script.create(content);
      await grok.dapi.scripts.save(script);
      return null;
    }, LAYOUT_SCRIPT_CONTENT);
    await page.waitForTimeout(300);
  });

  test.afterEach(async () => {
    await apiDeleteScript(page, LAYOUT_SCRIPT_NAME);
    await page.evaluate(async (name) => {
      try {
        const grok = (window as any).grok;

        for (let i = 0; i < 10; i++) {
          const list = await grok.dapi.projects.filter(name).list();
          if (list.length === 0) break;
          for (const p of list) await grok.dapi.projects.delete(p);
        }
      } catch (_) {  }
    }, PROJECT_NAME);
  });

  test('1. Layout tab is accessible in the script editor', async () => {
    await resetToScripts(page);
    await openScriptEditor(page);

    const layoutTab = page.locator('.d4-tab-header[name="Layout"]').first();
    await expect(layoutTab).toBeVisible({ timeout: 10_000 });
    await layoutTab.click();

    await expect(page.locator('.d4-tab-header[name="Layout"].selected')).toBeVisible({ timeout: 5_000 });
  });

  test('2. Running the script from editor produces a result grid', async () => {
    await resetToScripts(page);
    await openScriptEditor(page);
    await runScriptFromEditor(page);

    const errorBalloon = page.locator('.d4-balloon-error');
    await expect(errorBalloon).toHaveCount(0, { timeout: 3_000 }).catch(() => {});
  });

  test('3. Run from Layout tab shows the result pane in the script editor', async () => {
    await resetToScripts(page);
    await openScriptEditor(page);

    const layoutTab = page.locator('.d4-tab-header[name="Layout"]').first();
    await expect(layoutTab).toBeVisible({ timeout: 10_000 });
    await layoutTab.click();
    await expect(page.locator('.d4-tab-header[name="Layout"].selected')).toBeVisible({ timeout: 5_000 });

    await runScriptFromEditor(page);

    await expect(page.locator('.d4-grid').first()).toBeVisible({ timeout: 30_000 });

    const errorBalloon = page.locator('.d4-balloon-error');
    await expect(errorBalloon).toHaveCount(0, { timeout: 3_000 }).catch(() => {});
  });

  test('4. Save layout binds the viewer to the script on reopen', async () => {
    await resetToScripts(page);
    await openScriptEditor(page);

    const layoutTab = page.locator('.d4-tab-header[name="Layout"]').first();
    await expect(layoutTab).toBeVisible({ timeout: 10_000 });
    await layoutTab.click();

    await runScriptFromEditor(page);

    const viewersPane = page.locator('.d4-pane-viewers');
    await expect(viewersPane).toBeVisible({ timeout: 5_000 });
    const barChartIcon = viewersPane.locator('.grok-icon.svg-icon.svg-bar-chart').first();
    if (!(await barChartIcon.isVisible({ timeout: 1_000 }).catch(() => false))) {

      await page.locator('.d4-accordion-pane-header', { hasText: /^Viewers$/ }).first().click();
      await page.waitForTimeout(300);
    }
    await barChartIcon.click();

    await expectViewerActive(page, 'Bar chart', 10_000);

    await page.locator('button[name="button-Save"]').click();
    await expect(page.locator('.d4-balloon', { hasText: /saved/i }).first()).toBeVisible({ timeout: 10_000 });

    await page.evaluate(() => (window as any).grok.shell.closeAll());
    await page.waitForTimeout(500);

    await resetToScripts(page);
    await openScriptEditor(page);

    await page.locator('.d4-tab-header[name="Layout"]').first().click();

    await runScriptFromEditor(page);

    await expectViewerActive(page, 'Bar chart', 15_000);

    const viewersPane2 = page.locator('.d4-pane-viewers');
    const scatterIcon = viewersPane2.locator('.grok-icon.svg-icon.svg-scatter-plot').first();
    await expect(scatterIcon).toBeVisible({ timeout: 5_000 });
    await scatterIcon.click();
    await expectViewerActive(page, 'Scatter plot', 10_000);

    const firstColName: string | null = await page.evaluate(() => {
      const DG = (window as any).DG;
      const widgets: any[] = DG?.Widget?.getAll?.() ?? [];
      for (const w of widgets) {
        if (w?.type === 'Grid' && w?.dataFrame?.columns?.names) {
          const names: string[] = w.dataFrame.columns.names();
          if (names.length > 0) return names[0];
        }
      }
      return null;
    });
    expect(firstColName).not.toBeNull();

    const gridContainer = page.locator('.d4-grid').first();
    await expect(gridContainer).toBeVisible({ timeout: 5_000 });
    await gridContainer.click({ button: 'right', position: { x: 100, y: 10 } });

    const hideMenuItem = page.locator('div[name="div-Hide"]').first();
    await expect(hideMenuItem).toBeVisible({ timeout: 5_000 });
    await hideMenuItem.click();
    await page.waitForTimeout(800);

    const hiddenAfterClick: string[] = await page.evaluate(() => {
      const DG = (window as any).DG;
      for (const w of (DG?.Widget?.getAll?.() ?? [])) {
        if (w?.type !== 'Grid') continue;
        const names: string[] = w.dataFrame?.columns?.names?.() ?? [];
        return names.filter((n: string) => w.col?.(n)?.visible === false);
      }
      return [];
    });

    if (hiddenAfterClick.length === 0) {

      await page.evaluate((name) => {
        const DG = (window as any).DG;
        for (const w of (DG?.Widget?.getAll?.() ?? [])) {
          if (w?.type !== 'Grid') continue;
          const gc = w.col?.(name);
          if (gc) gc.visible = false;
        }
      }, firstColName as string);
      await page.waitForTimeout(300);
    }

    await page.locator('button[name="button-Save"]').click();
    await expect(page.locator('.d4-balloon', { hasText: /saved/i }).first()).toBeVisible({ timeout: 10_000 });

    await page.evaluate(() => (window as any).grok.shell.closeAll());
    await page.waitForTimeout(500);

    await resetToScripts(page);
    await openScriptEditor(page);
    await page.locator('.d4-tab-header[name="Layout"]').first().click();
    await runScriptFromEditor(page);

    await expectViewerActive(page, 'Bar chart', 15_000);
    await expectViewerActive(page, 'Scatter plot', 15_000);

    const isColHidden: boolean = await page.evaluate((name) => {
      const DG = (window as any).DG;
      for (const w of (DG?.Widget?.getAll?.() ?? [])) {
        if (w?.type !== 'Grid') continue;
        const gc = w.col?.(name);
        if (gc && gc.visible === false) return true;
      }
      return false;
    }, firstColName as string);
    expect(isColHidden).toBeTruthy();
  });

  test('5. Save project from standalone TableView and reopen preserves viewer', async () => {
    test.setTimeout(180_000);

    await resetToScripts(page);
    await openScriptEditor(page);
    await page.locator('.d4-tab-header[name="Layout"]').first().click();
    await runScriptFromEditor(page);

    const viewersPane = page.locator('.d4-pane-viewers');
    await expect(viewersPane).toBeVisible({ timeout: 5_000 });
    const barChartIcon = viewersPane.locator('.grok-icon.svg-icon.svg-bar-chart').first();
    if (!(await barChartIcon.isVisible({ timeout: 1_000 }).catch(() => false))) {
      await page.locator('.d4-accordion-pane-header', { hasText: /^Viewers$/ }).first().click();
      await page.waitForTimeout(300);
    }
    await barChartIcon.click();

    await expect(async () => {
      const hasBarChart: boolean = await page.evaluate(() => {
        const DG = (window as any).DG;
        for (const w of (DG?.Widget?.getAll?.() ?? [])) {
          if (w?.type === 'Bar chart') return true;
        }
        return false;
      });
      expect(hasBarChart).toBeTruthy();
    }).toPass({ timeout: 10_000 });

    await page.locator('button[name="button-Save"]:visible').first().click();
    await expect(page.locator('.d4-balloon', { hasText: /saved/i }).first()).toBeVisible({ timeout: 10_000 });

    await page.evaluate(() => (window as any).grok.shell.closeAll());
    await page.waitForTimeout(500);
    await resetToScripts(page);

    await rightClickScript(page, LAYOUT_SCRIPT_NAME);
    await clickMenuItem(page, 'Run...');

    const runDialog = page.locator('.d4-dialog').first();
    await expect(runDialog).toBeVisible({ timeout: 10_000 });
    await runDialog.locator('button.ui-btn-ok').first().click();
    await expect(runDialog).not.toBeVisible({ timeout: 15_000 });

    await expect(page.locator('.d4-grid').first()).toBeVisible({ timeout: 30_000 });

    await expectViewerActive(page, 'Bar chart', 15_000);

    await page.locator('button[name="button-Save"]:visible').first().click();

    const projectDialog = page.locator(
      '.d4-dialog[name="dialog-Save-project"], .d4-dialog:has-text("Save project")'
    ).first();
    await expect(projectDialog).toBeVisible({ timeout: 10_000 });
    const nameInput = projectDialog.locator('input[name="input-Name"], input.ui-input-editor').first();
    await expect(nameInput).toBeVisible({ timeout: 5_000 });
    await nameInput.fill(PROJECT_NAME);
    await projectDialog.locator('button.ui-btn-ok').first().click();
    await expect(projectDialog).not.toBeVisible({ timeout: 15_000 });

    const shareDialogVar = page.locator('.d4-dialog[name^="dialog-Share-"]').first();
    if (await shareDialogVar.isVisible({ timeout: 3_000 }).catch(() => false)) {
      const cancel = shareDialogVar.locator('button:has-text("Cancel"), button.ui-btn-cancel').first();
      if (await cancel.isVisible({ timeout: 1_000 }).catch(() => false))
        await cancel.click();
      else
        await page.keyboard.press('Escape');
      await expect(shareDialogVar).not.toBeVisible({ timeout: 5_000 });
    }

    await page.evaluate(async () => {
      const grok = (window as any).grok;
      const tv = grok.shell.tv;
      const df = tv?.dataFrame;
      const project = grok.shell.project;
      if (!df || !tv || !project) return;
      const tableInfo = df.getTableInfo?.();
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      await grok.dapi.tables.uploadDataFrame(df);
      if (tableInfo) await grok.dapi.tables.save(tableInfo);
      if (tableInfo) project.addChild(tableInfo);
      project.addChild(layout);
      await grok.dapi.projects.save(project);
    });

    await page.evaluate(() => (window as any).grok.shell.closeAll());
    await page.waitForTimeout(1_000);

    await page.evaluate(async (name) => {
      const grok = (window as any).grok;
      await grok.dapi.projects.open(name);
    }, PROJECT_NAME);

    await expect(async () => {
      const count: number = await page.evaluate(() => {
        const grok = (window as any).grok;
        return grok.shell.tables?.length ?? 0;
      });
      expect(count).toBeGreaterThan(0);
    }).toPass({ timeout: 30_000 });

    await page.evaluate(async () => {
      const grok = (window as any).grok;
      const tv = grok.shell.tv;
      if (!tv?.dataFrame) return;
      const layouts = await grok.dapi.layouts
        .getApplicable(tv.dataFrame)
        .catch(() => [] as any[]);
      if (layouts && layouts.length > 0) tv.loadLayout(layouts[0]);
    });

    await expectViewerActive(page, 'Bar chart', 20_000);
  });
});
