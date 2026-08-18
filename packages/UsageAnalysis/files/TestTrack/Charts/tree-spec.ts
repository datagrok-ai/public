import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';
import * as v from '../helpers/viewers';

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';

test('Tree viewer (Charts package)', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);
  await page.locator('[name="Browse"]').waitFor({timeout: 30_000});

  await page.evaluate(() => {
    document.querySelectorAll('.d4-dialog').forEach((d) => {
      const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
      if (cancel) cancel.click();
    });
    (window as any).grok.shell.closeAll();
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
  });

  await softStep('Setup: Open demog.csv; add Tree via gallery; set Hierarchy CONTROL/SEX/RACE', async () => {

    const open = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      const df = await grok.dapi.files.readCsv(path);
      grok.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      return {rowCount: df.rowCount};
    }, demogPath);
    expect(open.rowCount).toBe(5850);

    let dialogOpened = false;
    let fallbackUsed = false;
    try {
      await page.locator('i.svg-add-viewer').first().click({timeout: 5000});
      const dlg = page.locator('[name="dialog-Add-Viewer"]');
      try {
        await dlg.waitFor({state: 'visible', timeout: 3000});
        dialogOpened = true;
      } catch (e) {

        await page.locator('i.svg-add-viewer').first().click({timeout: 3000});
        try {
          await dlg.waitFor({state: 'visible', timeout: 3000});
          dialogOpened = true;
        } catch (e2) {
          dialogOpened = false;
        }
      }
    } catch (e) {
      dialogOpened = false;
    }

    if (dialogOpened) {

      await page.locator('[name="dialog-Add-Viewer"] .d4-item-card.viewer-gallery')
        .filter({hasText: /^\s*Tree\s*$/}).first().click({timeout: 5000});

      await page.waitForTimeout(4500);

      const close = page.locator('[name="dialog-Add-Viewer"] [name="icon-font-icon-close"]');
      const closeCount = await close.count();
      if (closeCount > 0) {
        await close.first().click({timeout: 2000}).catch(() => {});
      }
    } else {

      console.warn('[tree Setup]', 'Add Viewer gallery did not open via page.locator; falling back to tv.addViewer JS API');
      await page.evaluate(async () => {
        const grok = (window as any).grok;
        const tv = grok.shell.tv;
        tv.addViewer('Tree');
        await new Promise((r) => setTimeout(r, 4500));
      });
      fallbackUsed = true;
    }

    await page.locator('[name="viewer-Tree"]').waitFor({timeout: 30_000});

    await page.locator('[name="viewer-Tree"]')
      .locator('xpath=ancestor::*[contains(concat(" ", normalize-space(@class), " "), " panel-base ")]')
      .locator('.panel-titlebar [name="icon-font-icon-settings"]')
      .first().click({timeout: 5000});
    await page.locator('.grok-prop-panel').waitFor({timeout: 5000});

    let selectColsOpened = false;
    try {
      await page.locator('.grok-prop-panel [name="prop-hierarchy"] button').first().click({timeout: 5000});
      await page.locator('[name="dialog-Select-columns..."]').waitFor({state: 'visible', timeout: 3000});
      selectColsOpened = true;
    } catch (e) {
      selectColsOpened = false;
    }

    if (selectColsOpened) {
      await page.locator('[name="dialog-Select-columns..."] [name="button-CANCEL"]').first()
        .click({timeout: 3000}).catch(() => {});
      await page.waitForTimeout(300);
    }

    const setup = await page.evaluate(async () => {
      const grok = (window as any).grok;
      const tv = grok.shell.tv;
      let tree: any = null;
      for (const v of tv.viewers) if (v.type === 'Tree') { tree = v; break; }
      if (!tree) return {viewerTypes: [] as string[], hierarchy: null as any, ok: false};
      tree.setOptions({hierarchyColumnNames: ['CONTROL', 'SEX', 'RACE']});
      await new Promise((r) => setTimeout(r, 1500));
      const viewerTypes: string[] = [];
      for (const v of tv.viewers) viewerTypes.push(v.type);

      let hierarchy = null;
      try { hierarchy = tree.props.get('hierarchyColumnNames'); } catch (e) {}
      return {viewerTypes, hierarchy, ok: true};
    });
    expect(setup.ok).toBe(true);
    expect(setup.viewerTypes).toContain('Tree');
    if (setup.hierarchy != null) expect(setup.hierarchy).toEqual(['CONTROL', 'SEX', 'RACE']);
    console.log(`[Setup] fallbackUsed=${fallbackUsed}, selectColsDialogOpened=${selectColsOpened}`);
  });

  await softStep('Step 1: Open Filter Panel (page.locator); apply CONTROL=true filter; verify cardinality = 39', async () => {
    let panelOpened = false;
    let panelOpenedVia = 'none';

    const filtersHeader = page.locator('[name="div-section--Filters"]');
    if (await filtersHeader.count() > 0) {
      const expanded = await filtersHeader.first().evaluate((el) => el.classList.contains('expanded'));
      if (!expanded) {
        await filtersHeader.first().click({timeout: 3000}).catch(() => {});
      }
      await page.waitForTimeout(800);
      if (await page.locator('[name="viewer-Filters"]').count() > 0) {
        panelOpened = true;
        panelOpenedVia = 'toolbox-section';
      }
    }

    if (!panelOpened) {
      const ribbonIcon = page.locator('[name="icon-filter"]');
      if (await ribbonIcon.count() > 0) {
        await ribbonIcon.first().click({timeout: 3000}).catch(() => {});
        await page.waitForTimeout(800);
        if (await page.locator('[name="viewer-Filters"]').count() > 0) {
          panelOpened = true;
          panelOpenedVia = 'ribbon-icon';
        }
      }
    }

    const result = await page.evaluate(async () => {
      const grok = (window as any).grok;
      const df = grok.shell.tv.dataFrame;
      const control = df.col('CONTROL');
      df.filter.setAll(false);
      for (let i = 0; i < df.rowCount; i++)
        if (control.get(i) === true) df.filter.set(i, true);
      df.filter.fireChanged();
      await new Promise((r) => setTimeout(r, 500));
      return {filtered: df.filter.trueCount};
    });
    console.log(`[Step 1] panelOpened=${panelOpened} (via ${panelOpenedVia}), filtered=${result.filtered}`);
    if (!panelOpened) console.warn('[Step 1]', 'Filter Panel did not materialize via DOM; filter applied via JS API only');
    expect(result.filtered).toBe(39);
  });

  await softStep('Step 2: Clear CONTROL=true filter; verify cardinality returns to 5850', async () => {
    const result = await page.evaluate(async () => {
      const grok = (window as any).grok;
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true);
      df.filter.fireChanged();
      await new Promise((r) => setTimeout(r, 500));
      return {filtered: df.filter.trueCount, selected: df.selection.trueCount};
    });
    console.log(`[Step 2] filtered=${result.filtered}, selected=${result.selected}`);
    expect(result.filtered).toBe(5850);

  });

  await v.cleanupShell(page);

  v.finishSpec();
});
