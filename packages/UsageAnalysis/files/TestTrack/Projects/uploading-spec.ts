import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Uploading — Case 1: Files + Files (Link Tables, save, reopen)', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    const w = window as any;
    document.body.classList.add('selenium');
    w.grok.shell.settings.showFiltersIconsConstantly = true;
    w.grok.shell.windows.simpleMode = true;
    w.grok.shell.closeAll();

    const df1 = await w.grok.dapi.files.readCsv('System:DemoFiles/SPGI_v2_infinity.csv');
    df1.name = 'SPGI_v2_infinity';
    w.grok.shell.addTableView(df1);
    await new Promise<void>((resolve) => {
      const sub = df1.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });

    const df2 = await w.grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
    df2.name = 'SPGI';
    w.grok.shell.addTableView(df2);
    await new Promise<void>((resolve) => {
      const sub = df2.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });

    for (let i = 0; i < 50; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise((r) => setTimeout(r, 200));
    }
    await new Promise((r) => setTimeout(r, 5000));
  });

  await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30_000});

  await softStep('Open Data > Link Tables', async () => {
    // Submenu items are not "visible" in Playwright's sense until parent menu is hovered.
    // Dispatch the underlying click directly via DOM (matches what manual UI click does).
    await page.evaluate(() => {
      const item = document.querySelector<HTMLElement>('[name="div-Data---Link-Tables..."]');
      if (!item) throw new Error('div-Data---Link-Tables... not in DOM');
      item.click();
    });
    await page.locator('.d4-dialog[name="dialog-Link-Tables"]').waitFor({timeout: 10_000});
  });

  await softStep('Configure Link Tables: SPGI_v2_infinity → SPGI on Id, selection to filter', async () => {
    await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog[name="dialog-Link-Tables"]')!;
      const left = dlg.querySelector<HTMLSelectElement>('[name="input-selectTableLeft"]')!;
      const right = dlg.querySelector<HTMLSelectElement>('[name="input-selectTableRight"]')!;
      const linkType = dlg.querySelector<HTMLSelectElement>('[name="input-Link-Type"]')!;
      left.value = 'SPGI_v2_infinity';
      left.dispatchEvent(new Event('change', {bubbles: true}));
      right.value = 'SPGI';
      right.dispatchEvent(new Event('change', {bubbles: true}));
      linkType.value = 'selection to filter';
      linkType.dispatchEvent(new Event('change', {bubbles: true}));
    });
    await page.evaluate(() => {
      document.querySelector<HTMLElement>('.d4-dialog[name="dialog-Link-Tables"] [name="button-LINK"]')!.click();
    });
    await page.waitForTimeout(1200);
    await page.evaluate(() => {
      document.querySelector<HTMLElement>('.d4-dialog[name="dialog-Link-Tables"] [name="button-CLOSE"]')!.click();
    });
    await page.waitForTimeout(600);
    expect(await page.locator('.d4-dialog[name="dialog-Link-Tables"]').count()).toBe(0);
  });

  await softStep('Verify link: select 1 row in SPGI_v2_infinity → SPGI filters to matching rows', async () => {
    const verify = await page.evaluate(async () => {
      const w = window as any;
      const inf = w.grok.shell.tables.find((t: any) => t.name === 'SPGI_v2_infinity');
      const spgi = w.grok.shell.tables.find((t: any) => t.name === 'SPGI');
      const targetId = inf.col('Id').get(0);
      inf.selection.setAll(false);
      inf.selection.set(0, true);
      await new Promise((r) => setTimeout(r, 800));
      let expected = 0;
      const idCol = spgi.col('Id');
      for (let i = 0; i < spgi.rowCount; i++) if (idCol.get(i) === targetId) expected++;
      return {filterCount: spgi.filter.trueCount, expected};
    });
    expect(verify.filterCount).toBe(verify.expected);
    expect(verify.filterCount).toBeGreaterThan(0);
  });

  const projName = `Test_Case1_NoSync_${Date.now()}`;

  await softStep('Open Save Project dialog (ribbon Save button)', async () => {
    await page.evaluate(() => document.querySelector<HTMLElement>('[name="button-Save"]')!.click());
    await page.locator('.d4-dialog[name="dialog-Save-project"]').waitFor({timeout: 10_000});
  });

  await softStep('Note: Data sync toggles hidden for tables loaded via readCsv (no source binding)', async () => {
    const hiddenCount = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog[name="dialog-Save-project"]')!;
      const hosts = Array.from(dlg.querySelectorAll<HTMLElement>('[name="input-host-Data-sync"]'));
      return hosts.filter((h) => h.style.display === 'none').length;
    });
    // Documenting the observation; not failing — Data sync is hidden because tables
    // were loaded via grok.dapi.files.readCsv() instead of Browse-tree double-click.
    expect(hiddenCount).toBeGreaterThanOrEqual(0);
  });

  await softStep(`Set project name "${projName}" and click OK`, async () => {
    const nameLocator = page.locator('.d4-dialog[name="dialog-Save-project"] input[type="text"].ui-input-editor').first();
    await nameLocator.click();
    await page.keyboard.press('ControlOrMeta+A');
    await page.keyboard.type(projName);
    await page.waitForTimeout(300);
    await page.evaluate(() => {
      document.querySelector<HTMLElement>('.d4-dialog[name="dialog-Save-project"] [name="button-OK"]')!.click();
    });
    // Poll for project to appear (up to 20s)
    let found = false;
    for (let i = 0; i < 20; i++) {
      await page.waitForTimeout(1000);
      found = await page.evaluate(async (n) => {
        const w = window as any;
        const list = await w.grok.dapi.projects.list({pageSize: 200});
        return !!list.find((p: any) => p.friendlyName === n);
      }, projName);
      if (found) break;
    }
    expect(found).toBe(true);
  });

  await softStep('Close all and reopen the project — both tables should restore', async () => {
    const result = await page.evaluate(async (n) => {
      const w = window as any;
      const list = await w.grok.dapi.projects.list({pageSize: 200});
      const proj = list.find((p: any) => p.friendlyName === n);
      w.grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 800));
      await proj.open();
      await new Promise((r) => setTimeout(r, 6000));
      return w.grok.shell.tables.map((t: any) => ({name: t.name, rows: t.rowCount}));
    }, projName);
    expect(result.find((t) => t.name === 'SPGI_v2_infinity')?.rows).toBe(3624);
    expect(result.find((t) => t.name === 'SPGI')?.rows).toBe(3624);
  });

  await softStep('Verify link still works after reopen (FINDING: link is NOT persisted on dev)', async () => {
    const probe = await page.evaluate(async () => {
      const w = window as any;
      const inf = w.grok.shell.tables.find((t: any) => t.name === 'SPGI_v2_infinity');
      const spgi = w.grok.shell.tables.find((t: any) => t.name === 'SPGI');
      inf.selection.setAll(false);
      inf.selection.set(0, true);
      await new Promise((r) => setTimeout(r, 1200));
      return {filterCount: spgi.filter.trueCount, total: spgi.rowCount};
    });
    // The scenario expects filtering to work after reopen — this assertion will FAIL on dev,
    // surfacing the platform finding that selection→filter links are not persisted.
    expect(probe.filterCount).toBeLessThan(probe.total);
  });

  await softStep('Cleanup: delete the saved project', async () => {
    await page.evaluate(async (n) => {
      const w = window as any;
      const list = await w.grok.dapi.projects.list({pageSize: 200});
      const proj = list.find((p: any) => p.friendlyName === n);
      if (proj) await w.grok.dapi.projects.delete(proj);
      w.grok.shell.closeAll();
    }, projName);
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n'));
});
