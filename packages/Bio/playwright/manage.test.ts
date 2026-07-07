/* ---
sub_features_covered: [bio.detector, bio.manage.libraries-view]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';
test.use(specTestOptions);
test('Bio Manage Monomer Libraries (filter_HELM.csv)', async ({page}) => {
  test.setTimeout(120_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:AppData/Bio/tests/filter_HELM.csv');
    grok.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(() => resolve(), 4000);
    });
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    const hasMacromolecule = cols.some((c: any) => c.semType === 'Macromolecule');
    if (hasMacromolecule) {
      for (let i = 0; i < 60; i++) {
        const canvas = document.querySelector('[name="viewer-Grid"] canvas') as HTMLCanvasElement;
        if (canvas && canvas.width > 0) break;
        await new Promise((r) => setTimeout(r, 200));
      }
    }
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
  await page.waitForFunction(() => {
    return !!document.querySelector('[name="div-Bio"]');
  }, null, {timeout: 30_000});
  await softStep('Open filter_HELM.csv and detect Macromolecule semType (bio.detector)', async () => {
    const info = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macromolCol = cols.find((c: any) => c.semType === 'Macromolecule');
      return {
        rows: df.rowCount,
        hasMacromolecule: !!macromolCol,
        macromolColName: macromolCol?.name || null,
      };
    });
    expect(info.rows).toBeGreaterThan(0);
    expect(info.hasMacromolecule).toBe(true);
    expect(info.macromolColName).not.toBeNull();
  });
  await softStep('Dispatch Bio > Manage > Monomer Libraries top-menu (bio.manage.libraries-view)', async () => {
    await page.evaluate(() => (document.querySelector('[name="div-Bio"]') as HTMLElement).click());
    await page.waitForFunction(() => !!document.querySelector('[name="div-Bio---Manage"]'),
      null, {timeout: 15_000});
    await page.evaluate(() => {
      const manage = document.querySelector('[name="div-Bio---Manage"]')!;
      manage.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      manage.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
    });
    await page.waitForFunction(() => !!document.querySelector('[name="div-Bio---Manage---Monomer-Libraries"]'),
      null, {timeout: 15_000});
    await page.evaluate(() =>
      (document.querySelector('[name="div-Bio---Manage---Monomer-Libraries"]') as HTMLElement).click());
    await page.waitForFunction(() => (window as any).grok?.shell?.v?.name === 'Manage Monomer Libraries',
      null, {timeout: 30_000});
    const view = await page.evaluate(() => {
      const v = grok.shell.v;
      const root = v?.root;
      let hasDuplicateHeading = false;
      if (root) {
        const headings = root.querySelectorAll('h1, h2, h3, h4');
        for (const el of Array.from(headings)) {
          const t = (el.textContent || '').trim();
          if (t === 'Manage Duplicate Monomer Symbols') {
            hasDuplicateHeading = true;
            break;
          }
        }
      }
      return {
        name: v?.name || null,
        type: v?.type || null,
        rootClass: root?.className || null,
        hasGrokView: !!root?.closest('.grok-view') || (root?.className || '').includes('grok-view'),
        hasDuplicateHeading,
      };
    });
    expect(view.name).toBe('Manage Monomer Libraries');
    expect(view.type).toBe('view');
    expect(view.hasGrokView).toBe(true);
    expect(view.hasDuplicateHeading).toBe(true);
  });
  await softStep('View renders per-library checkbox listing (>=1 row) + checkbox element per row', async () => {
    const listing = await page.evaluate(() => {
      const root = grok.shell.v.root;
      const form = root.querySelector('.monomer-lib-controls-form');
      const libRows = form ? Array.from(form.querySelectorAll('.ui-input-bool')) : [];
      const checkboxes: HTMLInputElement[] = libRows
        .map((r) => r.querySelector('input[type="checkbox"]') as HTMLInputElement)
        .filter(Boolean);
      const labels = libRows.map((r) => {
        const span = r.querySelector('span');
        return span ? (span.textContent || '').trim() : '';
      });
      return {
        formPresent: !!form,
        rowCount: libRows.length,
        checkboxCount: checkboxes.length,
        labels,
        searchBoxPresent: !!root.querySelector('[name="input-host-Search"]'),
        eachRowHasCheckbox: libRows.length > 0 &&
          libRows.every((r) => !!r.querySelector('input[type="checkbox"]')),
      };
    });
    expect(listing.formPresent).toBe(true);
    expect(listing.rowCount).toBeGreaterThanOrEqual(1);
    expect(listing.checkboxCount).toBeGreaterThanOrEqual(1);
    expect(listing.searchBoxPresent).toBe(true);
    expect(listing.eachRowHasCheckbox).toBe(true);
  });
  await softStep('Toggle two libraries independently — settings persist, per-row isolation, no error', async () => {
    const rowCount = await page.evaluate(() => {
      const form = grok.shell.v.root.querySelector('.monomer-lib-controls-form');
      const cbs = (form ? Array.from(form.querySelectorAll('.ui-input-bool')) : [])
        .map((r) => r.querySelector('input[type="checkbox"]')).filter(Boolean);
      return cbs.length;
    });
    expect(rowCount, 'need >=2 library checkboxes to test per-row isolation').toBeGreaterThanOrEqual(2);

    const readState = () => page.evaluate(() => {
      const form = grok.shell.v.root.querySelector('.monomer-lib-controls-form')!;
      const cbs = Array.from(form.querySelectorAll('.ui-input-bool'))
        .map((r) => r.querySelector('input[type="checkbox"]') as HTMLInputElement).filter(Boolean);
      return {a: cbs[0].checked, b: cbs[1].checked};
    });
    const clickCb = (idx: number) => page.evaluate((i) => {
      document.querySelectorAll('.d4-balloon').forEach((b) => b.remove());
      const form = grok.shell.v.root.querySelector('.monomer-lib-controls-form')!;
      const cbs = Array.from(form.querySelectorAll('.ui-input-bool'))
        .map((r) => r.querySelector('input[type="checkbox"]') as HTMLInputElement).filter(Boolean);
      cbs[i].click();
    }, idx);
    // Anchor "toggle took effect" on the real downstream save, not just the native .checked flag.
    const waitSettingsSaved = () => page.waitForFunction(() =>
      Array.from(document.querySelectorAll('.d4-balloon'))
        .some((b) => /settings saved/i.test(b.textContent || '')), null, {timeout: 15_000});

    const s0 = await readState();
    await clickCb(0);
    await waitSettingsSaved();
    const s1 = await readState();
    expect(s1.a, 'checkbox A flips on click').toBe(!s0.a);
    expect(s1.b, 'checkbox B unchanged while A toggles (per-row isolation)').toBe(s0.b);

    await clickCb(1);
    await waitSettingsSaved();
    const s2 = await readState();
    expect(s2.b, 'checkbox B flips on click').toBe(!s1.b);
    expect(s2.a, 'checkbox A unchanged while B toggles (per-row isolation)').toBe(s1.a);

    await clickCb(0);
    await waitSettingsSaved();
    await clickCb(1);
    await waitSettingsSaved();
    const sR = await readState();
    expect(sR.a, 'checkbox A restored to initial state').toBe(s0.a);
    expect(sR.b, 'checkbox B restored to initial state').toBe(s0.b);

    const errorBalloons = await page.locator('.d4-balloon.error, .d4-balloon-error, .grok-balloon-error').count();
    expect(errorBalloons, 'no error balloon during library toggles (.md: toggle without error)').toBe(0);
  });
  finishSpec();
});
