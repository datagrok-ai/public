/* ---
sub_features_covered:
  - bio.manage.libraries-view
  - bio.detector
--- */
//     curated bug-library entries)
// view-name assertion → checkbox enumeration → toggle-and-restore.
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
test.use(specTestOptions);
test('Bio Manage Monomer Libraries (filter_HELM.csv)', async ({page}) => {
  test.setTimeout(300_000);
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
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 5000));
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
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 500));
      const manage = document.querySelector('[name="div-Bio---Manage"]')!;
      manage.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      manage.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 500));
      (document.querySelector('[name="div-Bio---Manage---Monomer-Libraries"]') as HTMLElement).click();
    });
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
    expect(view.hasDuplicateHeading).toBe(true);
  });
  // SCOPE_REDUCTION (within-slice): scenario step 3 asks to assert a
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
        eachRowHasCheckbox: libRows.length > 0 &&
          libRows.every((r) => !!r.querySelector('input[type="checkbox"]')),
      };
    });
    expect(listing.formPresent).toBe(true);
    expect(listing.rowCount).toBeGreaterThanOrEqual(1);
    expect(listing.checkboxCount).toBeGreaterThanOrEqual(1);
    expect(listing.eachRowHasCheckbox).toBe(true);
  });
  await softStep('Toggle two libraries independently — per-row state isolation holds', async () => {
    // the first). After the test, restore the original states so the
    const result = await page.evaluate(async () => {
      const root = grok.shell.v.root;
      const form = root.querySelector('.monomer-lib-controls-form');
      const rows = form ? Array.from(form.querySelectorAll('.ui-input-bool')) : [];
      const cbs = rows
        .map((r) => r.querySelector('input[type="checkbox"]') as HTMLInputElement)
        .filter(Boolean);
      if (cbs.length < 2) {
        return {
          haveTwoRows: false,
          rowCount: cbs.length,
        };
      }
      const idxA = 0;
      const idxB = 1;
      const a0 = cbs[idxA].checked;
      const b0 = cbs[idxB].checked;
      cbs[idxA].click();
      await new Promise((r) => setTimeout(r, 2500));
      const cbs2 = (Array.from(form!.querySelectorAll('.ui-input-bool')) as HTMLElement[])
        .map((r) => r.querySelector('input[type="checkbox"]') as HTMLInputElement)
        .filter(Boolean);
      const a1 = cbs2[idxA].checked;
      const b1 = cbs2[idxB].checked;
      cbs2[idxB].click();
      await new Promise((r) => setTimeout(r, 2500));
      const cbs3 = (Array.from(form!.querySelectorAll('.ui-input-bool')) as HTMLElement[])
        .map((r) => r.querySelector('input[type="checkbox"]') as HTMLInputElement)
        .filter(Boolean);
      const a2 = cbs3[idxA].checked;
      const b2 = cbs3[idxB].checked;
      // Restore: toggle A back, toggle B back. Order chosen so both end
      cbs3[idxA].click();
      await new Promise((r) => setTimeout(r, 2000));
      cbs3[idxB].click();
      await new Promise((r) => setTimeout(r, 2500));
      const cbs4 = (Array.from(form!.querySelectorAll('.ui-input-bool')) as HTMLElement[])
        .map((r) => r.querySelector('input[type="checkbox"]') as HTMLInputElement)
        .filter(Boolean);
      const aR = cbs4[idxA].checked;
      const bR = cbs4[idxB].checked;
      return {
        haveTwoRows: true,
        rowCount: cbs.length,
        a0, b0,
        a1, b1,
        a2, b2,
        aR, bR,
        aToggledOnFirstClick: a1 !== a0,
        bUnchangedAfterAToggle: b1 === b0,
        bToggledOnSecondClick: b2 !== b1,
        aUnchangedAfterBToggle: a2 === a1,
        aRestored: aR === a0,
        bRestored: bR === b0,
      };
    });
    expect(result.haveTwoRows).toBe(true);
    expect(result.aToggledOnFirstClick).toBe(true);
    expect(result.bUnchangedAfterAToggle).toBe(true);
    expect(result.bToggledOnSecondClick).toBe(true);
    expect(result.aUnchangedAfterBToggle).toBe(true);
    // Idempotent restore so subsequent runs/users don't inherit toggled state.
    expect(result.aRestored).toBe(true);
    expect(result.bRestored).toBe(true);
  });
  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
