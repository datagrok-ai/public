/* ---
sub_features_covered: [peptides.panels.peptides, peptides.rendering.weblogo-header, peptides.util.modify-selection, peptides.widgets.distribution]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';
test.use(specTestOptions);
const datasetPath = 'System:DemoFiles/bio/peptides.csv';
test('Peptides — SAR parameters and WebLogo', async ({page}) => {
  test.setTimeout(300_000);
  await loginToDatagrok(page);
  await softStep('Steps 1-2: Open dataset and select Macromolecule column', async () => {
    const result = await page.evaluate(async (path) => {
      document.querySelectorAll('.d4-dialog').forEach(d => {
        const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
        if (cancel) cancel.click();
      });
      grok.shell.closeAll();
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      grok.shell.windows.simpleMode = false;
      const df = await grok.dapi.files.readCsv(path);
      grok.shell.addTableView(df);
      await new Promise<void>(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 200));
      }
      await new Promise(r => setTimeout(r, 5000));
      const peptidesPanelFn = (DG.Func.find({name: 'peptidesPanel'}) || [])[0];
      const peptidesPkg = peptidesPanelFn?.package;
      if (peptidesPkg && typeof peptidesPkg.load === 'function')
        await peptidesPkg.load();
      const bioFns = DG.Func.find({package: 'Bio'}) || [];
      const bioPkg = bioFns[0]?.package;
      if (bioPkg && typeof bioPkg.load === 'function')
        await bioPkg.load();
      const col = df.col('AlignedSequence')!;
      df.currentCol = col;
      grok.shell.o = col;
      let paneMounted = false;
      const deadline = Date.now() + 40000;
      while (Date.now() < deadline) {
        const headers = document.querySelectorAll('.d4-accordion-pane-header').length;
        const pane = document.querySelector('[name="pane-Peptides"]');
        if (headers > 0 && pane && pane.querySelector('.d4-accordion-pane-header')) {
          paneMounted = true;
          break;
        }
        await new Promise(r => setTimeout(r, 500));
      }
      await new Promise(r => setTimeout(r, 2000));
      return {rows: df.rowCount, semType: df.col('AlignedSequence')?.semType, paneMounted};
    }, datasetPath);
    expect(result.rows).toBe(647);
    expect(result.semType).toBe('Macromolecule');
  });
  await softStep('Step 3: Expand Peptides context-panel pane', async () => {
    const pane = page.locator('[name="pane-Peptides"]');
    try {
      await pane.waitFor({state: 'attached', timeout: 15_000});
    } catch {
      await page.evaluate(async () => {
        const peptidesPanelFn = (DG.Func.find({name: 'peptidesPanel'}) || [])[0];
        const peptidesPkg = peptidesPanelFn?.package;
        if (peptidesPkg && typeof peptidesPkg.load === 'function')
          await peptidesPkg.load();
        const bioFns = DG.Func.find({package: 'Bio'}) || [];
        const bioPkg = bioFns[0]?.package;
        if (bioPkg && typeof bioPkg.load === 'function')
          await bioPkg.load();
        const df = grok.shell.tv?.dataFrame;
        const col = df?.col('AlignedSequence');
        if (!col) return;
        df.currentCol = col;
        grok.shell.o = col;
        const deadline = Date.now() + 40000;
        while (Date.now() < deadline) {
          if (document.querySelectorAll('.d4-accordion-pane-header').length > 0
            && document.querySelector('[name="pane-Peptides"] .d4-accordion-pane-header')) break;
          await new Promise(r => setTimeout(r, 500));
        }
        await new Promise(r => setTimeout(r, 2000));
      });
      await pane.waitFor({state: 'attached', timeout: 30_000});
    }
    await pane.waitFor({timeout: 30_000});
    const header = pane.locator('.d4-accordion-pane-header');
    const alreadyExpanded = await header.evaluate(el => el.classList.contains('expanded'));
    if (!alreadyExpanded) {
      await header.scrollIntoViewIfNeeded();
      await header.click();
    }
    await pane.locator('[name="input-Scaling"]').waitFor({timeout: 30_000});
    await pane.locator('.bio-wl-host canvas').waitFor({timeout: 30_000});
    const result = await pane.evaluate((el) => ({
      expanded: !!el.querySelector('.d4-accordion-pane-header')?.classList.contains('expanded'),
      hasWebLogo: !!el.querySelector('.bio-wl-host canvas'),
      hasParams: !!el.querySelector('[name="input-Scaling"]')
        && !!el.querySelector('[name="input-Activity"]')
        && !!el.querySelector('[name="input-Clusters"]'),
    }));
    expect(result.expanded).toBe(true);
    expect(result.hasWebLogo).toBe(true);
    expect(result.hasParams).toBe(true);
  });
  await softStep('Steps 4-5: Change Scaling parameter, panel re-renders', async () => {
    const pane = page.locator('[name="pane-Peptides"]');
    const scaling = pane.locator('[name="input-Scaling"]');
    await scaling.selectOption('lg');
    await page.waitForTimeout(1500);
    const result = await pane.evaluate((el) => {
      const sel = el.querySelector('[name="input-Scaling"]') as HTMLSelectElement | null;
      return {
        scalingChanged: sel?.value === 'lg',
        reRendered: !!el.querySelector('.bio-wl-host canvas') && !!el.querySelector('[name="viewer-Histogram"]'),
      };
    });
    expect(result.scalingChanged).toBe(true);
    expect(result.reRendered).toBe(true);
  });
  await softStep('Steps 6-9: WebLogo monomer click selects matching rows', async () => {
    await page.evaluate(async () => {
      const peptidesPanelFn = (DG.Func.find({name: 'peptidesPanel'}) || [])[0];
      const peptidesPkg = peptidesPanelFn?.package;
      if (peptidesPkg && typeof peptidesPkg.load === 'function')
        await peptidesPkg.load();
      const bioFns = DG.Func.find({package: 'Bio'}) || [];
      const bioPkg = bioFns[0]?.package;
      if (bioPkg && typeof bioPkg.load === 'function')
        await bioPkg.load();
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      const col = df.col('AlignedSequence')!;
      df.currentCol = col;
      grok.shell.o = col;
      const deadline = Date.now() + 40000;
      while (Date.now() < deadline) {
        if (document.querySelectorAll('.d4-accordion-pane-header').length > 0
          && document.querySelector('[name="pane-Peptides"] .d4-accordion-pane-header')) break;
        await new Promise(r => setTimeout(r, 500));
      }
      await new Promise(r => setTimeout(r, 2000));
    });
    const pane = page.locator('[name="pane-Peptides"]');
    await pane.waitFor({timeout: 30_000});
    const header = pane.locator('.d4-accordion-pane-header');
    const alreadyExpanded = await header.evaluate(el => el.classList.contains('expanded'));
    if (!alreadyExpanded)
      await header.click();
    const canvas = pane.locator('.bio-wl-host canvas');
    await canvas.waitFor({timeout: 30_000});
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const pane = document.querySelector('[name="pane-Peptides"]');
      const canvas = pane?.querySelector('.bio-wl-host canvas') as HTMLCanvasElement | null;
      if (!canvas) return {clicked: false, selectedRows: 0};
      const r = canvas.getBoundingClientRect();
      const cx = r.x + r.width * 0.40;
      const cy = r.y + r.height * 0.45;
      const selBefore = df.selection.trueCount;
      for (const type of ['pointerdown', 'mousedown', 'pointerup', 'mouseup', 'click']) {
        canvas.dispatchEvent(new MouseEvent(type, {
          bubbles: true, cancelable: true, clientX: cx, clientY: cy, button: 0, view: window,
        }));
      }
      await new Promise(r => setTimeout(r, 800));
      return {clicked: true, selBefore, selectedRows: df.selection.trueCount, total: df.rowCount};
    });
    expect(result.clicked).toBe(true);
    expect(result.selectedRows).toBeGreaterThan(0);
  });
  await page.evaluate(() => grok.shell.closeAll());
  finishSpec();
});
