import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
test.use(specTestOptions);
const datasetPath = 'System:DemoFiles/bio/peptides.csv';
test('Peptides — Info Panels', async ({page}) => {
  test.setTimeout(300_000);
  await loginToDatagrok(page);
  await softStep('Setup: open peptides.csv + ensure Context Panel wired', async () => {
    const result = await page.evaluate(async (path) => {
      document.querySelectorAll('.d4-dialog').forEach((d) => {
        const cancel = d.querySelector('[name="button-CANCEL"]');
        if (cancel) (cancel as HTMLElement).click();
      });
      grok.shell.closeAll();
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      grok.shell.windows.simpleMode = false;
      grok.shell.windows.showContextPanel = true;
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      try { (grok.shell as any).v = tv; } catch (_) {  }
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 5000));
      try { await grok.functions.call('Peptides:initPeptides'); }
      catch (e) { console.log('[note] Peptides:initPeptides pre-warm threw (non-fatal):', String(e)); }
      const wiredDeadline = Date.now() + 30_000;
      let detailsWired = false;
      while (Date.now() < wiredDeadline) {
        detailsWired = !!document.querySelector('.grok-prop-panel [name="pane-Details"]');
        if (detailsWired) break;
        await new Promise((r) => setTimeout(r, 250));
      }
      return {
        rows: df.rowCount,
        semType: df.col('AlignedSequence')?.semType,
        detailsWired,
        activeViewType: (grok.shell.v as any)?.type ?? (grok.shell.v as any)?.constructor?.name,
      };
    }, datasetPath);
    expect(result.rows).toBe(647);
    expect(result.semType).toBe('Macromolecule');
    expect(result.detailsWired).toBe(true);
    await page.locator('[name="viewer-Grid"]').waitFor({timeout: 30_000});
    await page.locator('.grok-prop-panel [name="pane-Details"]').waitFor({timeout: 30_000});
  });
  await softStep('Step 1: Verify amino acid coloring (cell.renderer === sequence)', async () => {
    const renderer = await page.evaluate(() => {
      const col = grok.shell.tv.dataFrame.col('AlignedSequence');
      return col?.getTag('cell.renderer');
    });
    expect(renderer).toBe('sequence');
  });
  await softStep('Step 2: Focus AlignedSequence column, wait for Context Panel rebuild', async () => {
    await page.locator('[name="viewer-Grid"]').waitFor({timeout: 15_000});
    await page.locator('.grok-prop-panel [name="pane-Details"]').waitFor({timeout: 15_000});
    const probe = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const col = df.col('AlignedSequence');
      df.currentCol = col;
      grok.shell.o = col;
      const deadline = Date.now() + 20_000;
      const retryAt = Date.now() + 10_000;
      let peptidesPresent = false;
      let detailsPresent = false;
      let retried = false;
      while (Date.now() < deadline) {
        detailsPresent = !!document.querySelector('.grok-prop-panel [name="pane-Details"]');
        peptidesPresent = !!document.querySelector('.grok-prop-panel [name="pane-Peptides"]');
        if (detailsPresent && peptidesPresent) break;
        if (!retried && Date.now() > retryAt) {
          grok.shell.o = col;
          retried = true;
        }
        await new Promise((r) => setTimeout(r, 250));
      }
      await new Promise((r) => setTimeout(r, 1500));
      return {
        currentCol: df.currentCol?.name,
        currentO: (grok.shell.o as any)?.name,
        detailsPresent,
        peptidesPresent,
        retried,
        paneCount: document.querySelectorAll('.grok-prop-panel .d4-accordion-pane-header').length,
      };
    });
    expect(probe.currentCol).toBe('AlignedSequence');
    expect(probe.detailsPresent).toBe(true);
    expect(probe.peptidesPresent).toBe(true);
    await page.locator('.grok-prop-panel [name="pane-Peptides"]').waitFor({timeout: 15_000});
  });
  await softStep('Step 3: Check Context Panel sections (Details + Peptides)', async () => {
    await page.locator('.grok-prop-panel [name="pane-Details"]').waitFor({timeout: 15_000});
    await page.locator('.grok-prop-panel [name="pane-Peptides"]').waitFor({timeout: 15_000});
    await expect(page.locator('.grok-prop-panel [name="pane-Details"]')).toHaveCount(1);
    await expect(page.locator('.grok-prop-panel [name="pane-Peptides"]')).toHaveCount(1);
  });
  await softStep('Step 4: Expand Details and Peptides panes', async () => {
    const detailsHeader = page.locator('.grok-prop-panel [name="pane-Details"] .d4-accordion-pane-header');
    const peptidesHeader = page.locator('.grok-prop-panel [name="pane-Peptides"] .d4-accordion-pane-header');
    await detailsHeader.waitFor({timeout: 15_000});
    await peptidesHeader.waitFor({timeout: 15_000});
    const detailsExpanded = await detailsHeader.evaluate((el) => el.classList.contains('expanded'));
    const peptidesExpanded = await peptidesHeader.evaluate((el) => el.classList.contains('expanded'));
    if (!detailsExpanded) await detailsHeader.click();
    if (!peptidesExpanded) await peptidesHeader.click();
    await page.waitForTimeout(2500);
    await expect(page.locator('.grok-prop-panel [name="pane-Details"] .d4-accordion-pane-header.expanded')).toHaveCount(1);
    await expect(page.locator('.grok-prop-panel [name="pane-Peptides"] .d4-accordion-pane-header.expanded')).toHaveCount(1);
  });
  await softStep('Step 5: Verify Details and Peptides expanded content', async () => {
    const detailsContent = page.locator('.grok-prop-panel [name="pane-Details"] .d4-accordion-pane-content');
    await detailsContent.waitFor({timeout: 15_000});
    const detailsText = (await detailsContent.innerText()) || '';
    expect(detailsText.length).toBeGreaterThan(0);
    expect(detailsText).toContain('Data type');
    expect(detailsText).toContain('Semantic type');
    const peptidesContent = page.locator('.grok-prop-panel [name="pane-Peptides"] .d4-accordion-pane-content');
    await peptidesContent.waitFor({timeout: 15_000});
    const peptidesText = (await peptidesContent.innerText()) || '';
    expect(peptidesText.length).toBeGreaterThan(0);
    await expect(page.locator('.grok-prop-panel [name="pane-Peptides"] [name="button-Launch-SAR"]')).toHaveCount(1);
    await expect(page.locator('.grok-prop-panel [name="pane-Peptides"] [name="input-host-Activity"]')).toHaveCount(1);
    expect(peptidesText.toUpperCase()).toContain('LAUNCH SAR');
    expect(peptidesText).toContain('Activity');
  });
  await page.evaluate(() => grok.shell.closeAll());
  finishSpec();
});
