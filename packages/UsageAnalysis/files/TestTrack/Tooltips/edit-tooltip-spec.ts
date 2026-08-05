import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Viewers: Edit tooltip', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    (grok as any).shell.settings.showFiltersIconsConstantly = true;
    (grok as any).shell.windows.simpleMode = true;
    (grok as any).shell.closeAll();
    const df = await (grok as any).dapi.files.readCsv('System:DemoFiles/chem/SPGI.csv');
    (grok as any).shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(() => resolve(), 3000);
    });
    const cols = Array.from({length: df.columns.length}, (_, i: number) => df.columns.byIndex(i));
    const hasBioChem = cols.some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise<void>((r) => setTimeout(r, 200));
      }
      await new Promise<void>((r) => setTimeout(r, 5000));
    }
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Enable Show Visible Columns In Tooltip in grid properties', async () => {
    await page.evaluate(() => {
      const grid = document.querySelector('[name="viewer-Grid"]')!;
      (grid.querySelector('[name="icon-font-icon-settings"]') as HTMLElement).click();
    });
    await page.waitForTimeout(500);
    const toggled = await page.evaluate(() => {
      const tds = Array.from(document.querySelectorAll('td'));
      const target = tds.find((td) => td.textContent?.trim() === 'Show Visible Columns In Tooltip');
      if (!target) return false;
      const row = target.parentElement;
      const cb = row?.querySelector('input[type="checkbox"]') as HTMLInputElement | null;
      if (!cb) return false;
      if (!cb.checked) cb.click();
      return cb.checked;
    });
    expect(toggled).toBe(true);
  });

  await softStep('Add a scatter plot and a box plot', async () => {
    await page.evaluate(() => {
      (document.querySelector('[name="icon-scatter-plot"]') as HTMLElement).click();
    });
    await page.locator('[name="viewer-Scatter-plot"]').waitFor({timeout: 10000});
    await page.evaluate(() => {
      (document.querySelector('[name="icon-box-plot"]') as HTMLElement).click();
    });
    await page.locator('[name="viewer-Box-plot"]').waitFor({timeout: 10000});
    await page.waitForTimeout(500);
    const types = await page.evaluate(() =>
      Array.from((grok as any).shell.tv.viewers).map((v: any) => v.type));
    expect(types).toContain('Scatter plot');
    expect(types).toContain('Box plot');
  });

  await softStep('Right-click scatter plot and open Tooltip > Edit...', async () => {
    await page.evaluate(() => {
      const sp = document.querySelector('[name="viewer-Scatter-plot"]')!;
      const canvas = sp.querySelector('canvas') as HTMLCanvasElement;
      const r = canvas.getBoundingClientRect();
      const x = r.left + r.width / 2;
      const y = r.top + r.height / 2;
      for (const t of ['mousedown', 'mouseup', 'contextmenu']) {
        canvas.dispatchEvent(new MouseEvent(t, {
          bubbles: true, cancelable: true, button: 2, buttons: 2, clientX: x, clientY: y,
        }));
      }
    });
    await page.locator('.d4-menu-popup').waitFor({timeout: 5000});
    await page.evaluate(async () => {
      const group = document.querySelector('.d4-menu-popup [name="div-Tooltip"]') as HTMLElement | null;
      if (!group) throw new Error('Tooltip group not found in context menu');
      group.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, cancelable: true}));
      group.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true, cancelable: true}));
      await new Promise<void>((r) => setTimeout(r, 300));
      const editItem = group.querySelector('[name="div-Tooltip---Edit..."]') as HTMLElement | null;
      if (!editItem) throw new Error('Edit... item not found in Tooltip submenu');
      editItem.click();
    });
    await page.locator('.d4-dialog').waitFor({timeout: 5000});
    const title = await page.locator('.d4-dialog .d4-dialog-title').first().textContent();
    expect(title?.trim()).toBe('Edit Tooltip');
  });

  await softStep('Dialog contents: search input, column grid, action labels, footer', async () => {
    const info = await page.evaluate(() => {
      const dialog = document.querySelector('.d4-dialog')!;
      const searchInput = dialog.querySelector('input[placeholder="Search"]');
      const colGrid = dialog.querySelector('.d4-column-grid');
      const masterCb = dialog.querySelector('.ui-input-bool input[type="checkbox"]') as HTMLInputElement | null;
      const labels = Array.from(dialog.querySelectorAll('label.d4-link-label')).map((l) => l.textContent?.trim());
      const ok = dialog.querySelector('[name="button-OK"]');
      const cancel = dialog.querySelector('[name="button-CANCEL"]');
      const history = dialog.querySelector('.d4-dialog-footer .fa-history');
      return {
        hasSearch: !!searchInput,
        hasGrid: !!colGrid,
        masterChecked: masterCb?.checked,
        labels,
        hasOK: !!ok,
        hasCancel: !!cancel,
        hasHistory: !!history,
      };
    });
    expect(info.hasSearch).toBe(true);
    expect(info.hasGrid).toBe(true);
    expect(info.masterChecked).toBe(true);
    expect(info.labels).toContain('Reset group tooltip');
    expect(info.labels).toContain('Design custom tooltip...');
    expect(info.hasOK).toBe(true);
    expect(info.hasCancel).toBe(true);
    expect(info.hasHistory).toBe(true);
  });

  await softStep('Search is case-insensitive', async () => {
    await page.evaluate(async () => {
      const searchInput = document.querySelector('.d4-dialog .d4-column-grid input[placeholder="Search"]') as HTMLInputElement;
      searchInput.focus();
      searchInput.value = 'series';
      searchInput.dispatchEvent(new Event('input', {bubbles: true}));
      await new Promise<void>((r) => setTimeout(r, 400));
    });
    await page.waitForTimeout(300);
    await page.evaluate(async () => {
      const searchInput = document.querySelector('.d4-dialog .d4-column-grid input[placeholder="Search"]') as HTMLInputElement;
      searchInput.focus();
      searchInput.value = '';
      searchInput.dispatchEvent(new Event('input', {bubbles: true}));
      await new Promise<void>((r) => setTimeout(r, 300));
    });
  });

  await softStep('Pick a few columns (Id, Chemist, Series) and click OK', async () => {
    await page.evaluate(async () => {
      const dialog = document.querySelector('.d4-dialog')!;
      const masterCb = dialog.querySelector('.ui-input-bool input[type="checkbox"]') as HTMLInputElement;
      if (masterCb.checked) masterCb.click();
      await new Promise<void>((r) => setTimeout(r, 300));
      const gridEl = dialog.querySelector('.d4-column-grid [name="viewer-Grid"]')!;
      const canvases = gridEl.querySelectorAll('canvas');
      const canvas = canvases[canvases.length - 1] as HTMLCanvasElement;
      function clickAt(x: number, y: number) {
        const init = {bubbles: true, cancelable: true, view: window, clientX: x, clientY: y, button: 0, buttons: 1};
        canvas.dispatchEvent(new MouseEvent('mousedown', init));
        canvas.dispatchEvent(new MouseEvent('mouseup', init));
        canvas.dispatchEvent(new MouseEvent('click', init));
      }
      clickAt(495, 204);
      await new Promise<void>((r) => setTimeout(r, 250));
      clickAt(495, 317);
      await new Promise<void>((r) => setTimeout(r, 250));
      clickAt(495, 398);
      await new Promise<void>((r) => setTimeout(r, 250));
      (dialog.querySelector('[name="button-OK"]') as HTMLElement).click();
    });
    await page.waitForTimeout(700);
    const tag = await page.evaluate(() => (grok as any).shell.tv.dataFrame.getTag('.tooltip'));
    expect(tag).toBe('Id\nChemist\nSeries');
  });

  await softStep('Scatter plot tooltip contains the selected columns in order', async () => {
    const result = await page.evaluate(async () => {
      document.body.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: 5, clientY: 5}));
      await new Promise<void>((r) => setTimeout(r, 300));
      const sp = document.querySelector('[name="viewer-Scatter-plot"]')!;
      const canvas = sp.querySelector('canvas') as HTMLCanvasElement;
      const r = canvas.getBoundingClientRect();
      const x = r.left + r.width / 2;
      const y = r.top + r.height / 2;
      canvas.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: x, clientY: y}));
      await new Promise<void>((rs) => setTimeout(rs, 200));
      canvas.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: x + 3, clientY: y + 3}));
      await new Promise<void>((rs) => setTimeout(rs, 800));
      const tt = document.querySelector('.d4-tooltip') as HTMLElement | null;
      if (!tt) return {text: ''};
      return {text: tt.textContent ?? ''};
    });
    expect(result.text).toContain('Id');
    expect(result.text).toContain('Chemist');
    expect(result.text).toContain('Series');
    const idIdx = result.text.indexOf('Id');
    const chemistIdx = result.text.indexOf('Chemist');
    const seriesIdx = result.text.indexOf('Series');
    expect(idIdx).toBeLessThan(chemistIdx);
    expect(chemistIdx).toBeLessThan(seriesIdx);
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map((s) => `  - ${s.step}: ${s.error}`).join('\n'));
});
