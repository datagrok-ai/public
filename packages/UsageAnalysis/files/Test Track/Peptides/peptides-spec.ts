import { test, expect, Page } from '@playwright/test';

const BASE_URL = 'https://public.datagrok.ai';

test.describe('Peptides / Peptides Panel', () => {
  let page: Page;

  test.beforeAll(async ({ browser }) => {
    page = await browser.newPage();
    await page.goto(BASE_URL, { waitUntil: 'networkidle', timeout: 30000 });
    await page.waitForTimeout(3000);
  });

  test.afterAll(async () => {
    await page.close();
  });

  test('Step 1: Open peptides.csv', async () => {
    await page.evaluate(async () => {
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:DemoFiles/bio/peptides.csv');
      grok.shell.addTableView(df);
      await new Promise(r => setTimeout(r, 3000));
    });
    await expect(page.locator('text=Rows: 647')).toBeVisible({ timeout: 10000 });
  });

  test('Step 2: Click AlignedSequence column title', async () => {
    await page.evaluate(async () => {
      grok.shell.windows.showContextPanel = true;
      grok.shell.windows.showHelp = false;
      await new Promise(r => setTimeout(r, 500));
      const col = grok.shell.tv.dataFrame.columns.byName('AlignedSequence');
      grok.shell.o = col;
      await new Promise(r => setTimeout(r, 2000));
    });
    await expect(page.locator('text=AlignedSequence').first()).toBeVisible({ timeout: 5000 });
  });

  test('Step 3: Expand Peptides panel', async () => {
    await page.evaluate(async () => {
      const panes = Array.from(document.querySelectorAll('.d4-accordion-pane-header, .grok-accordion-pane-header'));
      const peptidesPane = panes.find(p => p.textContent?.trim() === 'Peptides');
      if (peptidesPane) peptidesPane.click();
      await new Promise(r => setTimeout(r, 3000));
    });
    const content = await page.evaluate(() => {
      const panes = Array.from(document.querySelectorAll('.d4-accordion-pane-header, .grok-accordion-pane-header'));
      const peptidesPane = panes.find(p => p.textContent?.trim() === 'Peptides');
      return peptidesPane?.nextElementSibling?.textContent?.slice(0, 100) || '';
    });
    expect(content).toContain('Activity');
    expect(content).toContain('Scaling');
  });

  test('Step 4: Change Scaling parameter', async () => {
    const newVal = await page.evaluate(async () => {
      const panes = Array.from(document.querySelectorAll('.d4-accordion-pane-header, .grok-accordion-pane-header'));
      const peptidesPane = panes.find(p => p.textContent?.trim() === 'Peptides');
      const content = peptidesPane?.nextElementSibling;
      const scalingSelect = Array.from(content?.querySelectorAll('select') || [])
        .find(s => Array.from(s.options).some(o => o.value === 'lg'));
      if (scalingSelect) {
        const nativeSetter = Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value').set;
        nativeSetter.call(scalingSelect, 'lg');
        scalingSelect.dispatchEvent(new Event('change', { bubbles: true }));
        await new Promise(r => setTimeout(r, 2000));
        return scalingSelect.value;
      }
      return 'not found';
    });
    expect(newVal).toBe('lg');
  });

  test('Steps 5-6: Click amino acid in weblogo; rows selected', async () => {
    const result = await page.evaluate(async () => {
      const panes = Array.from(document.querySelectorAll('.d4-accordion-pane-header, .grok-accordion-pane-header'));
      const peptidesPane = panes.find(p => p.textContent?.trim() === 'Peptides');
      const content = peptidesPane?.nextElementSibling;
      const canvas = content?.querySelectorAll('canvas')[0];
      if (!canvas) return { error: 'no canvas' };

      const r = canvas.getBoundingClientRect();
      const clickX = Math.min(r.left + 100, 1440);
      const clickY = r.top + r.height / 2;

      ['mousedown', 'mouseup', 'click'].forEach(type => {
        canvas.dispatchEvent(new MouseEvent(type, { bubbles: true, clientX: clickX, clientY: clickY, cancelable: true }));
      });
      await new Promise(r => setTimeout(r, 1500));

      return { selectedRows: grok.shell.tv.dataFrame.selection.trueCount };
    });
    expect(result.selectedRows).toBeGreaterThan(0);
    await expect(page.locator('text=selected rows')).toBeVisible({ timeout: 5000 });
  });
});
