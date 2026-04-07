import { test, expect, Page } from '@playwright/test';

const BASE_URL = 'https://public.datagrok.ai';

test.describe('Peptides / SAR Analysis', () => {
  let page: Page;

  test.beforeAll(async ({ browser }) => {
    page = await browser.newPage();
    await page.goto(BASE_URL, { waitUntil: 'networkidle', timeout: 30000 });
    await page.waitForTimeout(3000);
  });

  test.afterAll(async () => {
    await page.close();
  });

  test('Steps 1-3: Load dataset, select column, expand Peptides panel', async () => {
    await page.evaluate(async () => {
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:DemoFiles/bio/peptides.csv');
      grok.shell.addTableView(df);
      await new Promise(r => setTimeout(r, 3000));
      grok.shell.windows.showContextPanel = true;
      grok.shell.windows.showHelp = false;
      const col = grok.shell.tv.dataFrame.columns.byName('AlignedSequence');
      grok.shell.o = col;
      await new Promise(r => setTimeout(r, 2000));
      const panes = Array.from(document.querySelectorAll('.d4-accordion-pane-header, .grok-accordion-pane-header'));
      const peptidesPane = panes.find(p => p.textContent?.trim() === 'Peptides');
      if (peptidesPane) peptidesPane.click();
      await new Promise(r => setTimeout(r, 2000));
    });
    await expect(page.locator('text=Launch SAR').first()).toBeVisible({ timeout: 10000 });
  });

  test('Step 4: Launch SAR', async () => {
    await page.evaluate(async () => {
      const btn = Array.from(document.querySelectorAll('button, .ui-btn'))
        .find(b => b.textContent?.trim() === 'Launch SAR');
      if (btn) btn.click();
      await new Promise(r => setTimeout(r, 12000));
    });
    await expect(page.locator('text=Mutation Cliffs').first()).toBeVisible({ timeout: 15000 });
  });

  test('Step 5: Four viewers visible', async () => {
    // MCL viewer
    await expect(page.locator('text=MCL')).toBeVisible({ timeout: 10000 });
    // Logo Summary Table
    await expect(page.locator('text=Logo Summary Table')).toBeVisible();
    // Most Potent Residues
    await expect(page.locator('text=Mos').first()).toBeVisible();
    // Sequence Variability Map
    await expect(page.locator('text=Sequence Variabi').first()).toBeVisible();
  });

  test('Steps 6-8: Open settings, change parameters, click OK', async () => {
    // Click wrench icon
    await page.evaluate(() => {
      const wrench = document.querySelector('.d4-ribbon-item .fa-wrench');
      if (wrench) wrench.click();
    });
    await expect(page.locator('text=Peptides settings')).toBeVisible({ timeout: 5000 });

    // Change scaling and similarity threshold
    await page.evaluate(() => {
      const dialog = document.querySelector('.d4-dialog');
      const scalingSelect = Array.from(dialog.querySelectorAll('select'))
        .find(s => Array.from(s.options).some(o => o.value === 'lg'));
      if (scalingSelect) {
        const ns = Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value').set;
        ns.call(scalingSelect, 'lg');
        scalingSelect.dispatchEvent(new Event('change', { bubbles: true }));
      }
      const simInput = Array.from(dialog.querySelectorAll('input')).find(i => i.value === '70');
      if (simInput) {
        const ns = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value').set;
        ns.call(simInput, '60');
        simInput.dispatchEvent(new Event('change', { bubbles: true }));
      }
    });

    // Click OK
    await page.evaluate(() => {
      const btns = Array.from(document.querySelectorAll('button, .ui-btn'));
      const ok = btns.find(b => b.textContent?.trim() === 'OK');
      if (ok) ok.click();
    });
    await page.waitForTimeout(12000);

    // Dialog should be gone, viewers reloaded
    const dialogGone = await page.evaluate(() => !document.querySelector('.d4-dialog'));
    expect(dialogGone).toBe(true);
  });

  test('Step 9: Viewers reload with new parameters', async () => {
    const updating = await page.evaluate(() =>
      Array.from(document.querySelectorAll('*'))
        .some(el => el.textContent?.trim() === 'Updating...' && el.children.length === 0)
    );
    expect(updating).toBe(false);
    await expect(page.locator('text=Columns: 26')).toBeVisible({ timeout: 5000 });
  });

  test('Step 10: Switch between Mutation Cliffs and Invariant Map', async () => {
    await page.evaluate(async () => {
      const label = Array.from(document.querySelectorAll('label, span'))
        .find(el => el.textContent?.trim() === 'Invariant Map');
      if (label) label.click();
      await new Promise(r => setTimeout(r, 1000));
    });
    await expect(page.locator('text=Invariant Map').first()).toBeVisible();

    // Switch back
    await page.evaluate(() => {
      const label = Array.from(document.querySelectorAll('label, span'))
        .find(el => el.textContent?.trim() === 'Mutation Cliffs');
      if (label) label.click();
    });
    await page.waitForTimeout(500);
  });

  test('Steps 11-12: Click cell; Mutation Cliffs pairs and Distribution panels appear', async () => {
    const panels = await page.evaluate(async () => {
      const svmCanvas = Array.from(document.querySelectorAll('canvas'))[4];
      const r = svmCanvas.getBoundingClientRect();
      ['mousedown', 'mouseup', 'click'].forEach(type => {
        svmCanvas.dispatchEvent(new MouseEvent(type, {
          bubbles: true, clientX: r.left + 95, clientY: r.top + 33, cancelable: true
        }));
      });
      await new Promise(r => setTimeout(r, 2000));
      return Array.from(document.querySelectorAll('.d4-accordion-pane-header, .grok-accordion-pane-header'))
        .map(el => el.textContent?.trim()).filter(t => t);
    });
    expect(panels).toContain('Mutation Cliffs pairs');
    expect(panels).toContain('Distribution');
  });
});
