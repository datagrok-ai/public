/* ---
realizes: [formsviewer.cp.fields-lifecycle-and-number-format, formsviewer.int.number-format-vs-grid, formsviewer.edge.empty-fields-column-names, formsviewer.edge.tilde-columns-excluded, formsviewer.edge.over-20-columns-capped-silently, formsviewer.edge.renamed-column-follows-rename, formsviewer.edge.number-format-float-only]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import {addViewerByIcon, cleanupShell, finishSpec, openTable, openViewerProperties, pollValue} from '../../helpers/viewers';
import {HOST, CURRENT, drawnLabelNames, balloonCount, withConsoleErrorCount} from '../../helpers/forms';

declare const grok: any;
declare const DG: any;

const datasetPath = 'System:DemoFiles/demog.csv';

async function cardFieldColumns(page: Page): Promise<string[]> {
  return page.evaluate((sel) => {
    const card = document.querySelector(sel) as HTMLElement | null;
    return card ? Array.from(card.querySelectorAll('[column]')).map((e) => e.getAttribute('column') ?? '') : [];
  }, CURRENT);
}

async function cardFieldText(page: Page, column: string): Promise<string | null> {
  return page.evaluate(({sel, col}) => {
    const card = document.querySelector(sel) as HTMLElement | null;
    const el = card?.querySelector(`[column="${col}"]`) as HTMLInputElement | null;
    return el ? (el.value ?? el.textContent ?? '') : null;
  }, {sel: CURRENT, col: column});
}

async function gridCellText(page: Page, column: string, row: number): Promise<string> {
  return page.evaluate(({col, r}) => {
    const c = grok.shell.t.col(col);
    if (c.isNone(r))
      return '';
    return DG.format(c.get(r), grok.shell.tv.grid.col(col).format);
  }, {col: column, r: row});
}

async function gridHeaderPoint(page: Page, column: string): Promise<{x: number; y: number}> {
  return page.evaluate((col) => {
    const grid = grok.shell.tv.grid;
    const gc = grid.columns.byName(col);
    const canvas = Array.from(document.querySelectorAll('[name="viewer-Grid"] canvas'))
      .find((c: any) => c.getBoundingClientRect().width > 100) as HTMLElement;
    const r = canvas.getBoundingClientRect();
    return {x: Math.round(r.x + gc.left + gc.width / 2), y: Math.round(r.y + grid.colHeaderHeight / 2)};
  }, column);
}

async function liveMenuItemCentre(page: Page, sel: string, ms: number): Promise<{x: number; y: number} | null> {
  const deadline = Date.now() + ms;
  for (;;) {
    const c = await page.evaluate((s) => {
      const e = document.querySelector(s);
      if (!e) return null;
      const r = e.getBoundingClientRect();
      if (r.width === 0 || r.height === 0) return null;
      return {x: r.x + r.width / 2, y: r.y + r.height / 2};
    }, sel);
    if (c || Date.now() > deadline) return c;
    await page.waitForTimeout(150);
  }
}

async function renameColumnViaGridUI(page: Page, oldName: string, newName: string): Promise<void> {
  await expect.poll(async () => {
    const names = await page.evaluate(() => grok.shell.t.columns.names() as string[]);
    if (names.includes(newName) && !names.includes(oldName))
      return true;
    if (!names.includes(oldName))
      return false; 
    await page.keyboard.press('Escape');
    const pt = await gridHeaderPoint(page, oldName);
    await page.mouse.click(pt.x, pt.y, {button: 'right'});
    const item = page.locator('[name="div-Column-Properties..."]').first();
    if (!await item.isVisible().catch(() => false))
      return false;
    await item.click();
    const nameInput = page.locator(`[name="dialog-${oldName}"] [name="input-New-name--"]`);
    if (!await nameInput.isVisible({timeout: 3000}).catch(() => false))
      return false;
    await nameInput.fill(newName);
    await page.locator(`[name="dialog-${oldName}"] [name="button-OK"]`).click();
    return false; 
  }, {timeout: 30_000, intervals: [300, 400, 600, 800, 1000]}).toBe(true);
}

async function removeColumnViaGridUI(page: Page, column: string): Promise<void> {
  await expect.poll(async () => {
    if (await page.evaluate((col) => !(grok.shell.t.columns.names() as string[]).includes(col), column))
      return true;
    await page.keyboard.press('Escape');
    const pt = await gridHeaderPoint(page, column);
    await page.mouse.click(pt.x, pt.y, {button: 'right'});
    const centre = await liveMenuItemCentre(page, '[name="div-Remove"]', 3000);
    if (centre)
      await page.mouse.click(centre.x, centre.y);
    return false; 
  }, {timeout: 20_000, intervals: [200, 300, 400, 600, 800]}).toBe(true);
}

async function setNumberFormatViaPanel(page: Page, value: string): Promise<void> {

  await page.evaluate(() => { grok.shell.windows.showContextPanel = false; });

  await page.waitForTimeout(300);
  await page.evaluate(() => {
    const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
    (vw?.root?.closest('.panel-base')
      ?.querySelector('.panel-titlebar [name="icon-font-icon-settings"]') as HTMLElement | null)?.click();
  });
  await page.locator('[name="prop-number-format"]').first().waitFor({state: 'visible', timeout: 8000}).catch(() => {});

  await expect.poll(async () => page.evaluate((v) => {
    const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
    if (vw.props.numberFormat === v) return true;
    const row = document.querySelector('[name="prop-number-format"]') as HTMLElement | null;
    if (!row) return false;

    if (row.offsetParent === null) {
      (document.querySelector('[name="prop-category-misc"]') as HTMLElement | null)?.click();
      return false;
    }
    (row.querySelector('[name="prop-view-number-format"]') as HTMLElement | null)?.click();
    const sel = row.querySelector('select') as HTMLSelectElement | null;
    if (!sel) return false;
    const opt = Array.from(sel.options).find((o) => o.value === v || (o.textContent ?? '').trim() === v);
    if (!opt) return false;
    sel.value = opt.value;
    sel.dispatchEvent(new Event('input', {bubbles: true}));
    sel.dispatchEvent(new Event('change', {bubbles: true}));
    return false;
  }, value), {timeout: 12_000, intervals: [250, 300, 400, 600, 800]}).toBe(true);

  await page.waitForTimeout(300);
}

async function removeFieldViaHeaderCloseIcon(page: Page, column: string): Promise<void> {
  await expect.poll(async () => {
    if (!await page.evaluate((col) => grok.shell.tv.viewers
      .find((x: any) => x.type === 'FormsViewer').fieldsColumnNames.includes(col), column))
      return true;

    const pts = await page.evaluate(({host, col}) => {
      const container = Array.from(document.querySelectorAll(
        `${host} .d4-multi-form-header .d4-multi-form-column-name`))
        .find((c) => c.querySelector(`div[name="div-${col}"]`)) as HTMLElement | undefined;
      if (!container)
        return null;
      const icon = container.querySelector('i.grok-icon.fal.fa-times') as HTMLElement | null;
      if (!icon)
        return null;
      const cr = container.getBoundingClientRect();
      const ir = icon.getBoundingClientRect();
      if (ir.width === 0 || ir.height === 0)
        return {zero: true as const};
      return {cx: cr.x + cr.width / 2, cy: cr.y + cr.height / 2, ix: ir.x + ir.width / 2, iy: ir.y + ir.height / 2};
    }, {host: HOST, col: column});
    if (pts === null)
      return false;

    if ('zero' in pts)
      throw new Error(`[removeFieldViaHeaderCloseIcon] '${column}' close icon has a zero-size box — a real mouse cannot target it; investigate before any synthetic fallback.`);

    await page.mouse.move(pts.cx, pts.cy);
    await page.mouse.click(pts.ix, pts.iy);
    return false; 
  }, {timeout: 20_000, intervals: [150, 150, 250, 400, 600]}).toBe(true);
}

test.use(specTestOptions);

test('Forms viewer — field lifecycle and number format', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  await softStep('Setup — open demog, add COMPUTED_H (no-format FLOAT), attach the Forms viewer', async () => {
    await openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

    const entry = await page.evaluate(() => {
      const df = grok.shell.t;
      const h = df.col('HEIGHT');
      const comp = df.columns.addNewFloat('COMPUTED_H');
      for (let i = 0; i < df.rowCount; i++) comp.set(i, h.get(i) * 1.0);
      df.currentRowIdx = 0;
      return {type: comp.type, format: comp.meta.format ?? null};
    });
    expect(entry.type).toBe('double');
    expect(entry.format).toBeNull();

    await addViewerByIcon(page, 'Forms', 'Forms', 30_000, 'FormsViewer');
    await page.locator(HOST).first().waitFor({timeout: 30_000});
    await page.locator('.d4-multi-form').first().waitFor({timeout: 30_000});
    await page.evaluate(() => { grok.shell.t.currentRowIdx = 0; });
  });

  await softStep('Scenario 1 Step 3 — Fields render in the picked order RACE, AGE, SEX (not table order)', async () => {

    await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      vw.setOptions({fieldsColumnNames: ['RACE', 'AGE', 'SEX']});
    });
    await expect.poll(() => drawnLabelNames(page), {timeout: 20_000}).toEqual(['RACE', 'AGE', 'SEX']);
    await expect.poll(() => cardFieldColumns(page), {timeout: 20_000}).toEqual(['RACE', 'AGE', 'SEX']);
  });

  await softStep('Scenario 1 Step 4 — Removing the AGE field drops its label and elements; order kept', async () => {
    await removeFieldViaHeaderCloseIcon(page, 'AGE');

    await expect.poll(() => drawnLabelNames(page), {timeout: 20_000}).toEqual(['RACE', 'SEX']);
    expect(await page.locator(`${HOST} [column="AGE"]`).count()).toBe(0);
    await expect.poll(() => cardFieldColumns(page), {timeout: 20_000}).toEqual(['RACE', 'SEX']);
  });

  await softStep('Scenario 1 Step 5 — Dropping RACE from the dataframe prunes its field with no error', async () => {
    const errCount = await withConsoleErrorCount(page, async () => {
      await removeColumnViaGridUI(page, 'RACE');
    });
    expect(await page.evaluate(() => grok.shell.t.columns.names().includes('RACE'))).toBe(false);
    await expect.poll(() => drawnLabelNames(page), {timeout: 20_000}).toEqual(['SEX']);
    expect(await page.locator(`${HOST} [column="RACE"]`).count()).toBe(0);
    expect(errCount).toBe(0);
    expect(await page.locator('.d4-balloon.error').count()).toBe(0);
  });

  await softStep('Scenario 1 Step 6 — Renaming SEX→GENDER via the grid UI: the field follows the new name', async () => {
    const errCount = await withConsoleErrorCount(page, async () => {
      await renameColumnViaGridUI(page, 'SEX', 'GENDER');
    });
    await expect.poll(() => drawnLabelNames(page), {timeout: 20_000}).toEqual(['GENDER']);
    await expect.poll(() => cardFieldColumns(page), {timeout: 20_000}).toEqual(['GENDER']);
    expect(await page.evaluate(() => grok.shell.t.columns.names().includes('GENDER'))).toBe(true);
    expect(errCount).toBe(0);
  });

  await softStep('Scenario 1 Step 7 — A >20-column table yields exactly 20 fields with no message', async () => {

    await cleanupShell(page);
    await page.evaluate(() => {
      const N = 30;
      const cols: any[] = [];
      for (let c = 0; c < 25; c++) {
        const name = 'C' + String(c).padStart(2, '0');
        cols.push(DG.Column.fromInt32Array(name, Int32Array.from(Array.from({length: N}, (_, i) => i + c * 100))));
      }
      cols.push(DG.Column.fromStrings('LABEL', Array.from({length: N}, (_, i) => 'r' + i)));
      cols.push(DG.Column.fromStrings('~hidden', Array.from({length: N}, () => 'h')));
      const df = DG.DataFrame.fromColumns(cols);
      df.name = 'wide26';
      grok.shell.addTableView(df);
      df.currentRowIdx = 0;
    });

    await pollValue(() => page.evaluate(() => grok.shell.t?.name ?? ''),
      (name) => name === 'wide26', 600, 100);
    await addViewerByIcon(page, 'Forms', 'Forms', 30_000, 'FormsViewer');
    await page.locator(HOST).first().waitFor({timeout: 30_000});

    await expect.poll(() => drawnLabelNames(page).then((l) => l.length), {timeout: 20_000}).toBe(20);
    expect(await cardFieldColumns(page)).toHaveLength(20);
    expect((await drawnLabelNames(page)).some((n) => n.startsWith('~'))).toBe(false);
    expect(await balloonCount(page)).toBe(0);
    expect(await page.evaluate((host) => (document.querySelector(host) as HTMLElement)?.textContent?.includes('Number of columns') ?? false, HOST)).toBe(false);
  });

  await softStep('Scenario 1 Step 8 — Renaming C05→~SERVICE drops it from the field set (no ~ label)', async () => {

    const beforeCount = (await drawnLabelNames(page)).length;
    expect(beforeCount).toBeGreaterThan(0);

    const errCount = await withConsoleErrorCount(page, async () => {
      await renameColumnViaGridUI(page, 'C05', '~SERVICE');
    });
    await expect.poll(() => drawnLabelNames(page).then((l) => l.length), {timeout: 20_000}).toBe(beforeCount - 1);
    const labels = await drawnLabelNames(page);
    expect(labels.some((n) => n.startsWith('~'))).toBe(false);
    expect(labels).not.toContain('C05');
    expect(await page.locator(`${HOST} [column="~SERVICE"]`).count()).toBe(0);
    expect(errCount).toBe(0);
    expect(await page.locator('.d4-balloon.error').count()).toBe(0);
  });

  await softStep('Scenario 1 Step 8b — Clearing every field draws nothing and raises no error', async () => {

    await openViewerProperties(page, 'Forms', '[name="prop-view-fields"]');
    const dlg = '[name="dialog-Select-columns..."]';
    const errCount = await withConsoleErrorCount(page, async () => {

      const fieldsButton = page.locator('[name="prop-view-fields"] button');
      await fieldsButton.waitFor({state: 'visible', timeout: 10_000});
      await fieldsButton.click();
      await page.locator(dlg).waitFor({state: 'visible', timeout: 10_000});
      await page.locator(`${dlg} [name="label-None"]`).click();
      await page.locator(`${dlg} [name="button-OK"]`).click();
      await expect.poll(() => drawnLabelNames(page).then((l) => l.length), {timeout: 20_000}).toBe(0);
    });
    expect(await drawnLabelNames(page)).toEqual([]);
    expect(await page.locator(`${HOST} [column]`).count()).toBe(0);
    expect(errCount).toBe(0);
    expect(await page.locator('.d4-balloon.error').count()).toBe(0);
  });

  await softStep('Scenario 2 Step 4 — COMPUTED_H field text equals the grid column format applied to the raw value (Same as grid)', async () => {
    await openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
    await page.evaluate(() => {
      const df = grok.shell.t;
      const h = df.col('HEIGHT');
      const comp = df.columns.addNewFloat('COMPUTED_H');
      for (let i = 0; i < df.rowCount; i++) comp.set(i, h.get(i) * 1.0);
      df.currentRowIdx = 0;
    });
    await addViewerByIcon(page, 'Forms', 'Forms', 30_000, 'FormsViewer');
    await page.locator(HOST).first().waitFor({timeout: 30_000});

    const defaultFmt = await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      grok.shell.t.currentRowIdx = 0;
      return vw.props.numberFormat;
    });
    expect(defaultFmt).toBe('Same as grid');

    await page.waitForTimeout(800);

    const gridText = await gridCellText(page, 'COMPUTED_H', 0);
    expect(gridText.length).toBeGreaterThan(0);
    await expect.poll(() => cardFieldText(page, 'COMPUTED_H'), {timeout: 15_000}).toBe(gridText);
  });

  await softStep('Scenario 2 Step 6 — An explicit float mask changes COMPUTED_H and differs from the grid', async () => {

    const before = await page.evaluate(() => grok.shell.t.col('COMPUTED_H').get(0));
    await setNumberFormatViaPanel(page, '3 digits after comma');
    await expect.poll(() => cardFieldText(page, 'COMPUTED_H'), {timeout: 15_000})
      .toMatch(/^-?\d+\.\d{3}$/);
    const gridText = await gridCellText(page, 'COMPUTED_H', 0);
    expect(await cardFieldText(page, 'COMPUTED_H')).not.toBe(gridText);
    expect(await page.evaluate(() => grok.shell.t.col('COMPUTED_H').get(0))).toBe(before);
  });

  await softStep('Scenario 2 Step 7 — INT and string fields render identically before/after the mask change', async () => {

    await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      vw.setOptions({fieldsColumnNames: ['COMPUTED_H', 'AGE', 'SEX']});
      grok.shell.t.currentRowIdx = 0;
    });
    await setNumberFormatViaPanel(page, 'Same as grid');

    await page.waitForTimeout(400);
    const ageBefore = await cardFieldText(page, 'AGE');
    const sexBefore = await cardFieldText(page, 'SEX');

    expect((ageBefore ?? '').length).toBeGreaterThan(0);
    expect((sexBefore ?? '').length).toBeGreaterThan(0);

    await setNumberFormatViaPanel(page, '3 digits after comma');

    await expect.poll(() => cardFieldText(page, 'COMPUTED_H'), {timeout: 15_000}).toMatch(/^-?\d+\.\d{3}$/);
    expect(await cardFieldText(page, 'AGE')).toBe(ageBefore);
    expect(await cardFieldText(page, 'SEX')).toBe(sexBefore);
  });

  await cleanupShell(page);
  finishSpec();
});
