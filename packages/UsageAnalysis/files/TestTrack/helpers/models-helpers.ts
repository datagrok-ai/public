import {Page, expect} from '@playwright/test';

export async function setPredict(page: Page, columnName: string): Promise<void> {
  const current = await page.evaluate(() => (window as any).grok.shell.v.root
    .querySelector('[name="input-host-Predict"] .d4-column-selector-column')?.textContent?.trim());
  if (current === columnName) return;
  await page.evaluate(() => {
    const root = (window as any).grok.shell.v.root;
    const sel = root.querySelector('[name="input-host-Predict"] .d4-column-selector') as HTMLElement;
    sel.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, cancelable: true, view: window, button: 0}));
  });
  await page.waitForFunction(() => !!document.querySelector('.d4-column-selector-backdrop'),
    null, {timeout: 5_000});
  await page.evaluate(() => (document.querySelector('.d4-column-selector-backdrop') as HTMLElement).focus());
  await page.keyboard.type(columnName);
  await page.waitForTimeout(300);
  await page.keyboard.press('Enter').catch(() => {});
  await page.waitForTimeout(500);
  await expect.poll(async () => await page.evaluate(() =>
    (window as any).grok.shell.v.root
      .querySelector('[name="input-host-Predict"] .d4-column-selector-column')?.textContent?.trim()),
    {timeout: 10_000}).toBe(columnName);
}

async function openFeaturesPicker(page: Page): Promise<void> {
  await page.evaluate(() => {
    const editor = (window as any).grok.shell.v.root
      .querySelector('[name="div-Features"]') as HTMLElement;
    const r = editor.getBoundingClientRect();
    const ev = (t: string) => new MouseEvent(t, {
      bubbles: true, cancelable: true, view: window, button: 0,
      clientX: r.left + 10, clientY: r.top + 5,
    });
    editor.dispatchEvent(ev('mousedown'));
    editor.dispatchEvent(ev('mouseup'));
    editor.dispatchEvent(ev('click'));
  });
  await page.locator('[name="dialog-Select-columns..."]').waitFor({timeout: 10_000});
}

async function readCheckedCount(page: Page): Promise<number> {
  return await page.locator('[name="dialog-Select-columns..."]').evaluate((d) => {
    const lbl = Array.from(d.querySelectorAll('label'))
      .find((l) => /\d+ checked/.test(l.textContent || ''));
    const m = lbl?.textContent?.match(/^(\d+)/);
    return m ? parseInt(m[1], 10) : 0;
  });
}

export async function selectFeaturesByName(page: Page, columnNames: string[]): Promise<void> {
  await openFeaturesPicker(page);
  const dlg = page.locator('[name="dialog-Select-columns..."]');
  await dlg.locator('[name="label-None"]').click();
  await page.waitForTimeout(300);

  const search = dlg.locator('input[placeholder="Search"]');
  const overlay = dlg.locator('[name="viewer-Grid"] [name="overlay"]');

  for (const name of columnNames) {
    const before = await readCheckedCount(page);
    await search.fill('');
    await page.waitForTimeout(200);
    await search.fill(name);
    await page.waitForTimeout(500);
    const box = await overlay.boundingBox();
    if (!box) throw new Error('column-picker overlay canvas not measurable');
    await page.mouse.click(box.x + box.width - 40, box.y + 36);
    await page.waitForTimeout(300);
    const after = await readCheckedCount(page);
    if (after !== before + 1) {

      await search.fill('');
      await dlg.locator('[name="label-All"]').click();
      await page.waitForTimeout(300);
      break;
    }
  }
  await search.fill('');
  await page.waitForTimeout(200);
  await dlg.locator('[name="button-OK"]').click();
  await page.waitForFunction(() => !document.querySelector('[name="dialog-Select-columns..."]'),
    null, {timeout: 10_000});
}

export async function uncheckFeatureRowByPosition(
  page: Page, rowIndex: number, _totalRows?: number): Promise<void> {
  await openFeaturesPicker(page);
  await page.waitForTimeout(1_200);
  const dlg = page.locator('[name="dialog-Select-columns..."]');
  const target = await page.evaluate((idx) => {
    const overlay = document.querySelector(
      '[name="dialog-Select-columns..."] canvas[name="overlay"]') as HTMLCanvasElement;
    const r = overlay.getBoundingClientRect();
    return {x: r.left + r.width - 40, y: r.top + 36 + 28 * idx};
  }, rowIndex);
  await page.mouse.click(target.x, target.y);
  await page.waitForTimeout(400);
  await dlg.locator('[name="button-OK"]').click();
  await page.waitForFunction(() => !document.querySelector('[name="dialog-Select-columns..."]'),
    null, {timeout: 10_000});
}
