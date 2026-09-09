import { test, expect } from '@playwright/test';
import * as H from './helpers';

test.describe.configure({ mode: 'serial' });

test('Sticky Meta: entity type & schema lifecycle (create, edit, delete)', async ({ page }) => {
  test.setTimeout(180_000);

  const suffix = H.uniqueSuffix();
  const typeName = `PW_SM_Type_${suffix}`;
  const schemaName = `PW_SM_Schema_${suffix}`;
  const typeCheckbox = `[name="prop-view-${typeName.toLowerCase()}"]`;

  try {
    await H.gotoHome(page);
    await H.setupEnv(page);
    await H.apiDeleteAllTestSchemas(page); 

    await H.gotoStickyMeta(page, 'Types');
    await page.locator('[name="button-New-Entity-Type..."]').click();
    const typeDialog = page.locator('.d4-dialog[name="dialog-Create-a-new-entity-type"]');
    await typeDialog.locator('[name="input-Name"]').waitFor({ timeout: 10_000 });

    expect(await typeDialog.locator('[name="button-OK"]').evaluate((b) => b.classList.contains('disabled'))).toBe(true);

    await typeDialog.locator('[name="input-Name"]').click();
    await page.keyboard.press('Control+A');
    await page.keyboard.type(typeName, { delay: 10 });

    expect(await typeDialog.locator('[name="button-OK"]').evaluate((b) => b.classList.contains('disabled'))).toBe(true);

    await typeDialog.locator('[name="input-Matching-expression"]').click();
    await page.keyboard.press('Control+A');
    await page.keyboard.type('semtype=molecule', { delay: 10 });
    await expect(typeDialog.locator('[name="button-OK"]')).not.toHaveClass(/disabled/);

    await H.clickDialogButton(page, 'dialog-Create-a-new-entity-type', 'button-OK');
    await page.waitForTimeout(2500);

    await H.searchList(page, typeName);
    expect(await H.listHasCard(page, typeName)).toBe(true);

    await H.gotoStickyMeta(page, 'Schemas');
    await page.locator('[name="button-New-Schema..."]').click();
    const schemaDialog = page.locator('.d4-dialog[name="dialog-Create-a-new-schema"]');
    await schemaDialog.waitFor({ timeout: 10_000 });

    await schemaDialog.locator('.ui-form > [name="input-host-Name"] > [name="input-Name"]').click();
    await page.keyboard.press('Control+A');
    await page.keyboard.type(schemaName, { delay: 10 });

    await schemaDialog.locator('[name="div-Associated-with-"] .d4-link-action').click();
    await page.waitForTimeout(800);
    await page.evaluate((cbSel) => {
      const cb = document.querySelector(cbSel) as HTMLInputElement | null;
      if (cb && !cb.checked) cb.click();
    }, typeCheckbox);
    await page.evaluate(() => {
      const picker = Array.from(document.querySelectorAll('.d4-dialog'))
        .find((d) => /Select types/i.test(d.querySelector('.d4-dialog-title')?.textContent || ''));
      (picker?.querySelector('[name="button-OK"]') as HTMLElement | null)?.click();
    });
    await page.waitForTimeout(600);
    expect(await schemaDialog.locator('[name="div-Associated-with-"]').textContent()).toContain(typeName);

    const props: [string, string][] = [
      ['rating', 'int'], ['notes', 'string'], ['verified', 'bool'], ['review_date', 'datetime'],
    ];

    const typeOptions = await schemaDialog.locator('select[name="input-Property-Type"]').first()
      .evaluate((s) => Array.from((s as HTMLSelectElement).options).map((o) => o.value));
    for (const t of ['string', 'int', 'bool', 'double', 'datetime', 'string_list'])
      expect(typeOptions).toContain(t);

    for (let i = 0; i < props.length; i++) {
      const [name, type] = props[i];
      if (i > 0) {
        await schemaDialog.locator('[name="button-Add-new-property-to-schema"]').click();
        await page.waitForTimeout(250);
      }
      const row = schemaDialog.locator('table.d4-item-table input[name="input-Name"]').nth(i);
      await row.click();
      await page.keyboard.press('Control+A');
      await page.keyboard.type(name, { delay: 10 });
      await page.evaluate(({ i, type }) => {
        const sel = document.querySelectorAll('.d4-dialog table.d4-item-table select[name="input-Property-Type"]')[i] as HTMLSelectElement;
        sel.value = type;
        sel.dispatchEvent(new Event('input', { bubbles: true }));
        sel.dispatchEvent(new Event('change', { bubbles: true }));
      }, { i, type });
      await page.waitForTimeout(150);
    }

    await schemaDialog.locator('[name="button-OK"].enabled').waitFor({ timeout: 5000 });
    await H.clickDialogButton(page, 'dialog-Create-a-new-schema', 'button-OK');
    await page.waitForTimeout(2500);

    await H.searchList(page, schemaName);
    expect(await H.listHasCard(page, schemaName)).toBe(true);

    await page.evaluate(() => document.querySelectorAll('.d4-dialog')
      .forEach((d) => (d.querySelector('[name="button-CANCEL"]') as HTMLElement | null)?.click()));
    await page.waitForTimeout(500);
    await H.cardContextMenu(page, schemaName, 'Edit');
    await page.waitForTimeout(1500);

    const content = await page.evaluate(() => {
      const dialog = Array.from(document.querySelectorAll('.d4-dialog'))
        .find((d) => d.querySelector('.d4-dialog-title')?.textContent?.trim() === 'Edit schema');
      if (!dialog) return null;
      const table = dialog.querySelector('table.d4-item-table');
      const names = Array.from(table?.querySelectorAll('input[name="input-Name"]') || []) as HTMLInputElement[];
      return {
        assoc: dialog.querySelector('[name="div-Associated-with-"]')?.textContent?.trim(),
        names: names.map((n) => n.value),
      };
    });
    expect(content).not.toBeNull();
    expect(content!.assoc).toContain(typeName);
    expect(content!.names).toEqual(['rating', 'notes', 'verified', 'review_date']);
    await page.evaluate(() => {
      const dialog = Array.from(document.querySelectorAll('.d4-dialog'))
        .find((d) => d.querySelector('.d4-dialog-title')?.textContent?.trim() === 'Edit schema');
      (dialog?.querySelector('[name="button-CANCEL"]') as HTMLElement | null)?.click();
    });
    await page.waitForTimeout(500);

    const apiRows = await page.evaluate(async (schemaName) => {
      const g = (window as any).grok;
      const schema = (await g.dapi.stickyMeta.getSchemas()).find((s: any) => s.name === schemaName);
      return schema ? schema.properties.map((p: any) => ({ name: p.name, type: p.type })) : null;
    }, schemaName);
    expect(apiRows).toEqual([
      { name: 'rating', type: 'int' },
      { name: 'notes', type: 'string' },
      { name: 'verified', type: 'bool' },
      { name: 'review_date', type: 'datetime' },
    ]);

    await H.searchList(page, schemaName);
    await H.cardContextMenu(page, schemaName, 'Delete');
    await page.waitForTimeout(600);
    await page.evaluate(() => (document.querySelector('[name="button-DELETE"]') as HTMLElement | null)?.click());
    await page.waitForTimeout(2000);
    await H.searchList(page, schemaName);
    expect(await H.listHasCard(page, schemaName)).toBe(false);

    await H.gotoStickyMeta(page, 'Types');
    await H.searchList(page, typeName);
    await H.cardContextMenu(page, typeName, 'Delete');
    await page.evaluate(() => (document.querySelector('[name="button-DELETE"]') as HTMLElement | null)?.click());
    await page.waitForTimeout(2000);
    await H.searchList(page, typeName);
    expect(await H.listHasCard(page, typeName)).toBe(false);
  } finally {

    await H.apiDeleteSchema(page, schemaName).catch(() => {});
    await page.evaluate(async (typeName) => {
      try {

        const g = (window as any).grok;

        void g; void typeName;
      } catch {  }
    }, typeName).catch(() => {});
  }
});
