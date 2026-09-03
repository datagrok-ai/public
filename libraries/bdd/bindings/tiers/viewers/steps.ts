/* The `viewers` tier: adding and checking viewers on the current table view. A project opts in
   with `{"tiers": ["viewers"]}` in its bdd.config.json. */
import {expect, Page} from '@playwright/test';
import {Then, When} from '../../../src/registry.js';

declare const grok: any;

export const openToolbox = When('user opens toolbox', async (page: Page) => {
  const shown = await page.evaluate(() => grok.shell.windows.showToolbox as boolean);
  if (!shown)
    await page.locator('[name="Toolbox"]').first().click();
  await page.locator('[name="div-section--Viewers"]').first().waitFor();
}, {tier: 'ui'});

export const viewerAdded = Then('{viewer} viewer should be added to the open tableview', async (page: Page, viewer: string) => {
  await expect(page.locator(`[name="viewer-${viewer.replace(/\s+/g, '-')}" i]`).first()).toBeVisible();
  const types: string[] = await page.evaluate(() => Array.from(grok.shell.tv.viewers).map((v: any) => String(v.type)));
  expect(types.map((t) => t.toLowerCase())).toContain(viewer.toLowerCase());
}, {tier: 'api', description: 'both the DOM and grok.shell.tv.viewers'});
