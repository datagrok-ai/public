/* The steps only this package can define: how to reach its app. Everything generic (clicks,
   typing, selects, assertions) is the library's — `grok-bdd list-steps` prints it all; the model
   steps are in diff-studio.ts. */
import type {Page} from '@playwright/test';
import {Given} from '@datagrok-libraries/bdd';
import {expect, takeErrors} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/** The app's route lands on its hub — templates and the library as cards — and never builds a
 * model's UI; a feature about a model asks for the model (diff-studio.ts). */
export const openApp = Given('user opens the Diff Studio app', async (page: Page) => {
  await page.goto('/apps/DiffStudio', {waitUntil: 'domcontentloaded', timeout: 180000});
  takeErrors(page);
  await page.locator('[name="Browse"]').waitFor({timeout: 180000});
  await expect.poll(() => page.evaluate(() => String(grok.shell.v?.name ?? '')),
    {message: 'the current view', timeout: 60000}).toBe('Diff Studio');
}, {tier: 'ui', description: 'needs this package published on the stand; done when the hub is the current view'});
