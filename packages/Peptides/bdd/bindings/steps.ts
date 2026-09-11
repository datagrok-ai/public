/* The steps only this package can define: how to reach its app. Everything generic (clicks,
   typing, selects, assertions) is the library's — `grok-bdd list-steps` prints it all. */
import type {Page} from '@playwright/test';
import {Given} from '@datagrok-libraries/bdd';

/** The app's route on the stand — adjust: one app is `/apps/<Package>`, several are
 * `/apps/<Package>/<App>`. */
export const APP_PATH = '/apps/peptides';

export const openApp = Given('user opens the Peptides app', async (page: Page) => {
  await page.goto(APP_PATH, {waitUntil: 'domcontentloaded'});
  await page.locator('[data-u2-name="Peptides"]').waitFor({timeout: 60000});
}, {tier: 'ui', enters: 'Peptides app', description: 'needs this package published on the stand'});
