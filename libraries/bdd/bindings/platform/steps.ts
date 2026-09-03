/* Platform base steps: setup through the JS API (the openers of @datagrok-libraries/test keep the
   provenance tags the UI would set). Viewer steps live in the `viewers` tier. */
import type {Page} from '@playwright/test';
import {openTableFromFile} from '@datagrok-libraries/test/src/playwright/openers.js';
import {DatasetEntry, Given} from '../../src/registry.js';

export const openDataset = Given('user opens {dataset} dataset', async (page: Page, dataset: DatasetEntry) => {
  await openTableFromFile(page, dataset.path);
  await page.locator('[name="viewer-Grid"]').first().waitFor();
}, {tier: 'api', description: 'OpenFile through the JS API — provenance as in the UI'});
