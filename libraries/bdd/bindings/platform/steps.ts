/* Platform base steps: setup through the JS API (the openers of @datagrok-libraries/test keep the
   provenance tags the UI would set). Viewer steps live in the `viewers` tier. */
import type {Page} from '@playwright/test';
import {openTableFromFile} from '@datagrok-libraries/test/src/playwright/openers.js';
import {DatasetEntry, Given} from '../../src/registry.js';

export const openDataset = Given('user opens {dataset} dataset', async (page: Page, dataset: DatasetEntry) => {
  await openTableFromFile(page, dataset.path);
  await page.locator('[name="viewer-Grid"]').first().waitFor();
}, {tier: 'api', description: 'OpenFile through the JS API — provenance as in the UI'});

export const switchTableView = Given('user switches to (the ){string} table view', async (page: Page, name: string) => {
  await page.evaluate((n) => {
    const grok = (window as any).grok;
    const views = Array.from(grok.shell.tableViews) as any[];
    const view = views.find((x) => String(x.dataFrame?.name).toLowerCase() === n.toLowerCase());
    if (!view)
      throw new Error(`no table view for "${n}"; open: ${views.map((x) => x.dataFrame?.name).join(', ')}`);
    grok.shell.v = view;
  }, name);
  await page.waitForFunction((n) => (window as any).grok.shell.tv?.dataFrame?.name?.toLowerCase() === n.toLowerCase(), name);
}, {tier: 'api', description: 'by table name — the view of a dataset opened earlier in the scenario'});
