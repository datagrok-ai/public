/* Platform base steps: setup through the JS API (the openers of @datagrok-libraries/test keep the
   provenance tags the UI would set). Viewer steps live in the `viewers` tier. */
import type {Page} from '@playwright/test';
import {openTableFromFile} from '@datagrok-libraries/test/src/playwright/openers.js';
import {DatasetEntry, Given} from '../../src/registry.js';

/** Opening a table starts semantic-type detection in the background (package detectors, a few
 * hundred ms on a molecule table); the platform reports its end on the global event bus, and the
 * step is over only then — otherwise that work lands on whatever step comes next. */
export const openDataset = Given('user opens {dataset} dataset', async (page: Page, dataset: DatasetEntry) => {
  await page.evaluate(() => {
    const w = window as any;
    if (w.__bddDetected)
      return;
    w.__bddDetected = [];
    w.grok.events.onEvent('ddt-semantic-type-detected').subscribe((a: any) => {
      w.__bddDetected = [...w.__bddDetected.slice(-19), a?.args?.dataFrame?.dart];
    });
  });
  await openTableFromFile(page, dataset.path);
  await page.locator('[name="viewer-Grid"]').first().waitFor();
  await page.waitForFunction(() => {
    const w = window as any;
    return w.__bddDetected.includes(w.grok.shell.tv?.dataFrame?.dart);
  }, undefined, {timeout: 15000}).catch(() => {
    throw new Error(`${dataset.name}: semantic types were not detected within 15 s (is auto-detection on?)`);
  });
}, {tier: 'api', description: 'OpenFile through the JS API — provenance as in the UI; done when semantic types are detected'});

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
