import type {Page} from '@playwright/test';
import {Given} from '@datagrok-libraries/bdd';
import {openSubDemo} from './demo.js';

/** The sub-demo (src/pages/msa-workbench.ts) by its route under the app. */
export const MSA_WORKBENCH_ROUTE = 'automation/msa-workbench';

export const openWorkbench = Given('user opens the MSA workbench', async (page: Page) => {
  await openSubDemo(page, MSA_WORKBENCH_ROUTE, '[data-u2-name="msaWorkbench"]');
}, {tier: 'ui', enters: 'MSA workbench', description: 'needs this package published on the stand'});
