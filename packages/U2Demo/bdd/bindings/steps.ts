import type {Page} from '@playwright/test';
import {Given} from '@datagrok-libraries/bdd';
import {openSubDemo} from './demo.js';

export const openWorkbench = Given('user opens the MSA workbench', async (page: Page) => {
  await openSubDemo(page, 'MSA workbench', '[data-u2-name="msaWorkbench"]');
}, {tier: 'ui', enters: 'MSA workbench', description: 'needs this package published on the stand'});
