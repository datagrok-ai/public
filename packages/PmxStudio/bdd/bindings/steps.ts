/* The steps only this package can define: how to reach its app. Everything generic (clicks,
   typing, selects, assertions) is the library's — `grok-bdd list-steps` prints it all. */
import type {Page} from '@playwright/test';
import {Given} from '@datagrok-libraries/bdd';

const open = async (page: Page, path?: string) => {
  await page.evaluate(async (p) => {
    const g = (window as any).grok;
    g.shell.addView(await g.functions.call('PmxStudio:pmxStudioApp', p ? {path: p} : {}));
  }, path);
  await page.locator('.pmx-stages').filter({visible: true}).waitFor({timeout: 60000});
};

export const openApp = Given('user opens PMX Studio', (page: Page) => open(page),
  {tier: 'api', description: 'the app\'s Plan view, added to the shell; needs this package published on the stand'});

export const openAnalyze = Given('user opens PMX Studio on Analyze', (page: Page) => open(page, '/analyze'),
  {tier: 'api', description: 'the app\'s Analyze view, added to the shell'});
