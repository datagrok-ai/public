/* The session is the storage state global setup writes (@datagrok-libraries/test); this step only
   lands in the shell — once per browser page, the scenarios of a feature share it — and applies the
   automation-friendly shell settings every suite applies. */
import type {Page} from '@playwright/test';
import {Given} from '../../src/registry.js';

declare const grok: any;

export const loggedIn = Given('user is logged in', async (page: Page) => {
  const inShell = await page.evaluate(() => typeof (window as any).grok?.shell?.closeAll === 'function').catch(() => false);
  if (!inShell) {
    await page.goto('/', {waitUntil: 'domcontentloaded'});
    await page.locator('[name="Browse"]').first().waitFor({timeout: 60000});
  }
  await page.evaluate(() => {
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = true;
  });
  // closeAll re-adds the Home view asynchronously; a table opened before it lands ends up behind it
  await page.waitForFunction(() => grok.shell.v?.type === 'datagrok', null, {timeout: 60000});
}, {tier: 'ui'});
