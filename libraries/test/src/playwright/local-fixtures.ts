import {test as base, Page, BrowserContext} from '@playwright/test';
import {openLocalDatagrok} from './spec-login';
import {installRenderWait} from './viewers';

/**
 * Worker-scoped local-mode client: the platform boots once per worker and every spec in
 * that worker runs on the booted page. Boot is the single largest cost in a spec (a cold
 * client is ~9 s, a cached one ~2.5 s, against ~3 s for a whole property suite), and a
 * per-file boot pays it again for every file. Local mode is what makes sharing safe —
 * there is no server session to leak between specs, so the only state to undo is the
 * shell's own, which `closeAll()` does.
 *
 * Use it exactly like the default fixture, taking `localPage` instead of `page`:
 *
 *     import {test} from '@datagrok-libraries/test/src/playwright/local-fixtures';
 *     test('my spec', async ({localPage}) => { ... });
 *
 * A spec that needs a server (dapi persistence, queries, file shares) must not use this —
 * take the standard `page` fixture and `loginToDatagrok` instead.
 */
export const test = base.extend<{localPage: Page}, {localContext: BrowserContext, localBoot: Page}>({
  localContext: [async ({browser}, use) => {
    const context = await browser.newContext({viewport: {width: 1920, height: 1080}});
    await use(context);
    await context.close();
  }, {scope: 'worker'}],

  localBoot: [async ({localContext}, use) => {
    const page = await localContext.newPage();
    await openLocalDatagrok(page);
    await use(page);
    await page.close();
  }, {scope: 'worker'}],

  localPage: async ({localBoot}, use) => {
    await use(localBoot);
    // Leave the shell as the next spec expects to find it. Views and tables are the only
    // state local mode can carry across a spec boundary.
    await localBoot.evaluate(() => (window as any).grok.shell.closeAll());
  },
});

export {expect} from '@playwright/test';
export {installRenderWait};
