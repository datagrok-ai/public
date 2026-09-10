/* The session is the storage state global setup writes (@datagrok-libraries/test); this step only
   lands in the shell — once per browser page, the scenarios of a feature share it — and applies the
   automation-friendly shell settings every suite applies, plus the viewer runtime: every viewer the
   page holds or adds later renders immediately (no debounce, no animation-frame wait). */
import type {Page} from '@playwright/test';
import {Given} from '../../src/registry.js';
import {takeErrors} from '../../src/runtime/harness.js';
import {installViewerRuntime, takeBalloons} from '../../src/runtime/viewers.js';

declare const grok: any;

export const loggedIn = Given('user is logged in', async (page: Page) => {
  const inShell = await page.evaluate(() => typeof (window as any).grok?.shell?.closeAll === 'function').catch(() => false);
  if (!inShell) {
    // a dev stand's pub serve can take minutes to hand out the bundle while it recompiles or is
    // starved: that is a delay once per page, not a failure of the feature
    const start = Date.now();
    await page.goto('/', {waitUntil: 'domcontentloaded', timeout: 180000});
    await page.locator('[name="Browse"]').first().waitFor({timeout: 180000});
    const seconds = Math.round((Date.now() - start) / 1000);
    if (seconds >= 30)
      console.warn(`bdd: the shell took ${seconds} s to load (a dev stand serving a bundle it is recompiling?)`);
  }
  await page.evaluate(() => {
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = true;
  });
  // closeAll re-adds the Home view asynchronously; a table opened before it lands ends up behind it
  await page.waitForFunction(() => grok.shell.v?.type === 'datagrok', null, {timeout: 60000});
  await installViewerRuntime(page);
  // what the stand logs or shows while booting (a broken package's autostart, "Debugging
  // packages") is not the scenario's
  takeErrors(page);
  await takeBalloons(page);
}, {tier: 'ui', description: 'the error floor starts here: "no errors should have been logged" counts from this step'});
