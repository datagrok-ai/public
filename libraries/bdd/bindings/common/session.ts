/* The session is the storage state global setup writes (@datagrok-libraries/test); this step only
   lands in the shell — once per browser page, the scenarios of a feature share it — and applies the
   automation-friendly shell settings every suite applies, plus the viewer runtime: every viewer the
   page holds or adds later renders immediately (no debounce, no animation-frame wait). */
import type {Page} from '@playwright/test';
import {Given} from '../../src/registry.js';
import {takeErrors} from '../../src/runtime/harness.js';
import {installViewerRuntime, takeBalloons} from '../../src/runtime/viewers.js';
import * as guide from '../../src/runtime/guide.js';

declare const grok: any;

const homeNotWaitedFor = new WeakSet<Page>();

/* What a PowerPack Home widget logs while it loads would land on whichever scenario runs by then. The
   widgets host shows up within a second of the shell; a stand without it is not waited for again. */
async function homeWidgetsSettled(page: Page): Promise<void> {
  if (homeNotWaitedFor.has(page))
    return;
  const hasHost = await page.waitForFunction(() => grok.shell.v?.root?.querySelector('.power-pack-widgets-host') != null,
    null, {timeout: 10000}).then(() => true, () => false);
  const settled = hasHost && await page.waitForFunction(() => {
    const contents = [...grok.shell.v.root.querySelectorAll('.power-pack-widgets-host .power-pack-widget-content')] as HTMLElement[];
    return contents.length > 0 && contents.every((c) => c.children.length > 0 && c.querySelector('.grok-loader') == null);
  }, null, {timeout: 30000}).then(() => true, () => false);
  if (!settled) {
    homeNotWaitedFor.add(page);
    console.warn(`bdd: ${hasHost ? 'the Home widgets did not finish loading in 30 s' : 'no PowerPack Home widgets'}; not waited for again on this page`);
  }
}

export const loggedIn = Given('user is logged in', async (page: Page) => {
  // a worker runs one spec after another on the same page, so what a feature leaves behind (an open
  // dialog, a docked panel, a sticky option) reaches the next one; BDD_FRESH_PAGE starts each
  // feature from a reload, at the cost of a shell load per feature
  guide.silent(page);
  const inShell = process.env.BDD_FRESH_PAGE === undefined &&
    await page.evaluate(() => typeof (window as any).grok?.shell?.closeAll === 'function').catch(() => false);
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
  await page.evaluate((simple) => {
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = simple;
  }, guide.shellSimpleMode());
  // closeAll re-adds the Home view asynchronously; a table opened before it lands ends up behind it
  await page.waitForFunction(() => grok.shell.v?.type === 'datagrok', null, {timeout: 60000});
  await homeWidgetsSettled(page);
  await installViewerRuntime(page);
  // what the stand logs or shows while booting (a broken package's autostart, "Debugging
  // packages") is not the scenario's
  takeErrors(page);
  await takeBalloons(page);
}, {tier: 'ui', description: 'the error floor starts here: "no errors should have been logged" counts from this step'});
