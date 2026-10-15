/* The demo shell (src/demo.ts, src/nav.ts) as a context: "user opens the … demo page" enters it,
   after which generic kinds resolve inside the sub-demo's content first — so "Tables tree node"
   is the page's tree and not the navigation's, which keeps its own name here. */
import type {Page} from '@playwright/test';
import {context, Given} from '@datagrok-libraries/bdd';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';

declare const grok: any;

export const demo = context('U2 Demo', {selector: '.u2demo-content', aliases: ['demo page', 'sub-demo']});
demo.element('demo navigation', {selector: '.u2demo-nav', aliases: ['navigation tree', 'nav tree', 'demo tree']});

/** Opens a sub-demo by its navigation label (nav.ts `leafForPath`) and waits for `ready`. The app view
 * opens in place — a reload would close the tables earlier steps opened — on a page in the shell. */
export async function openSubDemo(page: Page, label: string, ready: string): Promise<void> {
  if (!await page.evaluate(() => typeof (window as any).grok?.shell?.closeAll === 'function').catch(() => false))
    await loggedIn(page);
  await page.evaluate(async (path) => {
    grok.shell.addView(await grok.functions.call('U2Demo:u2DemoApp', {path}));
  }, '/' + encodeURIComponent(label.trim()));
  await page.locator(ready).first().waitFor({timeout: 60000});
}

export const openDemoPage = Given('user opens the {string} demo page', async (page: Page, label: string) => {
  await openSubDemo(page, label, '.u2demo-content .u2demo-page');
  // the shell's status bar names the sub-demo once the view is docked and laid out
  const status = new RegExp(`/ ${label.trim().replace(/[.*+?^${}()|[\]\\]/g, '\\$&')}$`, 'i');
  await page.locator('.d4-view-status-panel').filter({hasText: status}).first().waitFor({state: 'attached', timeout: 60000});
}, {tier: 'ui', enters: 'U2 Demo', description: 'a sub-demo by its navigation label; needs this package published'});
