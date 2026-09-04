/* The demo shell (src/demo.ts, src/nav.ts) as a context: "user opens the … demo page" enters it,
   after which generic kinds resolve inside the sub-demo's content first — so "Tables tree node"
   is the page's tree and not the navigation's, which keeps its own name here. The page's
   "name = value" lines (readout() in src/pages/common.ts) are a kind of this package alone. */
import type {Page} from '@playwright/test';
import {context, Given, kind} from '@datagrok-libraries/bdd';

declare const grok: any;

const APP_PATH = '/apps/U2demo/U2Demo';

/** Sub-demos by label, as nav.ts lists them; the route's last segment is the leaf id. */
const PAGES: Record<string, string> = {
  'overview': 'start/overview',
  'all inputs': 'inputs/all-inputs',
  'basic inputs': 'inputs/basic-inputs',
  'range slider': 'inputs/range-slider',
  'multi-select': 'inputs/multi-select',
  'async': 'inputs/async',
  'layout': 'containers/layout',
  'popups': 'containers/popups',
  'lists': 'collections/lists',
  'trees': 'collections/trees',
  'cards': 'display/cards',
  'feedback': 'display/feedback',
  'tables': 'display/tables',
  'sections & wizard': 'display/sections',
  'message input': 'display/messaging',
  'form': 'forms/form',
  'property grid': 'forms/property-grid',
  'object form': 'forms/objectform',
  'functions': 'forms/funcs',
  'funccalls': 'forms/func-convergence',
  'run history': 'forms/func-history',
  'dataframes': 'platform/dataframes',
  'files': 'platform/files',
  'entities': 'platform/entities',
  'spaces': 'platform/spaces',
  'molecules': 'platform/molecules',
  'bridge': 'platform/bridge',
  'msa workbench': 'automation/msa-workbench',
};

function routeOf(name: string): {route: string; label: string} {
  const key = name.trim().toLowerCase();
  const entry = Object.entries(PAGES).find(([label, r]) => label === key || r.endsWith('/' + key));
  if (!entry)
    throw new Error(`unknown sub-demo "${name}" — one of: ${Object.keys(PAGES).join(', ')}`);
  return {route: entry[1], label: entry[0]};
}

export const demo = context('U2 Demo', {selector: '.u2demo-content', aliases: ['demo page', 'sub-demo']});
demo.element('demo navigation', {selector: '.u2demo-nav', aliases: ['navigation tree', 'nav tree', 'demo tree']});
demo.element('user picker', {selector: '[data-u2="typeahead"]:has(input[placeholder="Search users…"])'});
demo.element('group picker', {selector: '[data-u2="typeahead"]:has(input[placeholder="Search groups…"])'});
demo.element('compound picker', {selector: '[data-u2="typeahead"]:has(input[placeholder="Find a compound…"])'});

kind('readout', {
  aliases: ['read-out'],
  selector: '.u2demo-status',
  match: ['label'],
  labelSelector: 'span:first-child',
  parts: {label: 'span:first-child', value: 'span:last-child'},
  description: 'a "name = value" line of a sub-demo; "value of <name> readout" is the value alone',
});

/** Opens a sub-demo by its route and waits for `ready`. A warm page opens the app view in place —
 * a reload would close the tables earlier steps opened — a cold one navigates. */
export async function openSubDemo(page: Page, route: string, ready: string): Promise<void> {
  const path = `${APP_PATH}/${route}`;
  const inShell = await page.evaluate(() => typeof (window as any).grok?.shell?.closeAll === 'function').catch(() => false);
  if (inShell) {
    await page.evaluate(async (p) => {
      grok.shell.addView(await grok.functions.call('U2Demo:u2DemoApp', {path: p}));
    }, path);
  }
  else
    await page.goto(path, {waitUntil: 'domcontentloaded'});
  await page.locator(ready).first().waitFor({timeout: 60000});
}

export const openDemoPage = Given('user opens the {string} demo page', async (page: Page, name: string) => {
  const {route, label} = routeOf(name);
  await openSubDemo(page, route, '.u2demo-content .u2demo-page');
  // the shell's status bar names the sub-demo once the view is docked and laid out
  const status = new RegExp(`/ ${label.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')}$`, 'i');
  await page.locator('.d4-view-status-panel').filter({hasText: status}).first().waitFor({state: 'attached', timeout: 60000});
}, {tier: 'ui', enters: 'U2 Demo', description: 'a sub-demo by its navigation label (or leaf id); needs this package published'});
