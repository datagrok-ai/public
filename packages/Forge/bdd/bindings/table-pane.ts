/* The table pane Forge models (`table-models-panel.ts`): a named table as the current object, the place of the pane
   in the context panel, and the feature's own lines in it. The pane is built in `ui.wait` each time it is shown, so
   every claim polls. */
import type {Page} from '@playwright/test';
import {Given, Then} from '@datagrok-libraries/bdd';
import {el, type ElementRef, expect, locate, viewers} from '@datagrok-libraries/bdd/runtime';
import {listOf} from './grok-s.js';
import {contextPanelTitle} from './results-state.js';

declare const grok: any;

/** 'current' when the named open table is the current object and the context panel's title is its name, or what is
 * shown instead. With [set], a forced set of the table is re-applied while another object is current. */
async function currentTable(page: Page, name: string, set = false): Promise<string> {
  const [object, title] = await Promise.all([viewers.evaluate(page, ([n, s]: [string, boolean]) => {
    const table = (window as any).__bdd.tableNamed(n);
    if (grok.shell.o?.dart === table.dart)
      return 'the table';
    // Not frozen: a frozen set drops the next unforced change of the current object for a second.
    if (s)
      grok.shell.setCurrentObject(table, false, true);
    return String(grok.shell.o?.name ?? grok.shell.o);
  }, [name, set]), contextPanelTitle(page)]);
  return object === 'the table' && title === name ? 'current' :
    `the current object is ${object}, the context panel's title "${title}"`;
}

/* LIB-PROMOTE: the manual's corner click does not make a table current again (a platform issue); a step that sets a
   named table as the current object stands for "the context panel shows the table". A library step "make X the
   current object" would re-apply the set until it stays, as this one does: a view opened right before it sets its own
   current object after the set. */
export const tableIsCurrentObject = Given('the table {string} is the current object', async (page: Page, name: string) => {
  await expect.poll(() => currentTable(page, name, true), {message: `the table "${name}" as the current object`})
    .toBe('current');
}, {tier: 'api', description: 'grok.shell.setCurrentObject(the open table, not frozen, forced), re-applied while ' +
  'another object is current; done when it is the current object and the context panel\'s title is its name'});

/* LIB-PROMOTE: the order of the context panel's panes, which the platform publishes only as the DOM order. */
export const paneComesAfter = Then('the {string} pane in context panel should come after the {string} pane', async (page: Page,
  pane: string, before: string) => {
  await expect.poll(async () => {
    const order = await page.locator('.grok-prop-panel .d4-accordion-pane-header').evaluateAll((all) =>
      all.map((h) => (h.textContent ?? '').trim()));
    const [i, j] = [order.indexOf(pane), order.indexOf(before)];
    return i >= 0 && j >= 0 && i > j ? 'after' : order.join(', ');
  }, {message: `"${pane}" after "${before}" in the context panel's panes`}).toBe('after');
}, {description: 'both pane headers are in the context panel and the first comes later in its order'});

/* The run's own lines: a fixture name ends in the run's {time}, so the lines of the user's models and of other runs
   (a killed run's model under an hour old) are ignored. */
export const modelLinesList = Then('the model lines of {element} should list exactly {string} among those ending in {string}',
  async (page: Page, target: ElementRef, list: string, suffix: string) => {
    const names = await locate(page, el(`name of model line in ${target.phrase}`));
    await expect.poll(() => names.evaluateAll((all, s) =>
      all.map((n) => (n.textContent ?? '').trim()).filter((n) => n.startsWith('forge-test-') && n.endsWith(s)).join(', '),
    suffix), {message: `the forge-test- model lines of ${target.phrase} ending in "${suffix}", in order`})
      .toBe(listOf(list).join(', '));
  }, {description: 'the names of the model lines that start with forge-test- and end in the suffix (the run\'s ' +
  '"-{time}"), in their order (comma-separated; "" for none); the lines of other models are ignored'});
