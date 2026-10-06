/* Multi-host controls arrive with one pane. What a dropped tabs used to be: an empty box with a
   'Drop here' hint, no tab strip, nothing to click — the user could not see what they had made. Now
   the pane is part of the insert, so it renders, it is selectable, the first 'Add tab' counts it,
   and one Ctrl+Z still takes the whole thing back. */
import {ok, shot} from '../local.mjs';
import {dropControl, expandTree, focusTree, newForm, nodeCount, row, selectRow, statusText, waitStatus}
  from '../lib.mjs';

export async function fixture(page) {
  await newForm(page);
}

async function checkSeededPanes(page) {
  const blank = await nodeCount(page);
  await dropControl(page, 'u2-tabs', 'tabs', 'form1');
  await waitStatus(page, 'tabs1');
  const tabLabels = () => page.evaluate(() => [...document.querySelectorAll(
    '.u2-designer-surface .u2-tabs-label')].map((l) => l.textContent));
  const seeded = await tabLabels();
  await shot(page, 'seeded-panes-1-seeded-tab');
  ok('seeded-panes/1a/a-dropped-tabs-renders-its-first-tab',
    JSON.stringify(seeded) === JSON.stringify(['Tab 1']) && await nodeCount(page) === blank + 2,
    `tabs=${JSON.stringify(seeded)} ${await nodeCount(page)} nodes`);

  await expandTree(page);
  const pane = await selectRow(page, 'panel1');
  ok('seeded-panes/1b/the-seeded-pane-is-a-node-like-any-other', pane.startsWith('panel1 ('), `"${pane}"`);

  await row(page, 'tabs1').click({button: 'right'});
  await page.waitForSelector('.u2-menu-item', {timeout: 5000});
  await page.locator('.u2-menu-item').filter({hasText: /^Add tab$/}).first().click();
  await page.waitForTimeout(300);
  const added = await tabLabels();
  ok('seeded-panes/1c/the-first-add-tab-counts-the-seeded-pane',
    JSON.stringify(added) === JSON.stringify(['Tab 1', 'Tab 2']) && /panel2 · \d+ nodes/
      .test(await statusText(page)), `tabs=${JSON.stringify(added)} status="${await statusText(page)}"`);

  await focusTree(page, 'tabs1');
  await page.keyboard.press('Control+z');
  await page.waitForTimeout(200);
  await page.keyboard.press('Control+z');
  await page.waitForTimeout(300);
  ok('seeded-panes/1d/one-undo-per-gesture-takes-the-host-and-its-panes-back',
    await nodeCount(page) === blank && (await tabLabels()).length === 0,
    `${await nodeCount(page)} nodes, tabs=${JSON.stringify(await tabLabels())}`);

  await dropControl(page, 'u2-accordion', 'accordion', 'form1');
  await waitStatus(page, 'accordion1');
  const panes = await page.evaluate(() => [...document.querySelectorAll(
    '.u2-designer-surface .u2-accordion-title')].map((t) => t.textContent));
  await shot(page, 'seeded-panes-1-seeded-pane');
  await focusTree(page, 'accordion1');
  await page.keyboard.press('Control+z');
  await page.waitForTimeout(300);
  ok('seeded-panes/1e/an-accordion-opens-on-pane-1-and-undoes-whole',
    JSON.stringify(panes) === JSON.stringify(['Pane 1']) && await nodeCount(page) === blank,
    `panes=${JSON.stringify(panes)} ${await nodeCount(page)} nodes after undo`);
}

export const checks = [
  {id: 'seeded-panes/1 seeded panes: a dropped tabs shows Tab 1, Add tab counts it, one undo takes the host whole',
    run: checkSeededPanes},
];
