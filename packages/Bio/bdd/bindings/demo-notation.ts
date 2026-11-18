// Steps of the notation demo clips (features/demos/01-03c), all set-up a filmed demo leaves out: the Browse panel
// narrowed so the sequences get the room, and the Renderer Settings of a sequence column (Multiline Rendering, Max
// Monomer Length) as the column's settings set them.
import type {Page} from '@playwright/test';
import {Given} from '@datagrok-libraries/bdd';
import {atFeatureEnd, expect, silent, viewers} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/** The splitter at the right edge of the docked Browse panel (alone, or in a column with the Dashboards panel). */
const BROWSE_SPLITTER = [
  '.panel-base:has(> .panel-content .grok-view-browse) + .splitbar-vertical',
  '.splitter-container-column:has(> .panel-base > .panel-content .grok-view-browse) + .splitbar-vertical',
].join(', ');

/** The Browse panel docks at its default width each time it is shown (`dockBrowse`): narrowed by its splitter for the
 * clip, it is docked afresh at feature end (hidden, then shown again), so a later feature meets the default width. */
export const browsePanelNarrower = Given('the browse panel is {int} px narrower', async (page: Page, px: number) => {
  silent(page);
  const tree = page.locator('.grok-view-browse [role="tree"], .layout-browse [role="tree"]').filter({visible: true}).first();
  await expect(tree, 'the browse tree').toBeVisible();
  const before = (await tree.boundingBox())!.width;
  // the splitbar right of the panel: by its place (at the tree's right edge), whichever container the dock put it in
  const treeBox = (await tree.boundingBox())!;
  let box: {x: number, y: number, width: number, height: number} | null = null;
  await expect.poll(async () => {
    const bars = page.locator(`${BROWSE_SPLITTER}, .splitbar-vertical, .splitbar-horizontal`).filter({visible: true});
    box = null;
    for (const b of await bars.all()) {
      const r = await b.boundingBox();
      if (r && r.height > treeBox.height / 2 && r.width < 20 && Math.abs(r.x - (treeBox.x + treeBox.width)) < 40 &&
        (!box || Math.abs(r.x - (treeBox.x + treeBox.width)) < Math.abs(box.x - (treeBox.x + treeBox.width))))
        box = r;
    }
    return box !== null;
  }, {message: 'the splitter right of the browse panel', timeout: 15000}).toBe(true);
  const at = box!;
  await viewers.dragFrom(page, {x: at.x + at.width / 2, y: at.y + at.height / 2}, {dx: -px, dy: 0});
  atFeatureEnd(page, async () => {
    await page.evaluate(() => {
      if (!grok.shell.windows.showBrowse)
        return;
      grok.shell.windows.showBrowse = false;
      grok.shell.windows.showBrowse = true;
    }).catch(() => undefined);
  });
  await expect.poll(async () => Math.round(before - ((await tree.boundingBox())?.width ?? before)),
    {message: 'how much narrower the browse tree is'}).toBeGreaterThan(px / 2);
  await viewers.settleAll(page);
}, {tier: 'ui', description: 'set-up a filmed demo leaves out: the splitter right of the Browse panel dragged left by that many px; the panel is docked afresh at its default width at feature end'});

/** Renderer Settings > Multiline Rendering of a sequence column (`representations.ts`): a sequence longer than its
 * cell wraps onto more lines, its bonds bracketed within a line and numbered where a line break cuts them. */
export const sequenceMultiline = Given('the {string} column is drawn on several lines', async (page: Page, column: string) => {
  silent(page);
  await page.evaluate((name) => {
    const col = grok.shell.t?.col(name);
    if (!col)
      throw new Error(`the current table has no "${name}" column`);
    col.setTag('renderMultiline', 'true');
    col.dataFrame.fireValuesChanged();
  }, column);
  await viewers.settleAll(page);
}, {tier: 'api', description: 'set-up a filmed demo leaves out: the column\'s Multiline Rendering, as its Renderer Settings set it'});

/** Renderer Settings > Max Monomer Length of a sequence column (`representations.ts`): a longer monomer name is cut
 * with an ellipsis; the package default keeps four letters, so Phe_34diCl shows as "Phe…". */
export const monomerNameLength = Given('the monomer names of the {string} column are cut at {int} characters',
  async (page: Page, column: string, length: number) => {
    silent(page);
    await page.evaluate(([name, n]) => {
      const col = grok.shell.t?.col(name);
      if (!col)
        throw new Error(`the current table has no "${name}" column`);
      col.temp['.mm.cellRenderer.maxMonomerLength'] = String(n);
      col.temp['.mm.cellRenderer.settingsChanged'] = '1';
      col.dataFrame.fireValuesChanged();
    }, [column, length] as const);
    await viewers.settleAll(page);
  }, {tier: 'api', description: 'set-up a filmed demo leaves out: the column\'s Max Monomer Length, as its Renderer Settings set it'});
