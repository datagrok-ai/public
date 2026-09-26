/* The steps only Dendrogram can define: the tree table its newick viewers show (TreeHelper.newickToDf:
   a row per node, with node, parent and leaf columns). */
import {expect, Page} from '@playwright/test';
import {Then, When} from '@datagrok-libraries/bdd';
import {type ElementRef, viewers} from '@datagrok-libraries/bdd/runtime';

/** Open and double-click in Browse show a newick file's preview; the view the file handler builds
 * (a tree table with a Dendrogram) is reached from the handler itself. */
export const openNewickFile = When('user opens the newick file {string} with its file handler', async (page: Page, path: string) => {
  await page.evaluate(async (p) => {
    const g = (window as any).grok;
    const text = await g.dapi.files.readAsText(p);
    await g.functions.call('Dendrogram:importNewick', {fileContent: text});
  }, path);
}, {tier: 'api', description: 'the file\'s text handed to Dendrogram:importNewick, the handler registered for .nwk and .newick'});

export const tableOfTreeLeaves = Then('the table of {widget} should hold a tree with leaves {string}', async (page: Page, widget: ElementRef, list: string) => {
  const names = await viewers.onViewer(page, widget, (e) => {
    const DG = (window as any).DG;
    const host = e.closest('[name^="viewer-"], .d4-viewer') ?? e;
    const df = DG.Widget.find(host)?.dataFrame;
    if (!df)
      throw new Error('the element is not a viewer bound to a table');
    const node = df.col('node');
    const leaf = df.col('leaf');
    if (!node || !leaf)
      throw new Error(`the viewer's table is not a tree table (node, leaf); it has: ${df.columns.names().join(', ')}`);
    const res: string[] = [];
    for (let i = 0; i < df.rowCount; i++)
      if (leaf.get(i))
        res.push(node.get(i));
    return res;
  });
  expect(names, 'the leaf rows of the viewer\'s tree table').toEqual(list.split(/\s*,\s*/).filter(Boolean));
}, {description: 'the leaf rows of the tree table a newick viewer is bound to, in table order'});
