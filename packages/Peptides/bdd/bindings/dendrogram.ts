/* The dendrogram the Dendrogram package attaches next to the analysis grid: not a viewer but a
   grid neighbour, which publishes its `dendrogram` status on the grid (`tree leaves`, `tree` area)
   while it is attached and withdraws it when it closes. Its own budget: the tree is computed after
   the analysis is ready, a distance matrix over every peptide. */
import type {Page} from '@playwright/test';
import {Then} from '@datagrok-libraries/bdd';
import {expect, pollMs} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

const treeLeaves = (page: Page): Promise<number | undefined> => page.evaluate(() => {
  const leaves = grok.shell.tv?.grid?.getWidgetStatus?.()?.values?.['tree leaves'];
  return typeof leaves === 'number' ? leaves : undefined;
});

export const dendrogramAttached = Then('the analysis grid should have a dendrogram', async (page: Page) => {
  await expect.poll(async () => (await treeLeaves(page) ?? 0) > 0,
    {message: 'no dendrogram is attached next to the grid', timeout: pollMs(120000)}).toBe(true);
}, {description: 'the grid reports the tree\'s "tree leaves" reading, above zero once the tree is drawn'});

export const dendrogramDetached = Then('the analysis grid should not have a dendrogram', async (page: Page) => {
  await expect.poll(async () => await treeLeaves(page) === undefined,
    {message: 'a dendrogram is still attached next to the grid', timeout: pollMs(10000)}).toBe(true);
}, {description: 'the grid reports no "tree leaves" reading: the tree withdrew its status when it closed'});
