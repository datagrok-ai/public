import type {Page} from '@playwright/test';
import {Then} from '@datagrok-libraries/bdd';
import {expect, pollMs} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

export const dendrogramAttached = Then('the analysis grid should have a dendrogram', async (page: Page) => {
  await expect.poll(() => page.evaluate(() => {
    const nb = grok.shell.tv?.grid?.temp?.['__dendrogram_neighbor_temp__'];
    return !!nb && !!nb.root?.isConnected && nb.root.getBoundingClientRect().width > 0;
  }), {message: 'no dendrogram is attached next to the grid', timeout: pollMs(120000)}).toBe(true);
}, {description: 'the Dendrogram package attaches its tree as a grid neighbour, not a viewer; reads the grid\'s neighbour and its visible root'});

export const dendrogramDetached = Then('the analysis grid should not have a dendrogram', async (page: Page) => {
  await expect.poll(() => page.evaluate(() => {
    const nb = grok.shell.tv?.grid?.temp?.['__dendrogram_neighbor_temp__'];
    return !nb || !nb.root?.isConnected;
  }), {message: 'a dendrogram is still attached next to the grid', timeout: pollMs(10000)}).toBe(true);
}, {description: 'no grid neighbour dendrogram, or its root detached'});
