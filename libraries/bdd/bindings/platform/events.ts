/* Custom platform events (`grok.events.fireCustomEvent`): the word a package gives about work it
   finished off-screen — Bio's `bio-monomer-lib-loaded` after the monomer libraries reload. A
   scenario listens by id, acts, and claims the event fired (or did not). */
import {Page} from '@playwright/test';
import {Given, Then} from '../../src/registry.js';
import {expectCustomEvent, listenCustomEvent} from '../../src/runtime/events.js';

export const listenCustom = Given('user listens for {string} custom event', (page: Page, id: string) => listenCustomEvent(page, id),
  {tier: 'api', description: 'grok.events.onCustomEvent(id), counted until the page resets; a "should have fired" read zeroes the count'});

export const customFired = Then('the {string} custom event should have fired', async (page: Page, id: string) => { await expectCustomEvent(page, id); },
  {description: 'at least once since "listens for" or the previous read (up to 30 s); reading zeroes the count'});
