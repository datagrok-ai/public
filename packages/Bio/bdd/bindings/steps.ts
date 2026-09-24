/* The steps only Bio can define: its readiness. Everything generic — the top menu, dialogs,
   columns, viewers — is the library's (`grok-bdd list-steps`). */
import type {Page} from '@playwright/test';
import {Given} from '@datagrok-libraries/bdd';

declare const grok: any;

/** Bio initializes on its first call (RDKit, the monomer libraries, the sequence helper — about
 * eight seconds on a fresh page) and the platform holds every call of the package until then;
 * a service getter is that hold made visible, so a feature's first command is not the one that
 * pays for it. Free once the package is up — provided the helper stays in the page: returned to
 * Node it is the RDKit module's 16 MB heap serialized to base64, ten seconds per call. */
export const bioInitialized = Given('the Bio package is initialized', async (page: Page) => {
  await page.evaluate(async () => { await grok.functions.call('Bio:getSeqHelper', {}); });
}, {tier: 'api', description: 'awaits Bio\'s init through the platform (a call of a Bio service getter returns once initBio has run)'});
