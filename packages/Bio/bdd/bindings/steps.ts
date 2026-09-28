/* The steps only Bio can define: its readiness and what it reads out of a sequence. Everything
   generic — the top menu, dialogs, columns, viewers — is the library's (`grok-bdd list-steps`). */
import type {Page} from '@playwright/test';
import {Given, Then} from '@datagrok-libraries/bdd';
import {expect, silent} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/** Bio initializes on its first call (RDKit, the monomer libraries, the sequence helper — about
 * eight seconds on a fresh page) and the platform holds every call of the package until then;
 * a service getter is that hold made visible, so a feature's first command is not the one that
 * pays for it. Free once the package is up — provided the helper stays in the page: returned to
 * Node it is the RDKit module's 16 MB heap serialized to base64, ten seconds per call. */
export const bioInitialized = Given('the Bio package is initialized', async (page: Page) => {
  silent(page);
  await page.evaluate(async () => { await grok.functions.call('Bio:getSeqHelper', {}); });
}, {tier: 'api', description: 'awaits Bio\'s init through the platform (a call of a Bio service getter returns once initBio has run); not in the video'});

/** The selection against the monomers Bio splits each sequence into, in any notation: the selected
 * rows are exactly the ones whose sequence has that monomer at that 1-based position. */
export const selectedByMonomer = Then('only the rows with {string} at position {int} of {string} column should be selected',
  async (page: Page, monomer: string, position: number, column: string) => {
    await expect.poll(() => page.evaluate(async ([m, pos, name]) => {
      const t = grok.shell.t;
      const col = t.col(name);
      if (!col)
        return `the table has no "${name}" column`;
      const sh = (await grok.functions.call('Bio:getSeqHelper', {})).getSeqHandler(col);
      let matching = 0;
      let wrong = 0;
      for (let i = 0; i < t.rowCount; i++) {
        const seq = sh.getSplitted(i);
        const has = seq.length >= pos && seq.getOriginal(pos - 1) === m;
        if (has)
          matching++;
        if (has !== t.selection.get(i))
          wrong++;
      }
      return matching === 0 ? `no row has ${m} at position ${pos}` : wrong === 0 ? 'exactly' : `${wrong} rows off`;
    }, [monomer, position, column] as const), {message: `the selection against the rows with ${monomer} at position ${position} of ${column}`})
      .toBe('exactly');
  }, {description: 'every selected row has the monomer there and no other row does; a monomer no row has there fails'});
