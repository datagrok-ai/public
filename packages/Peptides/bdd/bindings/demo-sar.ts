/* What the Peptides demo clips (features/demos) need and the library has not got: set-up a filmed demo leaves
   off camera. */
import type {Page} from '@playwright/test';
import {Given} from '@datagrok-libraries/bdd';
import {expect, silent, viewers} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;
declare const DG: any;

/** The SAR analysis docks a "Selection" grid (the selected rows) under the Logo Summary Table; with nothing selected
 * it shows only its headers, which a clip that ends before any selection closes off camera. The viewer is the
 * current view's one with that title (or `viewer-<title>` as its element's name), closed as its × closes it. */
export const viewerClosed = Given('the {string} viewer of the analysis is closed', async (page: Page, title: string) => {
  silent(page);
  const why = await page.evaluate((t) => {
    const tv = grok.shell.tv;
    const viewer = Array.from(tv?.viewers ?? []).find((v: any) => {
      let own = '';
      try { own = v.props?.title ?? ''; } catch (_) { own = ''; }
      return own === t || v.root?.getAttribute('name') === `viewer-${t}`;
    }) as any;
    if (!viewer)
      return `the current view has no "${t}" viewer`;
    viewer.close();
    return '';
  }, title);
  if (why)
    throw new Error(why);
  await expect.poll(() => page.evaluate((t) => document.querySelectorAll(`[name="viewer-${t}"]`).length, title),
    {message: `the "${title}" viewer closed`}).toBe(0);
  await viewers.settleAll(page);
}, {tier: 'api', description: 'set-up a filmed demo leaves out: the analysis view\'s viewer of that title closed'});
