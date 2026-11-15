/* What the Peptides demo clips (features/demos) need and the library has not got: set-up a filmed demo leaves
   off camera. */
import type {Page} from '@playwright/test';
import {Given} from '@datagrok-libraries/bdd';
import {expect, silent, viewers} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;
declare const DG: any;

/** A numbered column starts in the N-terminal flank: its header and cells show a window of positions from the
 * column's `.positionShift` tag (Bio's sequence-scrolling widget), named by its position names. */
export const headerStartsAt = Given('the header of the {string} column starts at position {string}',
  async (page: Page, column: string, position: string) => {
    silent(page);
    await page.evaluate(([c, p]) => {
      const col = grok.shell.tv?.dataFrame.col(c);
      if (!col)
        throw new Error(`the table has no "${c}" column`);
      const tag: string = col.getTag('.positionLabels') || col.getTag('.positionNames') || '';
      const names = tag ? tag.split(',').map((n: string) => n.trim()) : [];
      const index = names.length ? names.indexOf(p) : Number(p) - 1;
      if (!(index >= 0))
        throw new Error(`"${c}" has no position "${p}"; it starts: ${names.slice(0, 30).join(', ')}`);
      col.setTag('.positionShift', String(index));
      grok.shell.tv.grid.invalidate();
    }, [column, position] as [string, string]);
    await viewers.settleAll(page);
  }, {tier: 'api', description: 'the column\'s window of positions starts at the named one'});

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
