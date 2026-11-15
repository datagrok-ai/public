/* Steps of the Oligo Toolkit demo clip (features/demos/21): the app opened off camera, the Translator's structure
   box hidden off camera, and the structure picture of its Structure tool opened full size. */
import type {Page} from '@playwright/test';
import {Given, When} from '@datagrok-libraries/bdd';
import {expect, pollMs, silent} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;
declare const DG: any;

/** An app whose view is named after the tool it opens on, not after the app (Oligo Toolkit opens as the
 * "Oligo Translator" view): the library's `user opens the {string} app` waits for a view of the app's name. */
export const openAppOnView = Given('user opens the {string} app on its {string} view', async (page: Page, app: string, view: string) => {
  await page.evaluate(async (name) => {
    // an admin app (Usage Analysis) carries its role in its meta, not a tag
    const apps = [...DG.Func.find({tags: ['app']}), ...DG.Func.find({meta: {role: 'adminApp'}})];
    const fn = apps.find((f: any) => f.friendlyName === name || f.name === name);
    if (!fn)
      throw new Error(`no "${name}" app is registered`);
    const v = await fn.apply({});
    if (v?.root && !Array.from(grok.shell.views).some((x: any) => x.dart === v.dart))
      grok.shell.addView(v);
  }, app);
  await expect.poll(() => page.evaluate(() => String(grok.shell.v?.name ?? '')), {message: 'the current view'}).toBe(view);
}, {tier: 'api', description: 'runs the app function by name and shows the view it returns; done when the view of that name is current'});

/** The Translator draws the strand it reads in a box under its tables (`.mol-host`, through Chem's canvasMol), which a
 * Chem rendering regression since Nov 12 draws as two corner-to-corner lines and a tiny nucleoside: a clip about the
 * notations hides the box (its block) off camera, before the strand is typed, and leaves the rest of the view as it is. */
export const translatorStructureHidden = Given('the structure box of the {string} view is hidden', async (page: Page, view: string) => {
  silent(page);
  await expect.poll(() => page.evaluate(() => String(grok.shell.v?.name ?? '')), {message: 'the current view'}).toBe(view);
  const hidden = await page.evaluate(() => {
    const host = Array.from(document.querySelectorAll('.mol-host'))
      .find((e) => (e as HTMLElement).offsetParent !== null) as HTMLElement | undefined;
    if (!host)
      return 'no structure box is shown';
    const block = (host.closest('.ui-block') ?? host) as HTMLElement;
    block.style.display = 'none';
    return '';
  });
  if (hidden)
    throw new Error(`the ${view} view: ${hidden}`);
  await expect(page.locator('.mol-host').filter({visible: true}), 'the structure box').toHaveCount(0);
}, {tier: 'api', description: 'set-up a filmed demo leaves out: the Translator\'s structure box (its block) not displayed, the rest of the view untouched'});

/** The Structure tool draws the linked strands small and opens them full size on a click (`Click to zoom`). It
 * redraws the picture 300 ms after the last keystroke, on a new canvas, and a canvas made for an empty input
 * zooms into nothing: the click waits for a picture that has ink and has stayed on the page past that delay. */
export const zoomStructure = When('user clicks on the structure picture of the {string} view', async (page: Page, view: string) => {
  // the picture's canvas shows the zoom-in cursor (SequenceTranslator's MoleculeImage)
  const picture = page.locator('canvas[style*="zoom-in"]').filter({visible: true}).first();
  await expect.poll(() => page.evaluate(() => String(grok.shell.v?.name ?? '')), {message: 'the current view'}).toBe(view);
  await expect(picture, `the structure picture of the ${view} view`).toBeVisible();
  await expect.poll(() => picture.evaluate((c: HTMLCanvasElement) => {
    const since = Number(c.dataset.bddSince ?? (c.dataset.bddSince = String(Date.now())));
    const d = c.getContext('2d')?.getImageData(0, 0, c.width, c.height).data;
    let ink = 0;
    for (let i = 0; d && i < d.length; i += 16)
      if (d[i + 3] > 128 && (d[i] < 200 || d[i + 1] < 200 || d[i + 2] < 200))
        ink++;
    return ink < 200 ? 'not drawn yet' : Date.now() - since < 700 ? 'just redrawn' : 'drawn';
  }), {message: `the structure picture of the ${view} view`, timeout: pollMs(30000)}).toBe('drawn');
  await picture.click();
  await expect(page.locator('[name="dialog-Molecule"] canvas').filter({visible: true}), 'the full-size structure')
    .toHaveCount(1, {timeout: pollMs(30000)});
}, {tier: 'ui', description: 'a click on the drawing of the linked strands; done when the Molecule dialog shows it full size'});
