// The reaction centre's legend (src/widgets/reaction-centre-legend.ts): the legends under the pictures of the Reaction
// Sketcher, its full-screen view and the Reaction pane, read from the page (the line's words, each swatch's colour, the
// picture it is under); the app opened on a reaction of the scenario's; a table of the reaction the app draws; the
// centre's colours in a PNG the Export panel saved (reaction-export.ts). Which bonds each picture colours is Chem's
// package tests ('Reactions: centre').
import {readFileSync} from 'node:fs';
import {Page} from '@playwright/test';
import {Then, When} from '@datagrok-libraries/bdd';
import {el, expect, locate, pollMs} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;
declare const DG: any;

/** The file a page's last Save in the Export panel wrote, as reaction-export.ts keeps it (on the global object under
 * its symbol: importing that binding would load it twice). */
const savedFile = (page: Page): string | undefined =>
  ((globalThis as any)[Symbol.for('chem.bdd.export.saved')] as WeakMap<Page, string> | undefined)?.get(page);

/** The swatches' colours (reaction-centre-legend.ts CENTRE_COLOURS) as the browser gives them, by their words. */
const COLOUR_OF: {[words: string]: string} = {'bond made': 'rgb(0, 158, 76)', 'bond broken': 'rgb(219, 40, 40)',
  'bond order changed': 'rgb(30, 107, 234)'};

/** A legend as the page shows it: shown or hidden, its first words, each colour's words and swatch, and whether it is
 * right under the first picture (the route's: the first canvas, just before it, above it, centred on it). */
interface Legend {shown: boolean, title: string, items: string[], underPicture: boolean}

/** The legends in `root`, in order (run in the page). */
const readLegends = (root: Element): Legend[] =>
  Array.from(root.querySelectorAll('.chem-reaction-centre-legend')).map((l) => {
    const legend = l as HTMLElement;
    const before = legend.previousElementSibling;
    const items = Array.from(legend.children).slice(1).map((item) => {
      const [swatch, words] = Array.from(item.children) as HTMLElement[];
      return `${words?.textContent ?? ''} ${swatch ? getComputedStyle(swatch).backgroundColor : ''}`;
    });
    const [p, r] = [before?.getBoundingClientRect(), legend.getBoundingClientRect()];
    return {shown: legend.style.display !== 'none', title: legend.children[0]?.textContent ?? '', items,
      underPicture: before?.tagName === 'CANVAS' && before === root.querySelector('canvas') && !!p &&
        p.bottom <= r.top + 1 && Math.abs((p.left + p.right) / 2 - (r.left + r.right) / 2) <= 2};
  });

/** The legends the current view shows (the Reactions app). */
async function appLegends(page: Page): Promise<Legend[]> {
  const root = await page.evaluateHandle(() => grok.shell.v?.root as Element);
  return root.evaluate(readLegends);
}

/** Polls `read` until one legend, right under the route's picture and centred on it, names `words` (comma
 * separated), each in its colour, and no other legend is shown (a route's steps share the route's: a second one sat at
 * the end, under the steps). */
async function expectRouteLegend(read: () => Promise<Legend[]>, where: string, words: string): Promise<void> {
  const items = words.split(',').map((w) => w.trim()).map((w) => `${w} ${COLOUR_OF[w] ?? '(no such colour)'}`);
  const one: Legend = {shown: true, title: 'Reaction centre:', items, underPicture: true};
  await expect.poll(async () => JSON.stringify((await read()).filter((l) => l.shown)),
    {message: `${where}: one legend, right under the route's picture`, timeout: pollMs(20000)})
    .toBe(JSON.stringify([one]));
}

export const appRouteLegend = Then('the Reactions app should name the centre\'s colours once, right under its ' +
  'route: {string}', async (page: Page, words: string) => {
  await expectRouteLegend(() => appLegends(page), 'the Reactions app', words);
}, {tier: 'ui', description: 'the app\'s preview: one line right under the route\'s picture, centred on it, ' +
  '"Reaction centre:" and each colour the route and its steps have, its swatch in that colour, and no other'});

export const fullScreenRouteLegend = Then('the full-screen Reaction dialog should name the centre\'s colours once, ' +
  'right under its route: {string}', async (page: Page, words: string) => {
  const dialog = (await locate(page, el('"Reaction" dialog'))).first();
  await expect(dialog, 'the full-screen view').toBeVisible({timeout: pollMs(15000)});
  await expectRouteLegend(() => dialog.evaluate(readLegends), 'the full-screen view', words);
}, {tier: 'ui', description: 'the full-screen view\'s legend, as the app\'s: once, right under the route\'s picture'});

export const appFullScreen = When('user clicks the Full screen icon of the Reactions app\'s preview', async (page: Page) => {
  // the preview's own icon (the Reaction panel's, widgets/reaction-panel.ts), found in the app's view and marked
  await page.evaluate(() => grok.shell.v?.root?.querySelector('.fa-expand-alt')?.setAttribute('data-bdd-full-screen', ''));
  const icon = page.locator('[data-bdd-full-screen]').first();
  await expect(icon, 'the preview\'s Full screen icon').toBeVisible({timeout: pollMs(15000)});
  await icon.click();
  await expect((await locate(page, el('"Reaction" dialog'))).first(), 'the full-screen view')
    .toBeVisible({timeout: pollMs(15000)});
}, {tier: 'ui', description: 'the Full screen icon over the app\'s preview clicked; the full-screen Reaction dialog ' +
  'opens'});

export const previewLegend = Then('the Export panel\'s preview should carry the centre\'s legend', async (page: Page) => {
  const preview = page.locator('.d4-dialog .chem-reaction-export-preview').first();
  await expect(preview.locator('canvas'), 'the preview is drawn').toBeVisible({timeout: pollMs(15000)});
  await expect(preview, 'the legend under the scheme').toHaveAttribute('data-legend', 'true');
}, {tier: 'ui', description: 'the Export panel\'s preview drawn with the legend under its scheme (the colours on)'});

export const pngColours = Then('the saved PNG should have the centre\'s green and red in it', async (page: Page) => {
  const file = savedFile(page);
  if (!file)
    throw new Error('nothing saved: "user clicks Save in the Export panel" first');
  const counts = await page.evaluate(async (base64) => {
    const blob = await (await fetch(`data:image/png;base64,${base64}`)).blob();
    const bitmap = await createImageBitmap(blob);
    const canvas = new OffscreenCanvas(bitmap.width, bitmap.height);
    const g = canvas.getContext('2d')!;
    g.drawImage(bitmap, 0, 0);
    const d = g.getImageData(0, 0, canvas.width, canvas.height).data;
    // the centre's green and red (reaction-centre-legend.ts CENTRE_COLOURS), each channel within 24: a bond's core
    // pixels and the legend's swatches (an oxygen's red, #FF0000, is further)
    const near = (i: number, [r, gr, b]: number[]) =>
      d[i + 3] > 200 && Math.abs(d[i] - r) <= 24 && Math.abs(d[i + 1] - gr) <= 24 && Math.abs(d[i + 2] - b) <= 24;
    let green = 0;
    let red = 0;
    for (let i = 0; i < d.length; i += 4) {
      if (near(i, [0, 158, 76]))
        green++;
      else if (near(i, [219, 40, 40]))
        red++;
    }
    return {green, red};
  }, readFileSync(file).toString('base64'));
  expect(counts.green, 'green pixels (bonds made, their swatch)').toBeGreaterThan(50);
  expect(counts.red, 'red pixels (bonds broken, their swatch)').toBeGreaterThan(50);
}, {tier: 'api', description: 'the downloaded PNG has pixels of the centre\'s green and of its red: the bonds made and ' +
  'broken, and the legend\'s swatches'});
