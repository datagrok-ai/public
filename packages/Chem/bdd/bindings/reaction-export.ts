// The Export panel (src/widgets/reaction-export.ts): its dialog, what it offers, the
// choices a person makes in it, its preview at print size, and the PNG its Save writes, read from the downloaded file.
// What it draws for each choice (the sizes, the rows, the legend, the colours) is Chem's package tests
// ('Reactions: export').
import {readFileSync} from 'node:fs';
import {tmpdir} from 'node:os';
import {join} from 'node:path';
import {Page} from '@playwright/test';
import {Given, Then, When, element} from '@datagrok-libraries/bdd';
import {expect, pollMs} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;
declare const DG: any;

const PANEL = '.d4-dialog:has(.chem-reaction-export-preview)';

element('Export panel', {selector: PANEL,
  description: 'the Export panel: the dialog a reaction\'s Export opens (the Reaction pane\'s icon, a reaction cell\'s ' +
    'menu, the Reactions app\'s ribbon), with its preview'});

/** What the panel offers: its presets, widths and formats, and Copy as picture. */
const OFFERED = ['ACS Document 1996', 'Screen', 'single column (8.3 cm)', 'double column (17.8 cm)', 'free', 'SVG',
  'PNG 300 dpi', 'PNG 600 dpi', 'CDXML', 'Transparent background', 'Copy as picture'];

/** The panel's input captioned `caption`: its root. */
const inputOf = (page: Page, caption: string) =>
  page.locator(`${PANEL} .ui-input-root`).filter({has: page.locator('.ui-input-label', {hasText: caption})}).first();

export const openAppOnDefault = Given('user opens Chem\'s Reactions app on its default example', async (page: Page) => {
  await page.evaluate(async () => {
    // a dialog an earlier scenario left open on the shared page (the full-screen Reaction view fills the window: it
    // took the next scenario's clicks, in reaction-centre-legend.feature) is closed first
    DG.Dialog.getOpenDialogs().forEach((d) => d.close());
    window.localStorage.removeItem('Chem.ReactionsApp.document'); // reaction-document.ts DOCUMENT_KEY
    const view = await grok.functions.call('Chem:reactionsApp');
    grok.shell.addView(view);
  });
  await expect.poll(() => page.evaluate(() => String(grok.shell.v?.name ?? '')),
    {message: 'the current view', timeout: pollMs(30000)}).toBe('Reaction Sketcher');
}, {tier: 'api', description: 'the Reactions app run with nothing kept in the browser, so it opens on its default ' +
  'example (a route of five mapped steps)'});

export const exportFromRibbon = When('user clicks Export... in the Reactions app\'s ribbon', async (page: Page) => {
  await page.locator('.d4-ribbon button').filter({hasText: 'Export...'}).first().click();
}, {tier: 'ui', description: 'the ribbon\'s one Export button (the one-click Export SVG and Export PNG are gone)'});

export const panelOffers = Then('the Export panel should offer its presets, widths and formats, and no PDF',
  async (page: Page) => {
    const panel = page.locator(PANEL).first();
    await expect(panel, 'the Export panel').toBeVisible({timeout: pollMs(15000)});
    const text = await panel.textContent() ?? '';
    for (const offered of OFFERED)
      expect(text, `the panel offers "${offered}"`).toContain(offered);
    expect(text, 'no PDF (decision 17A)').not.toMatch(/PDF/i);
  }, {tier: 'ui', description: `the panel's text names ${OFFERED.join(', ')}, and never PDF`});

export const coloursPerPreset = Then('the Export panel\'s centre colours should be off for {string} and on for {string}',
  async (page: Page, off: string, on: string) => {
    const style = inputOf(page, 'Style').locator('select');
    const colours = inputOf(page, 'Reaction centre colours').locator('input[type="checkbox"]');
    await style.selectOption(on);
    await expect(colours, `${on}: the centre colours on`).toBeChecked();
    await style.selectOption(off);
    await expect(colours, `${off}: the centre colours off (P9)`).not.toBeChecked();
  }, {tier: 'ui', description: 'the colours switch follows the preset chosen: on for the one, off for the other'});

export const choose = When('user chooses {string}, {string} and {string} in the Export panel',
  async (page: Page, preset: string, width: string, format: string) => {
    await inputOf(page, 'Style').locator('select').selectOption(preset);
    await inputOf(page, 'Width').locator('select').selectOption(width);
    await inputOf(page, 'Format').locator('select').selectOption(format);
  }, {tier: 'ui', description: 'the preset, the width and the format picked in the panel\'s three choices'});

export const previewWidth = Then('the Export panel\'s preview should be {float} cm wide in the {string} look',
  async (page: Page, cm: number, look: string) => {
    const preview = page.locator(`${PANEL} .chem-reaction-export-preview`).first();
    await expect(preview.locator('canvas'), 'the preview is drawn').toBeVisible({timeout: pollMs(15000)});
    await expect(preview, 'the width on paper').toHaveAttribute('data-cm', cm.toFixed(1));
    await expect(preview, 'the look').toHaveAttribute('data-look', look);
    // on the screen: print size (96 CSS px an inch) times the zoom that fits it to the preview's box (the box fills the
    // dialog, the picture drawn anew to fit it)
    const zoom = Number(await preview.getAttribute('data-zoom'));
    expect(zoom > 0, `the preview's zoom: ${zoom}`).toBe(true);
    const shown = await preview.locator('canvas').evaluate((c) => (c as HTMLCanvasElement).getBoundingClientRect().width);
    expect(Math.abs(shown - cm / 2.54 * 96 * zoom), `the preview ${shown.toFixed(1)} CSS px wide at zoom ${zoom}`)
      .toBeLessThan(2);
  }, {tier: 'ui', description: 'the preview\'s width on paper and its look, and its canvas as wide on the screen as the ' +
    'width at 96 CSS px an inch times the zoom that fits it to the preview\'s box'});

export const previewAtPrintSize = Then('the Export panel\'s preview should be at print size', async (page: Page) => {
  const preview = page.locator(`${PANEL} .chem-reaction-export-preview`).first();
  await expect(preview, 'the bonds at print size (ACS 1996: 0.508 cm)').toHaveAttribute('data-percent', '100',
    {timeout: pollMs(15000)});
  await expect(page.locator(`${PANEL} .chem-reaction-export-note`).first(), 'the note').toContainText('at print size');
}, {tier: 'ui', description: 'the preview\'s bonds at print size (100%), as its note says'});

export const previewSmaller = Then('the Export panel\'s note should say how large the preview is against print size',
  async (page: Page) => {
    const preview = page.locator(`${PANEL} .chem-reaction-export-preview`).first();
    await expect(preview.locator('canvas'), 'the preview is drawn').toBeVisible({timeout: pollMs(15000)});
    const percent = Number(await preview.getAttribute('data-percent'));
    const note = page.locator(`${PANEL} .chem-reaction-export-note`).first();
    await expect(note, 'the note').toContainText(percent < 100 ? `drawn at ${percent}% of print size` : 'at print size');
  }, {tier: 'ui', description: 'the note under the preview: "at print size", or how much smaller it is drawn (a part of ' +
    'the scheme wider than the width)'});

/** The file each page's last Save in the Export panel wrote, kept on the global object under a symbol that
 * reaction-centre-legend.ts reads too (a binding importing another loads it a second time, and grok-bdd compile then
 * finds steps that are no export of the module it read). */
const saved: WeakMap<Page, string> = ((globalThis as any)[Symbol.for('chem.bdd.export.saved')] ??=
  new WeakMap<Page, string>());

export const save = When('user clicks Save in the Export panel', async (page: Page) => {
  const button = page.locator(PANEL).first().getByRole('button', {name: /^save$/i});
  const [download] = await Promise.all([page.waitForEvent('download', {timeout: pollMs(30000)}), button.click()]);
  const file = join(tmpdir(), `bdd-export-${process.pid}-${Date.now()}-${download.suggestedFilename()}`);
  await download.saveAs(file);
  saved.set(page, file);
  await expect(page.locator(PANEL), 'Save closes the panel').toHaveCount(0, {timeout: pollMs(10000)});
}, {tier: 'ui', description: 'the panel\'s Save clicked; the file it downloads kept for the checks after'});

export const savedPng = Then('the saved PNG should be {int} px wide at {int} dpi, with no colour in it',
  async (page: Page, width: number, dpi: number) => {
    const file = saved.get(page);
    if (!file)
      throw new Error('nothing saved: "user clicks Save in the Export panel" first');
    const bytes = readFileSync(file);
    expect(bytes.subarray(1, 4).toString('latin1'), 'a PNG').toBe('PNG');
    expect(bytes.readUInt32BE(16), 'its width (IHDR)').toBe(width);
    // the pHYs chunk right after IHDR (reaction-document.ts withResolution): pixels per metre
    expect(bytes.subarray(37, 41).toString('latin1'), 'a pHYs chunk after IHDR').toBe('pHYs');
    expect(Math.round(bytes.readUInt32BE(41) * 0.0254), 'its dpi (pHYs)').toBe(dpi);
    // no pixel with colour: no centre colour, so no legend's swatches either (ACS 1996 is black)
    const coloured = await page.evaluate(async (base64) => {
      const blob = await (await fetch(`data:image/png;base64,${base64}`)).blob();
      const bitmap = await createImageBitmap(blob);
      const canvas = new OffscreenCanvas(bitmap.width, bitmap.height);
      const g = canvas.getContext('2d')!;
      g.drawImage(bitmap, 0, 0);
      const d = g.getImageData(0, 0, canvas.width, canvas.height).data;
      let n = 0;
      for (let i = 0; i < d.length; i += 4) {
        if (d[i + 3] > 128 && Math.max(Math.abs(d[i] - d[i + 1]), Math.abs(d[i + 1] - d[i + 2]),
          Math.abs(d[i] - d[i + 2])) > 80)
          n++;
      }
      return n;
    }, bytes.toString('base64'));
    expect(coloured, 'coloured pixels in the PNG').toBe(0);
  }, {tier: 'api', description: 'the downloaded PNG\'s width (IHDR) and dpi (pHYs), and not one coloured pixel (the ' +
    'centre\'s green, red, blue or a legend swatch)'});
