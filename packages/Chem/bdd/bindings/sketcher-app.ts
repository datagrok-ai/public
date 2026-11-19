// Apps | Chem | Sketcher (src/utils/sketcher-app/sketcher-app.ts) and its Export dialog (src/widgets/molecule-export.ts):
// the link in the address bar read and compared as a molecule, the view measured for scrollbars, the window sized, a
// link followed in a fresh page load, and what the Export dialog's Save downloads, read from the file. What each format
// holds in each style, exactly, is Chem's package tests ('Sketcher: app').
import {readFileSync} from 'node:fs';
import {tmpdir} from 'node:os';
import {join} from 'node:path';
import {Page} from '@playwright/test';
import {Given, Then, When, element} from '@datagrok-libraries/bdd';
import {atFeatureEnd, el, expect, pollMs, viewers} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;
declare const DG: any;

const DIALOG = '.d4-dialog:has(.chem-molecule-export-preview)';

element('molecule Export dialog', {selector: DIALOG,
  description: 'the Sketcher app\'s Export dialog (its ribbon\'s Export...): Style, Format, Transparent background, ' +
    'the preview, Save and Copy'});

/** The dialog's input captioned `caption`: its root. */
const inputOf = (page: Page, caption: string) =>
  page.locator(`${DIALOG} .ui-input-root`).filter({has: page.locator('.ui-input-label', {hasText: caption})}).first();

/** RDKit's canonical SMILES of `smiles` (its CX part left out), in the page; '' when RDKit reads it as no molecule. */
function canonical(page: Page, smiles: string): Promise<string> {
  return page.evaluate(async (s) => {
    const rdkit = await grok.functions.call('Chem:getRdKitModule');
    const mol = rdkit.get_mol(s);
    try {
      return mol?.is_valid() ? mol.get_smiles() : '';
    } finally {
      mol?.delete();
    }
  }, smiles);
}

/** The `smiles` parameter of the address bar, decoded; null when the address has none. */
const linkOf = (page: Page) => page.evaluate(() => new URLSearchParams(location.search).get('smiles'));

export const linkIsMolecule = Then('the address\'s link should be the molecule {string}',
  async (page: Page, smiles: string) => {
    const want = await canonical(page, smiles);
    expect(want, `RDKit reads ${smiles}`).not.toBe('');
    let got = '';
    await expect.poll(async () => {
      got = (await linkOf(page)) ?? '(no smiles parameter)';
      return got.startsWith('(') ? got : await canonical(page, got);
    }, {message: `the molecule in the address (${page.url()})`, timeout: pollMs(15000)}).toBe(want);
  }, {tier: 'api', description: 'the address bar\'s smiles parameter, decoded, is that molecule (each read by RDKit, ' +
    'canonically; the CX part aside)'});

export const linkIs = Then('the address\'s link should be {string}', async (page: Page, text: string) => {
  await expect.poll(() => linkOf(page), {message: `the address's smiles parameter (${page.url()})`,
    timeout: pollMs(15000)}).toBe(text);
}, {tier: 'api', description: 'the address bar\'s smiles parameter, decoded, is that text exactly ("" for an empty one)'});

export const addressIsApp = Then('the address should be the Sketcher app\'s', async (page: Page) => {
  await expect.poll(() => page.evaluate(() => location.pathname), {message: 'the address\'s path',
    timeout: pollMs(15000)}).toBe('/apps/Chem/Sketcher');
}, {tier: 'api', description: 'the address bar\'s path is /apps/Chem/Sketcher'});

export const currentIsMolecule = Then('the current object should be the molecule {string}',
  async (page: Page, smiles: string) => {
    const want = await canonical(page, smiles);
    await expect.poll(async () => {
      const o = await page.evaluate(() => {
        const o = grok.shell.o;
        return o instanceof DG.SemanticValue && o.semType === DG.SEMTYPE.MOLECULE ? String(o.value) : null;
      });
      return o === null ? 'not a molecule' : await canonical(page, o);
    }, {message: 'the current object (grok.shell.o), which the context panel renders', timeout: pollMs(15000)})
      .toBe(want);
  }, {tier: 'api', description: 'grok.shell.o is a Molecule semantic value of that molecule (read by RDKit)'});

export const panelShowsPane = Then('the context panel should show the {string} pane', async (page: Page, pane: string) => {
  await expect(page.locator('.grok-prop-panel .d4-accordion-pane-header').filter({hasText: pane}).first(),
    `the ${pane} pane's header in the context panel`).toBeVisible({timeout: pollMs(30000)});
}, {tier: 'ui', description: 'a pane of that name in the context panel, its header shown'});

export const noScrollbar = Then('the Sketcher app should fill its view with Crux, no scrollbar and no field',
  async (page: Page) => {
    await viewers.settle(page, el('crux sketcher widget')).catch(() => undefined);
    const m = await page.evaluate(() => {
      const view = grok.shell.v?.root as HTMLElement | undefined;
      const crux = view?.querySelector('crux-sketch') as HTMLElement | null;
      if (!view || !crux)
        return null;
      const [v, c] = [view.getBoundingClientRect(), crux.getBoundingClientRect()];
      const chain: HTMLElement[] = [];
      for (let e: HTMLElement | null = crux; e && e !== view.parentElement; e = e.parentElement)
        chain.push(e);
      const scrolled = chain.filter((e) => e.scrollHeight > e.clientHeight + 1 || e.scrollWidth > e.clientWidth + 1)
        .map((e) => `${e.tagName.toLowerCase()}.${e.className}: ${e.scrollWidth}x${e.scrollHeight} in ` +
          `${e.clientWidth}x${e.clientHeight}`);
      const field = view.querySelector('.grok-sketcher-input') as HTMLElement | null;
      return {view: [v.width, v.height], crux: [c.width, c.height], scrolled,
        field: field === null ? 'none' : field.getClientRects().length > 0 ? 'shown' : 'hidden',
        page: [document.documentElement.scrollWidth > window.innerWidth, document.documentElement.scrollHeight > window.innerHeight]};
    });
    expect(m, 'the Sketcher app\'s view and its Crux').not.toBeNull();
    expect(m!.scrolled, 'what scrolls between Crux and the view').toEqual([]);
    expect(m!.page, 'the page scrolls neither way').toEqual([false, false]);
    expect(m!.field, 'the host\'s input row (left out by the host: showInput: false)').toBe('none');
    expect(Math.abs(m!.view[0] - m!.crux[0]) <= 2 && Math.abs(m!.view[1] - m!.crux[1]) <= 2,
      `Crux (${m!.crux}) fills the view (${m!.view})`).toBe(true);
  }, {tier: 'ui', description: 'Crux\'s box is the view\'s (2 px aside), nothing from Crux up to the view nor the page ' +
    'scrolls, and the sketcher host\'s input row (its molecule field) is not in the page'});

export const helpPanelOpen = Given('the help panel is open', async (page: Page) => {
  const was = await page.evaluate(() => grok.shell.windows.showHelp === true);
  await expect.poll(() => page.evaluate(() => {
    grok.shell.windows.showHelp = true;
    return (document.querySelector('.grok-help') as HTMLElement | null)?.offsetParent != null;
  }), {message: 'the help panel shown', timeout: pollMs(15000)}).toBe(true);
  atFeatureEnd(page, () => page.evaluate((shown) => {
    grok.shell.windows.showHelp = shown;
  }, was));
}, {tier: 'api', description: 'the shell setting, re-applied until the help panel is shown; at feature end the panel ' +
  'is as the step found it'});

export const windowSize = Given('the browser window is {int} by {int}', async (page: Page, width: number, height: number) => {
  const was = page.viewportSize();
  await page.setViewportSize({width, height});
  if (was)
    atFeatureEnd(page, () => page.setViewportSize(was));
}, {tier: 'api', description: 'the page\'s viewport, as a window of that size; the size it had comes back at feature end'});

export const sessionSketcher = Then('the session\'s sketcher should be {string}', async (page: Page, name: string) => {
  await expect.poll(() => page.evaluate(() => DG.chem.currentSketcherType), {message: 'DG.chem.currentSketcherType',
    timeout: pollMs(15000)}).toBe(name);
}, {tier: 'api', description: 'the sketcher every other molecule editor opens (DG.chem.currentSketcherType)'});

export const followLink = When('user follows the link {string} in a new page load', async (page: Page, path: string) => {
  await page.goto(new URL(path, page.url()).toString());
  // the shell's getters throw until the platform has started (grok exists before its Dart side answers)
  await page.waitForFunction(() => {
    try {
      return (window as any).grok?.shell?.v?.name === 'Sketcher';
    } catch {
      return false;
    }
  }, undefined, {timeout: pollMs(120000), polling: 250});
}, {tier: 'ui', description: 'the page loaded afresh at that address (a pasted link: the platform starts, its routing ' +
  'runs the app with the address\'s parameters); done when the Sketcher view is current'});

export const exportOpens = When('user clicks Export... in the Sketcher app\'s ribbon', async (page: Page) => {
  await page.locator('.d4-ribbon button').filter({hasText: 'Export...'}).first().click();
  await expect(page.locator(DIALOG).first(), 'the molecule Export dialog').toBeVisible({timeout: pollMs(15000)});
}, {tier: 'ui', description: 'the ribbon\'s Export... clicked; done when the dialog is up'});

export const exportChoose = When('user chooses the {string} style and the {string} format in the molecule Export dialog',
  async (page: Page, style: string, format: string) => {
    await inputOf(page, 'Format').locator('select').selectOption(format);
    const select = inputOf(page, 'Style').locator('select');
    if (await select.isEnabled())
      await select.selectOption(style);
    await expect.poll(() => page.locator(`${DIALOG} .chem-molecule-export-preview`).first()
      .getAttribute('data-ready'), {message: 'the preview made', timeout: pollMs(30000)}).toBe('true');
  }, {tier: 'ui', description: 'the format, then the style (off for a format it does not apply to) picked in the ' +
    'dialog\'s choices; done when its preview is made'});

export const styleOff = Then('the molecule Export dialog\'s style should be off', async (page: Page) => {
  await expect(inputOf(page, 'Style').locator('select'), 'the Style choice').toBeDisabled();
}, {tier: 'ui', description: 'the Style choice disabled: the format is the drawing\'s text whatever the style'});

/** The file the dialog's Save last downloaded, for the checks after (shared through the page's global, as the
 * reactions' panel's bindings do). */
const saved: WeakMap<Page, {file: string, name: string}> =
  ((globalThis as any)[Symbol.for('chem.bdd.molecule-export.saved')] ??= new WeakMap());

export const exportSave = When('user clicks Save in the molecule Export dialog', async (page: Page) => {
  const button = page.locator(DIALOG).first().getByRole('button', {name: /^save$/i});
  const [download] = await Promise.all([page.waitForEvent('download', {timeout: pollMs(30000)}), button.click()]);
  const file = join(tmpdir(), `bdd-molecule-export-${process.pid}-${Date.now()}-${download.suggestedFilename()}`);
  await download.saveAs(file);
  saved.set(page, {file, name: download.suggestedFilename()});
  await expect(page.locator(DIALOG), 'Save closes the dialog').toHaveCount(0, {timeout: pollMs(10000)});
}, {tier: 'ui', description: 'the dialog\'s Save clicked; the file it downloads kept for the checks after'});

export const savedFile = Then('the saved file should be {string}, holding {string}', async (page: Page, name: string,
  marker: string) => {
  const got = saved.get(page);
  if (!got)
    throw new Error('nothing saved: "user clicks Save in the molecule Export dialog" first');
  expect(got.name, 'the file\'s name').toBe(name);
  const bytes = readFileSync(got.file);
  const text = marker === 'PNG' ? bytes.subarray(1, 4).toString('latin1') : bytes.toString('utf8');
  expect(text, `${name} holds "${marker}"`).toContain(marker.replace(/\\n/g, '\n'));
}, {tier: 'api', description: 'the downloaded file\'s name, and a text it holds (a PNG: its signature)'});

export const exportCopy = When('user clicks Copy in the molecule Export dialog', async (page: Page) => {
  await page.context().grantPermissions(['clipboard-read', 'clipboard-write']);
  await page.locator(DIALOG).first().getByRole('button', {name: /^copy$/i}).click();
}, {tier: 'ui', description: 'the dialog\'s Copy clicked (the page allowed the clipboard)'});

export const clipboardIs = Then('the clipboard should hold the molecule {string}', async (page: Page, smiles: string) => {
  const want = await canonical(page, smiles);
  await expect.poll(async () => canonical(page, await page.evaluate(() => navigator.clipboard.readText())),
    {message: 'the clipboard\'s text read as a molecule', timeout: pollMs(10000)}).toBe(want);
}, {tier: 'api', description: 'the clipboard\'s text is that molecule (read by RDKit)'});

export const exportShowsText = Then('the molecule Export dialog should show the text {string}',
  async (page: Page, text: string) => {
    await expect(page.locator(`${DIALOG} .chem-molecule-export-text`).first(), 'the text shown').toHaveValue(text);
  }, {tier: 'ui', description: 'the preview\'s text box holds that text exactly'});

export const exportNote = Then('the molecule Export dialog\'s note should say {string}', async (page: Page, text: string) => {
  await expect(page.locator(`${DIALOG} .chem-molecule-export-note`).first(), 'the note').toContainText(text);
}, {tier: 'ui', description: 'the note under the preview contains the text'});

export const openSketcherApp = Given('user opens Chem\'s Sketcher app', async (page: Page) => {
  await page.evaluate(() => {
    DG.Dialog.getOpenDialogs().forEach((d: any) => d.close());
    grok.shell.route('/apps/Chem/Sketcher');
  });
  await expect.poll(() => page.evaluate(() => String(grok.shell.v?.name ?? '')),
    {message: 'the current view', timeout: pollMs(30000)}).toBe('Sketcher');
}, {tier: 'api', description: 'the app opened from its address (grok.shell.route, as a link inside the running ' +
  'platform), on its example; done when its view is current'});
