// The Reaction Sketcher's Text beside its drawing (src/utils/reactions/reactions-app.ts; one molecule drawn used to
// turn the Text into a MOL file and the preview into an error dump): the app opened on
// an empty reaction, the reaction arrow tool, a drag on the canvas, the Text read, the preview's line for a lone
// molecule and the words the view shows. What the Text holds for each drawing is Chem's package tests
// ('Reactions: app'). And the app's help (its dialogs said "Help unavailable", and the preview's buttons had no
// tooltips): the app on the default example, the preview's full-screen button hovered and
// clicked, a dialog's help shown.
import {Page} from '@playwright/test';
import {Given, Then, When, element} from '@datagrok-libraries/bdd';
import {el, expect, pollMs, viewers} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;
declare const DG: any;

const CRUX = '[data-u2-name="cruxSketch"] crux-sketch';

element('crux reaction arrow tool', {selector: `${CRUX} [data-testid="toolbar.reaction.arrows"]`,
  description: 'the reaction arrow tool, the arrow-styles palette\'s button in a reaction document\'s left bar (the ' +
    'plain arrow until another style is chosen): a drag on the canvas draws an arrow'});

/** The app's Text box's value (in the page). */
const appText = (page: Page) => page.evaluate(() =>
  (document.querySelector('.chem-reactions-app-text') as HTMLTextAreaElement | null)?.value ?? '');

export const openAppEmpty = Given('user opens Chem\'s Reactions app on an empty reaction', async (page: Page) => {
  await page.evaluate(async () => {
    window.localStorage.removeItem('Chem.ReactionsApp.document'); // reaction-document.ts DOCUMENT_KEY
    // one Reactions app on the page: the steps read the first one they find (a journey's scenarios share the page)
    for (const v of [...grok.shell.views].filter((v: any) => v.name === 'Reaction Sketcher'))
      v.close();
    const view = await grok.functions.call('Chem:reactionsApp');
    grok.shell.addView(view);
  });
  await expect.poll(() => page.evaluate(() => String(grok.shell.v?.name ?? '')),
    {message: 'the current view', timeout: pollMs(30000)}).toBe('Reaction Sketcher');
  // the app's sketcher host is the cruxSketch widget itself (reaction-cell-editor.ts ReactionSketchArea)
  const sketcher = '.chem-reactions-app-sketcher crux-sketch';
  await expect(page.locator(`${sketcher} [data-testid="toolbar.clear"]`).first(), 'the app\'s sketcher, started')
    .toBeVisible({timeout: pollMs(30000)});
  // emptied through the Text (crux-sketch 16868b81's Clear canvas left the reaction document as it was; the drawing
  // follows the Text when it is left)
  const text = page.locator('.chem-reactions-app-text');
  await text.fill('');
  await text.press('Tab');
  await expect.poll(() => appText(page), {message: 'the Text emptied'}).toBe('');
  await expect.poll(() => page.evaluate((sel) => (document.querySelector(sel) as any)?.isEmpty ?? false, sketcher),
    {message: 'the drawing emptied', timeout: pollMs(10000)}).toBe(true);
}, {tier: 'api', description: 'the Reactions app run with nothing kept in the browser, its Text emptied: the drawing ' +
  'and the Text empty'});

/** The app's sketcher (the cruxSketch widget is its host, reaction-cell-editor.ts ReactionSketchArea). */
const APP_SKETCHER = '.chem-reactions-app-sketcher crux-sketch';

export const clickAppClear = When('user clicks on the Reactions app\'s Clear canvas', async (page: Page) => {
  // as a person clicks: as soon as the button shows. Until the sketcher has read the reaction its box is under a cover
  // (reaction-cell-editor.ts SKETCH_LOADING_CLASS), which the click waits out as it waits out any box over its target
  await page.locator(`${APP_SKETCHER} [data-testid="toolbar.clear"]`).first().click({timeout: pollMs(30000)});
}, {tier: 'ui', description: 'a click on Clear canvas in the top bar of the app\'s sketcher, as soon as it shows and ' +
  'nothing covers it'});

export const appDrawingEmpty = Then('the Reactions app\'s drawing should be empty', async (page: Page) => {
  const empty = () => page.evaluate((sel) => {
    const el = document.querySelector(sel) as any;
    return !!el && el.isEmpty === true && el.shadowRoot.querySelectorAll('g[data-arrow]').length === 0;
  }, APP_SKETCHER);
  await expect.poll(empty, {message: 'the app\'s sketcher emptied', timeout: pollMs(10000)}).toBe(true);
  // and it stays so: nothing draws the reaction again (the Text's write once came back over the cleared canvas)
  await page.waitForTimeout(1500);
  expect(await empty(), 'the app\'s sketcher, 1.5 s later').toBe(true);
}, {tier: 'api', description: 'the app\'s sketcher holds nothing, no molecule and no arrow, and still holds nothing ' +
  '1.5 s later'});

export const previewAsks = Then('the Reactions app\'s preview should ask for a reaction', async (page: Page) => {
  const drawing = page.locator('.chem-reactions-app-drawing');
  await expect(drawing, 'the right pane').toContainText('Type or paste a reaction', {timeout: pollMs(10000)});
  await expect(drawing.locator('canvas'), 'no picture in the right pane').toHaveCount(0);
}, {tier: 'ui', description: 'the right pane draws no reaction and asks for one to be typed, pasted or chosen ' +
  '(reactions-app.ts draw())'});

export const dragOnCanvas = When('user drags on Crux canvas from {int}% to {int}% across, {int}% down',
  async (page: Page, from: number, to: number, down: number) => {
    const b = await viewers.hitArea(page, el('crux sketcher widget'), 'canvas', true);
    const y = b.y + b.height * down / 100;
    await page.mouse.move(b.x + b.width * from / 100, y);
    await page.mouse.down();
    await page.mouse.move(b.x + b.width * (from + to) / 200, y, {steps: 4});
    await page.mouse.move(b.x + b.width * to / 100, y, {steps: 4});
    await page.mouse.up();
  }, {tier: 'ui', description: 'a drag across Crux\'s canvas, from one point to another on a line, by their places ' +
    'across and down the canvas'});

export const textIs = Then('the Reactions app\'s Text should be {string}', async (page: Page, text: string) => {
  await expect.poll(() => appText(page), {message: 'the Reactions app\'s Text', timeout: pollMs(10000)}).toBe(text);
}, {tier: 'api', description: 'the text box under the app\'s sketcher holds exactly this text'});

export const asksForArrow = Then('the Reactions app\'s preview should say that a reaction needs an arrow',
  async (page: Page) => {
    const line = page.locator('.chem-reactions-app-drawing .chem-reactions-app-no-arrow');
    await expect(line, 'the preview\'s line for a lone molecule').toBeVisible({timeout: pollMs(10000)});
    await expect(line).toHaveText(/^A reaction needs an arrow/);
    await expect(page.locator('.chem-reactions-app-status')).toHaveText(/^A reaction needs an arrow/);
  }, {tier: 'ui', description: 'the right pane shows one short line asking for an arrow (reaction-document.ts NO_ARROW), ' +
    'and the status under the editing area says the same'});

export const viewWords = Then('no text of the Reactions app should contain {string}, {string}, {string} or {string}',
  async (page: Page, a: string, b: string, c: string, d: string) => {
    const shown = await page.evaluate(() => {
      const box = document.querySelector('.chem-reactions-app-text') as HTMLTextAreaElement | null;
      return `${(grok.shell.v?.root as HTMLElement | undefined)?.innerText ?? ''}\n${box?.value ?? ''}`;
    });
    for (const word of [a, b, c, d])
      expect(shown.includes(word), `"${word}" on the Reactions app's view`).toBe(false);
  }, {tier: 'api', description: 'the app\'s view, its Text included, shows none of these words (P10: no engine named, ' +
    'no file text)'});

/** The preview's full-screen button (widgets/reaction-panel.ts: `ui.iconFA('expand-alt', ..., 'Full screen')`). */
const FULL_SCREEN = '.chem-reactions-app-drawing i.fa-expand-alt';

export const openAppExample = Given('user opens Chem\'s Reactions app on the default example', async (page: Page) => {
  await page.evaluate(async () => {
    window.localStorage.removeItem('Chem.ReactionsApp.document'); // reaction-document.ts DOCUMENT_KEY
    // one Reactions app on the page: the steps read the first one they find (a journey's scenarios share the page)
    for (const v of [...grok.shell.views].filter((v: any) => v.name === 'Reaction Sketcher'))
      v.close();
    const view = await grok.functions.call('Chem:reactionsApp');
    grok.shell.addView(view);
  });
  await expect.poll(() => page.evaluate(() => String(grok.shell.v?.name ?? '')),
    {message: 'the current view', timeout: pollMs(30000)}).toBe('Reaction Sketcher');
  await expect(page.locator(FULL_SCREEN).first(), 'the preview drawn, with its full-screen button')
    .toBeVisible({timeout: pollMs(30000)});
}, {tier: 'api', description: 'the Reactions app run with nothing kept in the browser: the default example, drawn on ' +
  'the right with the Reaction panel\'s buttons'});

export const hoverFullScreen = When('user hovers over the Reactions app\'s full-screen button', async (page: Page) => {
  await page.locator(FULL_SCREEN).first().hover();
}, {tier: 'ui', description: 'the pointer over the full-screen button above the app\'s preview'});

export const tooltipSays = Then('a tooltip should say {string}', async (page: Page, text: string) => {
  const tooltip = page.locator('.d4-tooltip').filter({visible: true}).last();
  await expect(tooltip, 'a tooltip').toBeVisible({timeout: pollMs(5000)});
  await expect(tooltip).toHaveText(text);
}, {tier: 'ui', description: 'the tooltip shown now says exactly this'});

export const clickFullScreen = When('user clicks the Reactions app\'s full-screen button', async (page: Page) => {
  await page.locator(FULL_SCREEN).first().click();
  await expect.poll(() => page.evaluate(() => DG.Dialog.getOpenDialogs().some((d: any) => d.title === 'Reaction')),
    {message: 'the full-screen dialog', timeout: pollMs(10000)}).toBe(true);
}, {tier: 'ui', description: 'a click on the preview\'s full-screen button; done when the "Reaction" dialog is open'});

export const dialogHelp = When('user opens the help of the {string} dialog', async (page: Page, title: string) => {
  const url = await page.evaluate((t) => {
    const dialog = DG.Dialog.getOpenDialogs().find((d: any) => d.title === t);
    if (!dialog)
      return null;
    grok.shell.windows.showHelp = true;
    grok.shell.windows.help.showHelp(dialog.helpUrl);
    return String(dialog.helpUrl ?? '');
  }, title);
  expect(url, `an open dialog titled "${title}"`).not.toBeNull();
  expect(url, `the "${title}" dialog's help link`).not.toBe('');
}, {tier: 'api', description: 'the help pane shown on the open dialog\'s help link (its helpUrl), as its help does'});
