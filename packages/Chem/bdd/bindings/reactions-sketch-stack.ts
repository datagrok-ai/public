// Steps on Crux Sketch's text tool and Arrange, for reaction-cell-editor-lossy.feature and
// reactions-app-arrange.feature. The test ids are crux-sketch's, inside its open shadow root, which Playwright's CSS
// locators pierce:
// - the text tool `toolbar.reaction.text`: a click just above an arrow's middle (about a third of a bond over its line)
// opens `canvas.text-editor` over that arrow, focused; typing writes the text; Escape ends it, as one edit;
// - Clean Up and Layout's opener `toolbar.clean.open`, its menu `menu.clean-up`, the item `menu.clean-up.arrange`.
// Both features need a Crux Sketch with that text tool and Arrange, as the vendored one is.
import {Page} from '@playwright/test';
import {Then, When} from '@datagrok-libraries/bdd';
import {expect, pollMs} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/** The reaction cell editor's sketcher, and the Reactions app's. */
const EDITOR_SKETCH = '.d4-dialog:has(.chem-reaction-cell-editor) crux-sketch';
const APP_SKETCH = '.chem-reactions-app-sketcher crux-sketch';
/** How far above an arrow's middle the text tool is clicked (CSS px): about a third of a bond. */
const ABOVE_ARROW_PX = 12;
/** How long a drawing change takes to reach the app's Text (reactions-app.ts: at once on a sketcher change; the
 * picture after DRAWING_DELAY_MS), waited before the Text is read as unchanged (ms). */
const SETTLE_MS = 1500;

/** The app's Text box's value (in the page). */
const appText = (page: Page) => page.evaluate(() =>
  (document.querySelector('.chem-reactions-app-text') as HTMLTextAreaElement | null)?.value ?? '');

const rememberedTexts = new WeakMap<Page, string>();

export const writeOverArrow = When('user writes {string} over arrow {int} in the reaction editor',
  async (page: Page, text: string, n: number) => {
    await page.locator(`${EDITOR_SKETCH} [data-testid="toolbar.reaction.text"]`).first().click();
    const arrow = page.locator(`${EDITOR_SKETCH} g[data-arrow]`).nth(n - 1);
    await expect(arrow, `arrow ${n} drawn`).toBeVisible({timeout: pollMs(10000)});
    const b = (await arrow.boundingBox())!;
    await page.mouse.click(b.x + b.width / 2, b.y + b.height / 2 - ABOVE_ARROW_PX);
    const editor = page.locator(`${EDITOR_SKETCH} [data-testid="canvas.text-editor"]`).first();
    await expect(editor, 'the text editor over the arrow').toBeVisible({timeout: pollMs(5000)});
    await page.keyboard.type(text);
    await page.keyboard.press('Escape');
    await expect(editor, 'the text editor closed').toBeHidden({timeout: pollMs(5000)});
  }, {tier: 'ui', description: 'the text tool of the reaction editor\'s sketcher, a click just above the middle of that ' +
    'arrow (arrows counted from 1, in drawing order), the text typed, Escape: one edit'});

export const valueLacks = Then('the value in row {int} of {string} column does not contain {string}',
  async (page: Page, row: number, column: string, part: string) => {
    const value = await page.evaluate(([r, c]) => String(grok.shell.tv.dataFrame.get(c, r - 1)),
      [row, column] as [number, string]);
    expect(value.includes(part), `row ${row} of ${column}: ${value}`).toBe(false);
  }, {tier: 'api', description: 'the cell\'s value (row 1 the first) does not hold the text named'});

export const rememberAppText = When('user remembers the Reactions app\'s Text', async (page: Page) => {
  rememberedTexts.set(page, await appText(page));
}, {tier: 'api', description: 'the text box under the app\'s sketcher, kept for "as remembered"'});

export const appArrange = When('user chooses {string} in the Clean Up and Layout menu of the Reactions app\'s sketcher',
  async (page: Page, item: string) => {
    const id = ({'Arrange': 'arrange', 'Clean Up': 'clean', 'Layout': 'layout'} as {[k: string]: string})[item];
    expect(id, `an item of Clean Up and Layout: ${item}`).toBeTruthy();
    await page.locator(`${APP_SKETCH} [data-testid="toolbar.clean.open"]`).first().click();
    const entry = page.locator(`${APP_SKETCH} [data-testid="menu.clean-up.${id}"]`).first();
    await expect(entry, `"${item}" in the menu`).toBeVisible({timeout: pollMs(5000)});
    await entry.click();
  }, {tier: 'ui', description: 'the app\'s sketcher\'s Clean Up and Layout menu opened, the item clicked (Arrange, ' +
    'Clean Up or Layout)'});

export const appTextAsRemembered = Then('the Reactions app\'s Text should be as remembered', async (page: Page) => {
  const remembered = rememberedTexts.get(page);
  expect(remembered, '"user remembers the Reactions app\'s Text" first').toBeDefined();
  // a change would have reached the Text by now: read it once the drawing has settled
  await page.waitForTimeout(SETTLE_MS);
  expect(await appText(page), 'the Text').toBe(remembered);
}, {tier: 'api', description: 'once the drawing has settled, the text box holds exactly the remembered text'});
