import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import {knownOpenBug} from '../../helpers/known-open-bug';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

const VIEWER_NAME = 'Markup';
const VIEWER = `[name="viewer-${VIEWER_NAME}"]`;
const datasetPath = 'System:DemoFiles/demog.csv';

/** The rendered content — the viewer draws into the DOM, not a canvas. */
const content = (page: Page) => page.locator(`${VIEWER} .grok-help`);

/**
 * Bring a property-grid category of THIS viewer within reach. Clicking the grid
 * (or anything else in the view) replaces the Context Panel content, so the probe
 * has to be a row only the Markup viewer has — a bare `.property-grid` check
 * would happily accept the Grid's own properties.
 */
async function category(page: Page, cat: string, probe: string): Promise<void> {
  const own = page.locator('.property-grid tr[name="prop-markup-enabled"]');
  if (await own.count() === 0) {
    await v.clickViewerTitlebarIcon(page, VIEWER_NAME, 'icon-font-icon-settings');
    // Attached, not visible: the row sits in a category that may be collapsed.
    await own.first().waitFor({state: 'attached', timeout: 10_000});
  }
  await v.ensurePropertyCategory(page, VIEWER_NAME, cat, probe);
}

/** Value the table holds for the current row — what a `${COLUMN}` must render as. */
const currentRowValue = (page: Page, column: string) => page.evaluate((col) => {
  const t = (window as any).grok.shell.t;
  return String(t.col(col).get(t.currentRowIdx));
}, column);

/**
 * Replace the content through the Edit dialog, the viewer's own context-menu
 * path. Returns what the dialog was pre-filled with, so a caller can assert that
 * it offered the content currently on screen.
 */
async function editContentViaMenu(page: Page, text: string): Promise<string> {
  const box = (await page.locator(VIEWER).boundingBox())!;
  // Bottom of the viewer: the sample content is full of links, and a right-click
  // on one of them opens the browser's own menu instead.
  await page.mouse.click(box.x + box.width / 2, box.y + box.height - 40, {button: 'right'});
  await page.locator('.d4-menu-popup [name="div-Edit-content..."]').first().click();

  const editor = page.locator('.d4-dialog textarea.ui-input-editor').first();
  await editor.waitFor({timeout: 10_000});
  const before = await editor.inputValue();
  await editor.fill(text);
  await page.locator('.d4-dialog button[name="button-OK"]').click();
  await expect(page.locator('.d4-dialog')).toHaveCount(0);
  return before;
}

/**
 * Replace the content through the property grid's ellipsis editor, which opens a
 * *Content* dialog backed by CodeMirror — so the text goes in as real keystrokes,
 * not as a value assignment. Keep `text` on one line: the markdown editor
 * continues lists and indentation by itself.
 */
async function editContentViaPropertyGrid(page: Page, text: string): Promise<void> {
  await category(page, 'misc', 'content');
  await page.locator('.property-grid tr[name="prop-content"] button.property-grid-ellipsis-editor-ellipsis')
    .first().click();
  const editor = page.locator('.d4-dialog .CodeMirror').first();
  await editor.waitFor({timeout: 10_000});
  await editor.click();
  await page.keyboard.press('Control+a');
  await page.keyboard.type(text);
  await page.locator('.d4-dialog button[name="button-OK"]').click();
  await expect(page.locator('.d4-dialog')).toHaveCount(0);
}

test('Markup', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  // #### Add the viewer
  await softStep('Add Markup from the Viewers toolbox', async () => {
    await page.locator('[name="icon-markup"]').first().click();
    await page.locator(VIEWER).first().waitFor({timeout: 30_000});

    // The default content is a Markdown sample, and it arrives rendered: a real
    // heading, a real list and real links — not the markup source.
    await expect(content(page).locator('h1')).toContainText('Markdown');
    await expect(content(page).locator('ul li')).toHaveCount(4);
    expect(await content(page).locator('a').count()).toBeGreaterThan(3);
    await expect(content(page)).not.toContainText('#');
  });

  // #### Editing the content from the viewer's own menu
  await softStep('Edit content... replaces what the viewer renders', async () => {
    const before = await editContentViaMenu(page,
      '# Demographics\n\n* Age: ${AGE}\n* Sex: ${SEX}');
    // The dialog offers the content that is on screen, not an empty box.
    expect(before).toContain('Markdown');

    await expect(content(page).locator('h1')).toHaveText('Demographics');
    await expect(content(page).locator('ul li')).toHaveCount(2);
  });

  // #### Column references follow the current row
  await softStep('${COLUMN} renders the current row values', async () => {
    // Nothing is current on a freshly opened table, and the viewer substitutes
    // only for a current row — so the row is picked in the UI first, by clicking
    // a grid cell.
    // Well inside the data area: the leftmost strip is the row header, and a
    // click there selects rows instead of making one current.
    const grid = (await page.locator('.d4-grid[name="viewer-Grid"]').first().boundingBox())!;
    await page.mouse.click(grid.x + 150, grid.y + 140);
    await expect.poll(() => page.evaluate(() => (window as any).grok.shell.t.currentRowIdx),
      {timeout: 10_000}).toBeGreaterThan(-1);

    const age = await currentRowValue(page, 'AGE');
    const sex = await currentRowValue(page, 'SEX');
    await expect.poll(async () => await content(page).innerText(), {timeout: 10_000})
      .toContain(`Age: ${age}`);
    await expect(content(page)).toContainText(`Sex: ${sex}`);

    // Move to another row with the keyboard — the text follows it.
    const row = await page.evaluate(() => (window as any).grok.shell.t.currentRowIdx as number);
    await page.keyboard.press('ArrowDown');
    await expect.poll(() => page.evaluate(() => (window as any).grok.shell.t.currentRowIdx),
      {timeout: 10_000}).not.toBe(row);
    const newAge = await currentRowValue(page, 'AGE');
    await expect.poll(async () => await content(page).innerText(), {timeout: 10_000})
      .toContain(`Age: ${newAge}`);
  });

  // #### A column that does not exist stays literal
  await softStep('An unknown ${COLUMN} is left as written', async () => {
    await editContentViaPropertyGrid(page, 'Unknown: ${NO_SUCH_COLUMN}');
    await expect(content(page)).toContainText('Unknown: ${NO_SUCH_COLUMN}');
  });

  // #### Markup engine expressions
  await softStep('The Markup engine renders table expressions live', async () => {
    await editContentViaPropertyGrid(page, 'Rows: #{t.rowCount} Selected: #{t.selection.trueCount}');
    const rowCount = await page.evaluate(() => (window as any).grok.shell.t.rowCount as number);
    await expect(content(page)).toContainText(`Rows: ${rowCount}`);
    await expect(content(page)).toContainText('Selected: 0');

    // Select every row from the toolbar — the expression tracks the selection.
    await page.locator('[name="icon-select-all"]').first().click();
    await expect.poll(async () => (await content(page).innerText()).includes(`Selected: ${rowCount}`),
      {timeout: 10_000}).toBe(true);

    await page.locator('[name="icon-select-none"]').first().click();
    await expect.poll(async () => (await content(page).innerText()).includes('Selected: 0'),
      {timeout: 10_000}).toBe(true);
  });

  // #### Markup Enabled
  await softStep('Markup Enabled decides whether expressions are evaluated', async () => {
    await editContentViaPropertyGrid(page, 'Rows: #{t.rowCount}');
    const rowCount = await page.evaluate(() => (window as any).grok.shell.t.rowCount as number);
    await expect(content(page)).toContainText(`Rows: ${rowCount}`);

    await category(page, 'misc', 'markup-enabled');
    await v.setPropertyGridCheckbox(page, 'markup-enabled', false, 'misc');
    // The property is written — it is simply never read.
    expect(await v.propertyGridValue(page, 'markup-enabled', 'misc')).toBe('false');

    // GROK-20637 #1: with the engine off the expression should be left as
    // written. `markupEnabled` is declared in the look and read nowhere, so
    // `_refresh()` keeps calling Markup.render and the count stays resolved.
    await knownOpenBug('GROK-20637', async () => {
      await expect.poll(async () => await content(page).innerText(), {timeout: 8000})
        .toContain('#{t.rowCount}');
    });

    await v.setPropertyGridCheckbox(page, 'markup-enabled', true, 'misc');
  });

  // #### Interpretation mode
  await softStep('Mode decides how the content is interpreted', async () => {
    // One Markdown heading tells all four modes apart: it becomes an <h1> only
    // where the text is read as Markdown.
    await editContentViaPropertyGrid(page, '# Heading probe');

    // Auto: no leading angle bracket, so the content is taken as Markdown.
    expect(await v.propertyGridValue(page, 'mode', 'misc')).toBe('Auto');
    await expect(content(page).locator('h1')).toHaveText('Heading probe');

    // None: nothing is interpreted, the source is shown as preformatted text.
    await v.selectPropertyGridChoice(page, 'mode', 'None', 'misc');
    expect(await v.propertyGridValue(page, 'mode', 'misc')).toBe('None');
    await expect(content(page).locator('h1')).toHaveCount(0);
    await expect(content(page).locator('pre')).toContainText('# Heading probe');

    // Html: the content goes in as HTML, where a hash is just a hash.
    await v.selectPropertyGridChoice(page, 'mode', 'Html', 'misc');
    await expect(content(page).locator('h1')).toHaveCount(0);
    await expect(content(page)).toContainText('# Heading probe');

    // Markup: Markdown again, whatever the content looks like.
    await v.selectPropertyGridChoice(page, 'mode', 'Markup', 'misc');
    await expect(content(page).locator('h1')).toHaveText('Heading probe');

    await v.selectPropertyGridChoice(page, 'mode', 'Auto', 'misc');
    expect(await v.propertyGridValue(page, 'mode', 'misc')).toBe('Auto');
  });

  // #### None does not escape HTML
  await softStep('Mode = None shows HTML as source', async () => {
    await editContentViaPropertyGrid(page, '<b>bold probe</b> plain probe');
    await v.selectPropertyGridChoice(page, 'mode', 'None', 'misc');
    await expect(content(page).locator('pre')).toHaveCount(1);

    // GROK-20637 #2: None wraps the content in <pre> without escaping it, so the
    // browser parses the tag instead of showing it. The markers vanish and the
    // source can no longer be read — the element surviving in the DOM is the
    // proof, since the text itself does not look bold (that is #3 below).
    await knownOpenBug('GROK-20637', async () => {
      expect(await content(page).locator('b').count()).toBe(0);
      await expect(content(page).locator('pre')).toContainText('<b>bold probe</b>');
    });

    await v.selectPropertyGridChoice(page, 'mode', 'Auto', 'misc');
  });

  // #### Bold
  await softStep('Bold content is rendered bold', async () => {
    await v.selectPropertyGridChoice(page, 'mode', 'Markup', 'misc');
    await editContentViaPropertyGrid(page, '**md bold** plain tail');
    await expect(content(page).locator('strong')).toHaveText('md bold');

    // GROK-20637 #3: `.grok-help` sets the relative `font-weight: lighter`,
    // which computes to 100, and the browser's `b, strong { font-weight: bolder }`
    // is relative too — against 100 it resolves to 400. The two cancel out and
    // emphasis lands on the document's ordinary weight.
    const weight = await content(page).locator('strong').first()
      .evaluate((el) => Number(getComputedStyle(el).fontWeight));
    await knownOpenBug('GROK-20637', () => {
      expect(weight).toBeGreaterThanOrEqual(700);
    });

    await v.selectPropertyGridChoice(page, 'mode', 'Auto', 'misc');
  });

  // #### Title
  await softStep('Title and Show Title drive the viewer title bar', async () => {
    await category(page, 'description', 'title');
    await v.setPropertyGridValue(page, 'title', 'Patient card', 'description');
    expect(await v.propertyGridValue(page, 'title', 'description')).toBe('Patient card');

    const titleBar = page.locator(`${VIEWER} .d4-viewer-title textarea`).first();
    await v.setPropertyGridCheckbox(page, 'show-title', true, 'description');
    await expect(titleBar).toBeVisible();
    expect(await titleBar.inputValue()).toBe('Patient card');

    await v.setPropertyGridCheckbox(page, 'show-title', false, 'description');
    await expect(titleBar).toBeHidden();
  });

  // #### Closing the viewer
  await softStep('Close the viewer from its title bar', async () => {
    await v.clickViewerTitlebarIcon(page, VIEWER_NAME, 'Close');
    await expect(page.locator(VIEWER)).toHaveCount(0);
  });

  await v.cleanupShell(page);

  v.finishSpec();
});
