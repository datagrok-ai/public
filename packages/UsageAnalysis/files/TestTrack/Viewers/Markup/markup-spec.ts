/* ---
realizes: []
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import {knownOpenBug} from '../../helpers/known-open-bug';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

const VIEWER_NAME = 'Markup';
const VIEWER = `[name="viewer-${VIEWER_NAME}"]`;
const datasetPath = 'System:DemoFiles/demog.csv';

const content = (page: Page) => page.locator(`${VIEWER} .grok-help`);

async function category(page: Page, cat: string, probe: string): Promise<void> {
  const own = page.locator('.property-grid tr[name="prop-markup-enabled"]');
  if (await own.count() === 0) {
    await v.clickViewerTitlebarIcon(page, VIEWER_NAME, 'icon-font-icon-settings');

    await own.first().waitFor({state: 'attached', timeout: 10_000});
  }
  await v.ensurePropertyCategory(page, VIEWER_NAME, cat, probe);
}

const currentRowValue = (page: Page, column: string) => page.evaluate((col) => {
  const t = (window as any).grok.shell.t;
  return String(t.col(col).get(t.currentRowIdx));
}, column);

async function editContentViaMenu(page: Page, text: string): Promise<string> {
  const box = (await page.locator(VIEWER).boundingBox())!;

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

  await softStep('Add Markup from the Viewers toolbox', async () => {
    await page.locator('[name="icon-markup"]').first().click();
    await page.locator(VIEWER).first().waitFor({timeout: 30_000});

    await expect(content(page).locator('h1')).toContainText('Markdown');
    await expect(content(page).locator('ul li')).toHaveCount(4);
    expect(await content(page).locator('a').count()).toBeGreaterThan(3);
    await expect(content(page)).not.toContainText('#');
  });

  await softStep('Edit content... replaces what the viewer renders', async () => {
    const before = await editContentViaMenu(page,
      '# Demographics\n\n* Age: ${AGE}\n* Sex: ${SEX}');

    expect(before).toContain('Markdown');

    await expect(content(page).locator('h1')).toHaveText('Demographics');
    await expect(content(page).locator('ul li')).toHaveCount(2);
  });

  await softStep('${COLUMN} renders the current row values', async () => {

    const grid = (await page.locator('.d4-grid[name="viewer-Grid"]').first().boundingBox())!;
    await page.mouse.click(grid.x + 150, grid.y + 140);
    await expect.poll(() => page.evaluate(() => (window as any).grok.shell.t.currentRowIdx),
      {timeout: 10_000}).toBeGreaterThan(-1);

    const age = await currentRowValue(page, 'AGE');
    const sex = await currentRowValue(page, 'SEX');
    await expect.poll(async () => await content(page).innerText(), {timeout: 10_000})
      .toContain(`Age: ${age}`);
    await expect(content(page)).toContainText(`Sex: ${sex}`);

    const row = await page.evaluate(() => (window as any).grok.shell.t.currentRowIdx as number);
    await page.keyboard.press('ArrowDown');
    await expect.poll(() => page.evaluate(() => (window as any).grok.shell.t.currentRowIdx),
      {timeout: 10_000}).not.toBe(row);
    const newAge = await currentRowValue(page, 'AGE');
    await expect.poll(async () => await content(page).innerText(), {timeout: 10_000})
      .toContain(`Age: ${newAge}`);
  });

  await softStep('An unknown ${COLUMN} is left as written', async () => {
    await editContentViaPropertyGrid(page, 'Unknown: ${NO_SUCH_COLUMN}');
    await expect(content(page)).toContainText('Unknown: ${NO_SUCH_COLUMN}');
  });

  await softStep('The Markup engine renders table expressions live', async () => {
    await editContentViaPropertyGrid(page, 'Rows: #{t.rowCount} Selected: #{t.selection.trueCount}');
    const rowCount = await page.evaluate(() => (window as any).grok.shell.t.rowCount as number);
    await expect(content(page)).toContainText(`Rows: ${rowCount}`);
    await expect(content(page)).toContainText('Selected: 0');

    await page.locator('[name="icon-select-all"]').first().click();
    await expect.poll(async () => (await content(page).innerText()).includes(`Selected: ${rowCount}`),
      {timeout: 10_000}).toBe(true);

    await page.locator('[name="icon-select-none"]').first().click();
    await expect.poll(async () => (await content(page).innerText()).includes('Selected: 0'),
      {timeout: 10_000}).toBe(true);
  });

  await softStep('Markup Enabled decides whether expressions are evaluated', async () => {
    await editContentViaPropertyGrid(page, 'Rows: #{t.rowCount}');
    const rowCount = await page.evaluate(() => (window as any).grok.shell.t.rowCount as number);
    await expect(content(page)).toContainText(`Rows: ${rowCount}`);

    await category(page, 'misc', 'markup-enabled');
    await v.setPropertyGridCheckbox(page, 'markup-enabled', false, 'misc');

    expect(await v.propertyGridValue(page, 'markup-enabled', 'misc')).toBe('false');

    await knownOpenBug('GROK-20637', async () => {
      await expect.poll(async () => await content(page).innerText(), {timeout: 8000})
        .toContain('#{t.rowCount}');
    });

    await v.setPropertyGridCheckbox(page, 'markup-enabled', true, 'misc');
  });

  await softStep('Mode decides how the content is interpreted', async () => {

    await editContentViaPropertyGrid(page, '# Heading probe');

    expect(await v.propertyGridValue(page, 'mode', 'misc')).toBe('Auto');
    await expect(content(page).locator('h1')).toHaveText('Heading probe');

    await v.selectPropertyGridChoice(page, 'mode', 'None', 'misc');
    expect(await v.propertyGridValue(page, 'mode', 'misc')).toBe('None');
    await expect(content(page).locator('h1')).toHaveCount(0);
    await expect(content(page).locator('pre')).toContainText('# Heading probe');

    await v.selectPropertyGridChoice(page, 'mode', 'Html', 'misc');
    await expect(content(page).locator('h1')).toHaveCount(0);
    await expect(content(page)).toContainText('# Heading probe');

    await v.selectPropertyGridChoice(page, 'mode', 'Markup', 'misc');
    await expect(content(page).locator('h1')).toHaveText('Heading probe');

    await v.selectPropertyGridChoice(page, 'mode', 'Auto', 'misc');
    expect(await v.propertyGridValue(page, 'mode', 'misc')).toBe('Auto');
  });

  await softStep('Mode = None shows HTML as source', async () => {
    await editContentViaPropertyGrid(page, '<b>bold probe</b> plain probe');
    await v.selectPropertyGridChoice(page, 'mode', 'None', 'misc');
    await expect(content(page).locator('pre')).toHaveCount(1);

    await knownOpenBug('GROK-20637', async () => {
      expect(await content(page).locator('b').count()).toBe(0);
      await expect(content(page).locator('pre')).toContainText('<b>bold probe</b>');
    });

    await v.selectPropertyGridChoice(page, 'mode', 'Auto', 'misc');
  });

  await softStep('Bold content is rendered bold', async () => {
    await v.selectPropertyGridChoice(page, 'mode', 'Markup', 'misc');
    await editContentViaPropertyGrid(page, '**md bold** plain tail');
    await expect(content(page).locator('strong')).toHaveText('md bold');

    const weight = await content(page).locator('strong').first()
      .evaluate((el) => Number(getComputedStyle(el).fontWeight));
    await knownOpenBug('GROK-20637', () => {
      expect(weight).toBeGreaterThanOrEqual(700);
    });

    await v.selectPropertyGridChoice(page, 'mode', 'Auto', 'misc');
  });

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

  await softStep('Close the viewer from its title bar', async () => {
    await v.clickViewerTitlebarIcon(page, VIEWER_NAME, 'Close');
    await expect(page.locator(VIEWER)).toHaveCount(0);
  });

  await v.cleanupShell(page);

  v.finishSpec();
});
