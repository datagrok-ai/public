import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';

test.use(specTestOptions);

const BLUE = 'rgb(80, 169, 197)';

test('Add New Column — column name highlighting in expression editor', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    // @ts-ignore
    grok.shell.settings.showFiltersIconsConstantly = true;
    // @ts-ignore
    grok.shell.windows.simpleMode = true;
    // @ts-ignore
    grok.shell.closeAll();
    // @ts-ignore
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    // @ts-ignore
    grok.shell.addTableView(df);
    await new Promise(resolve => {
      // @ts-ignore
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
      setTimeout(() => resolve(null), 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('1. Open demog.csv dataset', async () => {
    const rows = await page.evaluate(() => (globalThis as any).grok.shell.tv.dataFrame.rowCount);
    expect(rows).toBe(5850);
  });

  await softStep('2. Open Add New Column dialog', async () => {
    await page.locator('[name="icon-add-new-column"]').click();
    await page.locator('[name="dialog-Add-New-Column"]').waitFor({timeout: 10_000});
    await page.locator('.d4-dialog .cm-content').waitFor({timeout: 5_000});
  });

  await softStep('3a. Paste Abs(${AGE}) — ${AGE} highlighted in blue', async () => {
    await page.locator('.d4-dialog .cm-content').click();
    await page.keyboard.press('Control+A');
    await page.keyboard.press('Delete');
    await page.keyboard.type('Abs(${AGE})');
    await page.waitForTimeout(300);
    const token = await page.evaluate(() => {
      const cm = document.querySelector('.d4-dialog .cm-content')!;
      const span = Array.from(cm.querySelectorAll('span'))
        .find(s => s.classList.contains('cm-column-name') && s.textContent === '${AGE}');
      return span
        ? {color: getComputedStyle(span as HTMLElement).color, hasClass: span.classList.contains('cm-column-name')}
        : null;
    });
    expect(token).not.toBeNull();
    expect(token!.hasClass).toBe(true);
    expect(token!.color).toBe(BLUE);
  });

  await softStep('3b. Paste Avg($[AGE]) — $[AGE] highlighted in blue', async () => {
    await page.locator('.d4-dialog .cm-content').click();
    await page.keyboard.press('Control+A');
    await page.keyboard.press('Delete');
    await page.keyboard.type('Avg($[AGE])');
    await page.waitForTimeout(300);
    const token = await page.evaluate(() => {
      const cm = document.querySelector('.d4-dialog .cm-content')!;
      const span = Array.from(cm.querySelectorAll('span'))
        .find(s => s.classList.contains('cm-column-name') && s.textContent === '$[AGE]');
      return span
        ? {color: getComputedStyle(span as HTMLElement).color, hasClass: span.classList.contains('cm-column-name')}
        : null;
    });
    expect(token).not.toBeNull();
    expect(token!.hasClass).toBe(true);
    expect(token!.color).toBe(BLUE);
  });

  await softStep('4. Add function + column via autocomplete — column highlighted in blue', async () => {
    await page.locator('.d4-dialog .cm-content').click();
    await page.keyboard.press('Control+A');
    await page.keyboard.press('Delete');
    await page.keyboard.type('Sqrt(');
    await page.keyboard.type('$');
    await page.locator('.cm-tooltip-autocomplete').waitFor({timeout: 5_000});
    await page.keyboard.type('HE');
    await page.waitForFunction(() => {
      const li = document.querySelector('.cm-tooltip-autocomplete li[aria-selected="true"]');
      return li?.textContent === 'HEIGHT';
    }, null, {timeout: 5_000});
    // CodeMirror autocomplete inserts on mousedown; Playwright's .click() dispatches
    // the full sequence and can race with the popup's stability check. Dispatch
    // mousedown directly so the completion is inserted reliably.
    await page.evaluate(() => {
      const popup = document.querySelector('.cm-tooltip-autocomplete')!;
      const li = Array.from(popup.querySelectorAll('li')).find(l => l.textContent === 'HEIGHT')!;
      li.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, cancelable: true, button: 0}));
    });
    await page.waitForFunction(() => {
      const cm = document.querySelector('.d4-dialog .cm-content');
      return cm?.textContent === 'Sqrt(${HEIGHT}';
    }, null, {timeout: 5_000});
    await page.keyboard.type(')');
    await page.waitForTimeout(500);
    const token = await page.evaluate(() => {
      const cm = document.querySelector('.d4-dialog .cm-content');
      if (!cm) return {dialogClosed: true} as any;
      // CM retokenizes after each keystroke. The HEIGHT reference may render as
      // a single span "${HEIGHT}" or as several spans ($, {, HEIGHT, }) — all
      // carrying the .cm-column-name class. Accept either shape: find any
      // cm-column-name span whose text contains "HEIGHT".
      const span = Array.from(cm.querySelectorAll('span.cm-column-name'))
        .find(s => s.textContent?.includes('HEIGHT'));
      return span
        ? {color: getComputedStyle(span as HTMLElement).color, hasClass: span.classList.contains('cm-column-name'), text: span.textContent}
        : null;
    });
    expect(token).not.toBeNull();
    expect((token as any).dialogClosed).toBeUndefined();
    expect((token as any).hasClass).toBe(true);
    expect((token as any).color).toBe(BLUE);
  });

  await page.locator('[name="button-CANCEL"]').click().catch(() => {});

  if (stepErrors.length > 0)
    throw new Error('Soft step failures:\n' + stepErrors.map(e => `- ${e.step}: ${e.error}`).join('\n'));
});
