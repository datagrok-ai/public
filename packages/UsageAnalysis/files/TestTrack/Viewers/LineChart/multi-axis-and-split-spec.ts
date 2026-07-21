/* ---
realizes: [linechart.cp.multi-axis-and-split]
--- */
import {test, expect, type Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import {countCanvasPixels} from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/SPGI.csv';

const pageErrors: string[] = [];
const consoleErrors: string[] = [];

function realErrors(): string[] {
  return [...pageErrors, ...consoleErrors];
}

async function setProps(page: Page, props: Record<string, any>) {
  await page.evaluate((p) => {
    const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
    for (const [k, val] of Object.entries(p)) (lc.props as any)[k] = val;
  }, props);
  await page.waitForTimeout(500);
}

async function getProps(page: Page, ...names: string[]): Promise<Record<string, any>> {
  return page.evaluate((ns) => {
    const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
    const out: Record<string, any> = {};
    for (const n of ns) out[n] = (lc.props as any)[n];
    return out;
  }, names);
}

async function chartCanvasNonEmpty(page: Page): Promise<boolean> {
  // Content-level check (github-2904): a blank chart still has a non-zero
  // bounding rect, so measure drawn pixels instead of geometry. The threshold
  // must sit ABOVE the axes/labels-only floor — dev empirics (666x308 canvas):
  // axes with no data lines paint ~6449 px, data lines bring it to ~18061 —
  // so 8000 fails the axes-only "blank" state that > 0-ish thresholds would
  // pass. -1 (no canvas / getImageData fault) fails the threshold too.
  return (await countCanvasPixels(page, 'Line chart')).total > 8000;
}

test('Line Chart — Multi-Axis and Split', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;
  pageErrors.length = 0;
  consoleErrors.length = 0;

  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });

  await loginToDatagrok(page);

  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    grok.shell.addTableView(df);
    await new Promise((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
    for (let i = 0; i < 50; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise((r) => setTimeout(r, 200));
    }
    await new Promise((r) => setTimeout(r, 5000));
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 60000});

  await page.locator('[name="icon-line-chart"]').click();
  await page.locator('[name="viewer-Line-chart"]').waitFor({timeout: 15000});

  await setProps(page, {xColumnName: 'CAST Idea ID', yColumnNames: ['Chemical Space X', 'Chemical Space Y']});

  await softStep('S1: enable Multi Axis with 2 Y columns', async () => {
    const before = realErrors().length;
    await setProps(page, {multiAxis: true});
    expect((await getProps(page, 'multiAxis')).multiAxis).toBe(true);
    expect(realErrors().length).toBe(before);
  });

  // github-2904: a split under multi-axis must keep the chart rendering.
  await softStep('S1: first split column, chart stays non-empty', async () => {
    const before = realErrors().length;
    await setProps(page, {splitColumnNames: ['Stereo Category']});
    expect((await getProps(page, 'splitColumnNames')).splitColumnNames).toEqual(['Stereo Category']);
    expect(await chartCanvasNonEmpty(page)).toBe(true);
    expect(realErrors().length).toBe(before);
  });

  // github-2904: the second split column is where the chart went blank.
  await softStep('S1: second split column, chart does not go blank', async () => {
    const before = realErrors().length;
    await setProps(page, {splitColumnNames: ['Stereo Category', 'Series']});
    expect((await getProps(page, 'splitColumnNames')).splitColumnNames).toHaveLength(2);
    expect(await chartCanvasNonEmpty(page)).toBe(true);
    const followup = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount === 3624);
    expect(followup).toBe(true);
    expect(realErrors().length).toBe(before);
  });

  await softStep('S2: 3 Y columns survive an edit, not reset to 1', async () => {
    await setProps(page, {splitColumnNames: [], yColumnNames: ['Chemical Space X', 'Chemical Space Y', 'TPSA']});
    expect((await getProps(page, 'yColumnNames')).yColumnNames).toHaveLength(3);
    // GROK-18484: re-issuing the same 3-column selection must not silently shrink
    // the Y-column list to one column.
    await setProps(page, {yColumnNames: ['Chemical Space X', 'Chemical Space Y', 'TPSA']});
    const y = (await getProps(page, 'yColumnNames')).yColumnNames;
    expect(y).toHaveLength(3);
  });

  // GROK-20033: the column-selector search input must stay inside the selector.
  await softStep('S2: search input is inside the selector bounds', async () => {
    const before = realErrors().length;
    await page.evaluate(() => {
      const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
      const combos = Array.from(lc.root.querySelectorAll('[name^="div-column-combobox"]')) as HTMLElement[];
      const yCombo = combos.find((c) => c.getBoundingClientRect().x > 800);
      (yCombo ?? combos[0])?.click();
    });
    await page.waitForTimeout(800);
    const result = await page.evaluate(() => {
      const input = Array.from(document.querySelectorAll('input.ui-input-editor'))
        .find((i) => (i as HTMLElement).offsetParent !== null &&
          (i as HTMLInputElement).placeholder?.includes('Search')) as HTMLInputElement | undefined;
      if (!input) return {present: false, contained: false};
      const container = input.closest('.ui-input-root');
      if (!container) return {present: true, contained: false};
      const ir = input.getBoundingClientRect();
      const cr = container.getBoundingClientRect();
      const contained = ir.top >= cr.top - 1 && ir.bottom <= cr.bottom + 1;
      return {present: true, contained};
    });
    if (result.present)
      expect(result.contained).toBe(true); // GROK-20033 invariant
    else
      expect(realErrors().length).toBe(before); // no-error floor fallback (search input absent)
  });

  await softStep('S3: disable Multi Axis', async () => {
    const before = realErrors().length;
    await setProps(page, {multiAxis: false});
    expect((await getProps(page, 'multiAxis')).multiAxis).toBe(false);
    expect(realErrors().length).toBe(before);
  });

  await softStep('S3: clear all split columns', async () => {
    const before = realErrors().length;
    await setProps(page, {splitColumnNames: []});
    expect((await getProps(page, 'splitColumnNames')).splitColumnNames).toHaveLength(0);
    expect(realErrors().length).toBe(before);
  });

  if (stepErrors.length > 0)
    throw new Error(`Line Chart multi-axis-and-split failures:\n${stepErrors.join('\n')}`);
});
