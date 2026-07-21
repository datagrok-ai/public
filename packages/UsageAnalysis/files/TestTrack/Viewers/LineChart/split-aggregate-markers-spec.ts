/* ---
realizes: [linechart.cp.setup-split-aggregate-markers]
--- */
import {test, expect, type Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';

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
  await page.waitForTimeout(400);
}

async function getProps(page: Page, ...names: string[]): Promise<Record<string, any>> {
  return page.evaluate((ns) => {
    const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
    const out: Record<string, any> = {};
    for (const n of ns) out[n] = (lc.props as any)[n];
    return out;
  }, names);
}

async function hoverChart(page: Page) {
  const box = await page.evaluate(() => {
    const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
    const canvases = lc.root.querySelectorAll('canvas');
    let mc: HTMLCanvasElement | null = null, ma = 0;
    for (const c of canvases) {
      const r = (c as HTMLCanvasElement).getBoundingClientRect();
      if (r.width * r.height > ma) { ma = r.width * r.height; mc = c as HTMLCanvasElement; }
    }
    const rect = mc!.getBoundingClientRect();
    return {cx: rect.left + rect.width * 0.5, cy: rect.top + rect.height * 0.5};
  });
  await page.mouse.move(box.cx, box.cy);
  await page.waitForTimeout(300);
}

test('Line Chart — Setup, Split, Aggregate, Markers', async ({page}) => {
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

  await softStep('S1: set X and Y columns (connected trend)', async () => {
    const before = realErrors().length;
    await setProps(page, {xColumnName: 'CAST Idea ID', yColumnNames: ['Chemical Space X']});
    const props = await getProps(page, 'xColumnName', 'yColumnNames');
    expect(props.xColumnName).toBe('CAST Idea ID');
    expect(props.yColumnNames).toEqual(['Chemical Space X']);
    expect(realErrors().length).toBe(before);
  });

  await softStep('S1: split by Stereo Category', async () => {
    const before = realErrors().length;
    await setProps(page, {splitColumnName: 'Stereo Category'});
    expect((await getProps(page, 'splitColumnName')).splitColumnName).toBe('Stereo Category');
    expect(realErrors().length).toBe(before);
  });

  await softStep('S1: aggregation avg + whiskers std err', async () => {
    const before = realErrors().length;
    await setProps(page, {aggrType: 'avg', whiskersType: 'Avg | ±StError'});
    const props = await getProps(page, 'aggrType', 'whiskersType');
    expect(props.aggrType).toBe('avg');
    expect(props.whiskersType).toBe('Avg | ±StError');
    expect(realErrors().length).toBe(before);
  });

  // GROK-19883: marker type together with a size-coding column.
  await softStep('S1: marker type + size-coding column', async () => {
    const before = realErrors().length;
    await setProps(page, {markerType: 'circle', markersSizeColumnName: 'Chemical Space Y'});
    const props = await getProps(page, 'markerType', 'markersSizeColumnName');
    expect(props.markerType).toBe('circle');
    expect(props.markersSizeColumnName).toBe('Chemical Space Y');
    expect(realErrors().length).toBe(before);
  });

  // GROK-20255: a second Y column combined with a split.
  await softStep('S1: second Y column with split', async () => {
    const before = realErrors().length;
    await setProps(page, {yColumnNames: ['Chemical Space X', 'Chemical Space Y']});
    expect((await getProps(page, 'yColumnNames')).yColumnNames).toHaveLength(2);
    await hoverChart(page);
    expect(realErrors().length).toBe(before);
  });

  // GROK-17519: many split columns must not hang the page.
  await softStep('S1: add R1/R2/R3 splits, page stays responsive', async () => {
    const before = realErrors().length;
    for (const cols of [
      ['Stereo Category', 'R1'],
      ['Stereo Category', 'R1', 'R2'],
      ['Stereo Category', 'R1', 'R2', 'R3'],
    ]) {
      await page.evaluate((c) => {
        const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
        lc.props.splitColumnNames = c;
      }, cols);
      await page.waitForTimeout(500);
      // Page still responds: a follow-up evaluate resolves within timeout.
      const responsive = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount === 3624);
      expect(responsive).toBe(true);
    }
    expect((await getProps(page, 'splitColumnNames')).splitColumnNames).toHaveLength(4);
    expect(realErrors().length).toBe(before);
  });

  await softStep('S2: X axis logarithmic', async () => {
    const before = realErrors().length;
    await setProps(page, {xAxisType: 'logarithmic'});
    expect((await getProps(page, 'xAxisType')).xAxisType).toBe('logarithmic');
    expect(realErrors().length).toBe(before);
  });

  await softStep('S2: X axis back to linear', async () => {
    const before = realErrors().length;
    await setProps(page, {xAxisType: 'linear'});
    expect((await getProps(page, 'xAxisType')).xAxisType).toBe('linear');
    expect(realErrors().length).toBe(before);
  });

  await softStep('S2: clear all split columns', async () => {
    const before = realErrors().length;
    await page.evaluate(() => {
      const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
      lc.props.splitColumnNames = [];
      lc.props.splitColumnName = '';
    });
    await page.waitForTimeout(400);
    expect((await getProps(page, 'splitColumnNames')).splitColumnNames).toHaveLength(0);
    expect(realErrors().length).toBe(before);
  });

  if (stepErrors.length > 0)
    throw new Error(`Line Chart split-aggregate-markers failures:\n${stepErrors.join('\n')}`);
});
