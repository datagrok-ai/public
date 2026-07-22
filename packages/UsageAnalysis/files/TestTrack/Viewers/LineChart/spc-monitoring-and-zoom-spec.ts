/* ---
realizes: [linechart.cp.spc-monitoring-and-zoom]
--- */
import {test, expect, type Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:AppData/Chem/tests/spgi-100.csv';

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

test('Line Chart — SPC Monitoring and Zoom', async ({page}) => {
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

  // Single Y, no split, no multi-axis — the SPC gating precondition.
  await setProps(page, {xColumnName: 'CAST Idea ID', yColumnNames: ['Chemical Space X']});

  await softStep('S1: SPC gating precondition — single Y, no split, no multi-axis', async () => {
    const cfg = await getProps(page, 'yColumnNames', 'splitColumnNames', 'multiAxis');
    expect(cfg.yColumnNames).toHaveLength(1);
    expect(cfg.splitColumnNames).toHaveLength(0);
    expect(cfg.multiAxis).toBe(false);
  });

  await softStep('S1: enable SPC — page stays responsive, no freeze (GROK-20126)', async () => {
    const before = realErrors().length;
    await setProps(page, {showStatisticalProcessControl: true});
    expect((await getProps(page, 'showStatisticalProcessControl')).showStatisticalProcessControl).toBe(true);
    // No-freeze guard: a follow-up JS roundtrip must resolve — a frozen page
    // would hang this evaluate until the test timeout.
    const alive = await page.evaluate(() => true);
    expect(alive).toBe(true);
    // grok.shell.warnings is undefined on this build, so page and console errors
    // are the error channel.
    expect(realErrors().length).toBe(before);
  });

  await softStep('S1: disable SPC — no-error teardown', async () => {
    const before = realErrors().length;
    await setProps(page, {showStatisticalProcessControl: false});
    expect((await getProps(page, 'showStatisticalProcessControl')).showStatisticalProcessControl).toBe(false);
    expect(realErrors().length).toBe(before);
  });

  await softStep('S2: Reset View clears the wheel-zoom but keeps explicit X Min/Max', async () => {
    const before = realErrors().length;
    const full = await page.evaluate(() => {
      const col = grok.shell.tv.dataFrame.columns.byName('CAST Idea ID');
      return {min: col.min, max: col.max};
    });
    const setMin = full.min + (full.max - full.min) * 0.1;
    const setMax = full.max - (full.max - full.min) * 0.1;
    await setProps(page, {xMin: setMin, xMax: setMax});
    const bounds = await getProps(page, 'xMin', 'xMax');
    expect(bounds.xMin).toBeGreaterThan(full.min);
    expect(bounds.xMax).toBeLessThan(full.max);

    // Zoom is canvas-only — no readable viewport prop — so the wheel-zoom and its
    // reset are checked via their events (d4-linechart-zoomed / -reset-view). The
    // explicit X Min/Max are the product-state signal: Reset View resets the zoom
    // but must NOT clear the configured axis bounds.
    const evt = await page.evaluate(async () => {
      const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
      let zoomFired = false;
      const subZ = lc.onEvent('d4-linechart-zoomed').subscribe(() => { zoomFired = true; });
      const cvs = Array.from(document.querySelectorAll('[name="viewer-Line-chart"] canvas')) as HTMLCanvasElement[];
      const t = cvs[cvs.length - 1];
      const r = t.getBoundingClientRect();
      for (let i = 0; i < 5; i++) {
        t.dispatchEvent(new WheelEvent('wheel', {bubbles: true, cancelable: true,
          clientX: r.left + r.width / 2, clientY: r.top + r.height / 2, deltaY: -200}));
        await new Promise((res) => setTimeout(res, 120));
      }
      await new Promise((res) => setTimeout(res, 400));
      subZ.unsubscribe();
      let resetFired = false;
      const subR = lc.onEvent('d4-linechart-reset-view').subscribe(() => { resetFired = true; });
      lc.resetView();
      await new Promise((res) => setTimeout(res, 700));
      subR.unsubscribe();
      return {zoomFired, resetFired};
    });
    expect(evt.zoomFired).toBe(true);
    expect(evt.resetFired).toBe(true);
    // Reset View cleared the zoom (event) but left the explicit bounds set.
    const after = await getProps(page, 'xMin', 'xMax');
    expect(after.xMin).toBe(setMin);
    expect(after.xMax).toBe(setMax);
    expect(realErrors().length).toBe(before);
  });

  if (stepErrors.length > 0)
    throw new Error(`Line Chart spc-monitoring-and-zoom failures:\n${stepErrors.join('\n')}`);
});
