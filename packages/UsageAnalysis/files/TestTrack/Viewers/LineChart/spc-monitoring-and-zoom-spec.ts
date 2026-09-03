/* ---
realizes: [linechart.cp.spc-monitoring-and-zoom]
--- */
import {expect, type Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, isLocalBootNoise} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:AppData/Chem/tests/spgi-100.csv';

async function setProps(page: Page, props: Record<string, any>) {
  await v.setViewerProps(page, 'Line chart', [{set: props, wait: 500}]);
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

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error' && !isLocalBootNoise(m.text())) consoleErrors.push(m.text()); });
  const errorCount = () => pageErrors.length + consoleErrors.length;

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'line-chart', 'Line-chart', 15_000, 'Line chart');

  await setProps(page, {xColumnName: 'CAST Idea ID', yColumnNames: ['Chemical Space X']});

  await softStep('S1: SPC gating precondition — single Y, no split, no multi-axis', async () => {
    const cfg = await getProps(page, 'yColumnNames', 'splitColumnNames', 'multiAxis');
    expect(cfg.yColumnNames).toHaveLength(1);
    expect(cfg.splitColumnNames).toHaveLength(0);
    expect(cfg.multiAxis).toBe(false);
  });

  await softStep('S1: enable SPC — page stays responsive, no freeze (GROK-20126)', async () => {
    const before = errorCount();
    await setProps(page, {showStatisticalProcessControl: true});
    expect((await getProps(page, 'showStatisticalProcessControl')).showStatisticalProcessControl).toBe(true);

    const alive = await page.evaluate(() => true);
    expect(alive).toBe(true);

    expect(errorCount()).toBe(before);
  });

  await softStep('S1: disable SPC — no-error teardown', async () => {
    const before = errorCount();
    await setProps(page, {showStatisticalProcessControl: false});
    expect((await getProps(page, 'showStatisticalProcessControl')).showStatisticalProcessControl).toBe(false);
    expect(errorCount()).toBe(before);
  });

  await softStep('S2: Reset View clears the wheel-zoom but keeps explicit X Min/Max', async () => {
    const before = errorCount();
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

    const evt = await page.evaluate(async () => {
      const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
      const until = async (fired: () => boolean, cap: number) => {
        const t0 = Date.now();
        while (!fired() && Date.now() - t0 < cap)
          await new Promise((res) => setTimeout(res, 20));
      };
      let zoomFired = false;
      const subZ = lc.onEvent('d4-linechart-zoomed').subscribe(() => { zoomFired = true; });
      const cvs = Array.from(document.querySelectorAll('[name="viewer-Line-chart"] canvas')) as HTMLCanvasElement[];
      const t = cvs[cvs.length - 1];
      const r = t.getBoundingClientRect();
      // the wheel steps are PACED: the zoom coalesces events that arrive inside one frame
      for (let i = 0; i < 5; i++) {
        t.dispatchEvent(new WheelEvent('wheel', {bubbles: true, cancelable: true,
          clientX: r.left + r.width / 2, clientY: r.top + r.height / 2, deltaY: -200}));
        await new Promise((res) => setTimeout(res, 120));
      }
      await until(() => zoomFired, 400);
      subZ.unsubscribe();
      let resetFired = false;
      const subR = lc.onEvent('d4-linechart-reset-view').subscribe(() => { resetFired = true; });
      lc.resetView();
      await until(() => resetFired, 700);
      subR.unsubscribe();
      return {zoomFired, resetFired};
    });
    expect(evt.zoomFired).toBe(true);
    expect(evt.resetFired).toBe(true);

    const after = await getProps(page, 'xMin', 'xMax');
    expect(after.xMin).toBe(setMin);
    expect(after.xMax).toBe(setMax);
    expect(errorCount()).toBe(before);
  });

  await v.cleanupShell(page);
  v.finishSpec();
});
