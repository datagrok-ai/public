/* ---
realizes: [pcplot.cp.panel-chart-filter-composition]
--- */
import {localTest as test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

test('PC Plot — Filter Panel + In-Chart Filter AND-Composition and Reset Scoping', async ({page}) => {
  test.setTimeout(300_000);

  await openDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await page.evaluate(() => {
    const icon = document.querySelector('[name="icon-pc-plot"]');
    if (icon)
      (icon as HTMLElement).click();
  });
  await page.locator('[name="viewer-PC-Plot"]').waitFor({timeout: 15000});

  await v.installEventWaits(page);

  await softStep('AND-composition ladder — panel filter drops the count, in-chart slider drops it further, Reset View returns to the panel-only value, re-drag drops it again, Reset filters restores the full count', async () => {
    await page.evaluate(() => grok.shell.tv.getFiltersGroup());
    await page.locator('.d4-filter-group-header').waitFor({timeout: 15000});
    const result = await page.evaluate(async () => {
      const w = window as any;
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      const df = grok.shell.tv.dataFrame;

      const settle = (stream: any, capMs: number, act?: () => void) => {
        const done = new Promise<void>((resolve) => {
          let sub: any = null;
          try { sub = stream.subscribe(() => { sub.unsubscribe(); resolve(); }); }
          catch (_) {  }
          setTimeout(() => { try { sub?.unsubscribe(); } catch (_) {} resolve(); }, capMs);
        });
        if (act) act();
        return done;
      };

      const settledFilterCount = async (capMs: number, act?: () => void) => {
        await settle(df.onRowsFiltered, capMs, act);
        let prev = -1;
        await w.__poll(() => df.filter.trueCount, (c: number) => {
          const stable = c === prev;
          prev = c;
          return stable;
        }, capMs, 50);
        return df.filter.trueCount;
      };
      await settle(pc.onViewerRendered, 800, () => {
        pc.props.columnNames = ['AGE', 'HEIGHT', 'WEIGHT'];
      });
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const vr = viewer.getBoundingClientRect();
      viewer.dispatchEvent(new MouseEvent('mousemove', {
        bubbles: true, clientX: vr.left + vr.width / 2, clientY: vr.top + vr.height / 2}));
      await w.__poll(() => document.querySelector('[name^="axis-slider-"]'),
        (el: Element | null) => el !== null, 400);
      const fullCount = df.filter.trueCount;

      const afterPanel = await settledFilterCount(700, () =>
        grok.shell.tv.getFiltersGroup().updateOrAdd({type: 'histogram', column: 'AGE', min: 30, max: 50}));

      const dragHeight = async () => {
        const svg = document.querySelector('[name="axis-slider-HEIGHT"]');
        if (!svg)
          return false;
        const maxHandle = svg.querySelector('[name="max-handle"]')!;
        const hr = maxHandle.getBoundingClientRect();
        const cx = hr.x + hr.width / 2;
        const cy = hr.y + hr.height / 2;
        const mk = (x: number, y: number) =>
          ({bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0});

        maxHandle.dispatchEvent(new MouseEvent('mousedown', mk(cx, cy)));
        await w.__drag(svg as HTMLElement, {x: cx, y: cy + 20}, {x: cx, y: cy + 200},
          {steps: 6, stepMs: 20, holdMs: 50});
        await settledFilterCount(500, () =>
          document.dispatchEvent(new MouseEvent('mouseup', mk(cx, cy + 200))));
        return true;
      };
      await dragHeight();
      const afterBoth = df.filter.trueCount;

      const canvas = viewer.querySelector('canvas[name="canvas"]')!;
      const cr = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: cr.left + cr.width / 2, clientY: cr.top + cr.height / 2}));

      const findResetView = () => Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find((el) => el.textContent!.trim() === 'Reset View');
      const rv: Element | undefined = await w.__poll(findResetView,
        (el: Element | undefined) => el !== undefined, 500, 25);
      const afterResetView = await settledFilterCount(700, () => {
        if (rv)
          (rv.closest('.d4-menu-item') as HTMLElement).click();
      });

      await dragHeight();
      const afterReDrag = df.filter.trueCount;
      const btn = document.querySelector('.d4-filter-group-header [name="icon-arrow-rotate-left"]') as HTMLElement | null;
      const afterPanelReset = await settledFilterCount(800, () => {
        if (btn)
          btn.click();
      });
      return {fullCount, afterPanel, afterBoth, afterResetView, afterReDrag, afterPanelReset};
    });

    expect(result.afterPanel).toBeLessThan(result.fullCount);
    expect(result.afterBoth).toBeLessThan(result.afterPanel);

    expect(result.afterResetView).toBe(result.afterPanel);

    expect(result.afterReDrag).toBeLessThan(result.afterPanel);
    expect(result.afterPanelReset).toBe(result.fullCount);
  });

  v.finishSpec();
});
