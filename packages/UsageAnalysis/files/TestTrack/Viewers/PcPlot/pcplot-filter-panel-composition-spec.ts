/* ---
realizes: [pcplot.cp.panel-chart-filter-composition]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';


declare const grok: any;


test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';


test('PC Plot — Filter Panel + In-Chart Filter AND-Composition and Reset Scoping', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await page.evaluate(() => {
    const icon = document.querySelector('[name="icon-pc-plot"]');
    if (icon)
      (icon as HTMLElement).click();
  });
  await page.locator('[name="viewer-PC-Plot"]').waitFor({timeout: 15000});

  // Both the Filter Panel and the in-chart range sliders write the shared
  // df.filter (AND-combined). Reset View clears ONLY the in-chart part; the
  // Filter Panel "Reset filters" button clears everything back to full.
  await softStep('AND-composition ladder — panel filter drops the count, in-chart slider drops it further, Reset View returns to the panel-only value, re-drag drops it again, Reset filters restores the full count', async () => {
    await page.evaluate(() => grok.shell.tv.getFiltersGroup());
    await page.locator('.d4-filter-group-header').waitFor({timeout: 15000});
    const result = await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      pc.props.columnNames = ['AGE', 'HEIGHT', 'WEIGHT'];
      await new Promise((r) => setTimeout(r, 800));
      const df = grok.shell.tv.dataFrame;
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const vr = viewer.getBoundingClientRect();
      viewer.dispatchEvent(new MouseEvent('mousemove', {
        bubbles: true, clientX: vr.left + vr.width / 2, clientY: vr.top + vr.height / 2}));
      await new Promise((r) => setTimeout(r, 400));
      const fullCount = df.filter.trueCount;

      // Filter Panel histogram narrows AGE.
      grok.shell.tv.getFiltersGroup().updateOrAdd({type: 'histogram', column: 'AGE', min: 30, max: 50});
      await new Promise((r) => setTimeout(r, 700));
      const afterPanel = df.filter.trueCount;

      // In-chart range slider narrows HEIGHT on top of the Filter Panel filter.
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
        await new Promise((r) => setTimeout(r, 50));
        for (let dy = 20; dy <= 200; dy += 30) {
          document.dispatchEvent(new MouseEvent('mousemove', mk(cx, cy + dy)));
          svg.dispatchEvent(new MouseEvent('mousemove', mk(cx, cy + dy)));
          await new Promise((r) => setTimeout(r, 20));
        }
        document.dispatchEvent(new MouseEvent('mouseup', mk(cx, cy + 200)));
        await new Promise((r) => setTimeout(r, 500));
        return true;
      };
      await dragHeight();
      const afterBoth = df.filter.trueCount;

      // Reset View clears only the in-chart slider; the Filter Panel filter survives.
      const canvas = viewer.querySelector('canvas[name="canvas"]')!;
      const cr = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: cr.left + cr.width / 2, clientY: cr.top + cr.height / 2}));
      await new Promise((r) => setTimeout(r, 500));
      const rv = Array.from(document.querySelectorAll('.d4-menu-item-label')).find((el) => el.textContent!.trim() === 'Reset View');
      if (rv)
        (rv.closest('.d4-menu-item') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 700));
      const afterResetView = df.filter.trueCount;

      // Re-narrow, then the Filter Panel "Reset filters" button clears everything.
      await dragHeight();
      const afterReDrag = df.filter.trueCount;
      const btn = document.querySelector('.d4-filter-group-header [name="icon-arrow-rotate-left"]') as HTMLElement | null;
      if (btn)
        btn.click();
      await new Promise((r) => setTimeout(r, 800));
      const afterPanelReset = df.filter.trueCount;
      return {fullCount, afterPanel, afterBoth, afterResetView, afterReDrag, afterPanelReset};
    });
    // Filter Panel filter takes effect, and the in-chart slider narrows further.
    expect(result.afterPanel).toBeLessThan(result.fullCount);
    expect(result.afterBoth).toBeLessThan(result.afterPanel);
    // Reset View resets ONLY the in-chart part: the Filter Panel filter remains.
    expect(result.afterResetView).toBe(result.afterPanel);
    // Re-narrowing works again, and the Filter Panel reset restores the full count.
    expect(result.afterReDrag).toBeLessThan(result.afterPanel);
    expect(result.afterPanelReset).toBe(result.fullCount);
  });

  v.finishSpec();
});
