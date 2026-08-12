/* ---
sub_features_covered: [legend.item.color-picker, legend.item.marker-picker, legend.mini-icon,
  legend.tooltip-mode, legend.corner.collapse, legend.size-limits, legend.section-resize]
--- */
// Regressions found in UI review of the frame-bound legend: the pickers, the mini icon /
// tooltip legend round trip, the corner collapse control, and the docked size ceiling.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

test('Legend chrome — pickers, mini icon, tooltip legend, collapse', async ({page}) => {
  test.setTimeout(900_000);

  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await v.openTable(page);
  await v.addLegendViewers(page, {column: 'Stereo Category', viewers: ['Scatter plot']});
  await v.setViewerProps(page, 'Scatter plot', [{set: {markersColumnName: 'Series'}, wait: 2500}]);
  await v.resizeViewer(page, 'Scatter plot', 900, 600);

  await softStep('A colour item offers the palette; a marker item offers only the shape', async () => {
    const pickers = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const out: Record<string, any> = {};
      for (const sid of ['main', 'extra']) {
        const item = x.root.querySelector(`[data-section="${sid}"] [data-item-key]`) as HTMLElement;
        if (!item) { out[sid] = null; continue; }
        item.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        const c = document.querySelector('[name="legend-icon-color-picker"]') as HTMLElement;
        const m = document.querySelector('[name="legend-icon-marker-picker"]') as HTMLElement;
        const ir = item.getBoundingClientRect();
        out[sid] = {
          color: !!c, marker: !!m,
          // distance from the item's left edge to the icon's right edge — a shape icon pushed
          // out past the colour slot ends up unreachable
          markerGap: m ? Math.round(ir.left - m.getBoundingClientRect().right) : null,
        };
        item.dispatchEvent(new MouseEvent('mouseleave', {bubbles: true}));
        document.querySelectorAll('[name="legend-icon-color-picker"],[name="legend-icon-marker-picker"]')
          .forEach((e) => e.remove());
      }
      return out;
    });

    expect(pickers.main.color, 'colour item has no palette icon').toBe(true);
    expect(pickers.extra.marker, 'marker item has no shape icon').toBe(true);
    expect(pickers.extra.color, 'marker items must not offer a colour picker').toBe(false);
    expect(pickers.extra.markerGap, 'the shape icon sits too far from its item').toBeLessThanOrEqual(24);
  });

  await softStep('The markers header is the selector alone, with no caption label', async () => {
    const text = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const h = x.root.querySelector('[name="legend-section-header"]') as HTMLElement;
      return h ? h.innerText.trim() : null;
    });
    expect(text).not.toContain('Markers');
  });

  await softStep('Right-clicking a colour item opens the picker', async () => {
    await page.locator('[name="viewer-Scatter-plot"] [data-section="main"] [data-item-key]').first()
      .click({button: 'right'});
    await page.waitForTimeout(1000);
    const opened = await page.evaluate(() =>
      !!document.querySelector('[class*="color-picker"], .d4-dialog'));
    expect(opened, 'right-click did not open the colour picker').toBe(true);
    await page.keyboard.press('Escape');
  });

  await softStep('Shrinking gives a mini icon whose tooltip holds the legend', async () => {
    await v.setViewerProps(page, 'Scatter plot', [{set: {legendVisibility: 'Auto'}, wait: 1200}]);
    await v.resizeViewer(page, 'Scatter plot', 200, 200);

    const mini = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const icon = x.root.querySelector('[name="mini-legend-icon"]') as HTMLElement;
      return {present: !!icon, display: icon?.style.display};
    });
    expect(mini.present).toBe(true);
    expect(mini.display).toBe('block');

    await page.locator('[name="viewer-Scatter-plot"] [name="mini-legend-icon"]').hover();
    await page.waitForTimeout(1500);
    const inTooltip = await page.evaluate(() => {
      const c = document.querySelector('.d4-tooltip, [name="tooltip"]') as HTMLElement;
      const l = c?.querySelector('[name="legend"]');
      return {found: !!l, mode: l?.getAttribute('data-legend-mode') ?? null};
    });
    expect(inTooltip.found, 'the tooltip legend never appears').toBe(true);
    expect(inTooltip.mode).toBe('tooltip');
  });

  await softStep('A corner legend can be collapsed and reopened', async () => {
    await v.resizeViewer(page, 'Scatter plot', 900, 600);
    await v.setViewerProps(page, 'Scatter plot', [{set: {markersColumnName: '', legendPosition: 'Auto'}, wait: 1500}]);

    const shown = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const root = x.root.querySelector('[name="legend"]') as HTMLElement;
      root?.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      const icon = x.root.querySelector('[name="icon-hide-corner-legend"]') as HTMLElement;
      return {mode: root?.getAttribute('data-legend-mode'), iconDisplay: icon?.style.display};
    });

    if (shown.mode !== 'corner')
      return; // the ladder docked it here; the collapse control only applies to corners

    expect(shown.iconDisplay, 'the collapse icon never becomes visible on hover').toBe('block');

    const collapsed = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      (x.root.querySelector('[name="icon-hide-corner-legend"]') as HTMLElement)?.click();
      await new Promise((r) => setTimeout(r, 1500));
      return (x.root.querySelector('[name="legend"]') as HTMLElement)?.getAttribute('data-legend-mode');
    });
    expect(collapsed).toBe('miniIcon');

    const restored = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      (x.root.querySelector('[name="mini-legend-icon"]') as HTMLElement)?.click();
      await new Promise((r) => setTimeout(r, 1500));
      return (x.root.querySelector('[name="legend"]') as HTMLElement)?.getAttribute('data-legend-mode');
    });
    expect(restored).toBe('corner');
  });

  await softStep('No uncaught errors', async () => {
    expect(errors, errors.join('\n')).toEqual([]);
  });

  await v.cleanupShell(page);
  v.finishSpec('Legend chrome failures');
});

test('Legend sizing — a docked legend never exceeds 40% and never scrolls sideways', async ({page}) => {
  test.setTimeout(900_000);

  await loginToDatagrok(page);
  await v.openTable(page);
  await v.addLegendViewers(page, {column: 'Series', viewers: ['Scatter plot']});

  const shapes: [number, number][] = [[420, 1100], [420, 700], [900, 500], [300, 900], [1200, 400]];
  for (const [w, h] of shapes) {
    await v.resizeViewer(page, 'Scatter plot', w, h);
    const info = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const root = x.root.querySelector('[name="legend"]') as HTMLElement;
      if (!root) return null;
      const list = root.querySelector('.d4-legend-list') as HTMLElement;
      return {
        slot: root.getAttribute('data-legend-slot'),
        mode: root.getAttribute('data-legend-mode'),
        fracW: root.offsetWidth / x.root.offsetWidth,
        fracH: root.offsetHeight / x.root.offsetHeight,
        rootScrollX: root.scrollWidth - root.clientWidth,
        listScrollX: list ? list.scrollWidth - list.clientWidth : 0,
      };
    });

    if (info == null || info.mode !== 'docked') continue;

    const label = `${w}x${h} slot=${info.slot}`;
    if (info.slot === 'left' || info.slot === 'right')
      expect(info.fracW, `${label}: legend wider than 40% of the viewer`).toBeLessThanOrEqual(0.42);
    else
      expect(info.fracH, `${label}: legend taller than 40% of the viewer`).toBeLessThanOrEqual(0.42);

    expect(info.rootScrollX, `${label}: legend root scrolls horizontally`).toBeLessThanOrEqual(0);
    expect(info.listScrollX, `${label}: legend list scrolls horizontally`).toBeLessThanOrEqual(0);
  }

  await v.cleanupShell(page);
  v.finishSpec('Legend sizing failures');
});
