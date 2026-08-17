/* ---
sub_features_covered: [legend.corner.close-icon, legend.corner.scroll, legend.corner.resize,
  legend.corner.dual-sizing, legend.tooltip.lifecycle]
--- */
// Corner and tooltip legend regressions of the frame-bound legend:
// the collapse icon visible only in the right-top corner, a phantom scrollbar under the
// last item, a resize that never re-evaluated the corner (or got it stuck), a dual corner
// legend whose second section shrank below one row, and a tooltip legend that rendered
// empty after a click until the user scrolled.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

async function legendGeom(page: any, viewerType: string) {
  return await page.evaluate((vt: string) => {
    const tv = (window as any).grok.shell.tv;
    const x = tv.viewers.find((q: any) => q.type === vt);
    const root = x.root.querySelector('[name="legend"]') as HTMLElement;
    const r = root.getBoundingClientRect();
    const lists = Array.from(root.querySelectorAll('.d4-legend-list')).map((l: any) => ({
      name: l.getAttribute('name'), clientH: l.clientHeight, scrollH: l.scrollHeight,
      items: l.querySelectorAll('[name="legend-item"]').length,
    }));
    const sections = Array.from(root.querySelectorAll('.d4-legend-section')).map((s: any) => ({
      name: s.getAttribute('name'), h: (s as HTMLElement).offsetHeight,
      items: Number(s.getAttribute('data-section-items') ?? -1),
      headerH: (s.querySelector('.d4-legend-header') as HTMLElement)?.offsetHeight ?? 0,
    }));
    return {mode: root.getAttribute('data-legend-mode'), slot: root.getAttribute('data-legend-slot'),
      box: {x: Math.round(r.x), y: Math.round(r.y), w: Math.round(r.width), h: Math.round(r.height)},
      lists, sections};
  }, viewerType);
}

test('Corner legend — close icon, no phantom scroll, resize re-evaluation', async ({page}) => {
  test.setTimeout(900_000);

  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await v.openTable(page);
  await v.addLegendViewers(page, {column: 'Stereo Category', viewers: [
    {type: 'Scatter plot', column: 'Stereo Category'}]});
  await v.resizeViewer(page, 'Scatter plot', 900, 600);
  await v.waitForLegendIdle(page, 'Scatter plot');

  await softStep('The collapse icon sits at the legend box in every corner', async () => {
    for (const slot of ['RightTop', 'RightBottom', 'LeftTop', 'LeftBottom']) {
      await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: slot}, wait: 1200}]);
      await v.waitForLegendIdle(page, 'Scatter plot');
      const icon = await page.evaluate(() => {
        const legend = document.querySelector('[name="legend"]') as HTMLElement;
        // dispatched: in headless a row tooltip can sit between the pointer and the legend
        legend.dispatchEvent(new MouseEvent('mouseenter', {bubbles: false}));
        const el = document.querySelector('[name="icon-hide-corner-legend"]') as HTMLElement;
        const r = el.getBoundingClientRect();
        const lr = legend.getBoundingClientRect();
        return {display: getComputedStyle(el).display, x: r.x, y: r.y, w: r.width, h: r.height,
          legend: {x: lr.x, y: lr.y, w: lr.width, h: lr.height}};
      });
      expect(icon.display, `collapse icon hidden on hover in ${slot}`).toBe('block');
      expect(icon.w, `collapse icon has no size in ${slot}`).toBeGreaterThan(0);
      // at the legend's top-right corner, allowing the -12/-7 offsets it always had
      expect(Math.abs(icon.y - icon.legend.y), `collapse icon vertically detached in ${slot}`)
        .toBeLessThanOrEqual(4);
      expect(Math.abs((icon.x + icon.w) - (icon.legend.x + icon.legend.w)),
        `collapse icon horizontally detached in ${slot}`).toBeLessThanOrEqual(16);
      await page.evaluate(() => {
        const el = document.querySelector('[name="icon-hide-corner-legend"]') as HTMLElement;
        el.style.display = 'none';
      });
    }
  });

  await softStep('A corner legend has no scrollable overflow below its last item', async () => {
    const g = (await legendGeom(page, 'Scatter plot'))!;
    expect(g.mode).toBe('corner');
    const main = g.lists.find((l: any) => l.name === 'legend-list-main')!;
    expect(main.scrollH, 'phantom scroll: content overflows the corner box')
      .toBeLessThanOrEqual(main.clientH);
  });

  await softStep('An explicit corner tracks the viewer through resizes', async () => {
    await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: 'RightBottom'}, wait: 1200}]);
    for (const [w, h] of [[500, 400], [1100, 500]]) {
      await v.resizeViewer(page, 'Scatter plot', w, h);
      await page.waitForTimeout(800);
      await v.waitForLegendIdle(page, 'Scatter plot');
      const {viewer, geom} = await page.evaluate(() => {
        const tv = (window as any).grok.shell.tv;
        const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
        const r = x.root.getBoundingClientRect();
        const root = x.root.querySelector('[name="legend"]') as HTMLElement;
        const lr = root.getBoundingClientRect();
        return {viewer: {x: r.x, y: r.y, w: r.width, h: r.height},
          geom: {mode: root.getAttribute('data-legend-mode'),
            x: lr.x, y: lr.y, w: lr.width, h: lr.height}};
      });
      expect(geom.mode, `left corner mode at ${w}x${h}`).toBe('corner');
      expect(geom.x + geom.w, `legend outside the viewer at ${w}x${h}`)
        .toBeLessThanOrEqual(viewer.x + viewer.w + 1);
      expect(geom.y + geom.h, `legend below the viewer at ${w}x${h}`)
        .toBeLessThanOrEqual(viewer.y + viewer.h + 1);
    }
  });

  await softStep('Auto placement re-evaluates on resize and a free corner is kept', async () => {
    // zoom out so the points cluster in the middle and the corners stay free
    await v.setViewerProps(page, 'Scatter plot', [
      {set: {zoomAndFilter: 'no action'}, wait: 800},
      {set: {legendPosition: 'Auto'}, wait: 1500}]);
    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const vp = x.viewport;
      x.viewport = {x: vp.x - vp.width, y: vp.y - vp.height, width: vp.width * 3, height: vp.height * 3};
      await new Promise((r) => setTimeout(r, 1000));
    });
    // a content change re-runs the ladder with the corners now free
    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      x.props.colorColumnName = 'Series';
    });
    await page.waitForTimeout(1500);
    await v.waitForLegendIdle(page, 'Scatter plot');
    const entered = (await legendGeom(page, 'Scatter plot'))!;
    expect(entered.mode, 'a free corner was not taken').toBe('corner');

    const modes: string[] = [entered.mode!];
    // sizes where the corner box (sized with the DOM-true item font + chrome) stays clear
    // of the clustered points — at 700x500 the honest estimate genuinely reaches them
    for (const [w, h] of [[900, 600], [1000, 600], [1100, 650]]) {
      await v.resizeViewer(page, 'Scatter plot', w, h);
      await page.waitForTimeout(800);
      await v.waitForLegendIdle(page, 'Scatter plot');
      modes.push((await legendGeom(page, 'Scatter plot'))!.mode!);
    }
    // the corners stay free through the sweep, so the settled answer stays corner —
    // this is what "resize re-evaluates placement without the legend getting stuck" means
    expect(modes.every((m) => m === 'corner'),
      `placement flapped or got stuck across resizes: ${modes.join(',')}`).toBe(true);
  });

  expect(errors, `page errors: ${errors.join(' | ')}`).toHaveLength(0);
  v.finishSpec('Corner legend failures');
  await v.cleanupShell(page);
});

test('Dual legend — the second section keeps at least its header and a row', async ({page}) => {
  test.setTimeout(900_000);

  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await v.openTable(page);
  await v.addLegendViewers(page, {column: 'Stereo Category', viewers: [
    {type: 'Scatter plot', column: 'Stereo Category'}]});
  await page.evaluate(() => {
    const tv = (window as any).grok.shell.tv;
    const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
    x.props.markersColumnName = 'Series';
  });
  await v.resizeViewer(page, 'Scatter plot', 900, 600);
  await page.waitForTimeout(1200);

  const checkSections = (g: any, where: string) => {
    const main = g.sections.find((s: any) => s.name === 'legend-section-main')!;
    const extra = g.sections.find((s: any) => s.name === 'legend-section-extra')!;
    const mainList = g.lists.find((l: any) => l.name === 'legend-list-main')!;
    expect(extra.items, `${where}: markers section lost its items`).toBeGreaterThan(0);
    expect(extra.headerH, `${where}: markers selector collapsed`).toBeGreaterThan(15);
    // header + at least one full row of items
    expect(extra.h, `${where}: markers section below header + one row`)
      .toBeGreaterThanOrEqual(extra.headerH + 18);
    expect(main.h, `${where}: main section squeezed out`).toBeGreaterThan(18);
    expect(mainList.scrollH, `${where}: main list clipped`)
      .toBeLessThanOrEqual(mainList.clientH + 1);
  };

  await softStep('Corner: both sections visible, main not clipped', async () => {
    await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: 'RightTop'}, wait: 1500}]);
    await v.waitForLegendIdle(page, 'Scatter plot');
    const g = (await legendGeom(page, 'Scatter plot'))!;
    expect(g.mode).toBe('corner');
    checkSections(g, 'corner');
    // the corner sizes to its entries when the 40% cap allows: header + all rows, or the cap
    const wanted = Math.min(600 * 0.4,
      g.sections.reduce((a: number, s: any) => a + s.h, 0));
    expect(g.box.h, 'corner box does not fit its entries').toBeGreaterThanOrEqual(wanted - 12);
  });

  await softStep('Top dock: both sections visible, main not clipped', async () => {
    await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: 'Top'}, wait: 1500}]);
    await v.waitForLegendIdle(page, 'Scatter plot');
    checkSections((await legendGeom(page, 'Scatter plot'))!, 'top');
  });

  expect(errors, `page errors: ${errors.join(' | ')}`).toHaveLength(0);
  v.finishSpec('Dual legend sizing failures');
  await v.cleanupShell(page);
});

test('Tooltip legend — items survive the click/hide/rehover cycle', async ({page}) => {
  test.setTimeout(900_000);

  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await v.openTable(page);
  await v.addLegendViewers(page, {column: 'Stereo Category', viewers: [
    {type: 'Scatter plot', column: 'Stereo Category'}]});
  // the context panel overlays the viewer's top-right corner where the mini icon lives,
  // and the spec harness force-hides `.d4-tooltip` — this test is about that tooltip
  await page.evaluate(() => {
    const w = (window as any).grok.shell.windows;
    w.showContextPanel = false;
    w.showProperties = false;
    w.showHelp = false;
    for (const s of Array.from(document.querySelectorAll('style')))
      if (s.textContent?.includes('.d4-tooltip'))
        s.textContent = s.textContent.replace(/\.d4-tooltip\s*{[^}]*}/g, '');
  });
  await v.resizeViewer(page, 'Scatter plot', 600, 500);
  await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: 'RightTop'}, wait: 1500}]);
  await v.waitForLegendIdle(page, 'Scatter plot');

  const collapsed = await page.evaluate(() => {
    const legend = document.querySelector('[name="legend"]') as HTMLElement;
    legend.dispatchEvent(new MouseEvent('mouseenter', {bubbles: false}));
    const icon = document.querySelector('[name="icon-hide-corner-legend"]') as HTMLElement;
    if (!icon || getComputedStyle(icon).display === 'none') return false;
    icon.click();
    return true;
  });
  expect(collapsed, 'the collapse icon never appeared').toBe(true);
  await page.waitForTimeout(1000);

  const hoverMini = async () => {
    const mini = await page.evaluate(() => {
      const el = document.querySelector('[name="mini-legend-icon"]') as HTMLElement;
      const r = el.getBoundingClientRect();
      return {x: r.x + r.width / 2, y: r.y + r.height / 2};
    });
    await page.mouse.move(mini.x - 40, mini.y + 40);
    await page.mouse.move(mini.x, mini.y, {steps: 5});
    await page.waitForTimeout(1200);
  };

  const state = async () => await page.evaluate(() => {
    const legend = document.querySelector('[name="legend"]') as HTMLElement;
    const tip = legend?.closest('.d4-tooltip') as HTMLElement;
    return {
      inTooltip: tip != null, tipShown: tip ? getComputedStyle(tip).display !== 'none' : false,
      visibleItems: legend ? Array.from(legend.querySelectorAll('[name="legend-item"]'))
        .filter((e: any) => e.getBoundingClientRect().height > 0).length : -1,
      rows: (window as any).grok.shell.tv.viewers
        .find((q: any) => q.type === 'Scatter plot').filter.trueCount,
    };
  });

  await softStep('Hovering the mini icon shows a populated legend', async () => {
    await hoverMini();
    const s = await state();
    expect(s.inTooltip, 'legend never moved into the tooltip').toBe(true);
    expect(s.tipShown, 'tooltip stayed hidden').toBe(true);
    expect(s.visibleItems, 'tooltip legend rendered empty').toBeGreaterThan(0);
  });

  await softStep('Clicking an item filters and the legend stays populated', async () => {
    const before = (await state()).rows;
    // dispatched: physical hit-testing of legend items is covered by legend-interaction;
    // this test is about the legend surviving the rebuild the click causes
    await page.evaluate(() => {
      const legend = document.querySelector('[name="legend"]') as HTMLElement;
      const it = legend.querySelector('[name="legend-item"]') as HTMLElement;
      const r = it.getBoundingClientRect();
      it.dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true,
        clientX: r.x + r.width / 2, clientY: r.y + r.height / 2}));
    });
    await page.waitForTimeout(1500);
    const s = await state();
    expect(s.rows, 'clicking the tooltip legend item did not filter').toBeLessThan(before);
    expect(s.visibleItems, 'the legend blanked after the click').toBeGreaterThan(0);
  });

  await softStep('After hiding and re-hovering the legend is populated without scrolling', async () => {
    await page.mouse.move(200, 500);
    await page.waitForTimeout(1500);
    await hoverMini();
    const s = await state();
    expect(s.tipShown, 'tooltip did not come back').toBe(true);
    expect(s.visibleItems, 'tooltip legend empty until scrolled').toBeGreaterThan(0);
  });

  expect(errors, `page errors: ${errors.join(' | ')}`).toHaveLength(0);
  v.finishSpec('Tooltip legend failures');
  await v.cleanupShell(page);
});

test('Markers selector — picking a column works, including after a tooltip cycle', async ({page}) => {
  test.setTimeout(900_000);

  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await v.openTable(page);
  await v.addLegendViewers(page, {column: 'Stereo Category', viewers: [
    {type: 'Scatter plot', column: 'Stereo Category'}]});
  await page.evaluate(() => {
    const w = (window as any).grok.shell.windows;
    w.showContextPanel = false; w.showProperties = false; w.showHelp = false;
    for (const s of Array.from(document.querySelectorAll('style')))
      if (s.textContent?.includes('.d4-tooltip'))
        s.textContent = s.textContent.replace(/\.d4-tooltip\s*{[^}]*}/g, '');
    const tv = (window as any).grok.shell.tv;
    const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
    x.props.markersColumnName = 'Series';
  });
  await page.waitForTimeout(1500);
  await v.resizeViewer(page, 'Scatter plot', 600, 500);
  await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: 'Right'}, wait: 1500}]);
  await v.waitForLegendIdle(page, 'Scatter plot');

  const pickViaMouse = async (where: string) => {
    const sel = await page.evaluate(() => {
      const legend = document.querySelector('[name="legend"]') as HTMLElement;
      const el = legend?.querySelector('.d4-legend-header-selector') as HTMLElement;
      if (!el) return null;
      const r = el.getBoundingClientRect();
      return {x: r.x + Math.min(20, r.width / 2), y: r.y + r.height / 2};
    });
    expect(sel, `${where}: no markers selector in the legend`).not.toBeNull();
    await page.mouse.click(sel!.x, sel!.y);
    await page.waitForTimeout(800);
    const popup = await page.evaluate(() => {
      const backdrop = document.querySelector('.d4-column-selector-backdrop') as HTMLElement;
      if (!backdrop) return null;
      const r = backdrop.getBoundingClientRect();
      return {x: r.x, y: r.y, h: r.height};
    });
    expect(popup, `${where}: the column popup did not open`).not.toBeNull();
    // the third row of the column grid — any row works, it just has to differ from 'Series'
    await page.mouse.click(popup!.x + 80, popup!.y + 70);
    await page.waitForTimeout(1500);
    const markers = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      return tv.viewers.find((q: any) => q.type === 'Scatter plot').props.markersColumnName;
    });
    expect(markers, `${where}: picking a column in the selector had no effect`)
      .not.toBe('Series');
    await page.keyboard.press('Escape');
    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      tv.viewers.find((q: any) => q.type === 'Scatter plot').props.markersColumnName = 'Series';
    });
    await page.waitForTimeout(1200);
  };

  await softStep('The selector applies a column pick while docked', async () =>
    await pickViaMouse('docked'));

  await softStep('The selector still works after a tooltip open/close cycle', async () => {
    await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: 'RightTop'}, wait: 1500}]);
    await v.waitForLegendIdle(page, 'Scatter plot');
    await page.evaluate(() => {
      const legend = document.querySelector('[name="legend"]') as HTMLElement;
      legend.dispatchEvent(new MouseEvent('mouseenter', {bubbles: false}));
      (document.querySelector('[name="icon-hide-corner-legend"]') as HTMLElement)?.click();
    });
    await page.waitForTimeout(800);
    const mini = await page.evaluate(() => {
      const el = document.querySelector('[name="mini-legend-icon"]') as HTMLElement;
      const r = el.getBoundingClientRect();
      return {x: r.x + r.width / 2, y: r.y + r.height / 2};
    });
    await page.mouse.move(mini.x - 40, mini.y + 40);
    await page.mouse.move(mini.x, mini.y, {steps: 5});
    await page.waitForTimeout(1200);
    await page.mouse.move(200, 500);
    await page.waitForTimeout(1500);
    await page.evaluate(() => (document.querySelector('[name="mini-legend-icon"]') as HTMLElement).click());
    await page.waitForTimeout(1200);
    await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: 'Right'}, wait: 1500}]);
    await v.waitForLegendIdle(page, 'Scatter plot');
    await pickViaMouse('after tooltip');
  });

  await softStep('Content added mid-hover resizes the tooltip legend', async () => {
    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      tv.viewers.find((q: any) => q.type === 'Scatter plot').props.markersColumnName = '';
    });
    await page.waitForTimeout(1200);
    await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: 'RightTop'}, wait: 1500}]);
    await v.waitForLegendIdle(page, 'Scatter plot');
    await page.evaluate(() => {
      const legend = document.querySelector('[name="legend"]') as HTMLElement;
      legend.dispatchEvent(new MouseEvent('mouseenter', {bubbles: false}));
      (document.querySelector('[name="icon-hide-corner-legend"]') as HTMLElement)?.click();
    });
    await page.waitForTimeout(800);
    const mini = await page.evaluate(() => {
      const el = document.querySelector('[name="mini-legend-icon"]') as HTMLElement;
      const r = el.getBoundingClientRect();
      return {x: r.x + r.width / 2, y: r.y + r.height / 2};
    });
    await page.mouse.move(mini.x - 40, mini.y + 40);
    await page.mouse.move(mini.x, mini.y, {steps: 5});
    await page.waitForTimeout(1200);

    const before = await page.evaluate(() => {
      const legend = document.querySelector('[name="legend"]') as HTMLElement;
      const r = legend.getBoundingClientRect();
      return {w: r.width, h: r.height};
    });
    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      tv.viewers.find((q: any) => q.type === 'Scatter plot').props.markersColumnName = 'Series';
    });
    await page.waitForTimeout(1500);
    const after = await page.evaluate(() => {
      const legend = document.querySelector('[name="legend"]') as HTMLElement;
      const r = legend.getBoundingClientRect();
      const sections = Array.from(legend.querySelectorAll('.d4-legend-section'))
        .map((s: any) => ({items: Number(s.getAttribute('data-section-items') ?? 0), h: s.offsetHeight}));
      return {w: r.width, h: r.height, sections};
    });
    expect(after.h, 'the tooltip legend kept its old height for twice the content')
      .toBeGreaterThan(before.h + 40);
    for (const s of after.sections.filter((q: any) => q.items > 0))
      expect(s.h, 'a tooltip legend section collapsed').toBeGreaterThan(30);
  });

  expect(errors, `page errors: ${errors.join(' | ')}`).toHaveLength(0);
  v.finishSpec('Markers selector failures');
  await v.cleanupShell(page);
});

test('Resize drag — no oscillation, and GROK-19041 equation avoidance', async ({page}) => {
  test.setTimeout(900_000);

  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await v.openTable(page);
  await v.addLegendViewers(page, {column: 'Stereo Category', viewers: [
    {type: 'Scatter plot', column: 'Stereo Category'}]});

  await softStep('A continuous drag out of the mini icon tracks the size without oscillating', async () => {
    await v.setViewerProps(page, 'Scatter plot',
      [{set: {legendVisibility: 'Auto', legendPosition: 'Auto'}, wait: 1200}]);
    await v.resizeViewer(page, 'Scatter plot', 240, 400);
    await v.waitForLegendIdle(page, 'Scatter plot');
    const small = await v.getLegendState(page, 'Scatter plot');
    expect(small?.mode).toBe('miniIcon');

    await page.evaluate(() => {
      const w = window as any;
      w.__modes = [];
      const legend = document.querySelector('[name="legend"]') as HTMLElement;
      new MutationObserver(() => {
        const m = `${legend.getAttribute('data-legend-mode')}/${legend.getAttribute('data-legend-slot')}`;
        if (w.__modes[w.__modes.length - 1] !== m) w.__modes.push(m);
      }).observe(legend, {attributes: true, attributeFilter: ['data-legend-mode', 'data-legend-slot']});
    });
    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const host = x.root.parentElement ?? x.root;
      for (let w = 270; w <= 900; w += 30) {
        host.style.width = `${w}px`;
        window.dispatchEvent(new Event('resize'));
        await new Promise((r) => setTimeout(r, 60));
      }
    });
    await page.waitForTimeout(1500);
    const modes = await page.evaluate(() => (window as any).__modes);
    // placement is decided per frame now, so a monotone drag may legitimately pass through
    // the states its sizes settle to (mini icon exit while still narrow, then the vertical
    // flip) — but each hysteresis band flips at most once and no state is ever revisited
    expect(modes.length, `too many moves mid-drag: ${modes.join(' -> ')}`)
      .toBeLessThanOrEqual(3);
    expect(new Set(modes).size, `the legend oscillated mid-drag: ${modes.join(' -> ')}`)
      .toBe(modes.length);
    const finalState = await v.getLegendState(page, 'Scatter plot');
    expect(finalState?.mode).toBe('docked');
    expect(finalState?.slot, 'a wide viewer settles on the right').toBe('right');
  });

  await softStep('GROK-19041 — the auto legend leaves the regression equation corner', async () => {
    await v.resizeViewer(page, 'Scatter plot', 900, 600);
    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      x.props.markersColumnName = 'Series';
    });
    await v.setViewerProps(page, 'Scatter plot', [
      {set: {zoomAndFilter: 'no action'}, wait: 500},
      {set: {colorColumnName: 'Average Mass'}, wait: 1000},
      {set: {legendPosition: 'Auto', legendVisibility: 'Auto'}, wait: 1000}]);
    // point cloud to the bottom-right so the top-left corner is the ladder's free pick
    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const vp = x.viewport;
      x.viewport = {x: vp.x - vp.width * 1.5, y: vp.y, width: vp.width * 2.5, height: vp.height * 2.5};
      await new Promise((r) => setTimeout(r, 1200));
    });
    await v.resizeViewer(page, 'Scatter plot', 901, 600);
    await v.waitForLegendIdle(page, 'Scatter plot');
    const noEq = await v.getLegendState(page, 'Scatter plot');
    expect(noEq?.mode, 'setup: the legend should sit in a corner').toBe('corner');
    expect(noEq?.slot, 'setup: without the equation the top-left corner is the pick')
      .toBe('leftTop');

    await v.setViewerProps(page, 'Scatter plot', [
      {set: {showRegressionLine: true}, wait: 500},
      {set: {showRegressionLineEquation: true}, wait: 2000}]);
    await v.waitForLegendIdle(page, 'Scatter plot');
    const withEq = await v.getLegendState(page, 'Scatter plot');
    expect(withEq?.slot, 'the legend sits on top of the regression equation')
      .not.toBe('leftTop');
  });

  expect(errors, `page errors: ${errors.join(' | ')}`).toHaveLength(0);
  v.finishSpec('Resize / overlay avoidance failures');
  await v.cleanupShell(page);
});

test('Collapse persistence, corner clearance, tooltip transparency (demog)', async ({page}) => {
  test.setTimeout(900_000);

  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await page.evaluate(async () => {
    const grok = (window as any).grok;
    grok.shell.addTableView(grok.data.demo.demog(3000));
    await new Promise((r) => setTimeout(r, 2000));
    grok.shell.tv.scatterPlot({x: 'weight', y: 'height', color: 'race', markers: 'sex'});
    await new Promise((r) => setTimeout(r, 2000));
    const w = grok.shell.windows;
    w.showContextPanel = false; w.showProperties = false; w.showHelp = false;
    for (const s of Array.from(document.querySelectorAll('style')))
      if (s.textContent?.includes('.d4-tooltip'))
        s.textContent = s.textContent.replace(/\.d4-tooltip\s*{[^}]*}/g, '');
  });
  await page.waitForTimeout(1500);
  await v.resizeViewer(page, 'Scatter plot', 700, 620);
  await v.setViewerProps(page, 'Scatter plot', [
    {set: {legendVisibility: 'Auto', legendPosition: 'RightBottom'}, wait: 1500}]);
  await v.waitForLegendIdle(page, 'Scatter plot');

  await softStep('A corner legend hugs its corner, clear of the selector strip only', async () => {
    const g = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const vr = x.root.getBoundingClientRect();
      const root = x.root.querySelector('[name="legend"]') as HTMLElement;
      const lr = root.getBoundingClientRect();
      return {mode: root.getAttribute('data-legend-mode'),
        rightGap: Math.round(vr.right - lr.right), bottomGap: Math.round(vr.bottom - lr.bottom)};
    });
    expect(g.mode).toBe('corner');
    // 15px base + the selector strip's thickness; a frozen long dimension pushed this past 100
    expect(g.rightGap, 'corner legend floats away from the right edge').toBeLessThanOrEqual(60);
    expect(g.bottomGap, 'corner legend floats away from the bottom').toBeLessThanOrEqual(70);
  });

  await softStep('A collapsed legend stays collapsed through any auto fallback', async () => {
    await page.evaluate(() => {
      const legend = document.querySelector('[name="legend"]') as HTMLElement;
      legend.dispatchEvent(new MouseEvent('mouseenter', {bubbles: false}));
      (document.querySelector('[name="icon-hide-corner-legend"]') as HTMLElement)?.click();
    });
    await page.waitForTimeout(1000);
    await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: 'Auto'}, wait: 1200}]);
    for (const [w, h] of [[380, 300], [950, 700], [700, 620]]) {
      await v.resizeViewer(page, 'Scatter plot', w, h);
      await v.waitForLegendIdle(page, 'Scatter plot');
      const mode = await page.evaluate(() =>
        document.querySelector('[name="legend"]')!.getAttribute('data-legend-mode'));
      expect(mode, `collapsed legend popped open at ${w}x${h}`).toBe('miniIcon');
    }
    await page.evaluate(() => (document.querySelector('[name="mini-legend-icon"]') as HTMLElement).click());
    await page.waitForTimeout(1200);
    const mode = await page.evaluate(() =>
      document.querySelector('[name="legend"]')!.getAttribute('data-legend-mode'));
    expect(mode, 'the mini icon click did not restore the legend').not.toBe('miniIcon');
  });

  await softStep('The tooltip legend has no background of its own', async () => {
    await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: 'RightTop'}, wait: 1200}]);
    await v.waitForLegendIdle(page, 'Scatter plot');
    await page.evaluate(() => {
      const legend = document.querySelector('[name="legend"]') as HTMLElement;
      legend.dispatchEvent(new MouseEvent('mouseenter', {bubbles: false}));
      (document.querySelector('[name="icon-hide-corner-legend"]') as HTMLElement)?.click();
    });
    await page.waitForTimeout(1000);
    const mini = await page.evaluate(() => {
      const el = document.querySelector('[name="mini-legend-icon"]') as HTMLElement;
      const r = el.getBoundingClientRect();
      return {x: r.x + r.width / 2, y: r.y + r.height / 2};
    });
    await page.mouse.move(mini.x - 40, mini.y + 40);
    await page.mouse.move(mini.x, mini.y, {steps: 5});
    await page.waitForTimeout(1500);
    const s = await page.evaluate(() => {
      const legend = document.querySelector('[name="legend"]') as HTMLElement;
      const tip = legend.closest('.d4-tooltip') as HTMLElement;
      return {inTip: tip != null,
        bg: getComputedStyle(legend).backgroundColor,
        corner: legend.classList.contains('d4-corner-legend'),
        shadow: getComputedStyle(legend).boxShadow};
    });
    expect(s.inTip, 'the tooltip legend never appeared').toBe(true);
    expect(s.bg, 'the tooltip legend paints its own background').toBe('rgba(0, 0, 0, 0)');
    expect(s.corner, 'corner styling leaked into the tooltip').toBe(false);
    expect(s.shadow, 'corner shadow leaked into the tooltip').toBe('none');
  });

  expect(errors, `page errors: ${errors.join(' | ')}`).toHaveLength(0);
  v.finishSpec('Collapse / clearance / tooltip-style failures');
  await v.cleanupShell(page);
});

test('Docked-to-corner offset, dual stacking, corner/tooltip backgrounds (SPGI)', async ({page}) => {
  test.setTimeout(900_000);
  // Docked/corner transition regressions: a docked→corner transition measured the flex part with the docked
  // legend still inside it, freezing a corner offset the size of the old legend; the dual
  // docked legend stretched its first section over the whole panel; the header between the
  // opaque white lists stayed see-through in a corner; and the tooltip legend's lists kept
  // .d4-root's white over the yellow tooltip.

  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await v.openTable(page);
  await v.addLegendViewers(page, {column: 'Primary Series Name', viewers: [
    {type: 'Scatter plot', column: 'Primary Series Name'}]});
  await v.setViewerProps(page, 'Scatter plot', [{set: {zoomAndFilter: 'no action'}, wait: 500}]);

  const legendGaps = async (viewerType: string) => await page.evaluate((vt: string) => {
    const tv = (window as any).grok.shell.tv;
    const x = tv.viewers.find((q: any) => q.type === vt);
    const vr = x.root.getBoundingClientRect();
    const root = x.root.querySelector('[name="legend"]') as HTMLElement;
    const r = root.getBoundingClientRect();
    return {mode: root.getAttribute('data-legend-mode'), slot: root.getAttribute('data-legend-slot'),
      right: Math.round(vr.right - r.right), top: Math.round(r.top - vr.top)};
  }, viewerType);

  await softStep('A docked legend leaving for a corner does not offset by its own width', async () => {
    // wide → docked right; pan the data toward the bottom-left; the resize then frees the
    // right-top corner and the ladder moves the legend there while it is still docked in
    // the DOM — the transition that used to freeze right≈180
    await v.resizeViewer(page, 'Scatter plot', 1100, 500);
    await v.waitForLegendIdle(page, 'Scatter plot');
    expect((await legendGaps('Scatter plot')).mode).toBe('docked');
    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const vp = x.viewport;
      x.viewport = {x: vp.x, y: vp.y - vp.height * 1.5, width: vp.width * 2.5, height: vp.height * 2.5};
    });
    await page.waitForTimeout(600);
    await v.resizeViewer(page, 'Scatter plot', 660, 800);
    await page.waitForTimeout(1000);
    await v.waitForLegendIdle(page, 'Scatter plot');
    const g = await legendGaps('Scatter plot');
    expect(g.mode, 'free right-top corner not taken').toBe('corner');
    expect(g.slot).toBe('rightTop');
    expect(g.right, 'corner offset by the departed docked legend').toBeLessThanOrEqual(60);
  });

  await softStep('Line chart: manual RightTop from docked lands in the corner', async () => {
    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      tv.addViewer('Line chart', {xColumnName: 'Competition assay Date',
        yColumnNames: ['CAST Idea ID'],
        splitColumnNames: ['Primary Series Name', 'Scaffold Names', 'Series', 'Stereo Category']});
    });
    await page.waitForTimeout(2000);
    await v.resizeViewer(page, 'Line chart', 1400, 480);
    await page.waitForTimeout(800);
    await v.setViewerProps(page, 'Line chart', [{set: {legendPosition: 'RightTop'}, wait: 1500}]);
    await v.waitForLegendIdle(page, 'Line chart');
    const g = await legendGaps('Line chart');
    expect(g.mode).toBe('corner');
    expect(g.right, 'corner offset by the departed docked legend').toBeLessThanOrEqual(60);
    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      tv.viewers.find((q: any) => q.type === 'Line chart')?.close();
    });
    await page.waitForTimeout(500);
  });

  await v.setViewerProps(page, 'Scatter plot', [
    {set: {colorColumnName: 'Scaffold Names'}, wait: 800},
    {set: {markersColumnName: 'Primary Series Name'}, wait: 1000}]);

  await softStep('Docked dual legend: the marker section follows the color section', async () => {
    await v.resizeViewer(page, 'Scatter plot', 900, 800);
    await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: 'Right'}, wait: 1200}]);
    await v.waitForLegendIdle(page, 'Scatter plot');
    const d = await page.evaluate(() => {
      const root = document.querySelector('[name="legend"]') as HTMLElement;
      const out: any = {};
      for (const s of Array.from(root.querySelectorAll('.d4-legend-section'))) {
        const r = (s as HTMLElement).getBoundingClientRect();
        out[s.getAttribute('data-section') as string] = {
          y: Math.round(r.y), h: Math.round(r.height),
          items: Number(s.getAttribute('data-section-items'))};
      }
      return out;
    });
    // main holds 9 items ≈ 9 rows; the extra section must start right after it, not at the
    // bottom of the 800px panel
    expect(d.extra.y - (d.main.y + d.main.h), 'gap between the sections').toBeLessThanOrEqual(20);
    expect(d.main.h, 'first section stretched over the panel').toBeLessThanOrEqual(d.main.items * 25 + 40);
  });

  await softStep('Corner dual legend: the selector header is not a see-through sliver', async () => {
    await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: 'RightTop'}, wait: 1200}]);
    await v.waitForLegendIdle(page, 'Scatter plot');
    const bg = await page.evaluate(() => {
      const root = document.querySelector('[name="legend"]') as HTMLElement;
      const header = root.querySelector('[name="legend-section-header"]') as HTMLElement;
      return getComputedStyle(header).backgroundColor;
    });
    expect(bg).toBe('rgb(255, 255, 255)');
  });

  await softStep('Corner→side: the lists rebuild for the new size without scrolling', async () => {
    await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: 'Right'}, wait: 1200}]);
    await v.waitForLegendIdle(page, 'Scatter plot');
    const cover = await page.evaluate(() => {
      const root = document.querySelector('[name="legend"]') as HTMLElement;
      const list = root.querySelector('[name="legend-list-main"]') as HTMLElement;
      const items = Array.from(list.querySelectorAll('[name="legend-item"]'))
        .map((el) => el.getBoundingClientRect());
      const lr = list.getBoundingClientRect();
      const lowest = items.reduce((m, r) => Math.max(m, r.bottom), 0);
      return {rendered: items.length, listBottom: Math.round(lr.bottom),
        contentBottom: Math.round(Math.min(lowest, lr.bottom)),
        listTop: Math.round(lr.top), listH: Math.round(lr.height),
        scrollH: list.scrollHeight};
    });
    expect(cover.rendered, 'no items rendered after the transition').toBeGreaterThan(0);
    // rendered items must reach the fold (or the end of the content if it is shorter)
    const visibleContent = Math.min(cover.scrollH, cover.listH);
    expect(cover.contentBottom - cover.listTop, 'stale window: legend empty below the fold')
      .toBeGreaterThanOrEqual(visibleContent - 25);
  });

  await softStep('Tooltip legend: the lists let the tooltip background through', async () => {
    await page.evaluate(() => {
      const w = (window as any).grok.shell.windows;
      w.showContextPanel = false; w.showProperties = false; w.showHelp = false;
      for (const s of Array.from(document.querySelectorAll('style')))
        if (s.textContent?.includes('.d4-tooltip'))
          s.textContent = s.textContent.replace(/\.d4-tooltip\s*{[^}]*}/g, '');
    });
    await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: 'RightTop'}, wait: 1000}]);
    await v.resizeViewer(page, 'Scatter plot', 700, 500);
    await page.waitForTimeout(1000);
    await page.evaluate(() => {
      const legend = document.querySelector('[name="legend"]') as HTMLElement;
      legend.dispatchEvent(new MouseEvent('mouseenter', {bubbles: false}));
      (document.querySelector('[name="icon-hide-corner-legend"]') as HTMLElement).click();
    });
    await page.waitForTimeout(1500);
    const mini = await page.evaluate(() => {
      const el = document.querySelector('[name="mini-legend-icon"]') as HTMLElement;
      const r = el.getBoundingClientRect();
      return {x: r.x + r.width / 2, y: r.y + r.height / 2};
    });
    await page.mouse.move(mini.x - 40, mini.y + 40);
    await page.mouse.move(mini.x, mini.y, {steps: 5});
    await page.waitForTimeout(1500);
    const bgs = await page.evaluate(() => {
      const legend = document.querySelector('[name="legend"]') as HTMLElement;
      const tip = legend?.closest('.d4-tooltip') as HTMLElement;
      if (!tip) return null;
      return Array.from(legend.querySelectorAll('.d4-legend-list'))
        .map((l) => getComputedStyle(l).backgroundColor);
    });
    expect(bgs, 'tooltip legend did not appear').not.toBeNull();
    for (const bg of bgs!)
      expect(bg, 'opaque list over the tooltip').toBe('rgba(0, 0, 0, 0)');
    await page.mouse.move(200, 900);
    await page.waitForTimeout(800);
  });

  expect(errors, `page errors: ${errors.join(' | ')}`).toHaveLength(0);
});
