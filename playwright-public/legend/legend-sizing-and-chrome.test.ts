/* ---
sub_features_covered: [legend.placement.markers-no-flap, legend.placement.bar-stack-no-flap,
  legend.sizing.whole-rows, legend.sizing.item-font-and-chrome, legend.header.column-priority,
  legend.mini-icon.axis-offsets, legend.tooltip.outside-viewer, legend.sections.splitter-scroll,
  legend.tooltip.no-stale-reclaim]
--- */
// Round-16 pins. A scatter with clustered data, a color column and a datetime marker
// column never flaps docked/top ↔ corner/rightBottom while widening (the undocked chart
// box is component-local: removing a docked strip grows it without moving its origin);
// a bounded custom-renderer legend is cut to whole rows and scrolls; a narrow legend
// header keeps the column label readable while the map selector shrinks, and the header's
// width is a floor on the legend estimate; the mini icon sits beside the axes in every
// corner; a collapsed corner's tooltip legend opens outside the viewer when there is
// room; and the inner splitter resizes continuously (whole-row snapping is for auto
// sizing only) while a squeezed section scrolls instead of clipping.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

async function closeViewer(page: any, vt: string) {
  await page.evaluate((t: string) => {
    (window as any).grok.shell.tv.viewers.find((q: any) => q.type === t)?.close();
  }, vt);
  await page.waitForTimeout(500);
}

test('Legend sizing, chrome and marker-placement stability (SPGI)', async ({page}) => {
  test.setTimeout(900_000);
  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await v.openTable(page);

  await softStep('Widening a marker+color scatter never flaps docked/top vs corner', async () => {
    await page.evaluate(() => {
      Object.defineProperty(window, 'devicePixelRatio', {get: () => 0.8999999, configurable: true});
      const w = (window as any).grok.shell.windows;
      w.showContextPanel = false; w.showProperties = false; w.showHelp = false;
      (window as any).grok.shell.tv.addViewer('Scatter plot',
        {xColumnName: 'Chemical Space X', yColumnName: 'Chemical Space Y',
         colorColumnName: 'Primary Series Name', markersColumnName: 'Competition assay Date', markersMap: 'year'});
    });
    await page.waitForTimeout(3000);
    await v.resizeViewer(page, 'Scatter plot', 575, 832);
    await page.waitForTimeout(1500);
    await page.evaluate(() => {
      const w = window as any;
      w.__moves = [];
      w.__sub = w.grok.events.onEvent('d4-legend-placement-changed').subscribe((e: any) => {
        const a = e.args ?? e;
        w.__moves.push(`${a.fromMode}/${a.fromSlot}->${a.toMode}/${a.toSlot}`);
      });
    });
    for (let w = 577; w <= 745; w += 2) {
      await page.evaluate(({width}) => {
        const tv = (window as any).grok.shell.tv;
        const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
        const host = x.root.parentElement ?? x.root;
        host.style.width = `${width}px`;
        host.style.height = `832px`;
        window.dispatchEvent(new Event('resize'));
      }, {width: w});
      await page.waitForTimeout(200);
    }
    await page.waitForTimeout(800);
    const moves = await page.evaluate(() => {
      const w = window as any; w.__sub?.unsubscribe(); return w.__moves;
    });
    expect(moves.length, `placement flapped while widening:\n${moves.join('\n')}`).toBeLessThanOrEqual(1);
    const slots = moves.map((s: string) => s.split('->')[1]);
    expect(slots.filter((s: string, i: number) => i >= 2 && slots[i - 2] === s),
      'returned to a slot it had just left').toEqual([]);
    await closeViewer(page, 'Scatter plot');
  });

  await softStep('Widening a stacked bar chart never flaps between top and a corner', async () => {
    await page.evaluate(() => {
      (window as any).grok.shell.tv.addViewer('Bar chart',
        {splitColumnName: 'Stereo Category', stackColumnName: 'Scaffold Names'});
    });
    await page.waitForTimeout(2500);
    await v.resizeViewer(page, 'Bar chart', 600, 700);
    await page.waitForTimeout(1500);
    await page.evaluate(() => {
      const w = window as any;
      w.__moves = [];
      w.__sub = w.grok.events.onEvent('d4-legend-placement-changed').subscribe((e: any) => {
        const a = e.args ?? e;
        w.__moves.push(`${a.fromMode}/${a.fromSlot}->${a.toMode}/${a.toSlot}`);
      });
    });
    for (let w = 604; w <= 800; w += 4) {
      await page.evaluate(({width}) => {
        const tv = (window as any).grok.shell.tv;
        const x = tv.viewers.find((q: any) => q.type === 'Bar chart');
        const host = x.root.parentElement ?? x.root;
        host.style.width = `${width}px`;
        host.style.height = `700px`;
        window.dispatchEvent(new Event('resize'));
      }, {width: w});
      await page.waitForTimeout(200);
    }
    await page.waitForTimeout(800);
    const moves = await page.evaluate(() => {
      const w = window as any; w.__sub?.unsubscribe(); return w.__moves;
    });
    expect(moves.length, `bar chart flapped while widening:\n${moves.join('\n')}`).toBeLessThanOrEqual(2);
    const slots = moves.map((s: string) => s.split('->')[1]);
    expect(slots.filter((s: string, i: number) => i >= 2 && slots[i - 2] === s),
      'returned to a slot it had just left').toEqual([]);
    await closeViewer(page, 'Bar chart');
  });

  await softStep('A long-label top legend sizes to its real columns — no phantom second column', async () => {
    await page.evaluate(() => {
      (window as any).grok.shell.tv.addViewer('Scatter plot', {colorColumnName: 'Scaffold Names'});
    });
    await page.waitForTimeout(2500);
    const read = () => page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const root = x.root.querySelector('[name="legend"]') as HTMLElement;
      const list = root.querySelector('.d4-legend-list') as HTMLElement;
      const items = Array.from(root.querySelectorAll('[name="legend-item"]')) as HTMLElement[];
      const lr = list.getBoundingClientRect();
      const visible = items.filter((i) => {
        const r = i.getBoundingClientRect();
        return r.bottom <= lr.bottom + 2 && r.top >= lr.top - 2;
      }).length;
      return {items: Number(root.getAttribute('data-legend-items')),
        state: `${root.getAttribute('data-legend-mode')}/${root.getAttribute('data-legend-slot')}`,
        visible, scroll: list.scrollHeight > list.clientHeight + 2};
    });
    for (const w of [606, 640]) {
      await v.resizeViewer(page, 'Scatter plot', w, 735);
      await page.waitForTimeout(1500);
      const s = await read();
      if (s.state !== 'docked/top') continue;
      expect(s.visible, `at ${w}px only ${s.visible}/${s.items} items visible: ${JSON.stringify(s)}`)
        .toBe(s.items);
      expect(s.scroll, `at ${w}px the estimate fit columns the DOM could not`).toBe(false);
    }
    await closeViewer(page, 'Scatter plot');
  });

  await softStep('A narrow legend header keeps the column label; the map selector shrinks', async () => {
    await page.evaluate(() => {
      (window as any).grok.shell.tv.addViewer('Scatter plot',
        {xColumnName: 'Chemical Space X', yColumnName: 'Chemical Space Y',
         colorColumnName: 'Stereo Category', markersColumnName: 'Competition assay Date', markersMap: 'year',
         legendPosition: 'RightBottom'});
    });
    await page.waitForTimeout(3000);
    await v.resizeViewer(page, 'Scatter plot', 800, 700);
    await page.waitForTimeout(1500);
    const h = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const root = x.root.querySelector('[name="legend"]') as HTMLElement;
      const sel = root.querySelector('.d4-legend-header-selector') as HTMLElement;
      const aggr = sel?.querySelector('.d4-legend-header-aggr') as HTMLElement;
      const colLabel = sel?.querySelector('.d4-column-selector-column') as HTMLElement;
      return {legendW: Math.round(root.getBoundingClientRect().width),
        selW: sel ? Math.round(sel.getBoundingClientRect().width) : null,
        aggrW: aggr ? Math.round(aggr.getBoundingClientRect().width) : null,
        colW: colLabel ? Math.round(colLabel.getBoundingClientRect().width) : null};
    });
    expect(h.colW, `column label crushed: ${JSON.stringify(h)}`).toBeGreaterThanOrEqual(55);
    expect(h.aggrW, 'map selector must shrink but stay usable').toBeGreaterThanOrEqual(20);
    expect(h.selW, 'the header must fit inside the legend').toBeLessThanOrEqual(h.legendW);
    expect(h.legendW, 'the header width is a floor on the legend estimate').toBeGreaterThanOrEqual(170);
    await closeViewer(page, 'Scatter plot');
  });

  await softStep('A bounded molecule legend on top shows whole rows and scrolls', async () => {
    await page.evaluate(() => {
      (window as any).grok.shell.tv.addViewer('Scatter plot', {colorColumnName: 'R2', legendPosition: 'Top'});
    });
    await page.waitForTimeout(4000);
    await v.resizeViewer(page, 'Scatter plot', 900, 700);
    await page.waitForTimeout(2500);
    const m = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const root = x.root.querySelector('[name="legend"]') as HTMLElement;
      const list = root.querySelector('.d4-legend-list') as HTMLElement;
      const items = Array.from(root.querySelectorAll('[name="legend-item"]')) as HTMLElement[];
      const lr = list.getBoundingClientRect();
      const cut = items.filter((i) => {
        const r = i.getBoundingClientRect();
        return r.top < lr.bottom - 2 && r.bottom > lr.bottom + 2;
      }).length;
      return {itemH: items.length ? Math.round(items[0].getBoundingClientRect().height) : 0,
        listH: Math.round(lr.height), cutItems: cut,
        canScroll: list.scrollHeight > list.clientHeight + 2};
    });
    expect(m.itemH, 'no molecule items rendered').toBeGreaterThan(40);
    expect(m.cutItems, 'a bounded legend must cut to whole rows').toBe(0);
    expect(m.listH % m.itemH, 'the list height must be a whole number of rows').toBe(0);
    expect(m.canScroll, 'the cut rows must be reachable by scrolling').toBe(true);
    await closeViewer(page, 'Scatter plot');
  });

  await softStep('The leftBottom mini icon sits beside the axes, not on the viewer edge', async () => {
    await page.evaluate(() => {
      (window as any).grok.shell.tv.addViewer('Scatter plot', {colorColumnName: 'Stereo Category'});
    });
    await page.waitForTimeout(2500);
    await v.resizeViewer(page, 'Scatter plot', 900, 600);
    await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: 'LeftBottom'}, wait: 1500}]);
    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const root = x.root.querySelector('[name="legend"]') as HTMLElement;
      root?.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 300));
      (x.root.querySelector('[name="icon-hide-corner-legend"]') as HTMLElement)?.click();
      await new Promise((r) => setTimeout(r, 1200));
    });
    const icon = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const el = x.root.querySelector('[name="mini-legend-icon"]') as HTMLElement;
      return {left: parseInt(el.style.left || ''), bottom: parseInt(el.style.bottom || '')};
    });
    expect(icon.left, 'icon must clear the y-axis').toBeGreaterThanOrEqual(20);
    expect(icon.left, 'icon must hug the data area, not drift inward').toBeLessThan(90);
    expect(icon.bottom, 'icon must clear the x-axis').toBeGreaterThanOrEqual(10);
    expect(icon.bottom).toBeLessThan(90);
    await closeViewer(page, 'Scatter plot');
  });

  await softStep('A collapsed corner\'s tooltip legend opens outside the viewer', async () => {
    await page.evaluate(() => {
      for (const s of Array.from(document.querySelectorAll('style')))
        if (s.textContent?.includes('.d4-tooltip'))
          s.textContent = s.textContent.replace(/\.d4-tooltip\s*{[^}]*}/g, '');
      (window as any).grok.shell.tv.addViewer('Scatter plot', {colorColumnName: 'Primary Series Name'});
    });
    await page.waitForTimeout(2500);
    await v.resizeViewer(page, 'Scatter plot', 700, 500);
    await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: 'RightBottom'}, wait: 1200}]);
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
    await page.mouse.move(mini.x - 40, mini.y - 40);
    await page.mouse.move(mini.x, mini.y, {steps: 5});
    await page.waitForTimeout(1500);
    const t = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const root = document.querySelector('[name="legend"]') as HTMLElement;
      const host = root?.closest('.d4-tooltip') as HTMLElement;
      if (!host) return null;
      const tr = host.getBoundingClientRect();
      const rr = x.root.getBoundingClientRect();
      const ov = Math.max(0, Math.min(tr.right, rr.right) - Math.max(tr.left, rr.left))
        * Math.max(0, Math.min(tr.bottom, rr.bottom) - Math.max(tr.top, rr.top));
      return {overlapFrac: ov / (tr.width * tr.height)};
    });
    expect(t, 'tooltip legend did not appear').not.toBeNull();
    expect(t!.overlapFrac, 'the tooltip legend must open outside the viewer when there is room')
      .toBeLessThanOrEqual(0.05);
    await page.mouse.move(200, 900);
    await page.waitForTimeout(600);
    await closeViewer(page, 'Scatter plot');
  });

  await softStep('The inner splitter resizes continuously; a squeezed section scrolls, never clips', async () => {
    await page.evaluate(() => {
      (window as any).grok.shell.tv.addViewer('Scatter plot',
        {colorColumnName: 'Primary Scaffold Name', markersColumnName: 'Series'});
    });
    await page.waitForTimeout(3000);
    await v.resizeViewer(page, 'Scatter plot', 1000, 620);
    await page.waitForTimeout(1500);
    const r = await page.evaluate(async () => {
      const root = document.querySelector('[name="legend"]') as HTMLElement;
      const splitter = root.querySelector('[name="legend-inner-splitter"]') as HTMLElement;
      if (!splitter || getComputedStyle(splitter).display === 'none') return {noSplitter: true};
      const state = (name: string) => {
        const sec = root.querySelector(`[name="legend-section-${name}"]`) as HTMLElement;
        const list = sec?.querySelector('.d4-legend-list') as HTMLElement;
        return {listH: Math.round(list.getBoundingClientRect().height),
          clipped: list.getBoundingClientRect().bottom > sec.getBoundingClientRect().bottom + 2};
      };
      const before = {main: state('main'), extra: state('extra')};
      const sr = splitter.getBoundingClientRect();
      const ev = (type: string, py: number, target: any) => target.dispatchEvent(new MouseEvent(type,
        {bubbles: true, cancelable: true, view: window, clientX: sr.x + 5, clientY: py}));
      ev('mousedown', sr.y + 2, splitter);
      for (let i = 1; i <= 10; i++) ev('mousemove', sr.y + 2 + 200 * i / 10, document);
      ev('mouseup', sr.y + 202, document);
      await new Promise((res) => setTimeout(res, 800));
      return {before, after: {main: state('main'), extra: state('extra')}};
    });
    expect((r as any).noSplitter, 'no dual legend with an inner splitter').toBeFalsy();
    const {before, after} = r as any;
    // continuous, not row-stepped: the drag lands exactly at the clamp (extra keeps 20px),
    // a value a whole-row scheme cannot produce
    const expected = Math.min(200, before.extra.listH - 20);
    expect(Math.abs(after.main.listH - before.main.listH - expected),
      `main went ${before.main.listH} -> ${after.main.listH}, expected +${expected}`).toBeLessThanOrEqual(2);
    expect(after.main.clipped, 'the main list must never overflow its section').toBe(false);
    expect(after.extra.clipped, 'the extra list must never overflow its section').toBe(false);

    // shrink the viewer: the fixed main share re-clamps and the extra list scrolls
    await v.resizeViewer(page, 'Scatter plot', 1000, 380);
    await page.waitForTimeout(1200);
    const s = await page.evaluate(() => {
      const root = document.querySelector('[name="legend"]') as HTMLElement;
      const list = root.querySelector('[name="legend-section-extra"] .d4-legend-list') as HTMLElement;
      const sec = root.querySelector('[name="legend-section-extra"]') as HTMLElement;
      const before = list.scrollTop;
      list.scrollTop = 60;
      return {canScroll: list.scrollHeight > list.clientHeight + 2,
        scrolled: list.scrollTop !== before || list.scrollTop > 0,
        clientH: list.clientHeight,
        clipped: list.getBoundingClientRect().bottom > sec.getBoundingClientRect().bottom + 2};
    });
    expect(s.clientH, 'the squeezed extra section must keep at least a row').toBeGreaterThan(15);
    expect(s.canScroll, 'the squeezed extra list must scroll').toBe(true);
    expect(s.scrolled, 'scrollTop must actually move').toBe(true);
    expect(s.clipped, 'no clipping without a scrollbar').toBe(false);
  });

  await softStep('Expanding from the mini icon never de-docks the legend on a later tooltip', async () => {
    // corner -> collapse -> Auto -> hover-expand leaves a tooltip-visible flag; the next
    // row tooltip must not reclaim the now-docked legend out of the flex layout (that
    // painted the canvas at full width under a shrunken CSS box: blur + broken hit tests)
    await closeViewer(page, 'Scatter plot');
    await page.evaluate(() => {
      (window as any).grok.shell.tv.addViewer('Scatter plot', {colorColumnName: 'Series'});
    });
    await page.waitForTimeout(2500);
    await v.resizeViewer(page, 'Scatter plot', 760, 600);
    await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: 'RightTop'}, wait: 1500}]);
    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const root = x.root.querySelector('[name="legend"]') as HTMLElement;
      root?.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 300));
      (x.root.querySelector('[name="icon-hide-corner-legend"]') as HTMLElement)?.click();
      await new Promise((r) => setTimeout(r, 1200));
    });
    await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: 'Auto'}, wait: 1500}]);

    const mini = await page.evaluate(() => {
      const el = document.querySelector('[name="mini-legend-icon"]') as HTMLElement;
      const r = el.getBoundingClientRect();
      return {x: r.x + r.width / 2, y: r.y + r.height / 2};
    });
    await page.mouse.move(mini.x - 40, mini.y + 40);
    await page.mouse.move(mini.x, mini.y, {steps: 5});
    await page.waitForTimeout(1500);
    await page.mouse.click(mini.x, mini.y);
    await page.waitForTimeout(2000);

    // hover data points (row tooltips fire), then force a repaint
    const c = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const r = (x.root.querySelector('canvas') as HTMLElement).getBoundingClientRect();
      return {x: r.x, y: r.y, w: r.width, h: r.height};
    });
    for (let i = 0; i < 10; i++)
      await page.mouse.move(c.x + c.w * (0.25 + 0.05 * i), c.y + c.h * 0.5);
    await page.waitForTimeout(600);
    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      x.props.showRegressionLine = true;
      await new Promise((r) => setTimeout(r, 800));
      x.props.showRegressionLine = false;
    });
    await page.waitForTimeout(1500);

    const s = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const root = x.root.querySelector('[name="legend"]') as HTMLElement;
      const dpr = window.devicePixelRatio;
      const canvases = Array.from(x.root.querySelectorAll('canvas')) as HTMLCanvasElement[];
      return {mode: root.getAttribute('data-legend-mode'),
        docked: x.root.contains(root) && root.parentElement !== x.root,
        canvases: canvases.map((cv) => {
          const r = cv.getBoundingClientRect();
          return {css: `${Math.round(r.width)}x${Math.round(r.height)}`, attr: `${cv.width}x${cv.height}`,
            ok: Math.abs(cv.width - r.width * dpr) <= 2 && Math.abs(cv.height - r.height * dpr) <= 2};
        })};
    });
    expect(s.mode, `legend must stay expanded: ${JSON.stringify(s)}`).toBe('docked');
    expect(s.docked, 'the docked legend was pulled out of the flex layout').toBe(true);
    for (const cv of s.canvases)
      expect(cv.ok, `canvas backing store out of sync: ${JSON.stringify(s.canvases)}`).toBe(true);
  });

  await softStep('A bar chart corner legend lands beside the bars on the very first placement', async () => {
    await closeViewer(page, 'Scatter plot');
    await page.evaluate(() => {
      (window as any).grok.shell.tv.addViewer('Bar chart',
        {splitColumnName: 'Stereo Category', stackColumnName: 'Series'});
    });
    await page.waitForTimeout(2500);
    await v.resizeViewer(page, 'Bar chart', 620, 760);
    await page.waitForTimeout(1500);
    const read = () => page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Bar chart');
      const root = x.root.querySelector('[name="legend"]') as HTMLElement;
      return {state: `${root.getAttribute('data-legend-mode')}/${root.getAttribute('data-legend-slot')}`,
        bottom: parseFloat(root.style.bottom || '0'), right: parseFloat(root.style.right || '0')};
    });
    const first = await read();
    await v.resizeViewer(page, 'Bar chart', 622, 760);
    await page.waitForTimeout(1500);
    const settled = await read();
    expect(first.state).toBe(settled.state);
    // the first anchor must already be the settled one — a stale-geometry anchor put the
    // legend on the axis until a resize re-resolved
    expect(Math.abs(first.bottom - settled.bottom),
      `first placement ${JSON.stringify(first)} vs settled ${JSON.stringify(settled)}`).toBeLessThanOrEqual(2);
    expect(Math.abs(first.right - settled.right)).toBeLessThanOrEqual(2);

    // the mini icon anchors where the legend sat (below the axis), not at the viewer edge
    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Bar chart');
      const root = x.root.querySelector('[name="legend"]') as HTMLElement;
      root?.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 300));
      (x.root.querySelector('[name="icon-hide-corner-legend"]') as HTMLElement)?.click();
      await new Promise((r) => setTimeout(r, 1200));
    });
    const icon = await page.evaluate(() => {
      const el = document.querySelector('[name="mini-legend-icon"]') as HTMLElement;
      return {bottom: parseFloat(el.style.bottom || '0'), right: parseFloat(el.style.right || '0')};
    });
    expect(Math.abs(icon.bottom - settled.bottom),
      `icon ${JSON.stringify(icon)} must anchor at the legend's corner ${JSON.stringify(settled)}`)
      .toBeLessThanOrEqual(20);
    expect(Math.abs(icon.right - settled.right)).toBeLessThanOrEqual(20);
    await closeViewer(page, 'Bar chart');
  });

  await softStep('A pie chart mini icon tucks beside the circle, not at the viewer corner', async () => {
    await page.evaluate(() => {
      (window as any).grok.shell.tv.addViewer('Pie chart', {categoryColumnName: 'Primary Series Name'});
    });
    await page.waitForTimeout(2500);
    await v.resizeViewer(page, 'Pie chart', 620, 760);
    await page.waitForTimeout(1500);
    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Pie chart');
      const root = x.root.querySelector('[name="legend"]') as HTMLElement;
      root?.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 300));
      (x.root.querySelector('[name="icon-hide-corner-legend"]') as HTMLElement)?.click();
      await new Promise((r) => setTimeout(r, 1200));
    });
    const icon = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Pie chart');
      const el = x.root.querySelector('[name="mini-legend-icon"]') as HTMLElement;
      const r = el.getBoundingClientRect();
      const vr = x.root.getBoundingClientRect();
      return {fromRight: Math.round(vr.right - r.right), fromTop: Math.round(r.top - vr.top),
        display: el.style.display};
    });
    expect(icon.display).toBe('block');
    expect(icon.fromRight, `icon must be pulled toward the circle: ${JSON.stringify(icon)}`)
      .toBeGreaterThan(40);
    expect(icon.fromTop, 'icon must sit beside the circle, not at the viewer top').toBeGreaterThan(40);
    await closeViewer(page, 'Pie chart');
  });

  expect(errors, `page errors:\n${errors.join('\n')}`).toEqual([]);
  await v.cleanupShell(page);
  v.finishSpec('Legend sizing / chrome failures');
});
