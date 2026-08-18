/* ---
sub_features_covered: [legend.item.hover, legend.item.click-filter, legend.item.color-picker,
  legend.sections.layout, legend.placement.stability, legend.item.text-shortening]
--- */
// Interaction regressions of the frame-bound legend. Everything here regressed against the legacy
// timer-driven path: hover highlighting, filtering from the extra section, picker lifetime,
// scroll position across a click, and placement stability while the viewport moves.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

test('Legend interaction — hover, click-to-filter, pickers', async ({page}) => {
  test.setTimeout(900_000);

  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await v.openTable(page);
  await v.addLegendViewers(page, {column: 'Stereo Category', viewers: ['Scatter plot']});
  await v.resizeViewer(page, 'Scatter plot', 900, 600);
  await v.waitForLegendIdle(page, 'Scatter plot');

  await softStep('Hovering a category highlights the matching rows', async () => {
    const hover = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const item = x.root.querySelector('[data-section="main"] [data-item-key]') as HTMLElement;
      item.dispatchEvent(new MouseEvent('mouseenter', {bubbles: false}));
      const during = !!tv.dataFrame.rows.mouseOverRowFunc;
      item.dispatchEvent(new MouseEvent('mouseleave', {bubbles: false}));
      return {during, after: !!tv.dataFrame.rows.mouseOverRowFunc};
    });

    expect(hover.during, 'hovering a legend item highlights no rows').toBe(true);
    expect(hover.after, 'the highlight outlives the hover').toBe(false);
  });

  await softStep('The picker icon appears left of the item and survives the trip to it', async () => {
    const picker = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const item = x.root.querySelector('[data-section="main"] [data-item-key]') as HTMLElement;
      item.dispatchEvent(new MouseEvent('mouseenter', {bubbles: false}));
      const icon = document.querySelector('[name="legend-icon-color-picker"]') as HTMLElement;
      if (!icon) return null;

      const ir = icon.getBoundingClientRect();
      const br = item.getBoundingClientRect();
      // leaving the item for the icon starts the hide timer; entering the icon must stop it
      item.dispatchEvent(new MouseEvent('mouseleave', {bubbles: false}));
      icon.dispatchEvent(new MouseEvent('mouseenter', {bubbles: false}));
      await new Promise((r) => setTimeout(r, 500));

      return {onLeft: ir.right <= br.left + 2, gap: Math.round(br.left - ir.right),
        aligned: Math.abs(ir.top - br.top) <= 6,
        stillThere: !!document.querySelector('[name="legend-icon-color-picker"]')};
    });

    expect(picker, 'no colour picker icon on hover').not.toBeNull();
    expect(picker!.onLeft, 'the picker icon is not to the left of its item').toBe(true);
    expect(picker!.gap, 'the picker icon is too far from its item').toBeLessThanOrEqual(24);
    expect(picker!.aligned, 'the picker icon is not on the item\'s row').toBe(true);
    expect(picker!.stillThere, 'the picker icon vanishes when the cursor reaches it').toBe(true);
  });

  await softStep('Moving between items leaves exactly one picker icon', async () => {
    const count = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const items = Array.from(
        x.root.querySelectorAll('[data-section="main"] [data-item-key]')) as HTMLElement[];
      for (const it of items.slice(0, 4)) {
        it.dispatchEvent(new MouseEvent('mouseenter', {bubbles: false}));
        it.dispatchEvent(new MouseEvent('mouseleave', {bubbles: false}));
      }
      items[0].dispatchEvent(new MouseEvent('mouseenter', {bubbles: false}));
      await new Promise((r) => setTimeout(r, 500));
      return document.querySelectorAll('[name="legend-icon-color-picker"]').length;
    });

    expect(count, 'picker icons accumulated across items').toBe(1);
  });

  await softStep('Right-clicking several items never stacks colour pickers', async () => {
    const items = page.locator('[name="viewer-Scatter-plot"] [data-section="main"] [data-item-key]');
    for (let i = 0; i < 3; i++) {
      await items.nth(i).click({button: 'right'});
      await page.waitForTimeout(500);
    }

    const open = await page.evaluate(() =>
      document.querySelectorAll('.d4-dialog, [class*="color-picker-host"]').length);
    expect(open, 'each right-click spawned another colour picker').toBeLessThanOrEqual(1);
    await page.keyboard.press('Escape');
  });

  await softStep('Clicking a marker item filters the viewer', async () => {
    await v.setViewerProps(page, 'Scatter plot', [{set: {markersColumnName: 'Series'}, wait: 2500}]);
    const filtered = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      // the legend filters the viewer, not the data frame — `Viewer.filter` is combinedFilter
      const before = x.filter.trueCount;
      const item = x.root.querySelector('[data-section="extra"] [data-item-key]') as HTMLElement;
      if (!item) return null;

      item.click();
      await new Promise((r) => setTimeout(r, 1500));
      return {before, after: x.filter.trueCount,
        selected: item.getAttribute('data-item-selected')};
    });

    expect(filtered, 'no markers section to click').not.toBeNull();
    expect(filtered!.selected, 'the clicked marker item is not marked selected').toBe('true');
    expect(filtered!.after, 'clicking a marker item filtered nothing').toBeLessThan(filtered!.before);

    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      (x.root.querySelector('[data-section="extra"] [data-item-key]') as HTMLElement)?.click();
      await new Promise((r) => setTimeout(r, 1000));
    });
  });

  await softStep('No uncaught errors', async () => {
    expect(errors, errors.join('\n')).toEqual([]);
  });

  await v.cleanupShell(page);
  v.finishSpec('Legend interaction failures');
});

test('Legend layout — sections, scroll, placement stability', async ({page}) => {
  test.setTimeout(900_000);

  await loginToDatagrok(page);
  await v.openTable(page);
  await v.addLegendViewers(page, {column: 'Stereo Category', viewers: ['Scatter plot']});
  await v.setViewerProps(page, 'Scatter plot', [
    {set: {markersColumnName: 'Series', legendPosition: 'Right'}, wait: 2500}]);
  await v.resizeViewer(page, 'Scatter plot', 900, 600);
  await v.waitForLegendIdle(page, 'Scatter plot');

  await softStep('Both sections are reachable without scrolling the legend', async () => {
    const box = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const root = x.root.querySelector('[name="legend"]') as HTMLElement;
      const main = root.querySelector('[data-section="main"]') as HTMLElement;
      const extra = root.querySelector('[data-section="extra"]') as HTMLElement;
      if (!main || !extra || extra.style.display === 'none') return null;

      const rr = root.getBoundingClientRect();
      return {
        rootScroll: root.scrollHeight - root.clientHeight,
        extraInside: extra.getBoundingClientRect().top < rr.bottom,
        extraNestedInMain: main.contains(extra),
        splitterShown: (root.querySelector('[name="legend-inner-splitter"]') as HTMLElement)
          ?.style.display,
      };
    });

    expect(box, 'the markers section never appeared').not.toBeNull();
    expect(box!.extraNestedInMain, 'the extra section is nested inside the main one').toBe(false);
    expect(box!.rootScroll, 'the legend root scrolls instead of its sections').toBeLessThanOrEqual(0);
    expect(box!.extraInside, 'the extra section is pushed below the legend').toBe(true);
    expect(box!.splitterShown, 'no splitter between two populated sections').toBe('block');
  });

  await softStep('Clicking a scrolled item does not jump the list back to the top', async () => {
    const scroll = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const list = x.root.querySelector('[data-section="main"] .d4-legend-list') as HTMLElement;
      if (list.scrollHeight - list.clientHeight < 20) return null;

      list.scrollTop = Math.floor((list.scrollHeight - list.clientHeight) / 2);
      const before = list.scrollTop;
      const items = Array.from(list.querySelectorAll('[data-item-key]')) as HTMLElement[];
      const visible = items.find((i) => i.offsetTop >= before);
      visible?.click();
      await new Promise((r) => setTimeout(r, 1500));
      const after = (x.root.querySelector('[data-section="main"] .d4-legend-list') as HTMLElement)
        .scrollTop;
      return {before, after};
    });

    if (scroll == null)
      return; // the list fits — nothing to scroll

    expect(Math.abs(scroll.after - scroll.before), 'the legend scrolled back on click')
      .toBeLessThanOrEqual(4);

    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      x.dataFrame.rows.filter(() => true);
      await new Promise((r) => setTimeout(r, 500));
    });
  });

  await softStep('A corner legend is bounded on both axes', async () => {
    await v.setViewerProps(page, 'Scatter plot',
      [{set: {legendPosition: 'RightTop'}, wait: 2000}]);
    const frac = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const root = x.root.querySelector('[name="legend"]') as HTMLElement;
      return {mode: root.getAttribute('data-legend-mode'),
        w: root.offsetWidth / x.root.offsetWidth, h: root.offsetHeight / x.root.offsetHeight};
    });

    if (frac.mode !== 'corner')
      return;

    expect(frac.w, 'the corner legend is wider than 40% of the viewer').toBeLessThanOrEqual(0.42);
    expect(frac.h, 'the corner legend is taller than 40% of the viewer').toBeLessThanOrEqual(0.42);
  });

  await softStep('The markers column selector opens outside the legend', async () => {
    const escaped = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const root = x.root.querySelector('[name="legend"]') as HTMLElement;
      const selector = root.querySelector('[name="legend-section-header"] .d4-column-selector-column')
        ?? root.querySelector('[name="legend-section-header"] div');
      (selector as HTMLElement)?.click();
      await new Promise((r) => setTimeout(r, 900));
      const popup = document.querySelector('.d4-column-selector-popup, .d4-column-selector-grid')
        ?? Array.from(document.querySelectorAll('div'))
          .find((d) => d.style.position === 'fixed' && d.style.zIndex === '4001');
      return popup == null ? null : {inLegend: root.contains(popup)};
    });

    if (escaped == null)
      return; // no selector popup on this build path

    expect(escaped.inLegend, 'the column selector popup opened inside the legend').toBe(false);
    await page.keyboard.press('Escape');
  });

  await v.cleanupShell(page);
  v.finishSpec('Legend layout failures');
});

test('Legend text — long categories are clipped, not laid out in full', async ({page}) => {
  test.setTimeout(900_000);

  await loginToDatagrok(page);
  await v.openTable(page);

  await page.evaluate(async () => {
    const df = (window as any).grok.shell.tv.dataFrame;
    df.columns.addNewString('LongCat').init((i: number) => 'category-'.repeat(40) + (i % 5));
    await new Promise((r) => setTimeout(r, 800));
  });

  await v.addLegendViewers(page, {column: 'LongCat', viewers: ['Scatter plot']});
  await v.resizeViewer(page, 'Scatter plot', 900, 600);
  await v.waitForLegendIdle(page, 'Scatter plot');

  const text = await page.evaluate(() => {
    const tv = (window as any).grok.shell.tv;
    const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
    const root = x.root.querySelector('[name="legend"]') as HTMLElement;
    const spans = Array.from(root.querySelectorAll('.d4-legend-item > span.d4-legend-value')) as HTMLElement[];
    const style = spans.length === 0 ? null : getComputedStyle(spans[0]);
    return {
      count: spans.length,
      longest: spans.reduce((m, s) => Math.max(m, s.textContent!.length), 0),
      ellipsis: style?.textOverflow,
      overflowsRoot: spans.some((s) =>
        s.getBoundingClientRect().right > root.getBoundingClientRect().right + 2),
    };
  });

  expect(text.count, 'the long-category legend rendered no items').toBeGreaterThan(0);
  expect(text.longest, 'a multi-hundred-character label reached the DOM verbatim')
    .toBeLessThanOrEqual(121);
  expect(text.ellipsis, 'labels are not ellipsised').toBe('ellipsis');
  expect(text.overflowsRoot, 'a label spills out of the legend').toBe(false);

  await v.cleanupShell(page);
});
