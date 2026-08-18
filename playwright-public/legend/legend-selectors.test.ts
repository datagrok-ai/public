/* ---
sub_features_covered: [legend.markers.datetime-map, legend.line-chart.split-eq-marker,
  legend.mini-icon.corner-slot, legend.hover.cleared-on-column-switch]
--- */
// Legend column selectors and column-switching behavior. Pins: a datetime marker column
// gets a categorization (year/quarter/…) selector in the markers header, like the axis and
// color selectors; a line chart whose single split column doubles as the marker column
// keeps its legend (the split items carry the marker icons); a legend collapsed from an
// explicitly chosen corner shows its mini icon in that corner; and switching a legend
// column clears the mouse-over highlight predicate a hovered item left on the frame.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

async function legendInfo(page: any, vt: string) {
  return await page.evaluate((t: string) => {
    const tv = (window as any).grok.shell.tv;
    const x = tv.viewers.find((q: any) => q.type === t);
    const root = x?.root.querySelector('[name="legend"]') as HTMLElement;
    if (!root) return {noLegend: true};
    const sel = root.querySelector('.d4-legend-header-selector') as HTMLElement;
    const aggr = sel?.querySelector('select') as HTMLSelectElement;
    const icon = x.root.querySelector('[name="mini-legend-icon"]') as HTMLElement;
    return {mode: root.getAttribute('data-legend-mode'), slot: root.getAttribute('data-legend-slot'),
      items: Number(root.getAttribute('data-legend-items')),
      visible: root.getAttribute('data-legend-visible'),
      selector: !!sel, aggrVisible: aggr ? getComputedStyle(aggr).display !== 'none' : null,
      aggrValue: aggr ? aggr.value : null,
      iconStyle: icon ? {display: icon.style.display, top: icon.style.top, bottom: icon.style.bottom,
        left: icon.style.left, right: icon.style.right} : null};
  }, vt);
}

async function closeViewer(page: any, vt: string) {
  await page.evaluate((t: string) => {
    (window as any).grok.shell.tv.viewers.find((q: any) => q.type === t)?.close();
  }, vt);
  await page.waitForTimeout(500);
}

test('Legend selectors and column switching (SPGI)', async ({page}) => {
  test.setTimeout(900_000);
  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await v.openTable(page);

  await softStep('A datetime marker column gets a map selector in the markers header', async () => {
    await page.evaluate(() => {
      (window as any).grok.shell.tv.addViewer('Scatter plot',
        {colorColumnName: 'Stereo Category', markersColumnName: 'Competition assay Date'});
    });
    await page.waitForTimeout(3000);
    await v.resizeViewer(page, 'Scatter plot', 900, 600);
    await page.waitForTimeout(1500);
    const s = await legendInfo(page, 'Scatter plot');
    expect(s.selector, `no markers selector: ${JSON.stringify(s)}`).toBe(true);
    expect(s.aggrVisible, 'the datetime map selector is not shown').toBe(true);
    expect(s.aggrValue).toBe('year');

    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const sel = x.root.querySelector('[name="legend"] .d4-legend-header-selector select') as HTMLSelectElement;
      sel.value = 'quarter';
      sel.dispatchEvent(new Event('input', {bubbles: true}));
    });
    await page.waitForTimeout(2000);
    const mapped = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      return tv.viewers.find((q: any) => q.type === 'Scatter plot').props.markersMap;
    });
    expect(mapped, 'picking a map in the legend did not write the property').toBe('quarter');

    // a non-datetime marker column hides the map selector again
    await v.setViewerProps(page, 'Scatter plot', [{set: {markersColumnName: 'Series'}, wait: 2000}]);
    const plain = await legendInfo(page, 'Scatter plot');
    expect(plain.aggrVisible, 'the map selector must hide for a categorical column').toBe(false);
    await closeViewer(page, 'Scatter plot');
  });

  await softStep('A line chart split column that doubles as the marker column keeps its legend', async () => {
    await page.evaluate(() => {
      (window as any).grok.shell.tv.addViewer('Line chart',
        {yColumnNames: ['Average Mass'], splitColumnNames: ['Series'], showMarkers: 'Always'});
    });
    await page.waitForTimeout(2500);
    await v.resizeViewer(page, 'Line chart', 900, 600);
    await page.waitForTimeout(1500);
    const before = await legendInfo(page, 'Line chart');
    expect(before.items).toBeGreaterThan(0);

    await v.setViewerProps(page, 'Line chart', [{set: {markersColumnName: 'Series'}, wait: 2500}]);
    const after = await legendInfo(page, 'Line chart');
    expect(after.visible, `legend disappeared: ${JSON.stringify(after)}`).toBe('true');
    expect(after.items, 'the split items must survive split == marker').toBe(before.items);
    await closeViewer(page, 'Line chart');
  });

  await softStep('A legend collapsed from an explicit corner keeps its mini icon there', async () => {
    await page.evaluate(() => {
      (window as any).grok.shell.tv.addViewer('Scatter plot', {colorColumnName: 'Stereo Category'});
    });
    await page.waitForTimeout(2500);
    await v.resizeViewer(page, 'Scatter plot', 900, 600);
    await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: 'RightBottom'}, wait: 1500}]);
    const corner = await legendInfo(page, 'Scatter plot');
    expect(`${corner.mode}/${corner.slot}`).toBe('corner/rightBottom');

    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const root = x.root.querySelector('[name="legend"]') as HTMLElement;
      root?.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 300));
      (x.root.querySelector('[name="icon-hide-corner-legend"]') as HTMLElement)?.click();
      await new Promise((r) => setTimeout(r, 1200));
    });
    const mini = await legendInfo(page, 'Scatter plot');
    expect(`${mini.mode}/${mini.slot}`).toBe('miniIcon/rightBottom');
    // the icon anchors where the corner legend sat: from the bottom, above the x-axis strip
    const bottomPx = parseInt(mini.iconStyle.bottom || '');
    expect(bottomPx, `icon not in the bottom corner: ${JSON.stringify(mini.iconStyle)}`).toBeGreaterThanOrEqual(2);
    expect(bottomPx, 'icon must sit just above the x-axis, not mid-viewer').toBeLessThan(80);
    const rightPx = parseInt(mini.iconStyle.right || '');
    expect(rightPx, `icon must hug the right edge: ${JSON.stringify(mini.iconStyle)}`).toBeGreaterThanOrEqual(5);
    expect(rightPx).toBeLessThan(60);
    expect(mini.iconStyle.top).toBe('');
  });

  await softStep('Switching the legend column clears the hovered item\'s highlight predicate', async () => {
    const r = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const df = x.dataFrame;
      x.props.legendPosition = 'Auto';
      await new Promise((res) => setTimeout(res, 1000));
      let events = 0;
      const sub = df.onMouseOverRowGroupChanged.subscribe(() => events++);
      const item = x.root.querySelector('[name="legend"] [data-item-key]') as HTMLElement;
      item?.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      await new Promise((res) => setTimeout(res, 400));
      const afterHover = events;
      x.props.colorColumnName = 'Series';
      await new Promise((res) => setTimeout(res, 1500));
      sub.unsubscribe();
      return {item: !!item, afterHover, afterSwitch: events};
    });
    expect(r.item, 'no legend item to hover').toBe(true);
    expect(r.afterHover, 'hovering an item must set the highlight predicate').toBeGreaterThan(0);
    expect(r.afterSwitch, 'the column switch must clear the predicate (no event fired)')
      .toBeGreaterThan(r.afterHover);
  });

  expect(errors, `page errors:\n${errors.join('\n')}`).toEqual([]);
  await v.cleanupShell(page);
  v.finishSpec('Legend selector / column-switching failures');
});
