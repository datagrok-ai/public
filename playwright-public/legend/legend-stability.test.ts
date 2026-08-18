/* ---
sub_features_covered: [legend.placement.stability, legend.placement.resize-sweep,
  legend.list.virtualization, legend.list.scroll]
--- */
// The legend must be boring: a slow resize through every aspect ratio should move it a
// handful of times at most, and panning or zooming should not move it at all. Transitions
// are counted from `d4-legend-placement-changed`, which fires only when the placement
// actually changes.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';

async function recordPlacements(page: any): Promise<void> {
  await page.evaluate(() => {
    (window as any).__legendMoves = [];
    if ((window as any).__legendSub)
      return;

    (window as any).__legendSub = (window as any).grok.events
      .onEvent('d4-legend-placement-changed')
      .subscribe((a: any) => (window as any).__legendMoves.push(
        `${a.args.width}x${a.args.height} ${a.args.fromMode}/${a.args.fromSlot}` +
        `->${a.args.toMode}/${a.args.toSlot}`));
  });
}

const moves = (page: any): Promise<string[]> =>
  page.evaluate(() => (window as any).__legendMoves);

test('Legend stability — a slow resize does not make the legend hop', async ({page}) => {
  test.setTimeout(900_000);

  await loginToDatagrok(page);
  await v.openTable(page);
  await v.addLegendViewers(page, {column: 'Stereo Category',
    viewers: ['Scatter plot', 'Box plot']});
  await page.waitForTimeout(2000);
  await recordPlacements(page);

  const heights: number[] = [];
  for (let h = 300; h <= 1200; h += 75) heights.push(h);

  for (const viewer of ['Scatter plot', 'Box plot']) {
    await v.setViewerProps(page, viewer, [{set: {legendPosition: 'Auto'}, wait: 1000}]);

    for (const direction of ['widening', 'narrowing']) {
      await page.evaluate(() => { (window as any).__legendMoves = []; });
      for (const h of (direction === 'widening' ? heights : heights.slice().reverse())) {
        await v.resizeViewer(page, viewer, 900, h);
        await page.waitForTimeout(260);
      }

      const seen = await moves(page);
      await softStep(`${viewer}: ${direction} from 3:1 to 1:3 settles instead of oscillating`,
        async () => {
          // one slot change per genuine threshold crossing is fine; a flip-flop is not
          expect(seen.length, `${viewer} ${direction}:\n${seen.join('\n')}`)
            .toBeLessThanOrEqual(3);

          const slots = seen.map((s) => s.split('->')[1]);
          const repeats = slots.filter((s, i) => i >= 2 && slots[i - 2] === s);
          expect(repeats, `${viewer} ${direction} returned to a slot it had just left:\n` +
            seen.join('\n')).toEqual([]);
        });
    }
  }

  await v.cleanupShell(page);
  v.finishSpec('Legend resize stability failures');
});

test('Legend stability — panning and zooming never move the legend', async ({page}) => {
  test.setTimeout(900_000);

  await loginToDatagrok(page);
  await v.openTable(page);
  await v.addLegendViewers(page, {column: 'Stereo Category', viewers: ['Scatter plot']});
  await v.resizeViewer(page, 'Scatter plot', 900, 600);
  await page.waitForTimeout(2000);
  await recordPlacements(page);

  const result = await page.evaluate(async () => {
    const tv = (window as any).grok.shell.tv;
    const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
    const root = x.root.querySelector('[name="legend"]') as HTMLElement;
    const r0 = root.getBoundingClientRect();
    const vp0 = JSON.stringify(x.viewport);
    const before = {items: root.getAttribute('data-legend-items'),
      mode: root.getAttribute('data-legend-mode'), slot: root.getAttribute('data-legend-slot')};

    // `viewport` is the world rect — shrinking it around its centre is a real zoom in
    for (let i = 0; i < 4; i++) {
      const vp = x.viewport;
      x.viewport = {x: vp.x + vp.width * 0.1, y: vp.y + vp.height * 0.1,
        width: vp.width * 0.8, height: vp.height * 0.8};
      await new Promise((r) => setTimeout(r, 700));
    }

    for (let i = 0; i < 3; i++) {
      const vp = x.viewport;
      x.viewport = {x: vp.x + vp.width * 0.1, y: vp.y, width: vp.width, height: vp.height};
      await new Promise((r) => setTimeout(r, 700));
    }

    const r1 = root.getBoundingClientRect();
    return {before, viewportChanged: vp0 !== JSON.stringify(x.viewport),
      after: {items: root.getAttribute('data-legend-items'),
        mode: root.getAttribute('data-legend-mode'), slot: root.getAttribute('data-legend-slot')},
      dx: Math.round(r1.left - r0.left), dy: Math.round(r1.top - r0.top),
      dw: Math.round(r1.width - r0.width), dh: Math.round(r1.height - r0.height)};
  });

  // without this the rest of the test proves nothing
  expect(result.viewportChanged, 'the zoom never happened').toBe(true);
  expect(await moves(page), 'the viewport change re-resolved the placement').toEqual([]);
  expect(result.after, 'zoom/pan changed what the legend shows').toEqual(result.before);
  expect([result.dx, result.dy, result.dw, result.dh],
    'the legend moved or resized during zoom/pan').toEqual([0, 0, 0, 0]);

  await v.cleanupShell(page);
});

test('Legend list — hundreds of categories stay virtualised and scrollable', async ({page}) => {
  test.setTimeout(900_000);

  await loginToDatagrok(page);
  await v.openTable(page);
  await page.evaluate(async () => {
    const df = (window as any).grok.shell.tv.dataFrame;
    df.columns.addNewString('ManyCats').init((i: number) => 'cat-' + (i % 300));
    await new Promise((r) => setTimeout(r, 800));
  });

  await v.addLegendViewers(page, {column: 'ManyCats', viewers: ['Scatter plot']});
  await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: 'Right'}, wait: 2000}]);
  await v.resizeViewer(page, 'Scatter plot', 900, 600);
  await page.waitForTimeout(2500);

  const list = await page.evaluate(async () => {
    const tv = (window as any).grok.shell.tv;
    const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
    const root = x.root.querySelector('[name="legend"]') as HTMLElement;
    const el = root.querySelector('.d4-legend-list') as HTMLElement;
    const rendered = el.querySelectorAll('[data-item-key]').length;
    const firstKey = el.querySelector('[data-item-key]')?.getAttribute('data-item-key');

    el.scrollTop = Math.floor((el.scrollHeight - el.clientHeight) / 2);
    el.dispatchEvent(new Event('scroll'));
    await new Promise((r) => setTimeout(r, 700));

    return {items: parseInt(root.getAttribute('data-legend-items')!, 10),
      rendered, clientHeight: el.clientHeight, scrollHeight: el.scrollHeight,
      scrolled: el.scrollTop, firstKey,
      firstKeyAfter: el.querySelector('[data-item-key]')?.getAttribute('data-item-key'),
      renderedAfter: el.querySelectorAll('[data-item-key]').length};
  });

  expect(list.items, 'the legend did not build all the categories').toBeGreaterThan(250);
  expect(list.clientHeight, 'the list collapsed to nothing').toBeGreaterThan(100);
  expect(list.scrollHeight, 'the list is not scrollable')
    .toBeGreaterThan(list.clientHeight + 1);
  expect(list.scrolled, 'the list did not scroll').toBeGreaterThan(0);
  expect(list.rendered, 'every category is in the DOM — the list is not virtualised')
    .toBeLessThan(list.items / 2);
  expect(list.renderedAfter, 'the window grew unbounded while scrolling')
    .toBeLessThan(list.items / 2);
  expect(list.firstKeyAfter, 'scrolling did not move the rendered window')
    .not.toBe(list.firstKey);

  await v.cleanupShell(page);
});
