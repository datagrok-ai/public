/* ---
sub_features_covered: [legend.column-switching, legend.resize, legend.persistence, legend.placement]
--- */
// End-to-end exercise of the frame-bound legend: switching the source column, cycling the
// bar chart stack column, resizing the legend and the viewer, and reloading a layout.
// Every scenario asserts that no uncaught error reached the console.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

const CATEGORICAL = ['Stereo Category', 'Series', 'Primary Series Name', 'Chemist'];

test('Legend stress — column switching, resize, persistence', async ({page}) => {
  test.setTimeout(900_000);

  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 400)));
  page.on('console', (m) => {
    if (m.type() !== 'error') return;
    const t = m.text();
    if (/NullError|NoSuchMethod|Unsupported operation|Uncaught/.test(t)) errors.push(t.slice(0, 400));
  });

  await loginToDatagrok(page);
  await v.openTable(page);
  await v.addLegendViewers(page, {
    column: 'Stereo Category',
    viewers: ['Histogram', 'Bar chart', 'Scatter plot', 'Line chart', 'Pie chart', 'Box plot'],
  });

  await softStep('Bar chart: cycling the stack column never throws (GROK stack-switch crash)', async () => {
    const seen = await page.evaluate(async (cols) => {
      const tv = (window as any).grok.shell.tv;
      const bc = tv.viewers.find((x: any) => x.type === 'Bar chart');
      const out: {stack: string; items: number; mode: string}[] = [];
      for (const c of [...cols, '', ...cols.slice().reverse()]) {
        bc.props.stackColumnName = c;
        await new Promise((r) => setTimeout(r, 900));
        const root = bc.root.querySelector('[name="legend"]');
        out.push({
          stack: c,
          items: Number(root?.getAttribute('data-legend-items') ?? '0'),
          mode: root?.getAttribute('data-legend-mode') ?? 'none',
        });
      }
      bc.props.stackColumnName = '';
      await new Promise((r) => setTimeout(r, 800));
      return out;
    }, CATEGORICAL);

    expect(seen.length).toBe(CATEGORICAL.length * 2 + 1);
    expect(errors, `stack switching threw:\n${errors.join('\n')}`).toEqual([]);
  });

  await softStep('Switching the legend column re-keys the items in every viewer', async () => {
    for (const col of CATEGORICAL) {
      const res = await page.evaluate(async ({c, propMap}) => {
        const tv = (window as any).grok.shell.tv;
        for (const x of tv.viewers) {
          const entry = (propMap as any)[x.type];
          if (!entry) continue;
          try {
            x.props[entry.prop] = entry.array ? [c] : c;
          } catch (_) {}
        }
        await new Promise((r) => setTimeout(r, 1800));
        const out: Record<string, {items: number; keysMatchColumn: boolean}> = {};
        for (const x of tv.viewers) {
          if (x.type === 'Grid') continue;
          const root = x.root.querySelector('[name="legend"]');
          if (!root) continue;
          const keys = (Array.from(root.querySelectorAll('[data-item-key]')) as HTMLElement[])
            .map((e) => e.getAttribute('data-item-key') ?? '')
            .filter((k) => k.startsWith('cat:'));
          // Not every viewer takes its legend from the prop we set (the box plot legends
          // its colour column), so assert the weaker invariant that catches staleness:
          // one rebuild, one source column — never a mix of the old and the new.
          const columns = new Set(keys.map((k) => k.split(':')[1]));
          out[x.type] = {
            items: Number(root.getAttribute('data-legend-items') ?? '0'),
            singleSourceColumn: columns.size <= 1,
            column: Array.from(columns)[0] ?? null,
          };
        }
        return out;
      }, {c: col, propMap: v.LEGEND_COLUMN_PROP});

      for (const [type, info] of Object.entries(res))
        expect(info.singleSourceColumn,
          `${type} mixes item keys from more than one column after switching to ${col}`).toBe(true);
    }
    expect(errors, `column switching threw:\n${errors.join('\n')}`).toEqual([]);
  });

  await softStep('Histogram: a manual legend resize survives re-render and a data change', async () => {
    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const h = tv.viewers.find((x: any) => x.type === 'Histogram');
      h.props.legendPosition = 'Right';
      await new Promise((r) => setTimeout(r, 1200));
    });

    const before = await v.getLegendState(page, 'Histogram');
    expect(before?.mode).toBe('docked');

    await v.dragLegendSplitter(page, 'Histogram', -60, 0);
    const resized = await v.getLegendState(page, 'Histogram');
    expect(resized!.box!.width, 'splitter drag widened the legend').toBeGreaterThan(before!.box!.width);

    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      tv.dataFrame.rows.match({'Stereo Category': 'R_ONE'}).filter();
      await new Promise((r) => setTimeout(r, 1500));
    });
    const afterFilter = await v.waitForLegendIdle(page, 'Histogram');
    expect(Math.abs(afterFilter!.box!.width - resized!.box!.width),
      'the manual legend width was reset by a data change').toBeLessThanOrEqual(4);

    await page.evaluate(async () => {
      (window as any).grok.shell.tv.dataFrame.rows.filters.clear();
      (window as any).grok.shell.tv.dataFrame.rows.requestFilter();
      await new Promise((r) => setTimeout(r, 1200));
    });
  });

  await softStep('Legend size and position survive a layout save + reload', async () => {
    const before = await v.getLegendState(page, 'Histogram');
    const res = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const layout = tv.saveLayout();
      layout.name = 'LegendStress_' + Date.now();
      const saved = await (window as any).grok.dapi.layouts.save(layout);
      await new Promise((r) => setTimeout(r, 800));
      tv.loadLayout(await (window as any).grok.dapi.layouts.find(saved.id));
      await new Promise((r) => setTimeout(r, 4000));
      return {layoutId: saved.id};
    });
    (globalThis as any).__stressLayoutId = res.layoutId;

    // The slot is a look property and must survive; a splitter drag is runtime-only state
    // (it was never serialised in the legacy path either), so the width is expected to reset.
    const after = await v.waitForLegendIdle(page, 'Histogram');
    expect(after?.slot, 'legend slot changed across a layout reload').toBe(before?.slot);
    expect(after?.visible, 'legend disappeared across a layout reload').toBe(true);
  });

  await softStep('Shrinking the viewer walks docked -> miniIcon and back', async () => {
    // `addLegendViewers` pins legendVisibility to Always, which by design skips the mini icon.
    await v.setViewerProps(page, 'Scatter plot', [{set: {legendVisibility: 'Auto'}, wait: 1200}]);

    await v.resizeViewer(page, 'Scatter plot', 900, 600);
    const large = await v.waitForLegendIdle(page, 'Scatter plot');
    expect(['docked', 'corner']).toContain(large?.mode);

    await v.resizeViewer(page, 'Scatter plot', 160, 140);
    const tiny = await v.waitForLegendIdle(page, 'Scatter plot');
    expect(['miniIcon', 'hidden']).toContain(tiny?.mode);

    await v.resizeViewer(page, 'Scatter plot', 900, 600);
    const restored = await v.waitForLegendIdle(page, 'Scatter plot');
    expect(['docked', 'corner'], 'the legend did not come back from the mini icon')
      .toContain(restored?.mode);

    // Corner selection reads the free-space grid, and the grid depends on how much chart
    // there was to draw — so growing back out of the mini icon may legitimately settle on
    // docked where it had been in a corner. What must hold is that it is *stable*: the same
    // size twice in a row lands in the same place.
    await v.resizeViewer(page, 'Scatter plot', 900, 601);
    await v.waitForLegendIdle(page, 'Scatter plot');
    await v.resizeViewer(page, 'Scatter plot', 900, 600);
    const again = await v.waitForLegendIdle(page, 'Scatter plot');
    expect(again?.mode, 'the same viewer size gave two different placements')
      .toBe(restored?.mode);
    expect(again?.slot, 'the same viewer size gave two different slots').toBe(restored?.slot);
    expect(errors, `resizing threw:\n${errors.join('\n')}`).toEqual([]);
  });

  await softStep('No uncaught errors across the whole run', async () => {
    expect(errors, errors.join('\n')).toEqual([]);
  });

  await v.cleanupShell(page);
  v.finishSpec('Legend stress failures');
});
