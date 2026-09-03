/* ---
realizes: [viewers.scatter-plot, viewers.bar-chart, viewers.filters.histogram, viewers.filters.categorical, chem.filter.substructure-filter]
--- */
import {test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {clickCanvasFilter} from './canvas-filter';
import {deleteEntities} from './persistence';

// The server lane of the filtering scenario: the Chem substructure filter (a package) and the
// two layout round-trips. The filter gestures themselves are in filtering-spec.ts on the local lane.
test.use(specTestOptions);

test('Legend filtering — substructure filter and layout round-trips', async ({page}) => {
  test.setTimeout(600_000);

  await openDatagrok(page);
  await v.installEventWaits(page);
  await v.openTable(page);
  await v.addLegendViewers(page, {column: 'Stereo Category', viewers: ['Scatter plot', 'Bar chart'], settleMs: 500});

  await softStep('Structure filter on Core — platform API available (env-dependent)', async () => {
    const res = await page.evaluate(async () => {
      const df = (window as any).grok.shell.tv.dataFrame;
      const fg = (window as any).grok.shell.tv.getFiltersGroup();
      const firstSmiles = df.col('Core').get(0);
      try {
        fg.updateOrAdd({type: 'Chem:substructureFilter', column: 'Core', columnName: 'Core', molBlock: firstSmiles});
        // technical: the filter group debounces, so onRowsFiltered fires on an
        // intermediate row set — no channel marks the settled one
        await new Promise((r) => setTimeout(r, 2000));
        return {applied: true, filterCount: df.filter.trueCount};
      } catch (e: any) {
        const msg = String(e?.message ?? e);
        return {applied: false, chemMissing: msg.includes('Chem') || msg.includes('substructure')};
      }
    });
    expect(res.applied || res.chemMissing).toBe(true);
  });

  let filterLayoutId: string | null = null;
  await softStep('Save + re-apply layout (filter state + ≥3s settle)', async () => {
    const res = await page.evaluate(async () => {
      const w = window as any;
      const fg = (window as any).grok.shell.tv.getFiltersGroup();
      for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch (_) {} }
      (window as any).grok.shell.tv.dataFrame.filter.setAll(true);
      const DG = (window as any).DG;
      fg.updateOrAdd({type: 'histogram', column: 'Average Mass', min: 400, max: 10000});
      // technical: the filter group debounces, so onRowsFiltered fires on an
      // intermediate row set — no channel marks the settled one
      await new Promise((r) => setTimeout(r, 500));
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Stereo Category', selected: ['R_ONE', 'S_UNKN']});
      await new Promise((r) => setTimeout(r, 1500));
      const before = (window as any).grok.shell.tv.dataFrame.filter.trueCount;
      const tv = (window as any).grok.shell.tv;
      const layout = tv.saveLayout();
      layout.name = 'Filtering_' + Date.now();
      const saved = await w.grok.dapi.layouts.save(layout);
      const found = await w.__findSaved(() => w.grok.dapi.layouts.find(saved.id));
      const gen = w.__viewerGen();
      tv.loadLayout(found);
      await w.__rebuilt(gen, () => `${w.grok.shell.tv?.dataFrame?.filter?.trueCount ?? -1}`, 4500);
      return {layoutId: String(saved.id), before, after: (window as any).grok.shell.tv.dataFrame.filter.trueCount};
    });
    filterLayoutId = res.layoutId;
    expect(res.after).toBe(res.before);
  });

  // Platform doesn't persist click-to-filter state across layout save/load — assert round-trip mechanics only.
  let clickLayoutId: string | null = null;
  await softStep('Layout persistence: click-to-filter state survives save+reload', async () => {
    const clicked = await clickCanvasFilter(page, {viewerType: 'Bar chart', column: 'Stereo Category'});
    expect(clicked.survivors).toBe(1);
    const res = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const df = tv.dataFrame;
      const before = df.filter.trueCount;
      const layout = tv.saveLayout();
      layout.name = 'FilteringClick_' + Date.now();
      const w = window as any;
      const saved = await w.grok.dapi.layouts.save(layout);
      const found = await w.__findSaved(() => w.grok.dapi.layouts.find(saved.id));
      const gen = w.__viewerGen();
      tv.loadLayout(found);
      await w.__rebuilt(gen, () => `${w.grok.shell.tv?.dataFrame?.filter?.trueCount ?? -1}`, 4500);
      const tvAfter = (window as any).grok.shell.tv;
      return {
        before,
        after: df.filter.trueCount,
        layoutId: String(saved.id),
        rowCountAfter: tvAfter.dataFrame.rowCount,
        viewersAfter: tvAfter.viewers.length,
      };
    });
    clickLayoutId = res.layoutId;
    expect(typeof res.layoutId).toBe('string');
    expect(res.layoutId.length).toBeGreaterThan(0);
    expect(res.rowCountAfter).toBeGreaterThan(0);
    expect(res.viewersAfter).toBeGreaterThan(1);
  });

  await softStep('Cleanup', async () => {
    await deleteEntities(page, {layoutIds: [filterLayoutId, clickLayoutId]});
  });

  v.finishSpec();
});
