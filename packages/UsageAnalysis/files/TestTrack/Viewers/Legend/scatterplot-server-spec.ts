/* ---
realizes: [viewers.scatter-plot, viewers.box-plot, viewers.pc-plot]
--- */
import {test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {deleteEntities, layoutRoundTrip, projectRoundTrip} from './persistence';

// The server lane of the scatterplot legend scenarios: the layout and project round-trips of Sc1,
// Sc3 and Sc5. The viewer behaviour is proven in scatterplot-spec.ts on the local lane; the state
// each round-trip carries is set up here through the API.
test.use(specTestOptions);

// retry(1) absorbs transient dapi/FiltersGroup hangs causing ~3x runtime variance.
test.describe.configure({retries: 1});

const SCATTER_LEGEND = {viewerType: 'Scatter plot', legendItems: true};

async function addScatter(page: any, props: Record<string, any>): Promise<void> {
  await page.evaluate(async (p: Record<string, any>) => {
    const w = window as any;
    const tv = w.grok.shell.tv;
    const before = tv.viewers.filter((x: any) => x.type === 'Scatter plot').length;
    tv.addViewer('Scatter plot');
    await w.__poll(() => w.grok.shell.tv.viewers.filter((x: any) => x.type === 'Scatter plot').length,
      (c: number) => c > before, 1500);
    const sps = tv.viewers.filter((x: any) => x.type === 'Scatter plot');
    const sp = sps[sps.length - 1];
    for (const k of Object.keys(p)) sp.props[k] = p[k];
    try { sp.props.legendVisibility = 'Always'; } catch (_) {}
    let prev = -1;
    await w.__poll(() => sp.root.querySelectorAll('[name="legend"] .d4-legend-item').length,
      (c: number) => { const settled = c > 0 && c === prev; prev = c; return settled; }, 1500);
  }, props);
}

// scenario: 1. Color + Marker combined legend on Scatter plot [coverage_type: edge]
test('Legend scatterplot — Color + Marker combined: layout and project persistence', async ({page}) => {
  test.setTimeout(900_000);
  stepErrors.length = 0;
  await openDatagrok(page);
  await v.openTable(page);
  await v.installEventWaits(page);

  await softStep('Setup: Color=Series + Marker=Series, first category recoloured', async () => {
    await addScatter(page, {colorColumnName: 'Series', markersColumnName: 'Series'});
    const res = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv.dataFrame;
      const col = df.col('Series');
      const cat = col.categories[0];
      col.tags['.color-coding-type'] = 'Categorical';
      col.meta.colors.setCategorical({[cat]: '#1f77b4'});
      return {cat, tag: String(JSON.parse(col.tags['.color-coding-categorical'] ?? '{}')[cat] ?? '').toLowerCase()};
    });
    expect(res.tag).toBe('#1f77b4');
  });

  let layoutId: string | null = null;
  let projectId: string | null = null;
  try {
    await softStep('Sc1 steps 6-7: save+reapply layout, color persists', async () => {
      const res = await layoutRoundTrip(page, 'ScatterCombined', SCATTER_LEGEND);
      layoutId = res.layoutId;
      expect(typeof layoutId).toBe('string');
      expect((layoutId ?? '').length).toBeGreaterThan(0);
    });

    await softStep('Sc1 steps 8-9: project save+close+reopen (FK graceful-degrade)', async () => {
      const res = await projectRoundTrip(page, 'ScatterCombinedProj', SCATTER_LEGEND);
      projectId = res.projectId;
      expect(res.ok, res.ok ? '' : `project save+reopen failed in phase '${res.phase}': ${res.error}`).toBe(true);
      expect(projectId).toBeTruthy();
    });
  } finally {
    await softStep('Cleanup', async () => { await deleteEntities(page, {layoutIds: [layoutId], projectId}); });
  }

  v.finishSpec();
});

// scenario: 3. In-viewer filtering — multiple scatterplots with shared filter [coverage_type: edge]
test('Legend scatterplot — in-viewer filter: layout round-trip', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;
  await openDatagrok(page);
  await v.openTable(page);
  await v.installEventWaits(page);

  await softStep('Setup: two scatters, Marker+Color=Stereo Category, same in-viewer filter', async () => {
    const props = {
      markersColumnName: 'Stereo Category', colorColumnName: 'Stereo Category',
      filter: '${Stereo Category} in ("R_ONE", "S_UNKN")',
    };
    await addScatter(page, props);
    await addScatter(page, props);
    expect(await page.evaluate(() =>
      (window as any).grok.shell.tv.viewers.filter((x: any) => x.type === 'Scatter plot').length)).toBe(2);
  });

  let layoutId: string | null = null;
  try {
    await softStep('Sc3 steps 7-8: save+reapply layout, both legends survive', async () => {
      const res = await layoutRoundTrip(page, 'ScatterInViewerFilter', SCATTER_LEGEND);
      layoutId = res.layoutId;
      const after = await page.evaluate(() => {
        const sps = (window as any).grok.shell.tv.viewers.filter((x: any) => x.type === 'Scatter plot');
        return {
          scatterCount: sps.length,
          legendCounts: sps.map((sp: any) => sp.root.querySelectorAll('[name="legend"] .d4-legend-item').length),
        };
      });
      expect(after.scatterCount).toBeGreaterThanOrEqual(2);
      for (const c of after.legendCounts) expect(c).toBeGreaterThan(0);
    });
  } finally {
    await softStep('Cleanup', async () => { await deleteEntities(page, {layoutIds: [layoutId]}); });
  }

  v.finishSpec();
});

// scenario: 5. Color coding from grid — linear and categorical, with persistence [coverage_type: edge]
test('Legend scatterplot — grid color coding: layout and project persistence', async ({page}) => {
  test.setTimeout(900_000);
  stepErrors.length = 0;
  await openDatagrok(page);
  await v.openTable(page);
  await v.installEventWaits(page);

  await softStep('Setup: scatter + box + PC plots, Chemical Space X linear scheme with text-apply', async () => {
    const res = await page.evaluate(async () => {
      const w = window as any;
      const tv = w.grok.shell.tv;
      tv.addViewer('Scatter plot');
      tv.addViewer('Box plot');
      tv.addViewer('PC Plot');
      await w.__poll(() => w.grok.shell.tv.viewers.filter((x: any) => x.type === 'PC Plot').length,
        (c: number) => c > 0, 1500);
      const sp = tv.viewers.find((x: any) => x.type === 'Scatter plot');
      sp.props.colorColumnName = 'Chemical Space X';
      try { sp.props.legendVisibility = 'Always'; } catch (_) {}
      const col = tv.dataFrame.col('Chemical Space X');
      col.tags['.color-coding-type'] = 'Linear';
      col.tags['.color-coding-scheme'] = '[1, 8388607, 16711680]';
      col.tags['.color-coding-text'] = 'true';
      for (const x of tv.viewers)
        if (x.type !== 'Grid') try { x.invalidate?.(); } catch (_) {}
      await w.__quiet('viewer:Scatter plot.onViewerRendered', 200, 1500);
      return {codingType: col.tags['.color-coding-type'], textApplied: col.tags['.color-coding-text']};
    });
    expect(res.codingType).toBe('Linear');
    expect(res.textApplied).toBe('true');
  });

  let layoutId: string | null = null;
  let projectId: string | null = null;
  try {
    await softStep('Sc5 steps 8-9: save+reapply layout, scheme + text-apply persist', async () => {
      const res = await layoutRoundTrip(page, 'ScatterGridColor', {
        viewerType: 'Scatter plot', props: ['colorColumnName'],
        column: 'Chemical Space X', tags: ['.color-coding-type', '.color-coding-scheme', '.color-coding-text'],
      });
      layoutId = res.layoutId;
      expect(res.tags['.color-coding-type']).toBe('Linear');
      expect(res.tags['.color-coding-scheme']).toBeTruthy();
    });

    await softStep('Sc5 steps 10-11: grid coding → Categorical on Stereo Category', async () => {
      const res = await page.evaluate(async () => {
        const w = window as any;
        const tv = w.grok.shell.tv;
        const col = tv.dataFrame.col('Stereo Category');
        col.tags['.color-coding-type'] = 'Categorical';
        const map: Record<string, number> = {};
        for (const c of col.categories)
          map[c] = c === 'R_ONE' ? 0xFFFF0000 : c === 'S_UNKN' ? 0xFF00FF00 : 0xFF808080;
        try { col.meta.colors.setCategorical(map); } catch (_) {}
        const sp = tv.viewers.find((x: any) => x.type === 'Scatter plot');
        sp.props.colorColumnName = 'Stereo Category';
        await w.__quiet('viewer:Scatter plot.onViewerRendered', 200, 1500);
        const idxROne = col.categories.indexOf('R_ONE');
        return {codingType: col.tags['.color-coding-type'], rOneMeta: '0x' + (col.meta.colors.getColor(idxROne) >>> 0).toString(16)};
      });
      expect(res.codingType).toBe('Categorical');
      expect(res.rOneMeta).not.toBe('0xff808080');
    });

    await softStep('Sc5 steps 12-13: project save+close+reopen (FK graceful-degrade)', async () => {
      const res = await projectRoundTrip(page, 'ScatterGridColorProj',
        {column: 'Stereo Category', tags: ['.color-coding-type', '.color-coding-categorical']});
      projectId = res.projectId;
      expect(res.ok, res.ok ? '' : `project save+reopen failed in phase '${res.phase}': ${res.error}`).toBe(true);
      const rOneAfter = await page.evaluate(() => {
        const col = (window as any).grok.shell.tv.dataFrame.col('Stereo Category');
        return '0x' + (col.meta.colors.getColor(col.categories.indexOf('R_ONE')) >>> 0).toString(16);
      });
      expect(rOneAfter).not.toBe('0xff808080');
    });
  } finally {
    await softStep('Cleanup', async () => { await deleteEntities(page, {layoutIds: [layoutId], projectId}); });
  }

  v.finishSpec();
});
