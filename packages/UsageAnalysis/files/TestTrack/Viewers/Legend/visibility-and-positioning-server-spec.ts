/* ---
realizes: [viewers.scatter-plot, viewers.histogram]
--- */
import {test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {deleteEntities, layoutRoundTrip, projectRoundTrip} from './persistence';

// The server lane of the visibility-and-positioning scenario: Sc7, Sc8 steps 5-6, Sc10 steps 7-8
// and Sc11 — legend column, custom colour, visibility and position surviving layout and project
// round-trips. The gestures are in visibility-and-positioning-spec.ts on the local lane.
test.use(specTestOptions);

const probe = {viewerType: 'Scatter plot', props: ['colorColumnName', 'legendVisibility', 'legendPosition']};

function setScatterLegend(page: any, props: Record<string, any>): Promise<Record<string, any>> {
  return page.evaluate(async (p: Record<string, any>) => {
    const w = window as any;
    const sp = w.grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot');
    await w.__settled('viewer:Scatter plot.onViewerRendered', () => {
      for (const k of Object.keys(p)) try { sp.props[k] = p[k]; } catch (_) {}
    }, 1500);
    return {vis: sp.props.legendVisibility, pos: sp.props.legendPosition, color: sp.props.colorColumnName};
  }, props);
}

test('Legend visibility and positioning — layout and project persistence', async ({page}) => {
  test.setTimeout(900_000);

  await openDatagrok(page);
  await v.openTable(page);
  await v.installEventWaits(page);
  await v.addLegendViewers(page, {column: 'Stereo Category', viewers: ['Scatter plot', 'Histogram'], settleMs: 500});

  await softStep('Setup: Stereo Category legend, R_ONE recoloured, Visibility=Always, Position=Auto', async () => {
    await page.evaluate(() => {
      const col = (window as any).grok.shell.tv.dataFrame.col('Stereo Category');
      col.tags['.color-coding-type'] = 'Categorical';
      col.meta.colors.setCategorical({'R_ONE': '#1f77b4'});
    });
    const res = await setScatterLegend(page, {legendVisibility: 'Always', legendPosition: 'Auto'});
    expect(res.vis).toBe('Always');
    expect(res.color).toBe('Stereo Category');
  });

  let layoutId1: string | null = null;
  let layoutId2: string | null = null;
  let layoutId3: string | null = null;
  let projectId: string | null = null;
  try {
    await softStep('Sc7 steps 1-3: save+reapply layout, state persists', async () => {
      const res = await layoutRoundTrip(page, 'LegendVP', probe);
      layoutId1 = res.layoutId;
      expect(typeof layoutId1).toBe('string');
      expect(res.props.colorColumnName).toBe('Stereo Category');
      expect(await page.evaluate(() =>
        !!(window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot'))).toBe(true);
    });

    await softStep('Sc8 steps 5-6: layout round-trip (Always + Auto persist)', async () => {
      const res = await layoutRoundTrip(page, 'LegendVP2', probe);
      layoutId2 = res.layoutId;
      expect(res.props.legendVisibility).toBe('Always');
    });

    await softStep('Sc10 steps 7-8: layout round-trip (corner position persists)', async () => {
      const set = await setScatterLegend(page, {legendPosition: 'RightTop'});
      expect(set.pos).toBe('RightTop');
      const res = await layoutRoundTrip(page, 'LegendVP3', probe);
      layoutId3 = res.layoutId;
      expect(res.props.legendPosition).toBeTruthy();
    });

    await softStep('Sc11 steps 1-3: project save+close+reopen (FK graceful-degrade)', async () => {
      const res = await projectRoundTrip(page, 'LegendVPProj', probe);
      projectId = res.projectId;
      expect(res.ok, res.ok ? '' : `project save+reopen failed in phase '${res.phase}': ${res.error}`).toBe(true);
      expect(res.props.legendVisibility).toBeTruthy();
      expect(res.props.legendPosition).toBeTruthy();
    });
  } finally {
    await softStep('Cleanup: drop layouts/projects + closeAll', async () => {
      await deleteEntities(page, {layoutIds: [layoutId1, layoutId2, layoutId3], projectId});
    });
  }

  v.finishSpec();
});
