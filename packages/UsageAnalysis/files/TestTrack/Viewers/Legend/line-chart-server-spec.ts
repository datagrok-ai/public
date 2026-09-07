/* ---
realizes: [viewers.line-chart]
--- */
import {test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {deleteEntities, layoutRoundTrip, projectRoundTrip} from './persistence';

// The server lane of the line-chart legend scenario: Sc3 steps 3-5 and Sc4 steps 4-7, the multi-axis
// split and the replaced Y column surviving layout and project round-trips. The viewer steps that
// build that state are in line-chart-spec.ts on the local lane; here it is set through the API.
test.use(specTestOptions);

const probe = {viewerType: 'Line chart', props: ['multiAxis', 'splitColumnName', 'yColumnNames']};

function setLineChart(page: any, options: Record<string, any>): Promise<Record<string, any>> {
  return page.evaluate(async (o: Record<string, any>) => {
    const w = window as any;
    const lc = w.grok.shell.tv.viewers.find((x: any) => x.type === 'Line chart');
    await w.__settled('viewer:Line chart.onViewerRendered', () => { for (const k of Object.keys(o)) lc.props[k] = o[k]; }, 1500);
    return {multiAxis: lc.props.multiAxis, split: lc.props.splitColumnName, yCols: lc.props.yColumnNames};
  }, options);
}

test('Line chart legend — layout and project persistence', async ({page}) => {
  test.setTimeout(900_000);

  await openDatagrok(page);
  await v.openTable(page);
  await v.installEventWaits(page);

  await softStep('Setup: Line chart, Split=Series, Multi Axis on', async () => {
    await page.evaluate(async () => {
      const w = window as any;
      w.grok.shell.tv.addViewer('Line chart');
      await w.__poll(() => w.grok.shell.tv.viewers.filter((x: any) => x.type === 'Line chart').length,
        (c: number) => c > 0, 1000);
    });
    const res = await setLineChart(page, {splitColumnName: 'Series', legendVisibility: 'Always', multiAxis: true});
    expect(res.multiAxis).toBe(true);
    expect(res.split).toBe('Series');
  });

  let layoutId1: string | null = null;
  let layoutId2: string | null = null;
  let projectId: string | null = null;
  try {
    await softStep('Sc3 steps 1-3: save+reapply layout (multiAxis+split persist)', async () => {
      const res = await layoutRoundTrip(page, 'LineChart', probe);
      layoutId1 = res.layoutId;
      expect(res.props.multiAxis).toBe(true);
      expect(res.props.splitColumnName).toBe('Series');
    });

    await softStep('Sc4 step 3: replace Y column → NIBR logP', async () => {
      const res = await setLineChart(page, {yColumnNames: ['Average Mass', 'NIBR logP']});
      expect(res.yCols).toEqual(['Average Mass', 'NIBR logP']);
    });

    await softStep('Sc4 steps 4-5: save+reapply layout (new Y persists)', async () => {
      const res = await layoutRoundTrip(page, 'LineChart2', probe);
      layoutId2 = res.layoutId;
      expect(res.props.yColumnNames).toEqual(['Average Mass', 'NIBR logP']);
    });

    await softStep('Sc3 steps 4-5 / Sc4 steps 6-7: project save+close+reopen (FK graceful-degrade)', async () => {
      const res = await projectRoundTrip(page, 'LineChartProj', probe);
      projectId = res.projectId;
      expect(res.ok, res.ok ? '' : `project save+reopen failed in phase '${res.phase}': ${res.error}`).toBe(true);
      expect(res.props.multiAxis).toBe(true);
      expect(res.props.splitColumnName).toBe('Series');
      expect(res.props.yColumnNames).toEqual(['Average Mass', 'NIBR logP']);
    });
  } finally {
    await softStep('Cleanup', async () => {
      await deleteEntities(page, {layoutIds: [layoutId1, layoutId2], projectId});
    });
  }

  v.finishSpec();
});
