/* ---
realizes: [piechart.cp.setup-aggregation-legend-persistence]
--- */
import {test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaApi, deleteProjectWithCleanup} from '../../helpers/projects';

declare const grok: any;

// The pie section's server lane: the three scenarios whose subject IS server persistence —
// a layout saved through dapi.layouts, and a project saved, closed and reopened. Everything
// else about the pie chart is client behaviour and runs in the local lane, which needs no
// session and no round-trips (core/docs/features/ui2/LOCAL_MODE.md).
//
// The configured pie these scenarios persist is set up directly rather than re-driven
// through the UI: the ladder that proves each control works (legend colour dialog, canvas
// repaints, aggregation switches) is piechart-category-legend-aggregation-spec's subject,
// and repeating it here would pay for it twice to reach the same state.
test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

// The state all three scenarios persist and then assert survived.
const PIE_STATE = {
  categoryColumnName: 'RACE',
  legendVisibility: 'Always',
  showValue: true,
  segmentAngleColumnName: 'AGE',
  segmentAngleAggrType: 'sum',
  showTitle: true,
  title: 'Pie Persistence Probe',
};

test('Pie Chart — Layout and Project Persistence', async ({page}) => {
  test.setTimeout(300_000);

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'pie-chart', 'Pie-chart');
  await v.installEventWaits(page);

  // Writes the colour-coding tag the legend dialog produces. The dialog itself is covered
  // in the local lane; what persistence has to prove is that the tag survives a round-trip.
  const configurePie = () => page.evaluate(async (state) => {
    const w = window as any;
    const pie = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pie chart') as any;
    grok.shell.tv.dataFrame.col('RACE').tags['.color-coding-categorical'] = '{"Asian":"#d62728"}';
    await w.__settled('viewer:Pie chart.onViewerRendered', () => {
      for (const k of Object.keys(state)) pie.props[k] = (state as any)[k];
    }, 2000);
  }, PIE_STATE);

  await softStep('Layout persistence — saved layout restores the pie geometry', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const pie = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pie chart') as any;

      await w.__settled('viewer:Pie chart.onViewerRendered', () => {
        pie.props.categoryColumnName = 'RACE';
        pie.props.segmentAngleColumnName = 'AGE';
        pie.props.startAngle = 45;
        pie.props.shift = 5;
      }, 2000);

      const before = {
        cat: pie.props.categoryColumnName,
        angle: pie.props.segmentAngleColumnName,
        startAngle: pie.props.startAngle,
        shift: pie.props.shift,
      };

      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;

      await w.__settled('grok.events.onViewerClosed', () => pie.close(), 2000);

      const saved = await grok.dapi.layouts.find(layoutId);
      await w.__settled('grok.events.onViewLayoutApplied', () => grok.shell.tv.loadLayout(saved), 3000);

      const pie2 = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pie chart') as any;
      const after = pie2 ? {
        cat: pie2.props.categoryColumnName,
        angle: pie2.props.segmentAngleColumnName,
        startAngle: pie2.props.startAngle,
        shift: pie2.props.shift,
      } : null;

      await grok.dapi.layouts.delete(saved);
      return {before, after};
    });
    expect(result.after).toEqual(result.before);
  });

  await softStep('Layout round-trip — saved layout restores the configured pie and viewer set', async () => {
    await configurePie();
    const layoutId = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      return layout.id as string;
    });
    try {
      const result = await page.evaluate(async (id) => {
        const w = window as any;
        const tv = grok.shell.tv;
        await w.__settled('grok.events.onViewerAdded', () => tv.addViewer('Scatter plot'), 500);
        const saved = await grok.dapi.layouts.find(id);
        await w.__settled('grok.events.onViewLayoutApplied', () => tv.loadLayout(saved), 3000);
        await w.__eventFired('viewer:Pie chart.onViewerRendered', 4000);
        const item: HTMLElement | null = await w.__poll(() => {
          const t = grok.shell.tv;
          const p = t?.viewers ? Array.from(t.viewers).find((x: any) => x.type === 'Pie chart') as any : null;
          return p
            ? Array.from(p.root.querySelectorAll('[name="legend"] .d4-legend-item'))
              .find((i: any) => (i.textContent ?? '').includes('Asian')) as HTMLElement ?? null
            : null;
        }, (val: any) => val !== null, 3000);
        const viewers = Array.from(tv.viewers) as any[];
        const pie = viewers.find((x: any) => x.type === 'Pie chart');
        return {
          hasScatter: viewers.some((x: any) => x.type === 'Scatter plot'),
          hasPie: !!pie,
          cat: pie?.props.categoryColumnName,
          angleCol: pie?.props.segmentAngleColumnName,
          aggr: pie?.props.segmentAngleAggrType,
          showValue: pie?.props.showValue,
          title: pie?.props.title,
          asianSwatch: item ? item.style.color : null,
        };
      }, layoutId);

      expect(result.hasScatter).toBe(false);
      expect(result.hasPie).toBe(true);
      expect(result.cat).toBe('RACE');
      expect(result.angleCol).toBe('AGE');
      expect(result.aggr).toBe('sum');
      expect(result.showValue).toBe(true);
      expect(result.title).toBe('Pie Persistence Probe');
      expect(result.asianSwatch).toBe('rgb(214, 39, 40)');
    } finally {
      await page.evaluate(async (id) => {
        try {
          const saved = await grok.dapi.layouts.find(id);
          if (saved) await grok.dapi.layouts.delete(saved);
        } catch (_) {}
      }, layoutId);
    }
  });

  await softStep('Project save / Close All / reopen — project restores the configured pie', async () => {
    await configurePie();
    const projName = 'zz-piechart-persistence-probe-' + Date.now();
    let projectId: string | undefined;
    try {
      const saved = await saveProjectViaApi(page, projName);
      projectId = saved.projectId;

      await v.closeAllAndWait(page);

      const result = await page.evaluate(async (id) => {
        const w = window as any;
        await w.__settled('grok.events.onProjectOpened',
          () => grok.dapi.projects.find(id).then((f: any) => f.open()), 1000);

        const pie: any = await w.__poll(() => {
          const t = grok.shell.tv;
          return t ? Array.from(t.viewers).find((x: any) => x.type === 'Pie chart') as any : null;
        }, (val: any) => val != null, 6000, 200);
        await w.__eventFired('viewer:Pie chart.onViewerRendered', 4000);
        const item: HTMLElement | null = await w.__poll(() => {
          const p = pie ?? (grok.shell.tv
            ? Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pie chart') as any : null);
          return p
            ? Array.from(p.root.querySelectorAll('[name="legend"] .d4-legend-item'))
              .find((i: any) => (i.textContent ?? '').includes('Asian')) as HTMLElement ?? null
            : null;
        }, (val: any) => val !== null, 3000);
        const tv = grok.shell.tv;
        return {
          pieRestored: !!pie,
          cat: pie?.props?.categoryColumnName,
          angleCol: pie?.props?.segmentAngleColumnName,
          aggr: pie?.props?.segmentAngleAggrType,
          showValue: pie?.props?.showValue,
          title: pie?.props?.title,
          asianSwatch: item ? item.style.color : null,
          tag: tv ? (tv.dataFrame.col('RACE').tags['.color-coding-categorical'] ?? null) : null,
        };
      }, projectId);

      expect(projectId).toBeTruthy();
      expect(result.pieRestored).toBe(true);
      expect(result.cat).toBe('RACE');
      expect(result.angleCol).toBe('AGE');
      expect(result.aggr).toBe('sum');
      expect(result.showValue).toBe(true);
      expect(result.title).toBe('Pie Persistence Probe');
      expect(result.asianSwatch).toBe('rgb(214, 39, 40)');
      expect(result.tag).toContain('"Asian":"#d62728"');
    } finally {
      await deleteProjectWithCleanup(page, {projectId});
    }
  });

  v.finishSpec();
});
