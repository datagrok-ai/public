/* ---
realizes: [piechart.cp.setup-aggregation-legend-persistence, piechart.int.category-column-drives-legend]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaApi, deleteProjectWithCleanup} from '../../helpers/projects';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

test('Pie Chart — Category, Legend, Aggregation, and Persistence', async ({page}) => {
  test.setTimeout(300_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  const consoleErrors: string[] = [];

  const isBenignError = (text: string) =>
    /Failed to load resource/.test(text) || /404 \(\)/.test(text) || /favicon/.test(text) ||
    /Unable to find element in cloned iframe/.test(text);
  page.on('console', (m) => {
    if (m.type() === 'error' && !isBenignError(m.text()))
      consoleErrors.push(m.text());
  });

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await v.addViewerByIcon(page, 'pie-chart', 'Pie-chart');

  await v.installEventWaits(page);

  const readRootInDom = () => page.evaluate(() => {
    const pie = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pie chart') as any;
    return !!pie && document.body.contains(pie.root);
  });

  await softStep('Configuration ladder — category, Show Value, aggregations, legend color, title', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;

    console.time('step: category + legendVisibility');
    await v.setViewerProps(page, 'Pie chart', [{set: {categoryColumnName: 'RACE', legendVisibility: 'Always'}}]);
    await v.waitForViewerRendered(page, 'Pie chart', 1000);
    console.timeEnd('step: category + legendVisibility');

    console.time('step: read legend labels');
    const legend = await page.evaluate(() => {
      const pie = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pie chart') as any;
      const df = grok.shell.tv.dataFrame;
      const labels = Array.from(pie.root.querySelectorAll('[name="legend"] .d4-legend-item .d4-legend-value'))
        .map((e: any) => (e.textContent ?? '').trim());
      return {labels, cats: df.col('RACE').categories.slice()};
    });
    console.timeEnd('step: read legend labels');
    expect(legend.labels.length).toBeGreaterThan(0);
    expect([...legend.labels].sort()).toEqual([...legend.cats].sort());

    console.time('step: showValue toggle');
    await v.snapshotCanvasColors(page, 'Pie chart');
    await v.setViewerProps(page, 'Pie chart', [{set: {showValue: true}}]);
    const valueOn = await v.waitForCanvasChange(page, 'Pie chart', {minDelta: 500, timeoutMs: 8000});
    console.timeEnd('step: showValue toggle');

    console.time('step: segmentAngle = count(AGE)');
    await v.setViewerProps(page, 'Pie chart', [{set: {segmentAngleColumnName: 'AGE', segmentAngleAggrType: 'count'}}]);
    await v.waitForViewerRendered(page, 'Pie chart', 1000);
    console.timeEnd('step: segmentAngle = count(AGE)');

    console.time('step: aggr count -> avg');
    await v.snapshotCanvasColors(page, 'Pie chart');
    await v.setViewerProps(page, 'Pie chart', [{set: {segmentAngleAggrType: 'avg'}}]);
    const countToAvg = await v.waitForCanvasChange(page, 'Pie chart', {minDelta: 10000, timeoutMs: 8000});
    console.timeEnd('step: aggr count -> avg');

    console.time('step: aggr avg -> sum');
    await v.snapshotCanvasColors(page, 'Pie chart');
    await v.setViewerProps(page, 'Pie chart', [{set: {segmentAngleAggrType: 'sum'}}]);
    const avgToSum = await v.waitForCanvasChange(page, 'Pie chart', {minDelta: 10000, timeoutMs: 8000});
    console.timeEnd('step: aggr avg -> sum');

    console.log(`Repaint px: valueOn=${valueOn} countToAvg=${countToAvg} avgToSum=${avgToSum}`);
    expect(valueOn).toBeGreaterThan(500);
    expect(countToAvg).toBeGreaterThan(10000);
    expect(avgToSum).toBeGreaterThan(10000);

    console.time('step: legend color change (dialog)');
    await v.changeLegendItemColor(page, {
      viewerType: 'Pie chart',
      category: 'Asian',
      rgb: [214, 39, 40],
      hex: '#d62728',
      column: 'RACE',
    });
    const colorResult = await page.evaluate(() => {
      const pie = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pie chart') as any;
      const after = Array.from(pie.root.querySelectorAll('[name="legend"] .d4-legend-item'))
        .find((i: any) => (i.textContent ?? '').includes('Asian')) as HTMLElement;
      return {
        legendItemFound: !!after,
        swatchColor: after ? after.style.color : null,
        tag: grok.shell.tv.dataFrame.col('RACE').tags['.color-coding-categorical'] ?? null,
        dialogGone: !document.querySelector('.d4-dialog[name="dialog-Asian"]'),
      };
    });
    console.timeEnd('step: legend color change (dialog)');
    expect(colorResult.legendItemFound).toBe(true);
    expect(colorResult.swatchColor).toBe('rgb(214, 39, 40)');
    expect(colorResult.tag).toContain('"Asian":"#d62728"');
    expect(colorResult.dialogGone).toBe(true);

    console.time('step: title set');
    await v.setViewerProps(page, 'Pie chart', [{set: {showTitle: true, title: 'Pie Persistence Probe'}}]);
    await v.waitForViewerRendered(page, 'Pie chart', 1000);
    console.timeEnd('step: title set');

    expect(await readRootInDom()).toBe(true);
    expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
  });

  await softStep('Layout round-trip — saved layout restores the configured pie and viewer set', async () => {
    console.time('step: saveLayout + dapi.layouts.save');
    const layoutId = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      return layout.id as string;
    });
    console.timeEnd('step: saveLayout + dapi.layouts.save');
    try {
      console.time('step: addViewer + loadLayout + verify');
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
        }, (v: any) => v !== null, 3000);
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
      console.timeEnd('step: addViewer + loadLayout + verify');

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
          if (saved)
            await grok.dapi.layouts.delete(saved);
        } catch (_) {}
      }, layoutId);
    }
  });

  await softStep('Project save / Close All / reopen — project restores the configured pie', async () => {
    const projName = 'zz-piechart-persistence-probe-' + Date.now();
    let projectId: string | undefined;
    try {
      console.time('step: save project via API');
      const saved = await saveProjectViaApi(page, projName);
      projectId = saved.projectId;
      console.timeEnd('step: save project via API');

      console.time('step: Close All');
      await v.closeAllAndWait(page);
      console.timeEnd('step: Close All');

      console.time('step: reopen project + verify');
      const result = await page.evaluate(async (id) => {
        const w = window as any;
        await w.__settled('grok.events.onProjectOpened',
          () => grok.dapi.projects.find(id).then((f: any) => f.open()), 1000);

        const pie: any = await w.__poll(() => {
          const t = grok.shell.tv;
          return t ? Array.from(t.viewers).find((x: any) => x.type === 'Pie chart') as any : null;
        }, (v: any) => v != null, 6000, 200);
        await w.__eventFired('viewer:Pie chart.onViewerRendered', 4000);
        const item: HTMLElement | null = await w.__poll(() => {
          const p = pie ?? (grok.shell.tv
            ? Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pie chart') as any : null);
          return p
            ? Array.from(p.root.querySelectorAll('[name="legend"] .d4-legend-item'))
              .find((i: any) => (i.textContent ?? '').includes('Asian')) as HTMLElement ?? null
            : null;
        }, (v: any) => v !== null, 3000);
        const tv = grok.shell.tv;
        const p = pie as any;
        return {
          pieRestored: !!pie,
          cat: p?.props?.categoryColumnName,
          angleCol: p?.props?.segmentAngleColumnName,
          aggr: p?.props?.segmentAngleAggrType,
          showValue: p?.props?.showValue,
          title: p?.props?.title,
          asianSwatch: item ? item.style.color : null,
          tag: tv ? (tv.dataFrame.col('RACE').tags['.color-coding-categorical'] ?? null) : null,
        };
      }, projectId);
      console.timeEnd('step: reopen project + verify');

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
