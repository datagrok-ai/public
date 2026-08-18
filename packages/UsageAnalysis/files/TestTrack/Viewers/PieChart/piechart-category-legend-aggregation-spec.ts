/* ---
realizes: [piechart.cp.setup-aggregation-legend-persistence, piechart.int.category-column-drives-legend]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaUI, deleteProjectWithCleanup} from '../../helpers/projects';

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

  const readRootInDom = () => page.evaluate(() => {
    const pie = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pie chart') as any;
    return !!pie && document.body.contains(pie.root);
  });

  await softStep('Configuration ladder — category, Show Value, aggregations, legend color, title', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;

    await v.setViewerProps(page, 'Pie chart', [{set: {categoryColumnName: 'RACE', legendVisibility: 'Always'}}]);
    await v.waitForViewerRendered(page, 'Pie chart', 1000);
    const legend = await page.evaluate(() => {
      const pie = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pie chart') as any;
      const df = grok.shell.tv.dataFrame;
      const labels = Array.from(pie.root.querySelectorAll('[name="legend"] .d4-legend-item .d4-legend-value'))
        .map((e: any) => (e.textContent ?? '').trim());
      return {labels, cats: df.col('RACE').categories.slice()};
    });
    expect(legend.labels.length).toBeGreaterThan(0);
    expect([...legend.labels].sort()).toEqual([...legend.cats].sort());

    await v.snapshotCanvasColors(page, 'Pie chart');
    await v.setViewerProps(page, 'Pie chart', [{set: {showValue: true}}]);
    const valueOn = await v.waitForCanvasChange(page, 'Pie chart', {minDelta: 500, timeoutMs: 8000});

    await v.setViewerProps(page, 'Pie chart', [{set: {segmentAngleColumnName: 'AGE', segmentAngleAggrType: 'count'}}]);
    await v.waitForViewerRendered(page, 'Pie chart', 1000);
    await v.snapshotCanvasColors(page, 'Pie chart');
    await v.setViewerProps(page, 'Pie chart', [{set: {segmentAngleAggrType: 'avg'}}]);
    const countToAvg = await v.waitForCanvasChange(page, 'Pie chart', {minDelta: 10000, timeoutMs: 8000});
    await v.snapshotCanvasColors(page, 'Pie chart');
    await v.setViewerProps(page, 'Pie chart', [{set: {segmentAngleAggrType: 'sum'}}]);
    const avgToSum = await v.waitForCanvasChange(page, 'Pie chart', {minDelta: 10000, timeoutMs: 8000});

    console.log(`Repaint px: valueOn=${valueOn} countToAvg=${countToAvg} avgToSum=${avgToSum}`);
    expect(valueOn).toBeGreaterThan(500);
    expect(countToAvg).toBeGreaterThan(10000);
    expect(avgToSum).toBeGreaterThan(10000);

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
    expect(colorResult.legendItemFound).toBe(true);
    expect(colorResult.swatchColor).toBe('rgb(214, 39, 40)');
    expect(colorResult.tag).toContain('"Asian":"#d62728"');
    expect(colorResult.dialogGone).toBe(true);

    await v.setViewerProps(page, 'Pie chart', [{set: {showTitle: true, title: 'Pie Persistence Probe'}}]);
    await v.waitForViewerRendered(page, 'Pie chart', 1000);

    expect(await readRootInDom()).toBe(true);
    expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
  });

  await softStep('Layout round-trip — saved layout restores the configured pie and viewer set', async () => {
    const layoutId = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      return layout.id as string;
    });
    try {
      const result = await page.evaluate(async (id) => {
        const tv = grok.shell.tv;
        tv.addViewer('Scatter plot');
        await new Promise<void>((resolve) => {
          const sub = grok.events.onViewerAdded.subscribe(() => { sub.unsubscribe(); resolve(); });
          setTimeout(resolve, 500);
        });
        const saved = await grok.dapi.layouts.find(id);
        const applied = new Promise<void>((resolve) => {
          const sub = grok.events.onViewLayoutApplied.subscribe(() => { sub.unsubscribe(); resolve(); });
          setTimeout(resolve, 3000);
        });
        tv.loadLayout(saved);
        await applied;

        const pieAfter = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pie chart') as any;
        await new Promise<void>((resolve) => {
          let sub: any = null;
          try { sub = pieAfter?.onViewerRendered.subscribe(() => { sub.unsubscribe(); resolve(); }); }
          catch (_) { resolve(); }
          setTimeout(() => { try { sub?.unsubscribe(); } catch (_) {} resolve(); }, 4000);
        });
        let item: HTMLElement | null = null;
        for (let i = 0; i < 30; i++) {
          const t = grok.shell.tv;
          const p = t?.viewers ? Array.from(t.viewers).find((x: any) => x.type === 'Pie chart') as any : null;
          item = p
            ? Array.from(p.root.querySelectorAll('[name="legend"] .d4-legend-item'))
              .find((i: any) => (i.textContent ?? '').includes('Asian')) as HTMLElement ?? null
            : null;
          if (item) break;
          await new Promise((r) => setTimeout(r, 100));
        }
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
      const saved = await saveProjectViaUI(page, projName);
      projectId = saved.projectId;

      const result = await page.evaluate(async (id) => {
        grok.shell.closeAll();
        await new Promise<void>((resolve) => {
          const sub = grok.events.onProjectOpened.subscribe(() => { sub.unsubscribe(); resolve(); });
          grok.dapi.projects.find(id).then((f: any) => f.open()).finally(() => setTimeout(resolve, 500));
        });

        let pie: any = null;
        for (let a = 0; a < 30 && !pie; a++) {
          await new Promise((r) => setTimeout(r, 200));
          const t = grok.shell.tv;
          pie = t ? Array.from(t.viewers).find((x: any) => x.type === 'Pie chart') as any : null;
        }
        await new Promise<void>((resolve) => {
          let sub: any = null;
          try { sub = pie?.onViewerRendered.subscribe(() => { sub.unsubscribe(); resolve(); }); }
          catch (_) { resolve(); }
          setTimeout(() => { try { sub?.unsubscribe(); } catch (_) {} resolve(); }, 4000);
        });
        let item: HTMLElement | null = null;
        for (let a = 0; a < 30 && !item; a++) {
          const p = pie ?? (grok.shell.tv
            ? Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pie chart') as any : null);
          item = p
            ? Array.from(p.root.querySelectorAll('[name="legend"] .d4-legend-item'))
              .find((i: any) => (i.textContent ?? '').includes('Asian')) as HTMLElement ?? null
            : null;
          if (item) break;
          await new Promise((r) => setTimeout(r, 100));
        }
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
