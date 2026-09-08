/* ---
realizes: [piechart.cp.setup-aggregation-legend-persistence, piechart.int.category-column-drives-legend]
--- */
import {localTest as test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, isLocalBootNoise} from '../../spec-login';
import * as v from '../../helpers/viewers';

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
    if (m.type() === 'error' && !isBenignError(m.text()) && !isLocalBootNoise(m.text()))
      consoleErrors.push(m.text());
  });

  await openDatagrok(page);

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

  v.finishSpec();
});
