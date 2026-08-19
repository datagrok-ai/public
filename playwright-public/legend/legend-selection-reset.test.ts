/* ---
sub_features_covered: [legend.selection-reset-on-source-change]
--- */
// Legend selection is stored as positional category indexes and restored across model
// rebuilds. Restoring it across a SOURCE change (another legend column, another split set)
// silently selected whatever sat at the same position — pie on DIS_POP with 'Indigestion'
// selected switched to RACE with 'Black' selected, and its filter followed. A selection now
// survives only rebuilds whose section source is unchanged (LegendSection.sourceKey).

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

test('Legend selection resets when the source changes, survives when it does not', async ({page}) => {
  test.setTimeout(900_000);
  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await page.evaluate(async () => {
    const grok = (window as any).grok;
    grok.shell.closeAll();
    // side panels overlay the viewer's right edge at forced sizes — legend clicks would
    // land on them instead of the items
    grok.shell.windows.showProperties = false;
    grok.shell.windows.showContextPanel = false;
    grok.shell.windows.showHelp = false;
    grok.shell.windows.showToolbox = false;
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    grok.shell.addTableView(df);
  });
  await page.waitForTimeout(2500);

  const countPixels = (viewerName: string) => page.evaluate((vn) => {
    const root = document.querySelector(`[name="viewer-${vn.replace(/ /g, '-')}"]`) as HTMLElement;
    const canvases = Array.from(root?.querySelectorAll('canvas') ?? []) as HTMLCanvasElement[];
    const c = canvases.sort((a, b) => b.width * b.height - a.width * a.height)[0];
    if (!c) return -1;
    const d = c.getContext('2d')!.getImageData(0, 0, c.width, c.height).data;
    let n = 0;
    for (let i = 0; i < d.length; i += 4)
      if (d[i + 3] > 100 && (d[i] < 235 || d[i + 1] < 235 || d[i + 2] < 235)) n++;
    return n;
  }, viewerName);

  const legendState = () => page.evaluate(() => {
    const root = document.querySelector('[name="legend"]') as HTMLElement;
    const items = Array.from(root?.querySelectorAll('[name="legend-item"]') ?? []) as HTMLElement[];
    return {labels: items.map((i) => (i.textContent || '').trim()),
      selected: items.filter((i) => i.getAttribute('data-item-selected') === 'true')
        .map((i) => (i.textContent || '').trim())};
  });

  const clickItem = async (needle: string) => {
    const pos = await page.evaluate((n) => {
      const root = document.querySelector('[name="legend"]') as HTMLElement;
      const items = Array.from(root.querySelectorAll('[name="legend-item"]')) as HTMLElement[];
      const it = items.find((i) => (i.textContent || '').includes(n));
      if (!it) return null;
      const r = it.getBoundingClientRect();
      return {x: r.x + r.width / 2, y: r.y + r.height / 2};
    }, needle);
    expect(pos, `legend item ${needle} not found`).not.toBeNull();
    await page.mouse.click(pos!.x, pos!.y);
    await page.waitForTimeout(1500);
  };

  await softStep('Pie chart: switching the category column drops the selection', async () => {
    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      tv.addViewer('Pie chart');
      await new Promise((r) => setTimeout(r, 1200));
      tv.viewers.find((x: any) => x.type === 'Pie chart').setOptions({categoryColumnName: 'DIS_POP'});
    });
    await page.waitForTimeout(2500);
    await clickItem('Indigestion');
    expect((await legendState()).selected.join()).toContain('Indigestion');

    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      tv.viewers.find((x: any) => x.type === 'Pie chart').setOptions({categoryColumnName: 'RACE'});
    });
    await page.waitForTimeout(2500);
    const after = await legendState();
    expect(after.labels.join()).toContain('Black');
    expect(after.selected, 'a positional ghost selection survived the column switch').toEqual([]);
    await page.evaluate(() => {
      (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Pie chart').close();
    });
    await page.waitForTimeout(1000);
  });

  await softStep('Line chart: selection survives a row filter, resets on a split switch', async () => {
    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      tv.addViewer('Line chart');
      await new Promise((r) => setTimeout(r, 1200));
      tv.viewers.find((x: any) => x.type === 'Line chart')
        .setOptions({yColumnNames: ['AGE'], splitColumnNames: ['RACE']});
    });
    await page.waitForTimeout(3000);
    await v.resizeViewer(page, 'Line chart', 900, 500);
    await page.waitForTimeout(1500);
    await clickItem('Black');
    expect((await legendState()).selected.join()).toContain('Black');

    // same source, different row mask: the selection must survive this rebuild
    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      tv.dataFrame.rows.filter((r: any) => r.idx < 3000);
      await new Promise((res) => setTimeout(res, 2000));
    });
    expect((await legendState()).selected.join(),
      'a row-filter rebuild must keep the selection').toContain('Black');

    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      tv.viewers.find((x: any) => x.type === 'Line chart').setOptions({splitColumnNames: ['SEX']});
    });
    await page.waitForTimeout(2500);
    const after = await legendState();
    expect(after.labels.join()).toContain('M');
    expect(after.selected, 'a positional ghost selection survived the split switch').toEqual([]);
    await page.evaluate(() => {
      (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Line chart').close();
    });
    await page.waitForTimeout(1000);
  });

  await softStep('Scatter plot: the DATA unfilters on a color column switch, not just the classes', async () => {
    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      tv.addViewer('Scatter plot');
      await new Promise((r) => setTimeout(r, 1200));
      tv.viewers.find((x: any) => x.type === 'Scatter plot')
        .setOptions({xColumnName: 'AGE', yColumnName: 'HEIGHT', colorColumnName: 'DIS_POP'});
    });
    await page.waitForTimeout(2500);
    const full = await countPixels('Scatter plot');
    await clickItem('Indigestion');
    const filtered = await countPixels('Scatter plot');
    expect(filtered, 'the legend click must filter the drawn points').toBeLessThan(full * 0.7);

    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      tv.viewers.find((x: any) => x.type === 'Scatter plot').setOptions({colorColumnName: 'RACE'});
    });
    await page.waitForTimeout(3000);
    const after = await countPixels('Scatter plot');
    expect((await legendState()).selected).toEqual([]);
    // same marker shapes, new palette: a correct repaint is within a few percent of the
    // unfiltered paint; the ghost-filter bug leaves it at the filtered count
    expect(after, 'the points stayed filtered by the previous color column')
      .toBeGreaterThan(full * 0.9);
    await page.evaluate(() => {
      (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot').close();
    });
    await page.waitForTimeout(1000);
  });

  await softStep('Box plot: marker-color switch resets the selection and its filter', async () => {
    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      tv.addViewer('Box plot');
      await new Promise((r) => setTimeout(r, 1200));
      const bp = tv.viewers.find((x: any) => x.type === 'Box plot');
      bp.setOptions({valueColumnName: 'HEIGHT', category1ColumnName: 'RACE'});
      await new Promise((r) => setTimeout(r, 1500));
      bp.setOptions({markerColorColumnName: 'DIS_POP'});
    });
    await page.waitForTimeout(2500);
    await v.resizeViewer(page, 'Box plot', 700, 450);
    await page.waitForTimeout(1500);
    for (let i = 0; i < 10 && (await legendState()).labels.length === 0; i++)
      await page.waitForTimeout(1000);

    const full = await countPixels('Box plot');
    await clickItem('Indigestion');
    expect((await legendState()).selected.join()).toContain('Indigestion');
    const filtered = await countPixels('Box plot');
    expect(filtered, 'the legend click must filter the drawn points').toBeLessThan(full);

    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      tv.viewers.find((x: any) => x.type === 'Box plot').setOptions({markerColorColumnName: 'SEX'});
    });
    await page.waitForTimeout(3000);
    const after = await legendState();
    expect(after.labels.join()).toContain('M');
    expect(after.selected, 'a positional ghost selection survived the color switch').toEqual([]);
    expect(await countPixels('Box plot'), 'the points stayed filtered by the previous color column')
      .toBeGreaterThan(filtered);
  });

  expect(errors, `page errors:\n${errors.join('\n')}`).toEqual([]);
  await v.cleanupShell(page);
  v.finishSpec('Legend selection-reset failures');
});
