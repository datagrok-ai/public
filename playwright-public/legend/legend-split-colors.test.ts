/* ---
sub_features_covered: [legend.line-chart.split-colors-after-filter, legend.line-chart.custom-item-color-persistence]
--- */
// A line chart split by two columns colors its lines by the PACKED category index — the
// position among combinations that survive the filter. The legend must use the same
// packing: the renderer's categoryToColor cache is keyed by the old packing and is reset
// whenever the packed mapping recomputes, or the legend keeps the previous filter's
// colors while the lines repaint (regression vs the pre-refactor legend).

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

test('Line chart split legend colors follow the filter (demog)', async ({page}) => {
  test.setTimeout(900_000);
  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await page.evaluate(async () => {
    const grok = (window as any).grok;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    grok.shell.addTableView(df);
  });
  await page.waitForTimeout(2500);

  await softStep('After filtering to one RACE the legend colors match the repainted lines', async () => {
    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      tv.addViewer('Line chart');
      await new Promise((r) => setTimeout(r, 1200));
      const lc = tv.viewers.find((x: any) => x.type === 'Line chart');
      lc.setOptions({yColumnNames: ['AGE'], splitColumnNames: ['RACE', 'SEVERITY']});
    });
    await page.waitForTimeout(3000);
    await v.resizeViewer(page, 'Line chart', 900, 600);
    await page.waitForTimeout(1500);

    const before = await page.evaluate(() => {
      const root = document.querySelector('[name="legend"]') as HTMLElement;
      return root?.querySelectorAll('[name="legend-item"]').length ?? 0;
    });
    expect(before, 'no split legend before the filter').toBeGreaterThan(10);

    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const race = tv.dataFrame.col('RACE');
      tv.dataFrame.rows.filter((r: any) => race.get(r.idx) === 'Other');
    });
    await page.waitForTimeout(2500);

    const after = await page.evaluate(() => {
      const DG = (window as any).DG;
      const root = document.querySelector('[name="legend"]') as HTMLElement;
      const items = Array.from(root.querySelectorAll('[name="legend-item"]')) as HTMLElement[];
      return {colors: items.map((i) => getComputedStyle(i).color),
        labels: items.map((i) => (i.textContent || '').slice(0, 20)),
        palette: items.map((_, i) => DG.Color.toRgb(DG.Color.getCategoricalColor(i))
          .replace(/,/g, ', '))};
    });
    expect(after.colors.length, 'the legend must show only the filtered combinations')
      .toBeLessThan(before);
    expect(after.colors.length).toBeGreaterThan(2);
    // the lines repaint by packed index (position among surviving combos) — the legend
    // must show the same palette sequence, not the pre-filter colors
    for (let i = 0; i < after.colors.length; i++)
      expect(after.colors[i],
        `item ${after.labels[i]} keeps a stale pre-filter color`).toBe(after.palette[i]);
  });

  expect(errors, `page errors:\n${errors.join('\n')}`).toEqual([]);
  await v.cleanupShell(page);
  v.finishSpec('Line chart split-color failures');
});

test('User-picked split item color survives legend re-render (demog)', async ({page}) => {
  test.setTimeout(900_000);
  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await page.evaluate(async () => {
    const grok = (window as any).grok;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    grok.shell.addTableView(df);
  });
  await page.waitForTimeout(2500);

  await page.evaluate(async () => {
    const tv = (window as any).grok.shell.tv;
    tv.addViewer('Line chart');
    await new Promise((r) => setTimeout(r, 1200));
    const lc = tv.viewers.find((x: any) => x.type === 'Line chart');
    lc.setOptions({yColumnNames: ['AGE'], splitColumnNames: ['RACE', 'SEVERITY', 'SEX']});
  });
  await page.waitForTimeout(3000);
  await v.resizeViewer(page, 'Line chart', 900, 420);
  await page.waitForTimeout(1500);

  const readItemColor = (label: string) => page.evaluate((l) => {
    const root = document.querySelector('[name="legend"]') as HTMLElement;
    const items = Array.from(root?.querySelectorAll('[name="legend-item"]') ?? []) as HTMLElement[];
    const it = items.find((i) => i.textContent === l);
    return it ? getComputedStyle(it).color : 'ITEM NOT FOUND';
  }, label);

  let label = '';
  let swatchColor = '';
  await softStep('The picker recolors the item through the legend UI', async () => {
    const item = await page.evaluate(() => {
      const root = document.querySelector('[name="legend"]') as HTMLElement;
      const items = Array.from(root.querySelectorAll('[name="legend-item"]')) as HTMLElement[];
      const it = items[1];
      const r = it.getBoundingClientRect();
      return {count: items.length, label: it.textContent!, color: getComputedStyle(it).color,
        x: r.x + r.width / 2, y: r.y + r.height / 2};
    });
    label = item.label;
    expect(item.count, 'the three-split legend must be long enough to scroll').toBeGreaterThan(20);

    await page.mouse.move(item.x, item.y);
    await page.waitForTimeout(600);
    const icon = page.locator('[name="legend-icon-color-picker"]');
    await expect(icon).toBeVisible({timeout: 5000});
    await icon.click();
    await page.waitForTimeout(800);

    const swatch = page.locator('.d4-color-bar').nth(3);
    await expect(swatch).toBeVisible({timeout: 5000});
    swatchColor = await swatch.evaluate((el) => getComputedStyle(el).backgroundColor);
    expect(swatchColor).not.toBe(item.color);
    await swatch.click();
    await page.waitForTimeout(500);
    await page.locator('.d4-dialog button', {hasText: /^OK$/}).first().click();
    await page.waitForTimeout(1000);
    expect(await readItemColor(label)).toBe(swatchColor);
  });

  await softStep('The color survives a viewer resize', async () => {
    await v.resizeViewer(page, 'Line chart', 850, 400);
    await page.waitForTimeout(2000);
    expect(await readItemColor(label)).toBe(swatchColor);
  });

  await softStep('The color survives the virtualized list re-rendering on scroll', async () => {
    const result = await page.evaluate(async () => {
      const root = document.querySelector('[name="legend"]') as HTMLElement;
      const scroller = Array.from(root.querySelectorAll('div'))
        .find((d) => d.scrollHeight > d.clientHeight + 5) as HTMLElement;
      if (!scroller) return 'NO SCROLLER';
      scroller.scrollTop = 200;
      scroller.dispatchEvent(new Event('scroll'));
      await new Promise((r) => setTimeout(r, 600));
      scroller.scrollTop = 0;
      scroller.dispatchEvent(new Event('scroll'));
      await new Promise((r) => setTimeout(r, 600));
      return 'scrolled';
    });
    expect(result, 'the legend list must actually be scrollable').toBe('scrolled');
    expect(await readItemColor(label)).toBe(swatchColor);
  });

  await softStep('The color survives a legend model rebuild', async () => {
    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const sex = tv.dataFrame.col('SEX');
      tv.dataFrame.rows.filter((r: any) => sex.get(r.idx) !== null);
      await new Promise((r) => setTimeout(r, 1500));
    });
    expect(await readItemColor(label)).toBe(swatchColor);
  });

  expect(errors, `page errors:\n${errors.join('\n')}`).toEqual([]);
  await v.cleanupShell(page);
  v.finishSpec('Custom split-color persistence failures');
});
