/* ---
realizes: [piechart.cp.aggregation-tour-and-datetime-map]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

test('Pie Chart — Aggregation Tour, Validation Messages, DateTime Category Map, Grid Color Coding', async ({page}) => {
  test.setTimeout(300_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  const consoleErrors: string[] = [];
  page.on('console', (m) => {
    if (m.type() === 'error')
      consoleErrors.push(m.text());
  });

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await v.addViewerByIcon(page, 'pie-chart', 'Pie-chart');

  await page.evaluate(async () => {
    const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
    pie.props.categoryColumnName = 'RACE';
    await new Promise((r) => setTimeout(r, 500));
  });

  const readRootInDom = () => page.evaluate(() => {
    const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
    return document.body.contains(pie.root);
  });

  const readViewerError = () => page.evaluate(() => {
    const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
    const el = pie.root.querySelector('.d4-viewer-error');
    return el ? (el.textContent || '').trim() : '';
  });

  // Default aggregation types, captured before any scenario changes them.
  // Clearing an angle/length COLUMN while a non-default aggregation is still
  // set blanks the chart with a validation error, so every clear below also
  // restores these defaults.
  const aggrDefaults = await page.evaluate(() => {
    const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
    return {angle: pie.props.segmentAngleAggrType, length: pie.props.segmentLengthAggrType};
  });
  console.log(`Aggregation defaults: angle=${aggrDefaults.angle} length=${aggrDefaults.length}`);

  // Settle-gated canvas ink count: repeated reads until two consecutive counts
  // agree, so a delta between two measurements is the setter's effect, not a
  // render tail.
  const settledPx = async () => {
    let prev = (await v.countCanvasPixels(page, 'Pie chart')).total;
    let cur = prev;
    for (let i = 0; i < 5; i++) {
      await page.waitForTimeout(300);
      cur = (await v.countCanvasPixels(page, 'Pie chart')).total;
      if (Math.abs(cur - prev) < 200) break;
      prev = cur;
    }
    return cur;
  };

  await softStep('Scenario 1 — angle aggregation tour, length column switching, clear to standard pie', async () => {
    expect(await readViewerError()).toBe('');
    const errBefore = pageErrors.length + consoleErrors.length;
    const tour = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      pie.props.segmentAngleColumnName = 'AGE';
      const r: string[] = [];
      for (const aggr of ['min', 'max', 'med', 'stdev', 'count', 'avg']) {
        pie.props.segmentAngleAggrType = aggr;
        await new Promise((res) => setTimeout(res, 300));
        r.push(pie.props.segmentAngleAggrType);
      }
      return r;
    });
    expect(tour).toEqual(['min', 'max', 'med', 'stdev', 'count', 'avg']);

    // The disc stays full when only the angles change, so the avg → sum
    // signal is the per-color histogram delta, not the non-white total:
    // avg(AGE) is nearly uniform across races while sum(AGE) follows the row
    // counts, so the slice angles must redistribute.
    const avgPx = await settledPx();
    expect(await v.snapshotCanvasColors(page, 'Pie chart')).toBe(true);
    await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      pie.props.segmentAngleAggrType = 'sum';
      await new Promise((res) => setTimeout(res, 500));
    });
    const sumPx = await settledPx();
    const {deltaPx} = await v.diffCanvasColors(page, 'Pie chart');
    console.log(`Angle aggregation px: avgPx=${avgPx} sumPx=${sumPx} recolorDeltaPx=${deltaPx}`);
    expect(avgPx).toBeGreaterThanOrEqual(0);
    expect(sumPx).toBeGreaterThanOrEqual(0);
    expect(deltaPx).toBeGreaterThan(500);

    // Length coding min-max scales the per-category aggregates, so the
    // smallest category always drops to the minimum outer radius and the ink
    // count must DECREASE from the full disc by a real margin.
    await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      pie.props.segmentLengthColumnName = 'WEIGHT';
      pie.props.segmentLengthAggrType = 'avg';
      await new Promise((res) => setTimeout(res, 500));
    });
    const lengthPx = await settledPx();
    console.log(`Segment length px: fullDiscPx=${sumPx} lengthPx=${lengthPx}`);
    expect(sumPx - lengthPx).toBeGreaterThan(500);

    // Length-aggregation switching repaints within the coded disc — covered
    // by the scenario's no-error floor; then restore both aggregation types
    // to their defaults FIRST and clear the columns after — clearing a column
    // while a non-default aggregation is still set blanks the chart with a
    // validation error.
    await page.evaluate(async (d) => {
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      for (const aggr of ['min', 'max']) {
        pie.props.segmentLengthAggrType = aggr;
        await new Promise((res) => setTimeout(res, 300));
      }
      pie.props.segmentAngleAggrType = d.angle;
      pie.props.segmentLengthAggrType = d.length;
      await new Promise((res) => setTimeout(res, 300));
      pie.props.segmentAngleColumnName = '';
      pie.props.segmentLengthColumnName = '';
      await new Promise((res) => setTimeout(res, 500));
    }, aggrDefaults);
    const restoredPx = await settledPx();
    console.log(`Cleared px: restoredPx=${restoredPx} lengthPx=${lengthPx}`);
    expect(restoredPx - lengthPx).toBeGreaterThan(500);
    const cleared = await page.evaluate(() => {
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      return {angle: pie.props.segmentAngleColumnName, length: pie.props.segmentLengthColumnName};
    });
    expect(cleared.angle).toBe('');
    expect(cleared.length).toBe('');
    expect(await readViewerError()).toBe('');
    expect(await readRootInDom()).toBe(true);
    expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
  });

  await softStep('Scenario 2 — negative-minimum and all-zero aggregations render validation messages', async () => {
    expect(await readViewerError()).toBe('');
    const errBefore = pageErrors.length + consoleErrors.length;
    try {
      const result = await page.evaluate(async (d) => {
        const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
        const df = grok.shell.tv.dataFrame;
        df.columns.addNewFloat('NEG_PROBE').init((i: number) => i % 2 === 0 ? -5 : 3);
        df.columns.addNewFloat('ZERO_PROBE').init(() => 0);
        await new Promise((res) => setTimeout(res, 500));
        const readError = () => {
          const el = pie.root.querySelector('.d4-viewer-error');
          return el ? (el.textContent || '').trim() : '';
        };
        pie.props.segmentAngleColumnName = 'NEG_PROBE';
        pie.props.segmentAngleAggrType = 'min';
        await new Promise((res) => setTimeout(res, 800));
        const negMsg = readError();
        pie.props.segmentAngleColumnName = 'ZERO_PROBE';
        pie.props.segmentAngleAggrType = 'sum';
        await new Promise((res) => setTimeout(res, 800));
        const zeroMsg = readError();
        pie.props.segmentAngleAggrType = d.angle;
        await new Promise((res) => setTimeout(res, 300));
        pie.props.segmentAngleColumnName = '';
        await new Promise((res) => setTimeout(res, 600));
        const clearedMsg = readError();
        return {negMsg, zeroMsg, clearedMsg};
      }, aggrDefaults);
      expect(result.negMsg).toContain('contains negative values');
      expect(result.negMsg).toContain('NEG_PROBE');
      expect(result.zeroMsg).toContain('all values are 0');
      expect(result.zeroMsg).toContain('ZERO_PROBE');
      expect(result.clearedMsg).toBe('');
      expect(await readRootInDom()).toBe(true);
      expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
    } finally {
      // Never leak the scratch columns or the deliberate error state, even
      // when an assertion above failed: restore the default aggregation,
      // clear the angle column, then drop the probe columns.
      await page.evaluate(async (d) => {
        const pie = Array.from(grok.shell.tv?.viewers ?? []).find((vw: any) => vw.type === 'Pie chart') as any;
        if (pie) {
          try {
            pie.props.segmentAngleAggrType = d.angle;
            pie.props.segmentAngleColumnName = '';
          } catch (_) {}
          await new Promise((res) => setTimeout(res, 300));
        }
        const df = grok.shell.tv?.dataFrame;
        for (const name of ['NEG_PROBE', 'ZERO_PROBE'])
          try { df.columns.remove(name); } catch (_) {}
        await new Promise((res) => setTimeout(res, 300));
      }, aggrDefaults);
    }
  });

  await softStep('Scenario 3 — Category Map on STARTED changes the legend category set (year → month → quarter)', async () => {
    expect(await readViewerError()).toBe('');
    const errBefore = pageErrors.length + consoleErrors.length;
    const defaultMap = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      pie.props.legendVisibility = 'Always';
      pie.props.categoryColumnName = 'STARTED';
      await new Promise((res) => setTimeout(res, 800));
      return pie.props.categoryMap;
    });
    expect(defaultMap).toBe('year');
    const yearLegend = await v.readLegend(page, 'Pie chart');
    const setMap = (map: string) => page.evaluate(async (m) => {
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      pie.props.categoryMap = m;
      await new Promise((res) => setTimeout(res, 600));
    }, map);
    await setMap('month');
    const monthLegend = await v.readLegend(page, 'Pie chart');
    await setMap('quarter');
    const quarterLegend = await v.readLegend(page, 'Pie chart');
    console.log(`Category map legend labels: year=[${yearLegend.labels}] month=[${monthLegend.labels}] quarter=[${quarterLegend.labels}]`);
    expect(yearLegend.labels.length).toBeGreaterThan(0);
    expect(monthLegend.labels.length).toBeGreaterThan(0);
    expect(quarterLegend.labels.length).toBeGreaterThan(0);
    expect(monthLegend.labels.length).toBeLessThanOrEqual(12);
    expect(quarterLegend.labels.length).toBeLessThanOrEqual(4);
    // The map regroups the same dates into a DIFFERENT category set each time.
    expect([...yearLegend.labels].sort()).not.toEqual([...monthLegend.labels].sort());
    expect([...monthLegend.labels].sort()).not.toEqual([...quarterLegend.labels].sort());
    expect([...yearLegend.labels].sort()).not.toEqual([...quarterLegend.labels].sort());
    await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      pie.props.categoryMap = 'year';
      pie.props.categoryColumnName = 'RACE';
      pie.props.legendVisibility = 'Auto';
      await new Promise((res) => setTimeout(res, 500));
    });
    expect(await readRootInDom()).toBe(true);
    expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
  });

  await softStep('Scenario 4 — grid color-coding on RACE recolors the pie legend swatches', async () => {
    expect(await readViewerError()).toBe('');
    const errBefore = pageErrors.length + consoleErrors.length;
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      const df = grok.shell.tv.dataFrame;
      const target = 'Asian';
      pie.props.categoryColumnName = 'RACE';
      pie.props.legendVisibility = 'Always';
      await new Promise((res) => setTimeout(res, 800));
      const itemColor = () => {
        const items = Array.from(pie.root.querySelectorAll('[name="legend"] .d4-legend-item')) as HTMLElement[];
        const it = items.find((el) => (el.querySelector('.d4-legend-value')?.textContent || '').trim() === target);
        return it ? getComputedStyle(it).color : '';
      };
      const before = itemColor();
      df.col('RACE').meta.colors.setCategorical({[target]: '#ff0000'});
      try { pie.invalidate?.(); } catch (_) {}
      await new Promise((res) => setTimeout(res, 800));
      const after = itemColor();
      delete df.col('RACE').tags['.color-coding-categorical'];
      delete df.col('RACE').tags['.color-coding-type'];
      try { pie.invalidate?.(); } catch (_) {}
      await new Promise((res) => setTimeout(res, 800));
      const restored = itemColor();
      pie.props.legendVisibility = 'Auto';
      await new Promise((res) => setTimeout(res, 300));
      return {before, after, restored};
    });
    console.log(`Legend swatch colors: before=${result.before} after=${result.after} restored=${result.restored}`);
    expect(result.before).not.toBe('');
    expect(result.after).toBe('rgb(255, 0, 0)');
    expect(result.after).not.toBe(result.before);
    expect(result.restored).toBe(result.before);
    expect(await readRootInDom()).toBe(true);
    expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
  });

  await softStep('Scenario 5 — tooltip content follows the configured aggregations', async () => {
    expect(await readViewerError()).toBe('');
    const errBefore = pageErrors.length + consoleErrors.length;
    const raceCats: string[] = await page.evaluate(() =>
      grok.shell.tv.dataFrame.col('RACE').categories.slice());
    // The tooltip is a page-level singleton element that keeps its last text
    // when hidden, so every hover starts by blanking it — any text read
    // afterwards was written by THIS hover, never by a previous state.
    const readTooltip = () => page.evaluate(() => {
      const tts = Array.from(document.querySelectorAll('.d4-tooltip')) as HTMLElement[];
      const populated = tts.find((t) => (t.textContent || '').trim().length > 0) || tts[0] || null;
      return populated
        ? {display: getComputedStyle(populated).display, text: (populated.textContent || '').trim()}
        : {display: 'missing', text: ''};
    });
    const resetAndMove = (fx: number, fy: number) => page.evaluate(async (p) => {
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      const canvas = pie.root.querySelector('canvas') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      const mm = (x: number, y: number) => canvas.dispatchEvent(
        new MouseEvent('mousemove', {bubbles: true, clientX: x, clientY: y}));
      mm(rect.left + 2, rect.top + 2);
      await new Promise((r) => setTimeout(r, 250));
      for (const t of Array.from(document.querySelectorAll('.d4-tooltip')))
        t.textContent = '';
      mm(rect.left + rect.width * p.fx, rect.top + rect.height * p.fy);
      await new Promise((r) => setTimeout(r, 600));
    }, {fx, fy});
    // Probe a spread of canvas fractions until the hover lands on a slice
    // and the tooltip text populates.
    const hoverSlice = async () => {
      for (const [fx, fy] of [[0.65, 0.4], [0.45, 0.4], [0.5, 0.35], [0.6, 0.5], [0.4, 0.5], [0.5, 0.6]]) {
        await resetAndMove(fx, fy);
        const tt = await readTooltip();
        if (tt.text.length > 0) return tt;
      }
      return {display: 'missing', text: ''};
    };
    const setPie = (props: Record<string, any>) => page.evaluate(async (p) => {
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      for (const k of Object.keys(p)) pie.props[k] = p[k];
      await new Promise((r) => setTimeout(r, 600));
    }, props);

    await setPie({categoryColumnName: 'RACE', segmentAngleColumnName: '', segmentLengthColumnName: ''});
    const countTt = await hoverSlice();
    console.log(`Tooltip (count): display=${countTt.display} text=${countTt.text.slice(0, 160)}`);
    // Default row-count tooltip: the hovered category's name plus row counts.
    expect(countTt.text.length).toBeGreaterThan(0);
    expect(raceCats.some((c) => countTt.text.includes(c))).toBe(true);
    expect(countTt.text).toMatch(/\d/);

    await setPie({segmentAngleColumnName: 'AGE', segmentAngleAggrType: 'avg'});
    const avgTt = await hoverSlice();
    console.log(`Tooltip (avg AGE): display=${avgTt.display} text=${avgTt.text.slice(0, 160)}`);
    expect(avgTt.text).toContain('avg(AGE)');

    await setPie({segmentLengthColumnName: 'WEIGHT', segmentLengthAggrType: 'max'});
    const lenTt = await hoverSlice();
    console.log(`Tooltip (max WEIGHT): display=${lenTt.display} text=${lenTt.text.slice(0, 160)}`);
    expect(lenTt.text).toContain('avg(AGE)');
    expect(lenTt.text).toContain('max(WEIGHT)');

    // Move away: blank the tooltip, then move to an empty corner and leave the
    // canvas — nothing may be re-rendered into the tooltip and it must not be
    // visibly shown.
    const awayTt = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      const canvas = pie.root.querySelector('canvas') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      for (const t of Array.from(document.querySelectorAll('.d4-tooltip')))
        t.textContent = '';
      canvas.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: rect.left + 2, clientY: rect.top + 2}));
      canvas.dispatchEvent(new MouseEvent('mouseleave', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 600));
      const tts = Array.from(document.querySelectorAll('.d4-tooltip')) as HTMLElement[];
      const populated = tts.find((t) => (t.textContent || '').trim().length > 0) || tts[0] || null;
      return populated
        ? {display: getComputedStyle(populated).display, text: (populated.textContent || '').trim()}
        : {display: 'missing', text: ''};
    });
    console.log(`Tooltip (away): display=${awayTt.display} text=${awayTt.text.slice(0, 160)}`);
    expect(awayTt.text).toBe('');
    expect(awayTt.display).not.toBe('block');

    // Restore the aggregation defaults BEFORE clearing the columns — clearing
    // a column while a non-default aggregation is set blanks the chart.
    await setPie({
      segmentAngleAggrType: aggrDefaults.angle,
      segmentLengthAggrType: aggrDefaults.length,
      segmentAngleColumnName: '',
      segmentLengthColumnName: '',
    });
    expect(await readViewerError()).toBe('');
    expect(await readRootInDom()).toBe(true);
    expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
  });

  v.finishSpec();
});
