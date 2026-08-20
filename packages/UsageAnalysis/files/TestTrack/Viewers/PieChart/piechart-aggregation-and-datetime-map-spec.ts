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

  await v.installEventWaits(page);

  await page.evaluate(async () => {
    const w = window as any;
    const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
    await w.__settled('viewer:Pie chart.onViewerRendered', () => {
      pie.props.categoryColumnName = 'RACE';
    }, 2000);
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

  const aggrDefaults = await page.evaluate(() => {
    const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
    return {angle: pie.props.segmentAngleAggrType, length: pie.props.segmentLengthAggrType};
  });
  console.log(`Aggregation defaults: angle=${aggrDefaults.angle} length=${aggrDefaults.length}`);

  const settledPx = async () => {
    await v.waitForCanvasQuiet(page, 'Pie chart', {optional: true});
    return (await v.countCanvasPixels(page, 'Pie chart')).total;
  };

  await softStep('Scenario 1 — angle aggregation tour, length column switching, clear to standard pie', async () => {
    expect(await readViewerError()).toBe('');
    const errBefore = pageErrors.length + consoleErrors.length;
    const tour = await page.evaluate(async () => {
      const w = window as any;
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      pie.props.segmentAngleColumnName = 'AGE';
      const r: string[] = [];
      for (const aggr of ['min', 'max', 'med', 'stdev', 'count', 'avg']) {
        await w.__settled('viewer:Pie chart.onViewerRendered', () => {
          pie.props.segmentAngleAggrType = aggr;
        }, 2000);
        r.push(pie.props.segmentAngleAggrType);
      }
      return r;
    });
    expect(tour).toEqual(['min', 'max', 'med', 'stdev', 'count', 'avg']);

    const avgPx = await settledPx();
    expect(await v.snapshotCanvasColors(page, 'Pie chart')).toBe(true);
    await page.evaluate(async () => {
      const w = window as any;
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      await w.__settled('viewer:Pie chart.onViewerRendered', () => {
        pie.props.segmentAngleAggrType = 'sum';
      }, 2000);
    });
    const sumPx = await settledPx();
    const {deltaPx} = await v.diffCanvasColors(page, 'Pie chart');
    console.log(`Angle aggregation px: avgPx=${avgPx} sumPx=${sumPx} recolorDeltaPx=${deltaPx}`);
    expect(avgPx).toBeGreaterThanOrEqual(0);
    expect(sumPx).toBeGreaterThanOrEqual(0);
    expect(deltaPx).toBeGreaterThan(500);

    await page.evaluate(async () => {
      const w = window as any;
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      await w.__settled('viewer:Pie chart.onViewerRendered', () => {
        pie.props.segmentLengthColumnName = 'WEIGHT';
        pie.props.segmentLengthAggrType = 'avg';
      }, 2000);
    });
    const lengthPx = await settledPx();
    console.log(`Segment length px: fullDiscPx=${sumPx} lengthPx=${lengthPx}`);
    expect(sumPx - lengthPx).toBeGreaterThan(500);

    await page.evaluate(async (d) => {
      const w = window as any;
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      for (const aggr of ['min', 'max']) {
        await w.__settled('viewer:Pie chart.onViewerRendered', () => {
          pie.props.segmentLengthAggrType = aggr;
        }, 2000);
      }
      await w.__settled('viewer:Pie chart.onViewerRendered', () => {
        pie.props.segmentAngleAggrType = d.angle;
        pie.props.segmentLengthAggrType = d.length;
      }, 2000);
      await w.__settled('viewer:Pie chart.onViewerRendered', () => {
        pie.props.segmentAngleColumnName = '';
        pie.props.segmentLengthColumnName = '';
      }, 2000);
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
        const w = window as any;
        const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
        const df = grok.shell.tv.dataFrame;
        await w.__settled('viewer:Pie chart.onViewerRendered', () => {
          df.columns.addNewFloat('NEG_PROBE').init((i: number) => i % 2 === 0 ? -5 : 3);
          df.columns.addNewFloat('ZERO_PROBE').init(() => 0);
        }, 2000);
        const readError = () => {
          const el = pie.root.querySelector('.d4-viewer-error');
          return el ? (el.textContent || '').trim() : '';
        };
        await w.__settled('viewer:Pie chart.onViewerRendered', () => {
          pie.props.segmentAngleColumnName = 'NEG_PROBE';
          pie.props.segmentAngleAggrType = 'min';
        }, 2000);
        const negMsg = readError();
        await w.__settled('viewer:Pie chart.onViewerRendered', () => {
          pie.props.segmentAngleColumnName = 'ZERO_PROBE';
          pie.props.segmentAngleAggrType = 'sum';
        }, 2000);
        const zeroMsg = readError();
        await w.__settled('viewer:Pie chart.onViewerRendered', () => {
          pie.props.segmentAngleAggrType = d.angle;
        }, 2000);
        await w.__settled('viewer:Pie chart.onViewerRendered', () => {
          pie.props.segmentAngleColumnName = '';
        }, 2000);
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

      await page.evaluate(async (d) => {
        const w = window as any;
        const pie = Array.from(grok.shell.tv?.viewers ?? []).find((vw: any) => vw.type === 'Pie chart') as any;
        if (pie) {
          await w.__settled('viewer:Pie chart.onViewerRendered', () => {
            try {
              pie.props.segmentAngleAggrType = d.angle;
              pie.props.segmentAngleColumnName = '';
            } catch (_) {}
          }, 2000);
        }
        const df = grok.shell.tv?.dataFrame;
        await w.__settled('viewer:Pie chart.onViewerRendered', () => {
          for (const name of ['NEG_PROBE', 'ZERO_PROBE'])
            try { df.columns.remove(name); } catch (_) {}
        }, 2000);
      }, aggrDefaults);
    }
  });

  await softStep('Scenario 3 — Category Map on STARTED changes the legend category set (year → month → quarter)', async () => {
    expect(await readViewerError()).toBe('');
    const errBefore = pageErrors.length + consoleErrors.length;
    const defaultMap = await page.evaluate(async () => {
      const w = window as any;
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      await w.__settled('viewer:Pie chart.onViewerRendered', () => {
        pie.props.legendVisibility = 'Always';
        pie.props.categoryColumnName = 'STARTED';
      }, 2000);
      return pie.props.categoryMap;
    });
    expect(defaultMap).toBe('year');
    const yearLegend = await v.readLegend(page, 'Pie chart');
    const setMap = (map: string) => page.evaluate(async (m) => {
      const w = window as any;
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      await w.__settled('viewer:Pie chart.onViewerRendered', () => {
        pie.props.categoryMap = m;
      }, 2000);
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

    expect([...yearLegend.labels].sort()).not.toEqual([...monthLegend.labels].sort());
    expect([...monthLegend.labels].sort()).not.toEqual([...quarterLegend.labels].sort());
    expect([...yearLegend.labels].sort()).not.toEqual([...quarterLegend.labels].sort());
    await page.evaluate(async () => {
      const w = window as any;
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      await w.__settled('viewer:Pie chart.onViewerRendered', () => {
        pie.props.categoryMap = 'year';
        pie.props.categoryColumnName = 'RACE';
        pie.props.legendVisibility = 'Auto';
      }, 2000);
    });
    expect(await readRootInDom()).toBe(true);
    expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
  });

  await softStep('Scenario 4 — grid color-coding on RACE recolors the pie legend swatches', async () => {
    expect(await readViewerError()).toBe('');
    const errBefore = pageErrors.length + consoleErrors.length;
    const result = await page.evaluate(async () => {
      const w = window as any;
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      const df = grok.shell.tv.dataFrame;
      const target = 'Asian';
      await w.__settled('viewer:Pie chart.onViewerRendered', () => {
        pie.props.categoryColumnName = 'RACE';
        pie.props.legendVisibility = 'Always';
      }, 2000);
      const itemColor = () => {
        const items = Array.from(pie.root.querySelectorAll('[name="legend"] .d4-legend-item')) as HTMLElement[];
        const it = items.find((el) => (el.querySelector('.d4-legend-value')?.textContent || '').trim() === target);
        return it ? getComputedStyle(it).color : '';
      };
      const before = itemColor();
      await w.__settled('viewer:Pie chart.onViewerRendered', () => {
        df.col('RACE').meta.colors.setCategorical({[target]: '#ff0000'});
        try { pie.invalidate?.(); } catch (_) {}
      }, 2000);
      const after = itemColor();
      await w.__settled('viewer:Pie chart.onViewerRendered', () => {
        delete df.col('RACE').tags['.color-coding-categorical'];
        delete df.col('RACE').tags['.color-coding-type'];
        try { pie.invalidate?.(); } catch (_) {}
      }, 2000);
      const restored = itemColor();
      await w.__settled('viewer:Pie chart.onViewerRendered', () => {
        pie.props.legendVisibility = 'Auto';
      }, 2000);
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

    const readTooltip = () => page.evaluate(() => {
      const tts = Array.from(document.querySelectorAll('.d4-tooltip')) as HTMLElement[];
      const populated = tts.find((t) => (t.textContent || '').trim().length > 0) || tts[0] || null;
      return populated
        ? {display: getComputedStyle(populated).display, text: (populated.textContent || '').trim()}
        : {display: 'missing', text: ''};
    });
    const resetAndMove = (fx: number, fy: number) => page.evaluate(async (p) => {
      const w = window as any;
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      const canvas = pie.root.querySelector('canvas') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      const mm = (x: number, y: number) => canvas.dispatchEvent(
        new MouseEvent('mousemove', {bubbles: true, clientX: x, clientY: y}));
      await w.__settled('grok.events.onTooltipClosed', () => mm(rect.left + 2, rect.top + 2), 250);
      for (const t of Array.from(document.querySelectorAll('.d4-tooltip')))
        t.textContent = '';
      await w.__settled('grok.events.onTooltipShown',
        () => mm(rect.left + rect.width * p.fx, rect.top + rect.height * p.fy), 2000);
    }, {fx, fy});

    const hoverSlice = async () => {
      for (const [fx, fy] of [[0.65, 0.4], [0.45, 0.4], [0.5, 0.35], [0.6, 0.5], [0.4, 0.5], [0.5, 0.6]]) {
        await resetAndMove(fx, fy);
        const tt = await readTooltip();
        if (tt.text.length > 0) return tt;
      }
      return {display: 'missing', text: ''};
    };
    const setPie = (props: Record<string, any>) => page.evaluate(async (p) => {
      const w = window as any;
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      await w.__settled('viewer:Pie chart.onViewerRendered', () => {
        for (const k of Object.keys(p)) pie.props[k] = p[k];
      }, 2000);
    }, props);

    await setPie({categoryColumnName: 'RACE', segmentAngleColumnName: '', segmentLengthColumnName: ''});
    const countTt = await hoverSlice();
    console.log(`Tooltip (count): display=${countTt.display} text=${countTt.text.slice(0, 160)}`);

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

    const awayTt = await page.evaluate(async () => {
      const w = window as any;
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      const canvas = pie.root.querySelector('canvas') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      for (const t of Array.from(document.querySelectorAll('.d4-tooltip')))
        t.textContent = '';
      await w.__settled('grok.events.onTooltipClosed', () => {
        canvas.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: rect.left + 2, clientY: rect.top + 2}));
        canvas.dispatchEvent(new MouseEvent('mouseleave', {bubbles: true}));
      }, 2000);
      const tts = Array.from(document.querySelectorAll('.d4-tooltip')) as HTMLElement[];
      const populated = tts.find((t) => (t.textContent || '').trim().length > 0) || tts[0] || null;
      return populated
        ? {display: getComputedStyle(populated).display, text: (populated.textContent || '').trim()}
        : {display: 'missing', text: ''};
    });
    console.log(`Tooltip (away): display=${awayTt.display} text=${awayTt.text.slice(0, 160)}`);
    expect(awayTt.text).toBe('');
    expect(awayTt.display).not.toBe('block');

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
