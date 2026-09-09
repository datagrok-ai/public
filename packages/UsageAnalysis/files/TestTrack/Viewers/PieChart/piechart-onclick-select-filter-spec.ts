/* ---
realizes: [piechart.cp.segment-click-select-filter, piechart.int.segment-click-select-syncs-dataframe, piechart.int.segment-click-filter-multi-category]
--- */
import {localTest as test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

test('Pie Chart — Segment Click Select and Filter Modes', async ({page}) => {
  test.setTimeout(300_000);

  await openDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await v.addViewerByIcon(page, 'pie-chart', 'Pie-chart');

  await v.installEventWaits(page);

  await page.evaluate(() => {
    const w = window as any;
    const pieCanvas = () => (document.querySelector('[name="viewer-Pie-chart"]') as HTMLElement)
      .querySelector('canvas') as HTMLCanvasElement;
    w.__pieRect = () => pieCanvas().getBoundingClientRect();
    w.__pieClickAt = async (x: number, y: number, ctrl: boolean, channel = 'df.onSelectionChanged') => {
      const canvas = pieCanvas();
      const o = {bubbles: true, clientX: x, clientY: y, ctrlKey: ctrl};
      await w.__settled(channel, () => {
        canvas.dispatchEvent(new MouseEvent('mousedown', o));
        canvas.dispatchEvent(new MouseEvent('mouseup', o));
        canvas.dispatchEvent(new MouseEvent('click', o));
      }, 2000);
    };
    w.__pieSpans = (colName: string) => {
      const pie = Array.from(w.grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      const df = w.grok.shell.tv.dataFrame;
      const c = df.col(colName);
      const counts = new Map<any, number>();
      for (let i = 0; i < df.rowCount; i++) {
        if (c.isNone(i) && !pie.props.includeNulls) continue;
        const val = c.isNone(i) ? null : c.get(i);
        counts.set(val, (counts.get(val) || 0) + 1);
      }
      const order = Array.from(counts.entries()).sort((a, b) => a[1] - b[1]);
      const total = order.reduce((s, e) => s + e[1], 0);
      let acc = ((pie.props.startAngle || 0) / 180) * Math.PI;
      const spans = new Map<any, number[]>();
      for (const [cat, n] of order) {
        const sw = (n / total) * 2 * Math.PI;
        spans.set(cat, [acc, acc + sw]);
        acc += sw;
      }
      return spans;
    };

    // Slice geometry is derived, not read back: the sweep is walked at several radii and
    // angles until a click lands rows of exactly one category, which pins the hit point.
    w.__pieCalibrate = async (colName: string, targets: any[]) => {
      const pie = Array.from(w.grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      const df = w.grok.shell.tv.dataFrame;
      const c = df.col(colName);
      const rect = w.__pieRect();
      const radius = Math.min(pie.props.maxRadius || 1e9, Math.min(rect.width, rect.height) / 2);
      const cx = rect.left + rect.width / 2;
      const cy = rect.top + rect.height / 2;
      const spans = w.__pieSpans(colName);
      const pts: any = {};
      for (const t of targets) {
        const [a0, a1] = spans.get(t);
        outer:
        for (const rf of [0.6, 0.75, 0.45]) {
          for (const af of [0.5, 0.35, 0.65, 0.2, 0.8]) {
            const a = a0 + (a1 - a0) * af;
            const x = cx + rf * radius * Math.cos(a);
            const y = cy + rf * radius * Math.sin(a);
            df.selection.setAll(false);
            await w.__pieClickAt(x, y, false);
            if (df.selection.trueCount === 0) continue;
            let ok = true;
            for (let i = 0; i < df.rowCount; i++) {
              if (!df.selection.get(i)) continue;
              if (t === null ? !c.isNone(i) : (c.isNone(i) || c.get(i) !== t)) {
                ok = false;
                break;
              }
            }
            if (ok) {
              pts[t === null ? '__null__' : t] = {x, y};
              break outer;
            }
          }
        }
      }
      df.selection.setAll(false);
      return pts;
    };
  });

  await page.evaluate(async () => {
    const w = window as any;
    const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
    await w.__settled('viewer:Pie chart.onViewerRendered', () => { pie.props.rowSource = 'All'; }, 2000);
  });

  await softStep('Select mode — slice click selects exactly the category row count, another slice switches the selection, clearing returns to zero', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      const df = grok.shell.tv.dataFrame;
      const race = df.col('RACE');
      await w.__settled('viewer:Pie chart.onViewerRendered', () => {
        pie.props.categoryColumnName = 'RACE';
        pie.props.onClick = 'Select';
      }, 2000);
      const countOf = (cat: string) => {
        let n = 0;
        for (let i = 0; i < df.rowCount; i++)
          if (!race.isNone(i) && race.get(i) === cat) n++;
        return n;
      };
      const selectedCats = () => {
        const s = new Set<string>();
        for (let i = 0; i < df.rowCount; i++)
          if (df.selection.get(i)) s.add(race.get(i));
        return Array.from(s);
      };

      const byCount = Array.from((w.__pieSpans('RACE') as Map<any, number[]>).keys());
      byCount.sort((a: string, b: string) => countOf(b) - countOf(a));
      const t1 = byCount[0];
      const t2 = byCount[1];
      const pts = await w.__pieCalibrate('RACE', [t1, t2]);
      const blank = {
        calibrated: false, t1, t2, sel1: -1, cats1: [] as string[], expected1: -2,
        sel2: -1, cats2: [] as string[], expected2: -2, cleared: -1,
      };
      if (!pts[t1] || !pts[t2])
        return blank;

      await w.__pieClickAt(pts[t1].x, pts[t1].y, false);
      const sel1 = df.selection.trueCount;
      const cats1 = selectedCats();
      const expected1 = countOf(t1);

      await w.__pieClickAt(pts[t2].x, pts[t2].y, false);
      const sel2 = df.selection.trueCount;
      const cats2 = selectedCats();
      const expected2 = countOf(t2);

      await w.__settled('df.onSelectionChanged', () => df.selection.setAll(false), 2000);
      const cleared = df.selection.trueCount;
      return {calibrated: true, t1, t2, sel1, cats1, expected1, sel2, cats2, expected2, cleared};
    });
    expect(result.calibrated).toBe(true);
    expect(result.cats1).toEqual([result.t1]);
    expect(result.sel1).toBeGreaterThan(0);
    expect(result.sel1).toBe(result.expected1);
    expect(result.cats2).toEqual([result.t2]);
    expect(result.sel2).toBe(result.expected2);
    expect(result.cleared).toBe(0);
  });

  await softStep('Missing-values slice — clicking it selects exactly the rows with missing values, Include Nulls off repaints the pie without the slice, Include Nulls on restores it', async () => {
    const prepared = await page.evaluate(async () => {
      const w = window as any;
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      const df = grok.shell.tv.dataFrame;
      const src = df.col('RACE');
      const gapsCol = df.columns.addNewString('RACE_GAPS');
      for (let i = 0; i < df.rowCount; i++)
        gapsCol.set(i, i % 10 === 0 ? null : src.get(i), false);
      await w.__settled('viewer:Pie chart.onViewerRendered', () => {
        pie.props.categoryColumnName = 'RACE_GAPS';
        pie.props.onClick = 'Select';
        pie.props.includeNulls = true;
      }, 3000);
      const missing = gapsCol.stats.missingValueCount;

      const pts = await w.__pieCalibrate('RACE_GAPS', [null]);
      if (!pts['__null__'])
        return {calibrated: false, missing, nullSel: -1, allNull: false};
      await w.__pieClickAt(pts['__null__'].x, pts['__null__'].y, false);
      const nullSel = df.selection.trueCount;
      let allNull = nullSel > 0;
      for (let i = 0; i < df.rowCount; i++) {
        if (df.selection.get(i) && !gapsCol.isNone(i)) {
          allNull = false;
          break;
        }
      }
      await w.__settled('df.onSelectionChanged', () => df.selection.setAll(false), 2000);
      return {calibrated: true, missing, nullSel, allNull};
    });

    let offDelta = -1;
    let reSel = -1;
    if (prepared.calibrated) {
      await v.snapshotCanvasColors(page, 'Pie chart');
      await page.evaluate(async () => {
        const w = window as any;
        const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
        await w.__settled('viewer:Pie chart.onViewerRendered', () => { pie.props.includeNulls = false; }, 3000);
      });
      offDelta = (await v.diffCanvasColors(page, 'Pie chart')).deltaPx;

      reSel = await page.evaluate(async () => {
        const w = window as any;
        const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
        const df = grok.shell.tv.dataFrame;
        await w.__settled('viewer:Pie chart.onViewerRendered', () => { pie.props.includeNulls = true; }, 3000);
        const pts2 = await w.__pieCalibrate('RACE_GAPS', [null]);
        let n = -1;
        if (pts2['__null__']) {
          await w.__pieClickAt(pts2['__null__'].x, pts2['__null__'].y, false);
          n = df.selection.trueCount;
        }
        await w.__settled('df.onSelectionChanged', () => df.selection.setAll(false), 2000);
        return n;
      });
    }

    await page.evaluate(async () => {
      const w = window as any;
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      await w.__settled('viewer:Pie chart.onViewerRendered', () => {
        pie.props.categoryColumnName = 'RACE';
        if (df.col('RACE_GAPS'))
          df.columns.remove('RACE_GAPS');
      }, 2000);
    });

    expect(prepared.calibrated).toBe(true);
    expect(prepared.missing).toBeGreaterThan(0);
    expect(prepared.allNull).toBe(true);
    expect(prepared.nullSel).toBe(prepared.missing);
    expect(offDelta).toBeGreaterThan(1000);
    expect(reSel).toBe(prepared.missing);
  });

  await softStep('Filter mode — slice click filters to the category row count, Ctrl+click adds a second category to the sum, empty-area click restores the full count', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      const df = grok.shell.tv.dataFrame;
      const race = df.col('RACE');
      await w.__settled('viewer:Pie chart.onViewerRendered', () => {
        pie.props.categoryColumnName = 'RACE';
        pie.props.onClick = 'Select';
      }, 2000);
      const countOf = (cat: string) => {
        let n = 0;
        for (let i = 0; i < df.rowCount; i++)
          if (!race.isNone(i) && race.get(i) === cat) n++;
        return n;
      };
      const filteredCats = () => {
        const s = new Set<string>();
        for (let i = 0; i < df.rowCount; i++)
          if (df.filter.get(i)) s.add(race.get(i));
        return Array.from(s);
      };
      const byCount = Array.from((w.__pieSpans('RACE') as Map<any, number[]>).keys());
      byCount.sort((a: string, b: string) => countOf(b) - countOf(a));
      const t1 = byCount[0];
      const t2 = byCount[1];
      const pts = await w.__pieCalibrate('RACE', [t1, t2]);
      const blank = {
        calibrated: false, t1, t2, full: -1, rowCount: df.rowCount, f1: -1,
        cats1: [] as string[], expected1: -2, f2: -1, cats2: [] as string[],
        expected2: -2, cleared: -1,
      };
      if (!pts[t1] || !pts[t2])
        return blank;

      await w.__settled('viewer:Pie chart.onViewerRendered', () => { pie.props.onClick = 'Filter'; }, 2000);
      const full = df.filter.trueCount;

      await w.__pieClickAt(pts[t1].x, pts[t1].y, false, 'df.onRowsFiltered');
      const f1 = df.filter.trueCount;
      const cats1 = filteredCats();
      const expected1 = countOf(t1);

      await w.__pieClickAt(pts[t2].x, pts[t2].y, true, 'df.onRowsFiltered');
      const f2 = df.filter.trueCount;
      const cats2 = filteredCats();
      const expected2 = countOf(t1) + countOf(t2);

      const rect = w.__pieRect();
      await w.__pieClickAt(rect.left + 5, rect.top + 5, false, 'df.onRowsFiltered');
      const cleared = df.filter.trueCount;
      pie.props.onClick = 'Select';
      return {calibrated: true, t1, t2, full, rowCount: df.rowCount, f1, cats1, expected1, f2, cats2, expected2, cleared};
    });
    expect(result.calibrated).toBe(true);
    expect(result.full).toBe(result.rowCount);
    expect(result.cats1).toEqual([result.t1]);
    expect(result.f1).toBe(result.expected1);
    expect(result.f1).toBeLessThan(result.full);
    expect(result.cats2.length).toBe(2);
    expect(result.f2).toBe(result.expected2);
    expect(result.cleared).toBe(result.full);
  });

  await softStep('Filter Panel composition — pie click narrows the panel-filtered count further, clearing the pie part returns to the panel-only value, Reset filters restores the full count', async () => {
    await page.evaluate(() => grok.shell.tv.getFiltersGroup());
    await page.locator('.d4-filter-group-header').waitFor({timeout: 15000});
    const calibration = await page.evaluate(async () => {
      const w = window as any;
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      const df = grok.shell.tv.dataFrame;
      const race = df.col('RACE');
      await w.__settled('viewer:Pie chart.onViewerRendered', () => {
        pie.props.categoryColumnName = 'RACE';
        pie.props.onClick = 'Select';
      }, 2000);
      const countOf = (cat: string) => {
        let n = 0;
        for (let i = 0; i < df.rowCount; i++)
          if (!race.isNone(i) && race.get(i) === cat) n++;
        return n;
      };

      const byCount = Array.from((w.__pieSpans('RACE') as Map<any, number[]>).keys());
      byCount.sort((a: string, b: string) => countOf(b) - countOf(a));
      const t1 = byCount[0];
      const pts = await w.__pieCalibrate('RACE', [t1]);
      if (!pts[t1])
        return {calibrated: false, full: -1, t1, pt: {x: -1, y: -1}};

      await w.__settled('viewer:Pie chart.onViewerRendered', () => { pie.props.onClick = 'Filter'; }, 2000);
      return {calibrated: true, full: df.filter.trueCount, t1, pt: pts[t1]};
    });

    const afterPanel = calibration.calibrated ? await v.applyNumericFilter(page, 'AGE', 30, 50) : -1;

    const result = await page.evaluate(async ({calibrated, t1, pt}) => {
      if (!calibrated)
        return {afterBoth: -1, expectedBoth: -2, backToPanel: -1, afterReset: -1};
      const w = window as any;
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      const df = grok.shell.tv.dataFrame;
      const race = df.col('RACE');
      const panelMask: boolean[] = [];
      for (let i = 0; i < df.rowCount; i++)
        panelMask.push(df.filter.get(i));

      await w.__pieClickAt(pt.x, pt.y, false, 'df.onRowsFiltered');
      const afterBoth = df.filter.trueCount;
      let expectedBoth = 0;
      for (let i = 0; i < df.rowCount; i++)
        if (panelMask[i] && race.get(i) === t1) expectedBoth++;

      const rect = w.__pieRect();
      await w.__pieClickAt(rect.left + 5, rect.top + 5, false, 'df.onRowsFiltered');
      const backToPanel = df.filter.trueCount;

      const btn = document.querySelector('.d4-filter-group-header [name="icon-arrow-rotate-left"]') as HTMLElement | null;
      await w.__settled('df.onRowsFiltered', () => btn?.click(), 3000);
      const afterReset = df.filter.trueCount;
      pie.props.onClick = 'Select';
      return {afterBoth, expectedBoth, backToPanel, afterReset};
    }, {calibrated: calibration.calibrated, t1: calibration.t1, pt: calibration.pt});

    expect(calibration.calibrated).toBe(true);
    expect(afterPanel).toBeLessThan(calibration.full);
    expect(result.afterBoth).toBeLessThan(afterPanel);
    expect(result.afterBoth).toBe(result.expectedBoth);
    expect(result.backToPanel).toBe(afterPanel);
    expect(result.afterReset).toBe(calibration.full);
  });

  await softStep('Viewer close releases the filter — closing the pie chart with an active click-filter restores the full count, re-adding the viewer succeeds', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      const df = grok.shell.tv.dataFrame;
      const race = df.col('RACE');
      await w.__settled('viewer:Pie chart.onViewerRendered', () => {
        pie.props.categoryColumnName = 'RACE';
        pie.props.onClick = 'Select';
      }, 2000);
      const countOf = (cat: string) => {
        let n = 0;
        for (let i = 0; i < df.rowCount; i++)
          if (!race.isNone(i) && race.get(i) === cat) n++;
        return n;
      };
      const byCount = Array.from((w.__pieSpans('RACE') as Map<any, number[]>).keys());
      byCount.sort((a: string, b: string) => countOf(b) - countOf(a));
      const t1 = byCount[0];
      const pts = await w.__pieCalibrate('RACE', [t1]);
      if (!pts[t1])
        return {calibrated: false, rowCount: df.rowCount, active: -1, restored: -1};

      await w.__settled('viewer:Pie chart.onViewerRendered', () => { pie.props.onClick = 'Filter'; }, 2000);
      await w.__pieClickAt(pts[t1].x, pts[t1].y, false, 'df.onRowsFiltered');
      const active = df.filter.trueCount;

      await w.__settled('df.onRowsFiltered', () => pie.close(), 3000);
      const restored = df.filter.trueCount;
      return {calibrated: true, rowCount: df.rowCount, active, restored};
    });
    expect(result.calibrated).toBe(true);
    expect(result.active).toBeLessThan(result.rowCount);
    expect(result.active).toBeGreaterThan(0);
    expect(result.restored).toBe(result.rowCount);

    await v.addViewerByIcon(page, 'pie-chart', 'Pie-chart');
    const after = await page.evaluate(() => ({
      present: !!Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart'),
      fullCount: grok.shell.tv.dataFrame.filter.trueCount,
      rowCount: grok.shell.tv.dataFrame.rowCount,
    }));
    expect(after.present).toBe(true);
    expect(after.fullCount).toBe(after.rowCount);
  });

  v.finishSpec();
});
