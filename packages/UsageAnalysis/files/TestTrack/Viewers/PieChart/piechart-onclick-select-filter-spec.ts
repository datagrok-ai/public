/* ---
realizes: [piechart.cp.segment-click-select-filter, piechart.int.segment-click-select-syncs-dataframe, piechart.int.segment-click-filter-multi-category]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

test('Pie Chart — Segment Click Select and Filter Modes', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await v.addViewerByIcon(page, 'pie-chart', 'Pie-chart');

  await page.evaluate(() => {
    const w = window as any;
    const pieCanvas = () => (document.querySelector('[name="viewer-Pie-chart"]') as HTMLElement)
      .querySelector('canvas') as HTMLCanvasElement;
    w.__pieRect = () => pieCanvas().getBoundingClientRect();
    w.__pieClickAt = async (x: number, y: number, ctrl: boolean) => {
      const canvas = pieCanvas();
      const o = {bubbles: true, clientX: x, clientY: y, ctrlKey: ctrl};
      canvas.dispatchEvent(new MouseEvent('mousedown', o));
      canvas.dispatchEvent(new MouseEvent('mouseup', o));
      canvas.dispatchEvent(new MouseEvent('click', o));
      await new Promise((r) => setTimeout(r, 400));
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
    const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
    pie.props.rowSource = 'All';
    await new Promise((r) => setTimeout(r, 400));
  });

  await softStep('Select mode — slice click selects exactly the category row count, another slice switches the selection, clearing returns to zero', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      const df = grok.shell.tv.dataFrame;
      const race = df.col('RACE');
      pie.props.categoryColumnName = 'RACE';
      pie.props.onClick = 'Select';
      await new Promise((r) => setTimeout(r, 500));
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
      if (!pts[t1] || !pts[t2])
        return {calibrated: false};

      await w.__pieClickAt(pts[t1].x, pts[t1].y, false);
      const sel1 = df.selection.trueCount;
      const cats1 = selectedCats();
      const expected1 = countOf(t1);

      await w.__pieClickAt(pts[t2].x, pts[t2].y, false);
      const sel2 = df.selection.trueCount;
      const cats2 = selectedCats();
      const expected2 = countOf(t2);

      df.selection.setAll(false);
      await new Promise((r) => setTimeout(r, 300));
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
    const result = await page.evaluate(async () => {
      const w = window as any;
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      const df = grok.shell.tv.dataFrame;
      try {

        const src = df.col('RACE');
        const gapsCol = df.columns.addNewString('RACE_GAPS');
        for (let i = 0; i < df.rowCount; i++)
          gapsCol.set(i, i % 10 === 0 ? null : src.get(i), false);
        pie.props.categoryColumnName = 'RACE_GAPS';
        pie.props.onClick = 'Select';
        pie.props.includeNulls = true;
        await new Promise((r) => setTimeout(r, 800));
        const missing = gapsCol.stats.missingValueCount;

        const pts = await w.__pieCalibrate('RACE_GAPS', [null]);
        if (!pts['__null__'])
          return {calibrated: false};
        await w.__pieClickAt(pts['__null__'].x, pts['__null__'].y, false);
        const nullSel = df.selection.trueCount;
        let allNull = nullSel > 0;
        for (let i = 0; i < df.rowCount; i++) {
          if (df.selection.get(i) && !gapsCol.isNone(i)) {
            allNull = false;
            break;
          }
        }
        df.selection.setAll(false);

        const canvasHist = () => {
          const canvas = (document.querySelector('[name="viewer-Pie-chart"]') as HTMLElement)
            .querySelector('canvas') as HTMLCanvasElement;
          const d = canvas.getContext('2d')!.getImageData(0, 0, canvas.width, canvas.height).data;
          const m = new Map<number, number>();
          for (let i = 0; i < d.length; i += 4) {
            const k = (d[i] << 16) | (d[i + 1] << 8) | d[i + 2];
            m.set(k, (m.get(k) || 0) + 1);
          }
          return m;
        };
        const histDiff = (a: Map<number, number>, b: Map<number, number>) => {
          let delta = 0;
          for (const [k, n] of b) delta += Math.abs(n - (a.get(k) || 0));
          for (const [k, n] of a) if (!b.has(k)) delta += n;
          return delta;
        };
        const withNulls = canvasHist();
        pie.props.includeNulls = false;
        await new Promise((r) => setTimeout(r, 800));
        const offDelta = histDiff(withNulls, canvasHist());

        pie.props.includeNulls = true;
        await new Promise((r) => setTimeout(r, 800));
        const pts2 = await w.__pieCalibrate('RACE_GAPS', [null]);
        let reSel = -1;
        if (pts2['__null__']) {
          await w.__pieClickAt(pts2['__null__'].x, pts2['__null__'].y, false);
          reSel = df.selection.trueCount;
        }
        df.selection.setAll(false);
        return {calibrated: true, missing, nullSel, allNull, offDelta, reSel};
      } finally {
        df.selection.setAll(false);
        pie.props.categoryColumnName = 'RACE';
        if (df.col('RACE_GAPS'))
          df.columns.remove('RACE_GAPS');
        await new Promise((r) => setTimeout(r, 400));
      }
    });
    expect(result.calibrated).toBe(true);
    expect(result.missing).toBeGreaterThan(0);
    expect(result.allNull).toBe(true);
    expect(result.nullSel).toBe(result.missing);
    expect(result.offDelta).toBeGreaterThan(1000);
    expect(result.reSel).toBe(result.missing);
  });

  await softStep('Filter mode — slice click filters to the category row count, Ctrl+click adds a second category to the sum, empty-area click restores the full count', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      const df = grok.shell.tv.dataFrame;
      const race = df.col('RACE');
      pie.props.categoryColumnName = 'RACE';
      pie.props.onClick = 'Select';
      await new Promise((r) => setTimeout(r, 500));
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
      if (!pts[t1] || !pts[t2])
        return {calibrated: false};

      pie.props.onClick = 'Filter';
      await new Promise((r) => setTimeout(r, 400));
      const full = df.filter.trueCount;

      await w.__pieClickAt(pts[t1].x, pts[t1].y, false);
      const f1 = df.filter.trueCount;
      const cats1 = filteredCats();
      const expected1 = countOf(t1);

      await w.__pieClickAt(pts[t2].x, pts[t2].y, true);
      const f2 = df.filter.trueCount;
      const cats2 = filteredCats();
      const expected2 = countOf(t1) + countOf(t2);

      const rect = w.__pieRect();
      await w.__pieClickAt(rect.left + 5, rect.top + 5, false);
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
    const result = await page.evaluate(async () => {
      const w = window as any;
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      const df = grok.shell.tv.dataFrame;
      const race = df.col('RACE');
      pie.props.categoryColumnName = 'RACE';
      pie.props.onClick = 'Select';
      await new Promise((r) => setTimeout(r, 400));
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
        return {calibrated: false};

      pie.props.onClick = 'Filter';
      await new Promise((r) => setTimeout(r, 300));
      const full = df.filter.trueCount;

      grok.shell.tv.getFiltersGroup().updateOrAdd({type: 'histogram', column: 'AGE', min: 30, max: 50});
      await new Promise((r) => setTimeout(r, 700));
      const afterPanel = df.filter.trueCount;
      const panelMask: boolean[] = [];
      for (let i = 0; i < df.rowCount; i++)
        panelMask.push(df.filter.get(i));

      await w.__pieClickAt(pts[t1].x, pts[t1].y, false);
      const afterBoth = df.filter.trueCount;
      let expectedBoth = 0;
      for (let i = 0; i < df.rowCount; i++)
        if (panelMask[i] && race.get(i) === t1) expectedBoth++;

      const rect = w.__pieRect();
      await w.__pieClickAt(rect.left + 5, rect.top + 5, false);
      const backToPanel = df.filter.trueCount;

      const btn = document.querySelector('.d4-filter-group-header [name="icon-arrow-rotate-left"]') as HTMLElement | null;
      if (btn)
        btn.click();
      await new Promise((r) => setTimeout(r, 800));
      const afterReset = df.filter.trueCount;
      pie.props.onClick = 'Select';
      return {calibrated: true, full, afterPanel, afterBoth, expectedBoth, backToPanel, afterReset};
    });
    expect(result.calibrated).toBe(true);
    expect(result.afterPanel).toBeLessThan(result.full);
    expect(result.afterBoth).toBeLessThan(result.afterPanel);
    expect(result.afterBoth).toBe(result.expectedBoth);
    expect(result.backToPanel).toBe(result.afterPanel);
    expect(result.afterReset).toBe(result.full);
  });

  await softStep('Viewer close releases the filter — closing the pie chart with an active click-filter restores the full count, re-adding the viewer succeeds', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const pie = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Pie chart') as any;
      const df = grok.shell.tv.dataFrame;
      const race = df.col('RACE');
      pie.props.categoryColumnName = 'RACE';
      pie.props.onClick = 'Select';
      await new Promise((r) => setTimeout(r, 300));
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
        return {calibrated: false};

      pie.props.onClick = 'Filter';
      await new Promise((r) => setTimeout(r, 300));
      await w.__pieClickAt(pts[t1].x, pts[t1].y, false);
      const active = df.filter.trueCount;

      pie.close();
      await new Promise((r) => setTimeout(r, 800));
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
