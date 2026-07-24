/* ---
realizes: [piechart.cp.setup-aggregation-legend-persistence, piechart.category-column-drives-legend]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';


declare const grok: any;


test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';


test('Pie Chart — Category, Legend, Aggregation, and Persistence', async ({page}) => {
  test.setTimeout(300_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  const consoleErrors: string[] = [];
  // Transient resource-load hiccups on the shared dev server are not product
  // errors; a real fault (TypeError/NullError/…) still fails the floor.
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

  // Per-color canvas histogram machinery. The pie redistributes slice angles
  // inside a constant disc and paints labels OVER already-inked slices, so the
  // total non-white pixel count stays flat across every transition this spec
  // drives — the honest repaint signal is the summed per-color pixel delta
  // between SETTLED frames (a frame is settled when two consecutive snapshots
  // are identical).
  const settledColorSnap = () => page.evaluate(async () => {
    const pie = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pie chart') as any;
    const snap = () => {
      const cv = pie.root.querySelector('canvas') as HTMLCanvasElement;
      const d = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
      const m: Record<number, number> = {};
      for (let i = 0; i < d.length; i += 4) {
        const k = (d[i] << 16) | (d[i + 1] << 8) | d[i + 2];
        m[k] = (m[k] ?? 0) + 1;
      }
      return m;
    };
    const diff = (a: Record<number, number>, b: Record<number, number>) => {
      let s = 0;
      for (const k of Object.keys(b)) s += Math.abs(b[+k] - (a[+k] ?? 0));
      for (const k of Object.keys(a)) if (!(k in b)) s += a[+k];
      return s;
    };
    let prev = snap();
    for (let i = 0; i < 8; i++) {
      await new Promise((r) => setTimeout(r, 400));
      const cur = snap();
      if (diff(prev, cur) === 0) return cur;
      prev = cur;
    }
    return prev;
  });

  const colorDiff = (a: Record<number, number>, b: Record<number, number>) => {
    let s = 0;
    for (const k of Object.keys(b)) s += Math.abs(b[+k] - (a[+k] ?? 0));
    for (const k of Object.keys(a)) if (!(k in b)) s += a[+k];
    return s;
  };

  const setPieProps = (props: Record<string, any>, settleMs = 600) =>
    page.evaluate(async ({p, ms}) => {
      const pie = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pie chart') as any;
      for (const k of Object.keys(p)) pie.props[k] = p[k];
      await new Promise((r) => setTimeout(r, ms));
    }, {p: props, ms: settleMs});

  // The single scenario is a configuration ladder: every setting below stacks
  // on the previous ones and is NEVER reverted — the final state is exactly
  // what the two persistence round-trips save and restore.
  await softStep('Configuration ladder — category, Show Value, aggregations, legend color, title', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;

    // Rung 1: Category = RACE, legend always visible → legend labels equal
    // the RACE categories exactly, no extras — set equality both ways.
    await setPieProps({categoryColumnName: 'RACE', legendVisibility: 'Always'}, 800);
    const legend = await page.evaluate(() => {
      const pie = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pie chart') as any;
      const df = grok.shell.tv.dataFrame;
      const labels = Array.from(pie.root.querySelectorAll('[name="legend"] .d4-legend-item .d4-legend-value'))
        .map((e: any) => (e.textContent ?? '').trim());
      return {labels, cats: df.col('RACE').categories.slice()};
    });
    expect(legend.labels.length).toBeGreaterThan(0);
    expect([...legend.labels].sort()).toEqual([...legend.cats].sort());

    // Rung 2: Show Value on — the value labels paint OVER the slices, so the
    // delta is text-sized, not slice-sized: a conservative floor above render
    // noise (noise is 0 by the settle gate) but far below the
    // aggregation-transition threshold. Show Value STAYS on.
    const preValueFrame = await settledColorSnap();
    await setPieProps({showValue: true});
    const valueOnFrame = await settledColorSnap();
    const valueOn = colorDiff(preValueFrame, valueOnFrame);

    // Rungs 3-5: angle column AGE, aggregation count → avg → sum, all with
    // Show Value already on. In demog, sum(AGE) per race is nearly
    // proportional to the row count, so the count/sum pair barely changes the
    // geometry — the asserted transitions are count→avg and avg→sum, where
    // the slice angles genuinely redistribute. The aggregation STAYS at sum.
    await setPieProps({segmentAngleColumnName: 'AGE', segmentAngleAggrType: 'count'});
    const countFrame = await settledColorSnap();
    await setPieProps({segmentAngleAggrType: 'avg'});
    const avgFrame = await settledColorSnap();
    const countToAvg = colorDiff(countFrame, avgFrame);
    await setPieProps({segmentAngleAggrType: 'sum'});
    const sumFrame = await settledColorSnap();
    const avgToSum = colorDiff(avgFrame, sumFrame);
    // Keep the measured deltas visible on green runs so the fixed thresholds
    // can be audited against live numbers.
    console.log(`Repaint px: valueOn=${valueOn} countToAvg=${countToAvg} avgToSum=${avgToSum}`);
    expect(valueOn).toBeGreaterThan(500);
    expect(countToAvg).toBeGreaterThan(10000);
    expect(avgToSum).toBeGreaterThan(10000);

    // Rung 6: the REAL legend color flow, no JS-API substitution: right-click
    // on the legend item opens the per-category color dialog; the same dialog
    // also opens from the [name="legend-icon-color-picker"] icon shown on
    // hover. The color STAYS.
    const colorResult = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pie chart') as any;
      const item = Array.from(pie.root.querySelectorAll('[name="legend"] .d4-legend-item'))
        .find((i: any) => (i.textContent ?? '').includes('Asian')) as HTMLElement;
      if (!item) return {stage: 'no-legend-item'};
      const r = item.getBoundingClientRect();
      item.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2, clientX: r.x + 5, clientY: r.y + 5,
      }));
      await new Promise((res) => setTimeout(res, 800));
      const dlg = document.querySelector('.d4-dialog[name="dialog-Asian"]');
      if (!dlg) return {stage: 'no-dialog'};
      const sw = (Array.from(dlg.querySelectorAll('.d4-color-bar')) as HTMLElement[])
        .find((s) => s.style.backgroundColor === 'rgb(214, 39, 40)');
      if (!sw) return {stage: 'no-swatch'};
      const o = {bubbles: true, cancelable: true, view: window, button: 0};
      sw.dispatchEvent(new MouseEvent('mousedown', o));
      sw.dispatchEvent(new MouseEvent('mouseup', o));
      sw.dispatchEvent(new MouseEvent('click', o));
      await new Promise((res) => setTimeout(res, 300));
      (dlg.querySelector('[name="button-OK"]') as HTMLElement).click();
      await new Promise((res) => setTimeout(res, 1000));
      const after = Array.from(pie.root.querySelectorAll('[name="legend"] .d4-legend-item'))
        .find((i: any) => (i.textContent ?? '').includes('Asian')) as HTMLElement;
      return {
        stage: 'done',
        swatchColor: after ? after.style.color : null,
        tag: grok.shell.tv.dataFrame.col('RACE').tags['.color-coding-categorical'] ?? null,
        dialogGone: !document.querySelector('.d4-dialog[name="dialog-Asian"]'),
      };
    });
    expect(colorResult.stage).toBe('done');
    expect(colorResult.swatchColor).toBe('rgb(214, 39, 40)');
    expect(colorResult.tag).toContain('"Asian":"#d62728"');
    expect(colorResult.dialogGone).toBe(true);

    // Rung 7: the probe title — the last setting the round-trips must carry.
    await setPieProps({showTitle: true, title: 'Pie Persistence Probe'}, 800);

    // Whole-ladder floor: the viewer root stayed attached and no new page or
    // console errors appeared across all seven rungs.
    expect(await readRootInDom()).toBe(true);
    expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
  });

  // Trailing persistence checks: the two round-trips below persist the state
  // the ladder just configured (category RACE, sum(AGE), Show Value on, the
  // custom Asian color from the legend flow, probe title) — nothing is
  // reconfigured here, the ladder state IS the saved state.
  await softStep('Layout round-trip — saved layout restores the configured pie and viewer set', async () => {
    // Save the layout first and capture its id on the node side, so the
    // finally teardown can delete the probe layout even when a later
    // evaluate or assertion fails.
    const layoutId = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      return layout.id as string;
    });
    try {
      const result = await page.evaluate(async (id) => {
        const tv = grok.shell.tv;
        await new Promise((r) => setTimeout(r, 1000));
        tv.addViewer('Scatter plot');
        await new Promise((r) => setTimeout(r, 500));
        const saved = await grok.dapi.layouts.find(id);
        tv.loadLayout(saved);
        await new Promise((r) => setTimeout(r, 3000));
        const viewers = Array.from(tv.viewers) as any[];
        const pie = viewers.find((x: any) => x.type === 'Pie chart');
        const item = pie
          ? Array.from(pie.root.querySelectorAll('[name="legend"] .d4-legend-item'))
            .find((i: any) => (i.textContent ?? '').includes('Asian')) as HTMLElement
          : null;
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
      // The restored viewer set equals the SAVED set…
      expect(result.hasScatter).toBe(false);
      expect(result.hasPie).toBe(true);
      // …and the restored pie carries the full ladder configuration plus the
      // custom color (the color lives in the dataframe column tag, which the
      // layout round-trip leaves untouched — the legend must still render it).
      expect(result.cat).toBe('RACE');
      expect(result.angleCol).toBe('AGE');
      expect(result.aggr).toBe('sum');
      expect(result.showValue).toBe(true);
      expect(result.title).toBe('Pie Persistence Probe');
      expect(result.asianSwatch).toBe('rgb(214, 39, 40)');
    } finally {
      // Never leak the probe layout, even when an assertion above failed.
      await page.evaluate(async (id) => {
        try {
          const saved = await grok.dapi.layouts.find(id);
          if (saved)
            await grok.dapi.layouts.delete(saved);
        } catch (_) {}
      }, layoutId);
    }
  });

  // The project is saved through the real ribbon Save button (only the UI Save
  // captures the view layout), then closeAll + reopen restores the pie. The
  // ladder settings are deliberately NOT reverted afterwards — only the probe
  // project is deleted.
  await softStep('Project save / Close All / reopen — project restores the configured pie', async () => {
    const projName = 'zz-piechart-persistence-probe-' + Date.now();
    try {
      await page.locator('[name="button-Save"]').first().click();
      await page.locator('.d4-dialog input[type="text"]').first().waitFor({timeout: 8000});
      await page.locator('.d4-dialog input[type="text"]').first().fill(projName);
      await page.locator('.d4-dialog .ui-btn-ok, .d4-dialog-footer button').filter({hasText: /^OK$/i}).first().click({force: true});
      await page.waitForTimeout(3000);
      // A "Share <project>" dialog pops up after a successful save — dismiss it.
      const cancel = page.locator('.d4-dialog .ui-btn, .d4-dialog button').filter({hasText: /^CANCEL$/i}).first();
      if (await cancel.count() > 0)
        await cancel.click({force: true});
      await page.waitForTimeout(800);

      const result = await page.evaluate(async (name) => {
        let proj = null;
        for (let a = 0; a < 6 && !proj; a++) {
          try {
            proj = await grok.dapi.projects.filter('name = "' + name + '"').first();
          } catch (e) {}
          if (!proj)
            await new Promise((r) => setTimeout(r, 1200));
        }
        if (!proj)
          return {found: false};
        grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1500));
        const full = await grok.dapi.projects.find(proj.id);
        await full.open();
        // Readiness poll: wait for the restored Pie chart instead of a fixed
        // sleep — the project open re-creates the view asynchronously.
        let pie = null;
        for (let a = 0; a < 15 && !pie; a++) {
          await new Promise((r) => setTimeout(r, 600));
          const tv = grok.shell.tv;
          pie = tv ? Array.from(tv.viewers).find((x: any) => x.type === 'Pie chart') as any : null;
        }
        const tv = grok.shell.tv;
        const item = pie
          ? Array.from((pie as any).root.querySelectorAll('[name="legend"] .d4-legend-item'))
            .find((i: any) => (i.textContent ?? '').includes('Asian')) as HTMLElement
          : null;
        const p = pie as any;
        return {
          found: true,
          pieRestored: !!pie,
          cat: p?.props?.categoryColumnName,
          angleCol: p?.props?.segmentAngleColumnName,
          aggr: p?.props?.segmentAngleAggrType,
          showValue: p?.props?.showValue,
          title: p?.props?.title,
          asianSwatch: item ? item.style.color : null,
          tag: tv ? (tv.dataFrame.col('RACE').tags['.color-coding-categorical'] ?? null) : null,
        };
      }, projName);

      expect(result.found).toBe(true);
      expect(result.pieRestored).toBe(true);
      // The reopened project restores the full ladder configuration AND the
      // custom color — a cross-session round-trip: the project serializes the
      // dataframe, so the color tag must come back with it.
      expect(result.cat).toBe('RACE');
      expect(result.angleCol).toBe('AGE');
      expect(result.aggr).toBe('sum');
      expect(result.showValue).toBe(true);
      expect(result.title).toBe('Pie Persistence Probe');
      expect(result.asianSwatch).toBe('rgb(214, 39, 40)');
      expect(result.tag).toContain('"Asian":"#d62728"');
    } finally {
      // Never leak the probe project, even when the save/reopen flow failed.
      await page.evaluate(async (name) => {
        try {
          const leftover = await grok.dapi.projects.filter('name = "' + name + '"').first();
          if (leftover)
            await grok.dapi.projects.delete(leftover);
        } catch (_) {}
      }, projName);
    }
  });

  v.finishSpec();
});
