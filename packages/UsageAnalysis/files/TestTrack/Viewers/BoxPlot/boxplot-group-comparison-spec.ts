/* ---
realizes: [boxplot.cp.group-comparison-ladder, boxplot.cp.covariate-adjust-baseline, boxplot.int.covariate-sets-adjustment]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

async function revealIcon(page: Page, iconName: string): Promise<{x: number; y: number}> {
  const origin = await page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Box-plot"]')!;
    const c = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
    return {x: c.x, y: c.y};
  });
  for (const [dx, dy] of [[35, 15], [40, 17], [30, 14], [45, 16]]) {
    const shown1 = await v.armEvent(page, 'grok.events.onTooltipShown', 150);

    await page.mouse.move(origin.x + dx, origin.y + dy);
    await shown1();
  }
  return page.evaluate((name) => {
    const el = document.querySelector(`[name="${name}"]`) as HTMLElement;
    const r = el.getBoundingClientRect();
    return {x: r.x + r.width / 2, y: r.y + r.height / 2};
  }, iconName);
}

async function bpProp(page: Page, prop: string): Promise<any> {
  return page.evaluate((p) => {
    const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
    return bp?.props?.[p];
  }, prop);
}

async function driveSelect(page: Page, which: string, value: string): Promise<void> {
  await page.evaluate(({which, value}) => {
    const root = document.querySelector('[name="viewer-Box-plot"]')!;
    let sel: HTMLSelectElement | undefined;
    if (which === 'adjust') {
      sel = Array.from(root.querySelectorAll('select')).find((s) =>
        Array.from(s.options).some((o) => o.value === 'ratio' || o.value === 'regressOut'));
    } else {
      const controls = document.querySelector('.d4-box-plot-group-comparison-controls')!;
      const selects = Array.from(controls.querySelectorAll('select'));
      if (which === 'method') sel = selects[0];
      else if (which === 'control') sel = selects[1];
      else if (which === 'baseline')
        sel = selects.find((s) => Array.from(s.options).some((o) => o.value === 'pooled'));
    }
    if (!sel) throw new Error(`group-comparison select "${which}" not found`);
    sel.value = value;
    sel.dispatchEvent(new Event('input', {bubbles: true}));
  }, {which, value});

  await v.waitForViewerRendered(page, 'Box plot', 1200);
}

async function addComparisonTable(page: Page): Promise<string> {
  const before: string[] = await page.evaluate(() => grok.shell.tables.map((t: any) => t.name));

  const homeDf: string = await page.evaluate(() => grok.shell.tv.dataFrame.name);
  let lastLabels: string[] = [];
  let itemName: string | null = null;
  for (const [dx, dy] of [[42, 16], [50, 16], [60, 14], [80, 14], [42, 30]]) {
    const o = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Box-plot"]')!;
      const c = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
      return {x: c.x, y: c.y};
    });
    await page.mouse.click(o.x + dx, o.y + dy, {button: 'right'});

    await v.pollValue(() => page.evaluate(() =>
      (Array.from(document.querySelectorAll('.d4-menu-item')) as HTMLElement[]).some((i) => {
        const r = i.getBoundingClientRect();
        return r.width > 0 && r.height > 0
          && /^Add .*Table$/i.test((i.getAttribute('d4-name') ?? i.textContent ?? '').trim());
      })), (ready) => ready, 700, 100);
    const found = await page.evaluate(() => {
      const items = (Array.from(document.querySelectorAll('.d4-menu-item')) as HTMLElement[])
        .filter((i) => {
          const r = i.getBoundingClientRect();
          return r.width > 0 && r.height > 0;
        });
      const label = (i: HTMLElement) => (i.getAttribute('d4-name') ?? i.textContent ?? '').trim();
      const add = items.find((i) => /^Add .*Table$/i.test(label(i)));
      if (!add) return {labels: items.map(label), itemName: null};
      add.click();
      return {labels: items.map(label), itemName: add.getAttribute('name') ?? label(add)};
    });
    lastLabels = found.labels;
    if (found.itemName) { itemName = found.itemName; break; }

    await page.keyboard.press('Escape');

    await v.pollStable(() => page.evaluate(() => document.querySelectorAll('.d4-menu-popup').length),
    (a, b) => a === b, 300, 100);
  }
  if (!itemName)
    throw new Error(`Add-table item not reached; last menu: [${lastLabels.join(' | ')}]`);
  console.log('addComparisonTable clicked item:', itemName);
  await page.waitForFunction((prev) => {
    const names = grok.shell.tables.map((t: any) => t.name);
    return names.some((n: string) => !prev.includes(n));
  }, before, {timeout: 15000}).catch(() => {
    throw new Error(`item ${itemName} clicked but no new workspace table within 15s`);
  });
  const name: string = await page.evaluate((prev) =>
    grok.shell.tables.find((tb: any) => !prev.includes(tb.name))!.name, before);

  await page.evaluate((dfName) => {
    const home = Array.from(grok.shell.tableViews)
      .find((vw: any) => vw.dataFrame?.name === dfName);
    if (!home) throw new Error(`home TableView for ${dfName} not found after add-table`);
    grok.shell.v = home;
  }, homeDf);
  await page.locator('[name="viewer-Box-plot"]').waitFor({timeout: 10000});
  await v.waitForViewerRendered(page, 'Box plot', 600);
  return name;
}

test('Box plot group comparison and covariate adjustment', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await page.evaluate(() => {
    const tv = grok.shell.tv;
    const bp = tv.addViewer('Box plot');
    bp.props.valueColumnName = 'AGE';
    bp.props.category1ColumnName = 'SEX';
  });
  await page.locator('[name="viewer-Box-plot"]').waitFor({timeout: 10000});
  await v.waitForViewerRendered(page, 'Box plot', 1500);

  await softStep('Scenario 1 Step 1: bare p overlay baseline', async () => {
    expect(await bpProp(page, 'showPValue')).toBe(true);
    expect(await bpProp(page, 'showGroupComparison')).toBe(false);
  });

  await softStep('Scenario 1 Step 3: enable group comparison via reveal icon', async () => {
    const pt = await revealIcon(page, 'show-group-stats');
    await page.mouse.click(pt.x, pt.y);

    await v.waitForViewerRendered(page, 'Box plot', 1500);
    expect(await bpProp(page, 'showGroupComparison')).toBe(true);
  });

  await softStep('Scenario 1 Step 4: overall-p tooltip test conclusion', async () => {
    const origin = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Box-plot"]')!;
      const c = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
      return {x: c.x, y: c.y};
    });
    let tip = '';
    for (const [dx, dy] of [[42, 16], [46, 16], [50, 16], [40, 15], [54, 16]]) {
      await page.mouse.move(origin.x + dx, origin.y + dy);

      tip = await v.pollValue(() => page.evaluate(() => {
        const tt = document.querySelector('.d4-tooltip');
        return (tt?.textContent ?? '').trim();
      }), (t) => /t-test/i.test(t), 350, 100);
      if (/t-test/i.test(tip)) break;
    }
    expect(tip).toMatch(/Welch's t-test/i);
  });

  await softStep('Scenario 1 Step 5: switch method via on-chart select', async () => {
    await driveSelect(page, 'method', 'Student');
    expect(await bpProp(page, 'method')).toBe('Student');
    await driveSelect(page, 'method', 'Welch');
  });

  await softStep('Scenario 1 Step 6: switch to RACE (one-way ANOVA)', async () => {
    await v.waitForCanvasQuiet(page, 'Box plot');
    await v.snapshotCanvasColors(page, 'Box plot');
    await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      bp.props.category1ColumnName = 'RACE';
    });

    const delta = await v.waitForCanvasChange(page, 'Box plot', {minDelta: 500, timeoutMs: 20000});
    console.log('Step 6 ANOVA recompute canvas delta:', delta,
      'category1:', await bpProp(page, 'category1ColumnName'));
  });

  await softStep('Scenario 1 Step 7: pick control group (control comparisons)', async () => {

    await driveSelect(page, 'control', 'Caucasian');
    expect(await bpProp(page, 'controlComparisons')).toBe(true);
    expect(await bpProp(page, 'controlGroup')).toBe('Caucasian');

    const px = await v.countCanvasPixels(page, 'Box plot');
    expect(px.total).toBeGreaterThan(0);
  });

  await softStep('Scenario 1 Step 8: add control comparisons table', async () => {
    const name = await addComparisonTable(page);
    console.log('Step 8 comparison table name:', name);

    expect(name.length).toBeGreaterThan(0);
    expect(name).toContain(':');
    expect(name).toContain('by RACE');
    const tableInfo = await page.evaluate((tname) => {
      const t = grok.shell.tables.find((tb: any) => tb.name === tname);
      if (!t) return null;
      const cols = t.columns.names();
      const pCol = cols.find((c: string) => /^p( |$|-)/i.test(c) || /p.?value|p \(/i.test(c));

      let firstP: number | null = null;
      if (pCol) {
        const col = t.columns.byName(pCol);
        for (let i = 0; i < t.rowCount; i++) {
          const raw = col.get(i);
          const num = typeof raw === 'number' ? raw
            : parseFloat(String(raw ?? '').replace(/[<>≈~\s]/g, ''));
          if (Number.isFinite(num)) {firstP = num; break;}
        }
      }
      return {cols, rowCount: t.rowCount, pCol, firstP};
    }, name);
    expect(tableInfo).not.toBeNull();
    expect(tableInfo!.rowCount).toBeGreaterThan(0);
    expect(tableInfo!.pCol).toBeTruthy();

    expect(tableInfo!.firstP).not.toBeNull();

    const origin = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Box-plot"]')!;
      const c = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
      return {x: c.x, y: c.y, w: c.width, h: c.height};
    });
    let accordionPanes: string[] = [];
    for (const [fx, fy] of [[0.10, 0.86], [0.14, 0.86], [0.18, 0.86], [0.10, 0.80], [0.22, 0.86]]) {
      await page.mouse.click(origin.x + origin.w * fx, origin.y + origin.h * fy);

      accordionPanes = await v.pollValue(() => page.evaluate(() =>
        Array.from(document.querySelectorAll('.grok-prop-panel .d4-accordion .d4-accordion-pane'))
          .map((p) => (p.getAttribute('d4-title') ?? '').trim())
          .filter((t) => t.length > 0)),
      (panes) => panes.some((t) => /Statistics|Results|Summary|Hypotheses|Distribution/i.test(t)), 600, 100);
      if (accordionPanes.some((t) => /Statistics|Results|Summary|Hypotheses|Distribution/i.test(t)))
        break;
    }
    console.log('Step 8 stats context-panel accordion panes:', accordionPanes);
    expect(accordionPanes.some((t) => /Statistics|Results|Summary|Hypotheses|Distribution/i.test(t)))
      .toBe(true);
  });

  await softStep('Scenario 1 Step 9: two-way ANOVA', async () => {
    await v.setViewerProps(page, 'Box plot', [{set: {controlComparisons: false, controlGroup: ''}, wait: 800}]);
    await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      bp.props.category2ColumnName = 'SEX';
    });

    await v.pollValue(() => page.evaluate(() => {
      const controls = document.querySelector('.d4-box-plot-group-comparison-controls');
      if (!controls) return false;
      const baseline = Array.from(controls.querySelectorAll('select')).find((s) =>
        Array.from((s as HTMLSelectElement).options).some((o) => o.value === 'pooled'));
      return !!baseline && getComputedStyle(baseline.closest('.ui-input-root')!).display !== 'none';
    }), (shown) => shown, 15000, 250);
    expect(await bpProp(page, 'category2ColumnName')).toBe('SEX');

    const px = await v.countCanvasPixels(page, 'Box plot');
    expect(px.total).toBeGreaterThan(0);
    const name = await addComparisonTable(page);
    console.log('Step 9 two-way table name:', name);
    expect(name.length).toBeGreaterThan(0);
    expect(name).toContain('Two-Way ANOVA:');
    expect(name).toContain('by RACE, SEX');
    const rowCount = await page.evaluate((tname) => {
      const t = grok.shell.tables.find((tb: any) => tb.name === tname);
      return t ? t.rowCount : 0;
    }, name);

    expect(rowCount).toBeGreaterThanOrEqual(3);
  });

  await softStep('Scenario 1 Step 10: close group comparison', async () => {
    const pt = await revealIcon(page, 'close-group-stats');
    await page.mouse.click(pt.x, pt.y);

    await v.pollValue(() => page.evaluate(() => {
      const controls = document.querySelector('.d4-box-plot-group-comparison-controls');
      return !controls || getComputedStyle(controls).visibility === 'hidden';
    }), (hidden) => hidden, 1200, 100);
    expect(await bpProp(page, 'showGroupComparison')).toBe(false);
    expect(await bpProp(page, 'showPValue')).toBe(true);
    const controlsHidden = await page.evaluate(() => {
      const controls = document.querySelector('.d4-box-plot-group-comparison-controls');
      return !controls || getComputedStyle(controls).visibility === 'hidden';
    });
    expect(controlsHidden).toBe(true);
  });

  await softStep('Scenario 2 Setup: WEIGHT / SEX, group comparison on', async () => {
    await page.evaluate(async () => {
      const w = window as any;
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      await w.__settled('viewer:Box plot.onViewerRendered', () => {
        bp.props.category2ColumnName = '';
        bp.props.controlComparisons = false;
        bp.props.controlGroup = '';
        bp.props.covariateColumnName = '';
        bp.props.category1ColumnName = 'SEX';
        bp.props.valueColumnName = 'WEIGHT';
      }, 1000);
    });
    if (await bpProp(page, 'showGroupComparison') !== true) {
      const pt = await revealIcon(page, 'show-group-stats');
      await page.mouse.click(pt.x, pt.y);

      await v.waitForViewerRendered(page, 'Box plot', 1200);
    }
    expect(await bpProp(page, 'showGroupComparison')).toBe(true);
  });

  await softStep('Scenario 2 Step 2: covariate HEIGHT (regress-out default)', async () => {
    await page.evaluate(async () => {
      const w = window as any;
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      await w.__settled('viewer:Box plot.onViewerRendered', () => {
        bp.props.covariateColumnName = 'HEIGHT';
      }, 1500);
    });
    expect(await bpProp(page, 'adjustmentMode')).toBe('regressOut');
    const dom = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Box-plot"]')!;
      const caption = Array.from(root.querySelectorAll('.d4-column-selector-caption'))
        .some((c) => /Adjust by:/i.test(c.textContent ?? ''));
      const toggle = Array.from(root.querySelectorAll('select')).find((s) =>
        Array.from(s.options).some((o) => o.value === 'regressOut'));
      const selectors = Array.from(root.querySelectorAll('.d4-column-selector'))
        .map((s) => (s.textContent ?? '').replace(/\s+/g, ' ').trim());
      return {caption, toggleValue: toggle?.value ?? '', selectors};
    });
    expect(dom.caption).toBe(true);
    expect(dom.toggleValue).toBe('regressOut');
    console.log('Step 2 axis selectors:', JSON.stringify(dom.selectors));
    expect(dom.selectors.some((t) => /Adjust by:\s*HEIGHT/i.test(t))).toBe(true);
    expect(dom.selectors.some((t) => t === 'WEIGHT')).toBe(true);
    expect(dom.selectors.some((t) => /WEIGHT\s*\(adj\)|WEIGHT\s*\/\s*HEIGHT/i.test(t))).toBe(false);
  });

  await softStep('Scenario 2 Step 3: flip adjustment to Ratio', async () => {
    await driveSelect(page, 'adjust', 'ratio');
    expect(await bpProp(page, 'adjustmentMode')).toBe('ratio');
    const toggleValue = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Box-plot"]')!;
      const t = Array.from(root.querySelectorAll('select')).find((s) =>
        Array.from(s.options).some((o) => o.value === 'ratio'));
      return t?.value ?? '';
    });
    expect(toggleValue).toBe('ratio');
  });

  await softStep('Scenario 2 Step 4: ANCOVA method', async () => {
    await page.evaluate(async () => {
      const w = window as any;
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      await w.__settled('viewer:Box plot.onViewerRendered', () => {
        bp.props.category1ColumnName = 'RACE';
      }, 1200);
    });
    await driveSelect(page, 'control', 'Caucasian');
    await driveSelect(page, 'method', 'ANCOVA');

    expect(await bpProp(page, 'method')).toBe('ANCOVA');
    const dom = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Box-plot"]')!;
      const toggle = Array.from(root.querySelectorAll('select')).find((s) =>
        Array.from(s.options).some((o) => o.value === 'regressOut' || o.value === 'ratio'));
      const box = toggle?.getBoundingClientRect();
      const toggleHidden = !toggle || !box || box.width === 0 || box.height === 0;
      const caption = Array.from(root.querySelectorAll('.d4-column-selector-caption'))
        .some((c) => /Adjust by:/i.test(c.textContent ?? ''));
      return {toggleHidden, caption};
    });
    expect(dom.toggleHidden).toBe(true);
    expect(dom.caption).toBe(true);
    expect(await bpProp(page, 'adjustmentMode')).toBe('ratio');
  });

  await softStep('Scenario 2 Step 5: add ANCOVA table', async () => {
    const name = await addComparisonTable(page);
    console.log('Step 5 ANCOVA table name:', name);
    expect(name.length).toBeGreaterThan(0);
    expect(name).toContain('ANCOVA');
    expect(name).toContain('by RACE');
    const tableInfo = await page.evaluate((tname) => {
      const t = grok.shell.tables.find((tb: any) => tb.name === tname);
      if (!t) return null;
      const cols: string[] = t.columns.names();
      const pCol = cols.find((c) => /^p( |$|-)/i.test(c) || /p.?value|p \(/i.test(c));
      const fCol = cols.find((c) => /^f( |$|-)|f.?stat|^f$/i.test(c));
      const statCol = pCol ?? fCol;

      let firstStat: number | null = null;
      if (statCol) {
        const col = t.columns.byName(statCol);
        for (let i = 0; i < t.rowCount; i++) {
          const raw = col.get(i);
          const num = typeof raw === 'number' ? raw
            : parseFloat(String(raw ?? '').replace(/[<>≈~\s]/g, ''));
          if (Number.isFinite(num)) {firstStat = num; break;}
        }
      }
      const sample = cols.map((c) => {
        const col = t.columns.byName(c);
        return `${c}[${col.type}]=${col.getString(0)}`;
      });
      return {cols, rowCount: t.rowCount, pCol, fCol, firstStat, sample};
    }, name);
    console.log('Step 5 ANCOVA table sample row:', JSON.stringify(tableInfo?.sample));
    expect(tableInfo).not.toBeNull();
    expect(tableInfo!.rowCount).toBeGreaterThan(0);
    expect(tableInfo!.pCol || tableInfo!.fCol).toBeTruthy();

    expect(tableInfo!.firstStat).not.toBeNull();
  });

  await softStep('Scenario 2 Step 6: matched-per-stratum baseline', async () => {
    await page.evaluate(async () => {
      const w = window as any;
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      await w.__settled('viewer:Box plot.onViewerRendered', () => {
        bp.props.category2ColumnName = 'SEX';
      }, 1500);
    });
    await driveSelect(page, 'baseline', 'matched');
    expect(await bpProp(page, 'baselineMode')).toBe('matched');

    const px = await v.countCanvasPixels(page, 'Box plot');
    expect(px.total).toBeGreaterThan(0);

    try {
      await page.evaluate(async () => {
        const df = grok.shell.tv.dataFrame;
        await df.columns.addNewCalculated('SIMPSON_STRAT',
          'if(Mod(Round(${HEIGHT} * 1000), 2) == 0, "A", "B")');
        await df.columns.addNewCalculated('SIMPSON_VAL',
          'if(${SEX} == "M", 0, if(Mod(Round(${HEIGHT} * 1000), 2) == 0, 2.5, -2.5)) + (Mod(Round(${WEIGHT} * 137), 600) / 30)');
      });
      const created = await page.evaluate(() => grok.shell.tv.dataFrame.columns.names());
      expect(created).toContain('SIMPSON_STRAT');
      expect(created).toContain('SIMPSON_VAL');

      await page.evaluate(async () => {
        const w = window as any;
        const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
        await w.__settled('viewer:Box plot.onViewerRendered', () => {
          bp.props.category1ColumnName = 'SEX';
          bp.props.category2ColumnName = 'SIMPSON_STRAT';
          bp.props.valueColumnName = 'SIMPSON_VAL';
        }, 1500);
      });

      await page.waitForFunction(() => {
        const controls = document.querySelector('.d4-box-plot-group-comparison-controls');
        const sels = Array.from(controls?.querySelectorAll('select') ?? []);
        return sels.some((s) =>
          Array.from((s as HTMLSelectElement).options).some((o) => o.value === 'M'));
      }, null, {timeout: 15000});
      await driveSelect(page, 'control', 'M');
      expect(await bpProp(page, 'baselineMode')).toBe('matched');

      await v.pollValue(() => page.evaluate(() => {
        const root = document.querySelector('[name="viewer-Box-plot"]');
        const spans = Array.from(root?.querySelectorAll('span') ?? []) as HTMLElement[];
        const cue = spans.find((el) => (el.textContent ?? '').trim() === '!'
          && getComputedStyle(el).color === 'rgb(235, 103, 103)');
        if (!cue) return false;
        const r = cue.getBoundingClientRect();
        return getComputedStyle(cue).display !== 'none' && r.width > 0 && r.height > 0;
      }), (shown) => shown, 20000, 250);
      const cue = await page.evaluate(() => {
        const root = document.querySelector('[name="viewer-Box-plot"]');
        const spans = Array.from(root?.querySelectorAll('span') ?? []) as HTMLElement[];
        const el = spans.find((s) => (s.textContent ?? '').trim() === '!'
          && getComputedStyle(s).color === 'rgb(235, 103, 103)');
        if (!el) return null;
        const r = el.getBoundingClientRect();
        return {display: getComputedStyle(el).display, w: r.width, h: r.height,
          x: r.x + r.width / 2, y: r.y + r.height / 2};
      });
      console.log('Step 6 Simpson cue:', JSON.stringify(cue));
      expect(cue).not.toBeNull();
      expect(cue!.display).not.toBe('none');
      expect(cue!.w).toBeGreaterThan(0);
      await page.mouse.move(cue!.x, cue!.y);
      const tip = await v.pollValue(() => page.evaluate(() =>
        (document.querySelector('.d4-tooltip')?.textContent ?? '').trim()),
      (t) => /Pooling cancels opposite within-stratum trends/i.test(t), 600, 100);
      console.log('Step 6 Simpson tooltip:', tip);
      expect(tip).toMatch(/Pooling cancels opposite within-stratum trends/i);
    } finally {

      const names = await page.evaluate(async () => {
        const w = window as any;
        const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
        await w.__settled('viewer:Box plot.onViewerRendered', () => {
          if (bp) {
            bp.props.valueColumnName = 'WEIGHT';
            bp.props.category1ColumnName = 'RACE';
            bp.props.category2ColumnName = 'SEX';
          }
        }, 800);
        const df = grok.shell.tv.dataFrame;
        for (const n of ['SIMPSON_STRAT', 'SIMPSON_VAL'])
          if (df.columns.names().includes(n)) df.columns.remove(n);
        return df.columns.names();
      });
      expect(names).not.toContain('SIMPSON_STRAT');
      expect(names).not.toContain('SIMPSON_VAL');
    }
  });

  await softStep('Scenario 2 Step 7: clear covariate', async () => {
    await page.evaluate(async () => {
      const w = window as any;
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      await w.__settled('viewer:Box plot.onViewerRendered', () => {
        bp.props.covariateColumnName = '';
      }, 1200);
    });
    expect(await bpProp(page, 'covariateColumnName')).toBe('');
    expect(await bpProp(page, 'adjustmentMode')).toBe('');
    expect(await bpProp(page, 'method')).toBe('');
  });

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});
