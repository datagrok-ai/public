/* ---
realizes: [boxplot.cp.group-comparison-ladder, boxplot.cp.covariate-adjust-baseline, boxplot.int.covariate-sets-adjustment]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

// Reveal a hover-gated group-comparison icon by moving the trusted pointer over
// the bare-p reveal zone (canvas top-left), then return its center for a click.
async function revealIcon(page: Page, iconName: string): Promise<{x: number; y: number}> {
  const origin = await page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Box-plot"]')!;
    const c = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
    return {x: c.x, y: c.y};
  });
  for (const [dx, dy] of [[35, 15], [40, 17], [30, 14], [45, 16]]) {
    await page.mouse.move(origin.x + dx, origin.y + dy);
    await page.waitForTimeout(150);
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

// The on-chart group-comparison selects commit viewer props only on the
// `input` event (not `change`).
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
  await page.waitForTimeout(1200);
}

// The exclusive strip menu with the "Add ... Table" item opens only when the
// right-click lands inside the comparison strip (~dx 40 from the canvas origin);
// smaller offsets open the FULL viewer menu, which has no Add items. Probe with
// a FRESH canvas rect per attempt; a synthetic .click() on the item suffices.
async function addComparisonTable(page: Page): Promise<string> {
  const before: string[] = await page.evaluate(() => grok.shell.tables.map((t: any) => t.name));
  // The added table OPENS ITS OWN TableView and steals the current view —
  // remember the home dataframe so the helper can switch back afterwards.
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
    await page.waitForTimeout(700);
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
    // Wrong (non-strip) menu opened — close it and re-aim.
    await page.keyboard.press('Escape');
    await page.waitForTimeout(300);
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
  // Switch back to the home TableView so the box plot stays the active surface.
  await page.evaluate((dfName) => {
    const home = Array.from(grok.shell.tableViews)
      .find((vw: any) => vw.dataFrame?.name === dfName);
    if (!home) throw new Error(`home TableView for ${dfName} not found after add-table`);
    grok.shell.v = home;
  }, homeDf);
  await page.locator('[name="viewer-Box-plot"]').waitFor({timeout: 10000});
  await page.waitForTimeout(600);
  return name;
}

test('Box plot group comparison and covariate adjustment', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  // #### Setup: Box plot with a large viewport, Value=AGE, Category1=SEX
  await page.evaluate(() => {
    const tv = grok.shell.tv;
    const bp = tv.addViewer('Box plot');
    bp.props.valueColumnName = 'AGE';
    bp.props.category1ColumnName = 'SEX';
  });
  await page.locator('[name="viewer-Box-plot"]').waitFor({timeout: 10000});
  await page.waitForTimeout(1500);

  // ==================================================================
  // Scenario 1: Group-comparison ladder
  // ==================================================================

  // #### Scenario 1 Step 1: plain t-test p overlay is the pre-comparison baseline
  await softStep('Scenario 1 Step 1: bare p overlay baseline', async () => {
    expect(await bpProp(page, 'showPValue')).toBe(true);
    expect(await bpProp(page, 'showGroupComparison')).toBe(false);
  });

  // #### Scenario 1 Step 3: reveal + click show-group-stats → showGroupComparison true
  await softStep('Scenario 1 Step 3: enable group comparison via reveal icon', async () => {
    const pt = await revealIcon(page, 'show-group-stats');
    await page.mouse.click(pt.x, pt.y);
    await page.waitForTimeout(1500);
    expect(await bpProp(page, 'showGroupComparison')).toBe(true);
  });

  // #### Scenario 1 Step 4: overall-p tooltip carries the test conclusion
  await softStep('Scenario 1 Step 4: overall-p tooltip test conclusion', async () => {
    const origin = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Box-plot"]')!;
      const c = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
      return {x: c.x, y: c.y};
    });
    let tip = '';
    for (const [dx, dy] of [[42, 16], [46, 16], [50, 16], [40, 15], [54, 16]]) {
      await page.mouse.move(origin.x + dx, origin.y + dy);
      await page.waitForTimeout(350);
      tip = await page.evaluate(() => {
        const tt = document.querySelector('.d4-tooltip');
        return (tt?.textContent ?? '').trim();
      });
      if (/t-test/i.test(tip)) break;
    }
    expect(tip).toMatch(/Welch's t-test/i);
  });

  // #### Scenario 1 Step 5: switch method Welch→Student via on-chart select (input event)
  await softStep('Scenario 1 Step 5: switch method via on-chart select', async () => {
    await driveSelect(page, 'method', 'Student');
    expect(await bpProp(page, 'method')).toBe('Student');
    await driveSelect(page, 'method', 'Welch');
  });

  // #### Scenario 1 Step 6: Category1=RACE → one-way ANOVA (canvas floor)
  await softStep('Scenario 1 Step 6: switch to RACE (one-way ANOVA)', async () => {
    await v.waitForCanvasQuiet(page, 'Box plot');
    await v.snapshotCanvasColors(page, 'Box plot');
    await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      bp.props.category1ColumnName = 'RACE';
    });
    // waitForCanvasChange IS the signal — it throws unless the recompute repaints
    // >= minDelta pixels; the prop read is logging only (set-then-read proves nothing).
    const delta = await v.waitForCanvasChange(page, 'Box plot', {minDelta: 500, timeoutMs: 20000});
    console.log('Step 6 ANOVA recompute canvas delta:', delta,
      'category1:', await bpProp(page, 'category1ColumnName'));
  });

  // #### Scenario 1 Step 7: pick control group via on-chart control select → controlComparisons true
  await softStep('Scenario 1 Step 7: pick control group (control comparisons)', async () => {
    // Driven through the on-chart select — setting props.controlGroup alone does NOT flip the flag.
    await driveSelect(page, 'control', 'Caucasian');
    expect(await bpProp(page, 'controlComparisons')).toBe(true);
    expect(await bpProp(page, 'controlGroup')).toBe('Caucasian');
    // Control band + per-group p row are canvas floors.
    const px = await v.countCanvasPixels(page, 'Box plot');
    expect(px.total).toBeGreaterThan(0);
  });

  // #### Scenario 1 Step 8: Add Control Comparisons Table → workspace DataFrame + BoxPlotStatsMeta panel
  await softStep('Scenario 1 Step 8: add control comparisons table', async () => {
    const name = await addComparisonTable(page);
    console.log('Step 8 comparison table name:', name);
    // Two separate contains — an adjusted value may insert an infix between the
    // test prefix and 'by ...'; never full-name equality.
    expect(name.length).toBeGreaterThan(0);
    expect(name).toContain(':');
    expect(name).toContain('by RACE');
    const tableInfo = await page.evaluate((tname) => {
      const t = grok.shell.tables.find((tb: any) => tb.name === tname);
      if (!t) return null;
      const cols = t.columns.names();
      const pCol = cols.find((c: string) => /^p( |$|-)/i.test(c) || /p.?value|p \(/i.test(c));
      // p columns may be STRING-typed ('.0038', '<.0001') — parse the numeric.
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
    // A p-value is READ from the table, not just the column structure.
    expect(tableInfo!.firstP).not.toBeNull();

    // The per-group p row is canvas-drawn under the category labels; clicking a
    // stats region opens the BoxPlotStatsMeta accordion in the Context Panel —
    // sweep the lower-left stats-strip band.
    const origin = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Box-plot"]')!;
      const c = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
      return {x: c.x, y: c.y, w: c.width, h: c.height};
    });
    let accordionPanes: string[] = [];
    for (const [fx, fy] of [[0.10, 0.86], [0.14, 0.86], [0.18, 0.86], [0.10, 0.80], [0.22, 0.86]]) {
      await page.mouse.click(origin.x + origin.w * fx, origin.y + origin.h * fy);
      await page.waitForTimeout(600);
      accordionPanes = await page.evaluate(() =>
        Array.from(document.querySelectorAll('.grok-prop-panel .d4-accordion .d4-accordion-pane'))
          .map((p) => (p.getAttribute('d4-title') ?? '').trim())
          .filter((t) => t.length > 0));
      if (accordionPanes.some((t) => /Statistics|Results|Summary|Hypotheses|Distribution/i.test(t)))
        break;
    }
    console.log('Step 8 stats context-panel accordion panes:', accordionPanes);
    expect(accordionPanes.some((t) => /Statistics|Results|Summary|Hypotheses|Distribution/i.test(t)))
      .toBe(true);
  });

  // #### Scenario 1 Step 9: Category2=SEX → two-way ANOVA (three effect rows) + Add Two-Way ANOVA Table
  await softStep('Scenario 1 Step 9: two-way ANOVA', async () => {
    await page.evaluate(async () => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      bp.props.controlComparisons = false;
      bp.props.controlGroup = '';
      await new Promise((r) => setTimeout(r, 800));
      bp.props.category2ColumnName = 'SEX';
    });
    // ACTUATION HAZARD: the overlay may briefly show the stale one-way after a
    // props-set Category2 — wait for the recompute (baseline select appears only
    // in two-way; re-toggling group comparison is the recovery path).
    await page.waitForFunction(() => {
      const controls = document.querySelector('.d4-box-plot-group-comparison-controls');
      if (!controls) return false;
      const baseline = Array.from(controls.querySelectorAll('select')).find((s) =>
        Array.from((s as HTMLSelectElement).options).some((o) => o.value === 'pooled'));
      return !!baseline && getComputedStyle(baseline.closest('.ui-input-root')!).display !== 'none';
    }, {timeout: 15000}).catch(() => {});
    expect(await bpProp(page, 'category2ColumnName')).toBe('SEX');
    // THREE effect rows for RACE, SEX, RACE × SEX are a canvas floor.
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
    // The two-way ANOVA table carries the three effect rows.
    expect(rowCount).toBeGreaterThanOrEqual(3);
  });

  // #### Scenario 1 Step 10: close-group-stats → all group-comparison UI hidden, bare p returns
  await softStep('Scenario 1 Step 10: close group comparison', async () => {
    const pt = await revealIcon(page, 'close-group-stats');
    await page.mouse.click(pt.x, pt.y);
    await page.waitForTimeout(1200);
    expect(await bpProp(page, 'showGroupComparison')).toBe(false);
    expect(await bpProp(page, 'showPValue')).toBe(true);
    const controlsHidden = await page.evaluate(() => {
      const controls = document.querySelector('.d4-box-plot-group-comparison-controls');
      return !controls || getComputedStyle(controls).visibility === 'hidden';
    });
    expect(controlsHidden).toBe(true);
  });

  // ==================================================================
  // Scenario 2: Covariate adjustment and control-baseline modes
  // ==================================================================

  // #### Scenario 2 Setup: fresh ladder — Value=WEIGHT, Category1=SEX, group comparison on
  await softStep('Scenario 2 Setup: WEIGHT / SEX, group comparison on', async () => {
    await page.evaluate(async () => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      bp.props.category2ColumnName = '';
      bp.props.controlComparisons = false;
      bp.props.controlGroup = '';
      bp.props.covariateColumnName = '';
      bp.props.category1ColumnName = 'SEX';
      bp.props.valueColumnName = 'WEIGHT';
      await new Promise((r) => setTimeout(r, 1000));
    });
    if (await bpProp(page, 'showGroupComparison') !== true) {
      const pt = await revealIcon(page, 'show-group-stats');
      await page.mouse.click(pt.x, pt.y);
      await page.waitForTimeout(1200);
    }
    expect(await bpProp(page, 'showGroupComparison')).toBe(true);
  });

  // #### Scenario 2 Step 2: set covariate HEIGHT → adjustmentMode regressOut + DOM caption/toggle
  await softStep('Scenario 2 Step 2: covariate HEIGHT (regress-out default)', async () => {
    await page.evaluate(async () => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      bp.props.covariateColumnName = 'HEIGHT';
      await new Promise((r) => setTimeout(r, 1500));
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

  // #### Scenario 2 Step 3: flip the on-axis toggle to Ratio (input event)
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

  // #### Scenario 2 Step 4: RACE + control + method=ANCOVA → method 'ANCOVA', toggle hidden, mode preserved
  await softStep('Scenario 2 Step 4: ANCOVA method', async () => {
    await page.evaluate(async () => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      bp.props.category1ColumnName = 'RACE';
      await new Promise((r) => setTimeout(r, 1200));
    });
    await driveSelect(page, 'control', 'Caucasian');
    await driveSelect(page, 'method', 'ANCOVA');
    // Uppercase literal 'ANCOVA' (a lowercase 'ancova' is invalid and resets to '').
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

  // #### Scenario 2 Step 5: Add ANCOVA Table → workspace DataFrame named '<TestName>: ... by RACE'
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
      // The p-value column is a STRING column (formatted like '.0038' or
      // '<.0001') — parse the numeric out of it instead of expecting a double.
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
    // A statistic value is read from the table (not vacuously true).
    expect(tableInfo!.firstStat).not.toBeNull();
  });

  // #### Scenario 2 Step 6: Category2 + Matched baseline → per-stratum bands + Simpson cue
  await softStep('Scenario 2 Step 6: matched-per-stratum baseline', async () => {
    await page.evaluate(async () => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      bp.props.category2ColumnName = 'SEX';
      await new Promise((r) => setTimeout(r, 1500));
    });
    await driveSelect(page, 'baseline', 'matched');
    expect(await bpProp(page, 'baselineMode')).toBe('matched');
    // Per-stratum control bands are a canvas floor (no-error / non-empty canvas).
    const px = await v.countCanvasPixels(page, 'Box plot');
    expect(px.total).toBeGreaterThan(0);

    // Simpson's-paradox fixture: stratum = HEIGHT sub-mm digit parity (balanced,
    // sex-independent); value = +2.5/-2.5 by stratum for F, 0 for M, plus wide
    // sex-decorrelated noise from WEIGHT — opposite significant within-stratum
    // arms whose pooled effect cancels (the Simpson trigger).
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
      // baselineMode 'matched' persists across the re-bind.
      await page.evaluate(async () => {
        const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
        bp.props.category1ColumnName = 'SEX';
        bp.props.category2ColumnName = 'SIMPSON_STRAT';
        bp.props.valueColumnName = 'SIMPSON_VAL';
        await new Promise((r) => setTimeout(r, 1500));
      });
      // Wait until the control select rebuilt its options for SEX, then pick 'M'.
      await page.waitForFunction(() => {
        const controls = document.querySelector('.d4-box-plot-group-comparison-controls');
        const sels = Array.from(controls?.querySelectorAll('select') ?? []);
        return sels.some((s) =>
          Array.from((s as HTMLSelectElement).options).some((o) => o.value === 'M'));
      }, null, {timeout: 15000});
      await driveSelect(page, 'control', 'M');
      expect(await bpProp(page, 'baselineMode')).toBe('matched');
      // The red '!' cue sits just left of the baseline toggle and is
      // display-gated — poll for the recompute before asserting.
      await page.waitForFunction(() => {
        const root = document.querySelector('[name="viewer-Box-plot"]');
        const spans = Array.from(root?.querySelectorAll('span') ?? []) as HTMLElement[];
        const cue = spans.find((el) => (el.textContent ?? '').trim() === '!'
          && getComputedStyle(el).color === 'rgb(235, 103, 103)');
        if (!cue) return false;
        const r = cue.getBoundingClientRect();
        return getComputedStyle(cue).display !== 'none' && r.width > 0 && r.height > 0;
      }, null, {timeout: 20000}).catch(() => {});
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
      await page.waitForTimeout(600);
      const tip = await page.evaluate(() =>
        (document.querySelector('.d4-tooltip')?.textContent ?? '').trim());
      console.log('Step 6 Simpson tooltip:', tip);
      expect(tip).toMatch(/Pooling cancels opposite within-stratum trends/i);
    } finally {
      // Re-point the viewer off the fixture columns before removing them.
      const names = await page.evaluate(async () => {
        const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
        if (bp) {
          bp.props.valueColumnName = 'WEIGHT';
          bp.props.category1ColumnName = 'RACE';
          bp.props.category2ColumnName = 'SEX';
        }
        await new Promise((r) => setTimeout(r, 800));
        const df = grok.shell.tv.dataFrame;
        for (const n of ['SIMPSON_STRAT', 'SIMPSON_VAL'])
          if (df.columns.names().includes(n)) df.columns.remove(n);
        return df.columns.names();
      });
      expect(names).not.toContain('SIMPSON_STRAT');
      expect(names).not.toContain('SIMPSON_VAL');
    }
  });

  // #### Scenario 2 Step 7: clear the covariate → adjustmentMode '' AND method cleared
  await softStep('Scenario 2 Step 7: clear covariate', async () => {
    await page.evaluate(async () => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      bp.props.covariateColumnName = '';
      await new Promise((r) => setTimeout(r, 1200));
    });
    expect(await bpProp(page, 'covariateColumnName')).toBe('');
    expect(await bpProp(page, 'adjustmentMode')).toBe('');
    expect(await bpProp(page, 'method')).toBe('');
  });

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});
