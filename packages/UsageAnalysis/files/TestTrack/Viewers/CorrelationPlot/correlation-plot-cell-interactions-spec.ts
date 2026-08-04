/* ---
realizes: [correlationplot.cp.cell-interactions, correlationplot.int.ignore-double-click]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;
declare const DG: any;


test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const TOL = 1e-3;

// CSS-pixel center of a correlation cell from the root rect + grid arithmetic (refdoc §Cell
// Geometry and Mouse Clicks); cellW is 40 with showPearsonR, else 20. pinnedW/headerH are
// auto-sized per dataset — the Setup probe click calibrates them.
interface Geometry {rootX: number; rootY: number; pinnedW: number; headerH: number; cellW: number; rowH: number}

async function readGeometry(page: Page): Promise<{rootX: number; rootY: number; cellW: number; xCols: string[]; yCols: string[]}> {
  return await page.evaluate(() => {
    const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
    const root = document.querySelector('[name="viewer-Correlation-plot"]')!;
    const R = root.getBoundingClientRect();
    return {
      rootX: R.x, rootY: R.y,
      cellW: cp.props.showPearsonR ? 40 : 20,
      xCols: cp.props.xColumnNames.slice(),
      yCols: cp.props.yColumnNames.slice(),
    };
  });
}

function cellCenter(g: Geometry, xi: number, yi: number): {x: number; y: number} {
  return {
    x: g.rootX + g.pinnedW + (xi + 0.5) * g.cellW,
    y: g.rootY + g.headerH + (yi + 0.5) * g.rowH,
  };
}

// Dock relayout (any viewer add/remove) SHIFTS the CP root — refresh rootX/rootY before every
// hover/click group that follows a viewer-set change; the cell metrics stay valid.
async function refreshRoot(page: Page, g: Geometry): Promise<void> {
  const r = await page.evaluate(() => {
    const R = document.querySelector('[name="viewer-Correlation-plot"]')!.getBoundingClientRect();
    return {x: R.x, y: R.y};
  });
  g.rootX = r.x;
  g.rootY = r.y;
}

// The cell tooltip fires on onCellMouseEnter DEBOUNCED 200ms and only while the mouse still rests
// on the SAME cell (grid_tooltip.dart) — the gesture must walk IN from an adjacent cell with a
// trusted mousemove stream, then REST past the debounce. The dismiss target must be OUTSIDE the
// viewer root (an inner point lands on another cell; the reusable .d4-tooltip retains stale text
// on mouse-away).
async function dismissTooltip(page: Page, g: Geometry): Promise<void> {
  await page.mouse.move(Math.max(4, g.rootX - 60), g.rootY + g.headerH + 100);
  await page.waitForTimeout(400);
}

async function hoverCell(page: Page, g: Geometry, xi: number, yi: number): Promise<void> {
  const c = cellCenter(g, xi, yi);
  // Walk in from the left neighbour (interpolated trusted moves cross a cell boundary INTO the
  // target), then rest past the 200ms debounce.
  await page.mouse.move(c.x - g.cellW, c.y, {steps: 3});
  await page.mouse.move(c.x, c.y, {steps: 4});
  await page.waitForTimeout(300);
}

// The platform can hold MORE THAN ONE .d4-tooltip node (singleton content + per-TooltipInfo hosts
// + _superTooltipContent) — a first-node read can land on a stale/empty sibling, so poll the
// COMBINED text of every node. Require the strict 'Pearson R: <d.ddd>' format (toStringAsFixed(3))
// so a stale prior tooltip cannot satisfy the wait.
async function waitForRTooltip(page: Page, timeoutMs = 6000): Promise<void> {
  await page.waitForFunction(() => {
    const text = Array.from(document.querySelectorAll('.d4-tooltip'))
      .map((t) => t.textContent ?? '').join(' ');
    return /Pearson R:\s*-?\d+\.\d+/.test(text);
  }, undefined, {timeout: timeoutMs});
}

// Record the last d4-correlation-plot-corr-cell-click event; args column1/column2 are column
// NAME strings (live recon), value is the coefficient.
async function armClickListener(page: Page): Promise<void> {
  await page.evaluate(() => {
    const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
    (window as any).__corrClicks = [];
    cp.onEvent('d4-correlation-plot-corr-cell-click').subscribe((e: any) => {
      const a = e.args ?? e;
      (window as any).__corrClicks.push({col1: a.column1, col2: a.column2, value: a.value});
    });
  });
}

async function lastClick(page: Page): Promise<{col1: string; col2: string; value: number} | null> {
  return await page.evaluate(() => {
    const c = (window as any).__corrClicks;
    return c && c.length ? c[c.length - 1] : null;
  });
}

test('Correlation Plot — Cell Interactions', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  // ## Setup — open demog, add the Correlation plot via its Toolbox icon (DOM-driven),
  // arm the cell-click listener, and calibrate the cell geometry with one probe click.
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'correlation-plot', 'Correlation-plot', 10000);
  await armClickListener(page);

  const base = await readGeometry(page);
  const xiHeight = base.xCols.indexOf('HEIGHT');
  const yiAge = base.yCols.indexOf('AGE');
  const xiWeight = base.xCols.indexOf('WEIGHT');

  // Calibrate pinnedW/headerH: click the computed HEIGHT/AGE center, read the reported pair,
  // and if it is off by a column/row shift pinnedW/headerH by one cell until it lands.
  const geom: Geometry = {rootX: base.rootX, rootY: base.rootY, pinnedW: 104, headerH: 60, cellW: base.cellW, rowH: 20};
  await softStep('Setup Step 5 — calibrate cell geometry via a probe click', async () => {
    let landed = false;
    for (let attempt = 0; attempt < 5 && !landed; attempt++) {
      await page.evaluate(() => { (window as any).__corrClicks = []; });
      const c = cellCenter(geom, xiHeight, yiAge);
      await page.mouse.click(c.x, c.y);
      await page.waitForTimeout(400);
      const ev = await lastClick(page);
      console.log(`[Setup] probe attempt ${attempt} at (${Math.round(c.x)},${Math.round(c.y)}) pinnedW=${geom.pinnedW} headerH=${geom.headerH} -> ${JSON.stringify(ev)}`);
      if (ev && ((ev.col1 === 'HEIGHT' && ev.col2 === 'AGE') || (ev.col1 === 'AGE' && ev.col2 === 'HEIGHT'))) {
        landed = true;
        break;
      }
      // nudge geometry: a miss to the right means pinnedW too small; a miss low means headerH too small.
      if (ev && ev.col1 && ev.col1 !== 'HEIGHT') {
        const idx = base.xCols.indexOf(ev.col1);
        if (idx >= 0) geom.pinnedW += (xiHeight - idx) * geom.cellW;
      } else if (!ev) {
        geom.headerH += 4;
      }
    }
    // A calibrated probe MUST have produced a valid cell-click event with the expected pair —
    // this both proves the coordinate maths and that real trusted mouse reaches the Dart handler.
    const ev = await lastClick(page);
    expect(ev).not.toBeNull();
    expect([ev!.col1, ev!.col2].sort()).toEqual(['AGE', 'HEIGHT']);
  });

  const baselineViewerCount: number = await page.evaluate(() =>
    grok.shell.tv.viewers.filter((x: any) => x.type === 'Scatter plot').length);

  // ### Scenario 1: Single-click fires event and updates context panel
  await softStep('Scenario 1 Step 3 — click HEIGHT/AGE fires event with correct pair and Pearson value', async () => {
    // Expected Pearson coefficient over the full demog set, computed at runtime (never hardcoded).
    const expected: number = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return DG.Stats.fromColumn(df.col('HEIGHT')).corr(df.col('AGE'));
    });
    await page.evaluate(() => { (window as any).__corrClicks = []; });
    const c = cellCenter(geom, xiHeight, yiAge);
    await page.mouse.click(c.x, c.y);
    await page.waitForTimeout(500);
    const ev = await lastClick(page);
    console.log(`[S1] click event=${JSON.stringify(ev)} expectedPearson=${expected}`);
    expect(ev).not.toBeNull();
    // The event names the clicked pair (order is X=column1, Y=column2 per source).
    expect([ev!.col1, ev!.col2].sort()).toEqual(['AGE', 'HEIGHT']);
    // and carries the runtime Pearson coefficient within tolerance.
    expect(Number.isFinite(ev!.value)).toBe(true);
    expect(Math.abs(ev!.value - expected)).toBeLessThanOrEqual(TOL);
  });

  await softStep('Scenario 1 Step 4 — context panel shows columnsCorrelation "<c1> vs <c2>" with a Scatter plot pane', async () => {
    // The click sets AppEvents.currentObject to a columnsCorrelation SemanticValue; the context
    // panel re-renders ASYNC — poll (~3s) for the "<col1> vs <col2>" header and the Scatter plot
    // pane, re-clicking if the panel is still cold.
    const readPanel = () => page.evaluate(() => {
      const vsHeaders = Array.from(document.querySelectorAll('*'))
        .filter((e) => e.children.length === 0 && /\bvs\b/.test(e.textContent ?? '') && (e.textContent ?? '').length < 40)
        .map((e) => (e.textContent ?? '').trim());
      const accordions = Array.from(document.querySelectorAll('.d4-accordion-pane-header, .d4-accordion-pane-title'))
        .map((e) => (e.textContent ?? '').trim());
      return {vsHeaders: [...new Set(vsHeaders)], hasScatterPane: accordions.includes('Scatter plot')};
    });
    const isPair = (hs: string[]) => hs.some((h) => /^(HEIGHT vs AGE|AGE vs HEIGHT)$/.test(h));
    let panel = await readPanel();
    for (let i = 0; i < 12 && !(isPair(panel.vsHeaders) && panel.hasScatterPane); i++) {
      await page.waitForTimeout(250);
      panel = await readPanel();
      if (i === 5 && !isPair(panel.vsHeaders)) {
        // panel still cold after ~1.5s — re-fire the cell click to re-push the semantic value.
        const c = cellCenter(geom, xiHeight, yiAge);
        await page.mouse.click(c.x, c.y);
      }
    }
    console.log(`[S1] panel vsHeaders=${JSON.stringify(panel.vsHeaders)} hasScatterPane=${panel.hasScatterPane}`);
    // The context panel carries a "HEIGHT vs AGE" (either order) correlation header.
    expect(isPair(panel.vsHeaders)).toBe(true);
    // and an expandable Scatter plot pane.
    expect(panel.hasScatterPane).toBe(true);
  });

  // ### Scenario 2: Double-click opens scatter plot; ignoreDoubleClick suppresses it
  await softStep('Scenario 2 Step 3-4 — double-click WEIGHT/AGE adds exactly one Scatter plot for the pair', async () => {
    const before: number = await page.evaluate(() =>
      grok.shell.tv.viewers.filter((x: any) => x.type === 'Scatter plot').length);
    const c = cellCenter(geom, xiWeight, yiAge);
    await page.mouse.dblclick(c.x, c.y);
    await page.waitForTimeout(1500);
    const r = await page.evaluate(() => {
      const scatters = grok.shell.tv.viewers.filter((x: any) => x.type === 'Scatter plot');
      const newest = scatters[scatters.length - 1];
      return {count: scatters.length, x: newest?.props.xColumnName ?? null, y: newest?.props.yColumnName ?? null};
    });
    console.log(`[S2] scatter before=${before} after=${r.count} x=${r.x} y=${r.y}`);
    // Exactly one new Scatter plot appeared, bound to the clicked pair.
    expect(r.count).toBe(before + 1);
    expect(r.x).toBe('WEIGHT');
    expect(r.y).toBe('AGE');
  });

  await softStep('Scenario 2 Step 5-6 — close the opened Scatter plot; count returns to baseline', async () => {
    const after: number = await page.evaluate(() => {
      const scatters = grok.shell.tv.viewers.filter((x: any) => x.type === 'Scatter plot');
      for (const sc of scatters) sc.close();
      return grok.shell.tv.viewers.filter((x: any) => x.type === 'Scatter plot').length;
    });
    await page.waitForTimeout(500);
    console.log(`[S2] scatter after close=${after} baseline=${baselineViewerCount}`);
    expect(after).toBe(baselineViewerCount);
  });

  await softStep('Scenario 2 Step 7-9 — Ignore Double Click suppresses the scatter (correlationplot.int.ignore-double-click)', async () => {
    // Set the property via the Context Panel property object (the scenario's Step 7 action).
    await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      cp.props.ignoreDoubleClick = true;
    });
    await page.waitForTimeout(400);
    const before: number = await page.evaluate(() =>
      grok.shell.tv.viewers.filter((x: any) => x.type === 'Scatter plot').length);
    const c = cellCenter(geom, xiWeight, yiAge);
    await page.mouse.dblclick(c.x, c.y);
    await page.waitForTimeout(1500);
    const after: number = await page.evaluate(() =>
      grok.shell.tv.viewers.filter((x: any) => x.type === 'Scatter plot').length);
    // Restore the default for later steps.
    await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      cp.props.ignoreDoubleClick = false;
    });
    console.log(`[S2] ignoreDoubleClick suppression before=${before} after=${after}`);
    // The suppression signal: the viewer count is UNCHANGED — no scatter opened.
    expect(after).toBe(before);
  });

  // ### Scenario 3: Hover tooltip — deduplication and per-cell scatter
  await softStep('Scenario 3 Step 3 — hover Cell A shows exactly one tooltip with "Pearson R: <value>"', async () => {
    const expected: number = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return DG.Stats.fromColumn(df.col('HEIGHT')).corr(df.col('AGE'));
    });
    await refreshRoot(page, geom);
    // Leave the viewer to reset hover state, then dwell-hover Cell A (HEIGHT/AGE).
    await dismissTooltip(page, geom);
    await hoverCell(page, geom, xiHeight, yiAge);
    await waitForRTooltip(page);
    const r = await page.evaluate(() => {
      const tips = Array.from(document.querySelectorAll('.d4-tooltip'));
      // Combined read across all .d4-tooltip nodes (see waitForRTooltip); an entirely EMPTY set
      // means the hover missed the cell, not that the text lives on a hidden layer.
      const combined = tips.map((t) => t.textContent ?? '').join(' ').replace(/\s+/g, ' ').trim();
      return {
        totalCount: tips.length,
        text: combined.slice(0, 400),
        hasCanvas: tips.some((t) => !!t.querySelector('canvas')),
      };
    });
    console.log(`[S3] tooltip total=${r.totalCount} text="${r.text}" expected=${expected.toFixed(3)}`);
    // GROK-19482 duplicate-tooltip guard: exactly one .d4-tooltip element in the DOM.
    expect(r.totalCount).toBe(1);
    // Its text carries "Pearson R: <value>" matching the runtime coefficient to 3 dp.
    const m = r.text.match(/Pearson R:\s*(-?\d+\.\d+)/);
    expect(m).not.toBeNull();
    expect(Math.abs(parseFloat(m![1]) - expected)).toBeLessThanOrEqual(TOL);
    // and embeds a scatter-plot canvas.
    expect(r.hasCanvas).toBe(true);
  });

  await softStep('Scenario 3 Step 5 — the embedded scatter canvas differs between Cell A and Cell B (GROK-20125)', async () => {
    // Snapshot Cell A's embedded scatter canvas (settle-gated), then hover Cell B and diff.
    await refreshRoot(page, geom);
    await dismissTooltip(page, geom);
    await hoverCell(page, geom, xiHeight, yiAge);
    await waitForRTooltip(page);
    const snapA = await page.evaluate(() => {
      // The embedded scatter <canvas> may sit in any of the .d4-tooltip nodes, not necessarily
      // the first — pick the first node that actually holds a sized canvas.
      const cv = Array.from(document.querySelectorAll('.d4-tooltip'))
        .map((t) => t.querySelector('canvas') as HTMLCanvasElement | null)
        .find((c) => !!c && !!c.width && !!c.height) ?? null;
      if (!cv || !cv.width || !cv.height) return null;
      try {
        const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
        const colors: Record<number, number> = {};
        for (let i = 0; i < img.length; i += 4) {
          const key = (img[i] << 16) | (img[i + 1] << 8) | img[i + 2];
          colors[key] = (colors[key] ?? 0) + 1;
        }
        (window as any).__tipSnapA = colors;
        return {w: cv.width, h: cv.height};
      } catch { return null; }
    });
    // Cell B = WEIGHT (X) vs STARTED (Y): a distinct column pair.
    const xiWeightB = base.xCols.indexOf('WEIGHT');
    const yiStarted = base.yCols.indexOf('STARTED');
    await dismissTooltip(page, geom);
    await hoverCell(page, geom, xiWeightB, yiStarted);
    await waitForRTooltip(page);
    const r = await page.evaluate(() => {
      const tips = Array.from(document.querySelectorAll('.d4-tooltip'));
      const cv = tips
        .map((t) => t.querySelector('canvas') as HTMLCanvasElement | null)
        .find((c) => !!c && !!c.width && !!c.height) ?? null;
      const prev = (window as any).__tipSnapA as Record<number, number> | undefined;
      if (!cv || !cv.width || !cv.height || !prev) return {tooltipCount: tips.length, deltaPx: -1};
      try {
        const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
        const colors: Record<number, number> = {};
        for (let i = 0; i < img.length; i += 4) {
          const key = (img[i] << 16) | (img[i + 1] << 8) | img[i + 2];
          colors[key] = (colors[key] ?? 0) + 1;
        }
        let deltaPx = 0;
        const keys = new Set([...Object.keys(colors), ...Object.keys(prev)].map(Number));
        for (const k of keys) deltaPx += Math.abs((colors[k] ?? 0) - (prev[k] ?? 0));
        return {tooltipCount: tips.length, deltaPx};
      } catch { return {tooltipCount: tips.length, deltaPx: -1}; }
    });
    console.log(`[S3] Cell A snap=${JSON.stringify(snapA)} Cell B tooltipCount=${r.tooltipCount} deltaPx=${r.deltaPx}`);
    // Still exactly one tooltip while hovering Cell B (GROK-19482 holds across cells).
    expect(r.tooltipCount).toBe(1);
    // Fault guard before the non-empty check: -1 is a canvas fault, not a valid diff.
    expect(r.deltaPx).toBeGreaterThanOrEqual(0);
    // GROK-20125 per-cell scatter guard: the embedded scatter differs between the two pairs.
    expect(r.deltaPx).toBeGreaterThan(0);
  });

  await softStep('Scenario 3 Step 6 — hover the pinned row-header name column shows a column-statistics tooltip (not an R tooltip)', async () => {
    // Dwell-hover the pinned NAME column at a data row: the tooltip there is the column-statistics
    // tooltip (min/avg/max/std/nulls), distinguished by the ABSENCE of an '<Type> R:' line.
    await refreshRoot(page, geom);
    await dismissTooltip(page, geom);
    const nameX = geom.rootX + geom.pinnedW - 45; // inside the pinned name column
    const nameY = geom.rootY + geom.headerH + geom.rowH * 1.5; // second data row
    // Same cell-enter debounce contract as hoverCell: walk in from the row below, rest past 200ms.
    await page.mouse.move(nameX, nameY + geom.rowH, {steps: 3});
    await page.mouse.move(nameX, nameY, {steps: 4});
    await page.waitForTimeout(300);
    // Poll for a stats marker (a stale corr tooltip cannot satisfy the wait).
    await page.waitForFunction(() => {
      const t = Array.from(document.querySelectorAll('.d4-tooltip'))
        .map((tip) => tip.textContent ?? '').join(' ');
      return /(min:|avg:|max:|std:|nulls:)/.test(t);
    }, undefined, {timeout: 6000});
    const r = await page.evaluate(() => {
      // Combined .d4-tooltip read (nodes stay display:none / zero-rect — no rect filter). The
      // stats tooltip has no '<Type> R:' line.
      const text = Array.from(document.querySelectorAll('.d4-tooltip'))
        .map((tip) => tip.textContent ?? '').join(' ').replace(/\s+/g, ' ').trim();
      return {hasText: text.length > 0, text: text.slice(0, 200), hasRLine: /\bR:\s*-?\d/.test(text)};
    });
    console.log(`[S3] row-header tooltip hasText=${r.hasText} hasRLine=${r.hasRLine} text="${r.text}"`);
    // A tooltip is shown for the row-header name cell ...
    expect(r.hasText).toBe(true);
    // ... and it is a column-statistics tooltip, NOT a correlation '<Type> R: <value>' tooltip.
    expect(r.hasRLine).toBe(false);
  });

  // ### Scenario 4: Right-click context menu items and Tooltip subgroup actions
  await softStep('Scenario 4 Step 2-4 — right-click menu has Show Pearson R, Columns > X/Y, and a Tooltip submenu with Edit.../Visible/Properties...', async () => {
    // Dismiss any tooltip, then real right-click the HEIGHT/AGE cell.
    await page.mouse.move(geom.rootX + 5, geom.rootY + geom.headerH + 200);
    await page.waitForTimeout(400);
    const c = cellCenter(geom, xiHeight, yiAge);
    await page.mouse.click(c.x, c.y, {button: 'right'});
    await page.waitForTimeout(700);
    const menu = await page.evaluate(() => {
      const labels = Array.from(document.querySelectorAll('.d4-menu-item-label')).map((l) => (l.textContent ?? '').trim());
      return {
        showPearsonR: labels.includes('Show Pearson R'),
        columns: labels.includes('Columns'),
        xColumns: labels.includes('X Columns'),
        yColumns: labels.includes('Y Columns'),
        tooltip: labels.includes('Tooltip'),
      };
    });
    console.log(`[S4] menu=${JSON.stringify(menu)}`);
    expect(menu.showPearsonR).toBe(true);
    expect(menu.columns).toBe(true);
    expect(menu.xColumns).toBe(true);
    expect(menu.yColumns).toBe(true);
    expect(menu.tooltip).toBe(true);

    // Step 3-4: the Tooltip submenu is pre-rendered but ZERO-SIZE until the group is hovered.
    // Its items are disambiguated by the sibling 'Use as Group Tooltip' (a second
    // Edit.../Properties... pair lives under the Grid Color Coding group).
    const sub = await page.evaluate(async () => {
      const items = Array.from(document.querySelectorAll('.d4-menu-item'));
      let ttNode: Element | null = null;
      for (const node of items) {
        const own = Array.from(node.querySelectorAll('.d4-menu-item-label')).find((l) => l.closest('.d4-menu-item') === node);
        if (own && (own.textContent ?? '').trim() === 'Tooltip') { ttNode = node; break; }
      }
      if (!ttNode) return {err: 'no Tooltip group'};
      const rect = ttNode.getBoundingClientRect();
      ttNode.dispatchEvent(new MouseEvent('mouseover', {bubbles: true, clientX: rect.x + rect.width / 2, clientY: rect.y + rect.height / 2}));
      ttNode.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: rect.x + rect.width / 2, clientY: rect.y + rect.height / 2}));
      await new Promise((r) => setTimeout(r, 800));
      // The tooltip submenu is the label container that also holds 'Use as Group Tooltip'.
      const anchor = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find((l) => (l.textContent ?? '').trim() === 'Use as Group Tooltip');
      const container = anchor?.closest('.d4-menu-item')?.parentElement ?? null;
      const sibs = container
        ? Array.from(container.querySelectorAll('.d4-menu-item-label')).map((s) => (s.textContent ?? '').trim())
        : [];
      return {edit: sibs.includes('Edit...'), visible: sibs.includes('Visible'), properties: sibs.includes('Properties...')};
    });
    console.log(`[S4] tooltip submenu=${JSON.stringify(sub)}`);
    // The Tooltip submenu carries at minimum Edit..., Visible, and Properties... items.
    expect((sub as any).edit).toBe(true);
    expect((sub as any).visible).toBe(true);
    expect((sub as any).properties).toBe(true);
  });

  await softStep('Scenario 4 Step 4 — Tooltip > Edit... opens the "Edit Tooltip" modal; CANCEL closes it cleanly (GROK-19054)', async () => {
    // The context menu is still open from the previous step. Hover the Tooltip group, then click
    // its own Edit... (disambiguated by the 'Use as Group Tooltip' sibling).
    const opened = await page.evaluate(async () => {
      const items = Array.from(document.querySelectorAll('.d4-menu-item'));
      let ttNode: Element | null = null;
      for (const node of items) {
        const own = Array.from(node.querySelectorAll('.d4-menu-item-label')).find((l) => l.closest('.d4-menu-item') === node);
        if (own && (own.textContent ?? '').trim() === 'Tooltip') { ttNode = node; break; }
      }
      if (!ttNode) return {err: 'no Tooltip group'};
      const rect = ttNode.getBoundingClientRect();
      ttNode.dispatchEvent(new MouseEvent('mouseover', {bubbles: true, clientX: rect.x + rect.width / 2, clientY: rect.y + rect.height / 2}));
      ttNode.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: rect.x + rect.width / 2, clientY: rect.y + rect.height / 2}));
      await new Promise((r) => setTimeout(r, 800));
      // find Edit... whose sibling set includes 'Use as Group Tooltip'
      const edits = Array.from(document.querySelectorAll('.d4-menu-item-label')).filter((l) => (l.textContent ?? '').trim() === 'Edit...');
      let target: {item: Element; container: Element | null} | null = null;
      for (const l of edits) {
        const parent = l.closest('.d4-menu-item')?.parentElement ?? null;
        const sibs = parent ? Array.from(parent.querySelectorAll('.d4-menu-item-label')).map((s) => (s.textContent ?? '').trim()) : [];
        if (sibs.includes('Use as Group Tooltip')) { target = {item: l.closest('.d4-menu-item')!, container: parent}; break; }
      }
      if (!target) return {err: 'no tooltip Edit...'};
      let cc: Element | null = target.container;
      while (cc && cc !== document.body) {
        if (getComputedStyle(cc as Element).display === 'none') (cc as HTMLElement).style.display = 'block';
        cc = cc.parentElement;
      }
      await new Promise((r) => setTimeout(r, 200));
      const r2 = target.item.getBoundingClientRect();
      for (const t of ['mousedown', 'mouseup', 'click'])
        target.item.dispatchEvent(new MouseEvent(t, {bubbles: true, cancelable: true, clientX: (r2.x || 10) + 5, clientY: (r2.y || 10) + 5, view: window}));
      await new Promise((r) => setTimeout(r, 1000));
      const dlg = document.querySelector('.d4-dialog[name="dialog-Edit-Tooltip"]');
      return {
        title: (dlg?.querySelector('.d4-dialog-header, .d4-dialog-title')?.textContent ?? '').trim(),
        present: !!dlg,
        hasCancel: !!dlg?.querySelector('[name="button-CANCEL"]'),
      };
    });
    console.log(`[S4] Edit... -> ${JSON.stringify(opened)}`);
    expect((opened as any).present).toBe(true);
    expect((opened as any).title).toBe('Edit Tooltip');
    // Close via CANCEL and confirm the dialog is gone (GROK-19054 dialog save-guard: opens and
    // closes cleanly without errors).
    const closed = await page.evaluate(async () => {
      const dlg = document.querySelector('.d4-dialog[name="dialog-Edit-Tooltip"]');
      const cancel = dlg?.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
      if (cancel) cancel.click();
      await new Promise((r) => setTimeout(r, 500));
      return document.querySelectorAll('.d4-dialog[name="dialog-Edit-Tooltip"]').length;
    });
    expect(closed).toBe(0);
  });

  await softStep('Scenario 4 Step 6-7 — Tooltip > Properties... pushes scatter props to the context panel, NO dialog', async () => {
    // Re-open the menu on the same cell and navigate Tooltip > Properties...
    await page.mouse.move(geom.rootX + 5, geom.rootY + geom.headerH + 200);
    await page.waitForTimeout(300);
    const c = cellCenter(geom, xiHeight, yiAge);
    await page.mouse.click(c.x, c.y, {button: 'right'});
    await page.waitForTimeout(700);
    const r = await page.evaluate(async () => {
      const dialogsBefore = document.querySelectorAll('.d4-dialog').length;
      const items = Array.from(document.querySelectorAll('.d4-menu-item'));
      let ttNode: Element | null = null;
      for (const node of items) {
        const own = Array.from(node.querySelectorAll('.d4-menu-item-label')).find((l) => l.closest('.d4-menu-item') === node);
        if (own && (own.textContent ?? '').trim() === 'Tooltip') { ttNode = node; break; }
      }
      if (!ttNode) return {err: 'no Tooltip group'};
      const rect = ttNode.getBoundingClientRect();
      ttNode.dispatchEvent(new MouseEvent('mouseover', {bubbles: true, clientX: rect.x + rect.width / 2, clientY: rect.y + rect.height / 2}));
      ttNode.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: rect.x + rect.width / 2, clientY: rect.y + rect.height / 2}));
      await new Promise((r) => setTimeout(r, 800));
      const props = Array.from(document.querySelectorAll('.d4-menu-item-label')).filter((l) => (l.textContent ?? '').trim() === 'Properties...');
      let target: {item: Element; container: Element | null} | null = null;
      for (const l of props) {
        const parent = l.closest('.d4-menu-item')?.parentElement ?? null;
        const sibs = parent ? Array.from(parent.querySelectorAll('.d4-menu-item-label')).map((s) => (s.textContent ?? '').trim()) : [];
        if (sibs.includes('Use as Group Tooltip')) { target = {item: l.closest('.d4-menu-item')!, container: parent}; break; }
      }
      if (!target) return {err: 'no tooltip Properties...'};
      let cc: Element | null = target.container;
      while (cc && cc !== document.body) {
        if (getComputedStyle(cc as Element).display === 'none') (cc as HTMLElement).style.display = 'block';
        cc = cc.parentElement;
      }
      await new Promise((r) => setTimeout(r, 200));
      const r2 = target.item.getBoundingClientRect();
      for (const t of ['mousedown', 'mouseup', 'click'])
        target.item.dispatchEvent(new MouseEvent(t, {bubbles: true, cancelable: true, clientX: (r2.x || 10) + 5, clientY: (r2.y || 10) + 5, view: window}));
      await new Promise((r) => setTimeout(r, 1000));
      const dialogsAfter = document.querySelectorAll('.d4-dialog').length;
      const panelText = document.querySelector('.grok-prop-panel')?.textContent ?? '';
      return {dialogsBefore, dialogsAfter, hasScatterProps: /Zoom and Filter|Marker|Row Source/.test(panelText)};
    });
    console.log(`[S4] Properties... -> ${JSON.stringify(r)}`);
    // No modal opens (dialog count does not increase) ...
    expect((r as any).dialogsAfter).toBe((r as any).dialogsBefore);
    // ... instead the embedded tooltip scatter's property panel appears in the context panel.
    expect((r as any).hasScatterProps).toBe(true);
    // Close the menu if still open.
    await page.keyboard.press('Escape');
    await page.evaluate(() => document.body.click());
    await page.waitForTimeout(300);
  });

  // ### Scenario 5: Open as table — column pairs and values match matrix
  await softStep('Scenario 5 Step 3 — Open as table yields a table whose HEIGHT/AGE value matches the matrix (GROK-19053)', async () => {
    const expected: number = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      const df = grok.shell.tv.dataFrame;
      return cp.getCorrelation(df.col('HEIGHT'), df.col('AGE'));
    });
    const beforeTables: string[] = await page.evaluate(() =>
      Array.from(grok.shell.tableViews).map((v: any) => v.table?.name));
    // Real right-click the HEIGHT/AGE cell and activate "Open as table" (exact lowercase text).
    const c = cellCenter(geom, xiHeight, yiAge);
    await page.mouse.click(c.x, c.y, {button: 'right'});
    await page.waitForTimeout(700);
    await page.evaluate(async () => {
      const items = Array.from(document.querySelectorAll('.d4-menu-item'));
      for (const node of items) {
        const own = Array.from(node.querySelectorAll('.d4-menu-item-label')).find((l) => l.closest('.d4-menu-item') === node);
        if (own && (own.textContent ?? '').trim() === 'Open as table') {
          const r = own.getBoundingClientRect();
          for (const t of ['mousedown', 'mouseup', 'click'])
            own.dispatchEvent(new MouseEvent(t, {bubbles: true, cancelable: true, clientX: r.x + 5, clientY: r.y + 5, view: window}));
          break;
        }
      }
      await new Promise((r) => setTimeout(r, 1800));
    });
    const r = await page.evaluate(({before, exp, tol}) => {
      const tvs = Array.from(grok.shell.tableViews);
      const newTv = tvs.find((v: any) => !before.includes(v.table?.name)) as any;
      if (!newTv || !newTv.table) return {opened: false};
      const t = newTv.table;
      const cols: string[] = Array.from({length: t.columns.length}, (_: any, i: number) => t.columns.byIndex(i).name);
      // The matrix table is symmetric: a 'name' column plus one value column per X column.
      // Find the row whose name === 'HEIGHT' and read its 'AGE' value (== the matrix cell value).
      const nameCol = t.col('name');
      let hAge: number | null = null;
      if (nameCol && t.col('AGE')) {
        for (let i = 0; i < t.rowCount; i++) {
          if (nameCol.get(i) === 'HEIGHT') { hAge = t.col('AGE').get(i); break; }
        }
      }
      return {
        opened: true,
        tableName: t.table?.name ?? t.name,
        cols,
        hasHeight: cols.includes('HEIGHT') || (nameCol ? Array.from({length: t.rowCount}, (_: any, i: number) => nameCol.get(i)).includes('HEIGHT') : false),
        hasAgeCol: cols.includes('AGE'),
        hAge,
        matches: hAge != null && Math.abs(hAge - exp) <= tol,
      };
    }, {before: beforeTables, exp: expected, tol: TOL});
    console.log(`[S5] Open as table -> ${JSON.stringify(r)} expected=${expected}`);
    // A new table view opened ...
    expect((r as any).opened).toBe(true);
    // ... whose fields identify the HEIGHT/AGE pair ...
    expect((r as any).hasHeight).toBe(true);
    expect((r as any).hasAgeCol).toBe(true);
    // ... and whose correlation value matches the matrix within tolerance (GROK-19053 content guard).
    expect(Number.isFinite((r as any).hAge)).toBe(true);
    expect((r as any).matches).toBe(true);
    // Close the produced table view.
    await page.evaluate(({before}) => {
      const tv = Array.from(grok.shell.tableViews).find((v: any) => !before.includes(v.table?.name)) as any;
      if (tv) tv.close();
    }, {before: beforeTables});
  });

  v.finishSpec();
});
