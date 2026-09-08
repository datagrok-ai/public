/* ---
realizes: [correlationplot.cp.cell-interactions, correlationplot.int.ignore-double-click]
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const TOL = 1e-3;

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

async function refreshRoot(page: Page, g: Geometry): Promise<void> {
  const r = await page.evaluate(() => {
    const R = document.querySelector('[name="viewer-Correlation-plot"]')!.getBoundingClientRect();
    return {x: R.x, y: R.y};
  });
  g.rootX = r.x;
  g.rootY = r.y;
}

async function tooltipShown(page: Page): Promise<boolean> {
  return page.evaluate(() => Array.from(document.querySelectorAll('.d4-tooltip'))
    .some((t) => getComputedStyle(t).display === 'block'));
}

async function dismissTooltip(page: Page, g: Geometry): Promise<void> {
  await page.mouse.move(Math.max(4, g.rootX - 60), g.rootY + g.headerH + 100);
  await v.pollValue(() => tooltipShown(page), (shown) => !shown, 400, 50);
}

async function hoverCell(page: Page, g: Geometry, xi: number, yi: number): Promise<void> {
  const c = cellCenter(g, xi, yi);
  await page.mouse.move(c.x - g.cellW, c.y, {steps: 3});
  await page.mouse.move(c.x, c.y, {steps: 4});
}

async function waitMenuOpen(page: Page, capMs: number): Promise<boolean> {
  return v.pollValue(() => page.evaluate(() => Array.from(document.querySelectorAll('.d4-menu-popup'))
    .some((el) => { const b = (el as HTMLElement).getBoundingClientRect(); return b.width > 0 && b.height > 0; })),
  (open) => open, capMs, 50);
}

async function waitMenuClosed(page: Page, capMs: number): Promise<void> {
  await v.pollValue(() => page.evaluate(() => Array.from(document.querySelectorAll('.d4-menu-popup'))
    .some((el) => { const b = (el as HTMLElement).getBoundingClientRect(); return b.width > 0 && b.height > 0; })),
  (open) => !open, capMs, 50);
}

async function waitForRTooltip(page: Page, timeoutMs = 6000): Promise<void> {
  await page.waitForFunction(() => {
    const text = Array.from(document.querySelectorAll('.d4-tooltip'))
      .map((t) => t.textContent ?? '').join(' ');
    return /Pearson R:\s*-?\d+\.\d+/.test(text);
  }, undefined, {timeout: timeoutMs});
}

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

async function scatterCount(page: Page): Promise<number> {
  return page.evaluate(() => grok.shell.tv.viewers.filter((x: any) => x.type === 'Scatter plot').length);
}

test('Correlation Plot — Cell Interactions', async ({page}) => {
  test.setTimeout(600_000);

  await openDatagrok(page);

  try {
    await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
    await v.addViewerByIcon(page, 'correlation-plot', 'Correlation-plot', 10000);
    await armClickListener(page);

    const base = await readGeometry(page);
    const xiHeight = base.xCols.indexOf('HEIGHT');
    const yiAge = base.yCols.indexOf('AGE');
    const xiWeight = base.xCols.indexOf('WEIGHT');

    const geom: Geometry = {rootX: base.rootX, rootY: base.rootY, pinnedW: 104, headerH: 60, cellW: base.cellW, rowH: 20};
    await softStep('Setup Step 5 — calibrate cell geometry via a probe click', async () => {
      let landed = false;
      for (let attempt = 0; attempt < 5 && !landed; attempt++) {
        await page.evaluate(() => { (window as any).__corrClicks = []; });
        const c = cellCenter(geom, xiHeight, yiAge);
        await page.mouse.click(c.x, c.y);

        const ev = await v.pollValue(() => lastClick(page), (e) => e !== null, 400, 50);
        console.log(`[Setup] probe attempt ${attempt} at (${Math.round(c.x)},${Math.round(c.y)}) pinnedW=${geom.pinnedW} headerH=${geom.headerH} -> ${JSON.stringify(ev)}`);
        if (ev && ((ev.col1 === 'HEIGHT' && ev.col2 === 'AGE') || (ev.col1 === 'AGE' && ev.col2 === 'HEIGHT'))) {
          landed = true;
          break;
        }

        if (ev && ev.col1 && ev.col1 !== 'HEIGHT') {
          const idx = base.xCols.indexOf(ev.col1);
          if (idx >= 0) geom.pinnedW += (xiHeight - idx) * geom.cellW;
        } else if (!ev) {
          geom.headerH += 4;
        }
      }

      const ev = await lastClick(page);
      expect(ev).not.toBeNull();
      expect([ev!.col1, ev!.col2].sort()).toEqual(['AGE', 'HEIGHT']);
    });

    const baselineViewerCount = await scatterCount(page);

    await softStep('Scenario 1 Step 3 — click HEIGHT/AGE fires event with correct pair and Pearson value', async () => {
      const expected: number = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        return DG.Stats.fromColumn(df.col('HEIGHT')).corr(df.col('AGE'));
      });
      await page.evaluate(() => { (window as any).__corrClicks = []; });
      const c = cellCenter(geom, xiHeight, yiAge);
      await page.mouse.click(c.x, c.y);
      const ev = await v.pollValue(() => lastClick(page), (e) => e !== null, 500, 50);
      console.log(`[S1] click event=${JSON.stringify(ev)} expectedPearson=${expected}`);
      expect(ev).not.toBeNull();

      expect([ev!.col1, ev!.col2].sort()).toEqual(['AGE', 'HEIGHT']);

      expect(Number.isFinite(ev!.value)).toBe(true);
      expect(Math.abs(ev!.value - expected)).toBeLessThanOrEqual(TOL);
    });

    await softStep('Scenario 1 Step 4 — context panel shows columnsCorrelation "<c1> vs <c2>" with a Scatter plot pane', async () => {
      const readPanel = () => page.evaluate(() => {
        const vsHeaders = Array.from(document.querySelectorAll('*'))
          .filter((e) => e.children.length === 0 && /\bvs\b/.test(e.textContent ?? '') && (e.textContent ?? '').length < 40)
          .map((e) => (e.textContent ?? '').trim());
        const accordions = Array.from(document.querySelectorAll('.d4-accordion-pane-header, .d4-accordion-pane-title'))
          .map((e) => (e.textContent ?? '').trim());
        return {vsHeaders: [...new Set(vsHeaders)], hasScatterPane: accordions.includes('Scatter plot')};
      });
      const isPair = (hs: string[]) => hs.some((h) => /^(HEIGHT vs AGE|AGE vs HEIGHT)$/.test(h));
      const shown = (p: {vsHeaders: string[]; hasScatterPane: boolean}) => isPair(p.vsHeaders) && p.hasScatterPane;
      let panel = await v.pollValue(readPanel, shown, 1500, 100);
      // the click sets the current object, but the panel does not always take the first one
      // (seen 1 run in 2 on the local lane); one more click on the same cell, then the same wait
      if (!shown(panel)) {
        console.log('[S1] panel did not show the pair after the first click; clicking once more');
        const c = cellCenter(geom, xiHeight, yiAge);
        await page.mouse.click(c.x, c.y);
        panel = await v.pollValue(readPanel, shown, 1500, 100);
      }
      console.log(`[S1] panel vsHeaders=${JSON.stringify(panel.vsHeaders)} hasScatterPane=${panel.hasScatterPane}`);

      expect(isPair(panel.vsHeaders)).toBe(true);
      expect(panel.hasScatterPane).toBe(true);
    });

    await softStep('Scenario 2 Step 3-4 — double-click WEIGHT/AGE adds exactly one Scatter plot for the pair', async () => {
      const before = await scatterCount(page);
      const c = cellCenter(geom, xiWeight, yiAge);
      await page.mouse.dblclick(c.x, c.y);
      const r = await v.pollValue(() => page.evaluate(() => {
        const scatters = grok.shell.tv.viewers.filter((x: any) => x.type === 'Scatter plot');
        const newest = scatters[scatters.length - 1];
        return {count: scatters.length, x: newest?.props.xColumnName ?? null, y: newest?.props.yColumnName ?? null};
      }), (res) => res.count === before + 1 && res.x === 'WEIGHT' && res.y === 'AGE', 1500, 50);
      console.log(`[S2] scatter before=${before} after=${r.count} x=${r.x} y=${r.y}`);

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
      await v.pollValue(() => page.evaluate(() =>
        document.querySelectorAll('[name="viewer-Scatter-plot"]').length), (n) => n === 0, 500, 50);
      console.log(`[S2] scatter after close=${after} baseline=${baselineViewerCount}`);
      expect(after).toBe(baselineViewerCount);
    });

    await softStep('Scenario 2 Step 7-9 — Ignore Double Click suppresses the scatter (correlationplot.int.ignore-double-click)', async () => {
      await v.setViewerProps(page, 'Correlation plot', [{set: {ignoreDoubleClick: true}, wait: 400}]);
      const before = await scatterCount(page);
      const c = cellCenter(geom, xiWeight, yiAge);
      await page.mouse.dblclick(c.x, c.y);
      // a hold for a scatter that must not appear; the un-suppressed add above lands well inside it
      const after = await v.pollValue(() => scatterCount(page), (n) => n !== before, 1000, 100);

      await page.evaluate(() => {
        const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
        cp.props.ignoreDoubleClick = false;
      });
      console.log(`[S2] ignoreDoubleClick suppression before=${before} after=${after}`);

      expect(after).toBe(before);
    });

    await softStep('Scenario 3 Step 3 — hover Cell A shows exactly one tooltip with "Pearson R: <value>"', async () => {
      const expected: number = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        return DG.Stats.fromColumn(df.col('HEIGHT')).corr(df.col('AGE'));
      });
      await refreshRoot(page, geom);

      await dismissTooltip(page, geom);
      await hoverCell(page, geom, xiHeight, yiAge);
      await waitForRTooltip(page);
      const r = await page.evaluate(() => {
        const tips = Array.from(document.querySelectorAll('.d4-tooltip'));

        const combined = tips.map((t) => t.textContent ?? '').join(' ').replace(/\s+/g, ' ').trim();
        return {
          totalCount: tips.length,
          text: combined.slice(0, 400),
          hasCanvas: tips.some((t) => !!t.querySelector('canvas')),
        };
      });
      console.log(`[S3] tooltip total=${r.totalCount} text="${r.text}" expected=${expected.toFixed(3)}`);

      expect(r.totalCount).toBe(1);

      const m = r.text.match(/Pearson R:\s*(-?\d+\.\d+)/);
      expect(m).not.toBeNull();
      expect(Math.abs(parseFloat(m![1]) - expected)).toBeLessThanOrEqual(TOL);

      expect(r.hasCanvas).toBe(true);
    });

    await softStep('Scenario 3 Step 5 — the embedded scatter canvas differs between Cell A and Cell B (GROK-20125)', async () => {
      await refreshRoot(page, geom);
      await dismissTooltip(page, geom);
      await hoverCell(page, geom, xiHeight, yiAge);
      await waitForRTooltip(page);
      // the tooltip's scatter is one reused instance repainted on a later frame than the text, so
      // Cell A is read once its canvas has ink, and Cell B once that ink has moved
      const snapA = await page.evaluate(async () => {
        const w = window as any;
        const tipCanvas = () => Array.from(document.querySelectorAll('.d4-tooltip'))
          .map((t) => t.querySelector('canvas') as HTMLCanvasElement | null)
          .find((c) => !!c && !!c.width && !!c.height) ?? null;
        const histogram = (cv: HTMLCanvasElement) => {
          const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
          const colors: Record<number, number> = {};
          for (let i = 0; i < img.length; i += 4) {
            const key = (img[i] << 16) | (img[i + 1] << 8) | img[i + 2];
            colors[key] = (colors[key] ?? 0) + 1;
          }
          return colors;
        };
        w.__tipHistogram = histogram;
        const inked = (colors: Record<number, number>) => Object.keys(colors).some((k) => Number(k) !== 0xffffff);
        try {
          const colors = await w.__poll(() => { const cv = tipCanvas(); return cv ? histogram(cv) : {}; },
            inked, 1500, 50);
          const cv = tipCanvas();
          if (!cv) return null;
          w.__tipSnapA = colors;
          return {w: cv.width, h: cv.height, inked: inked(colors)};
        } catch { return null; }
      });

      const xiWeightB = base.xCols.indexOf('WEIGHT');
      const yiStarted = base.yCols.indexOf('STARTED');
      await dismissTooltip(page, geom);
      await hoverCell(page, geom, xiWeightB, yiStarted);
      await waitForRTooltip(page);
      const r = await page.evaluate(async () => {
        const w = window as any;
        const tips = () => Array.from(document.querySelectorAll('.d4-tooltip'));
        const tipCanvas = () => tips()
          .map((t) => t.querySelector('canvas') as HTMLCanvasElement | null)
          .find((c) => !!c && !!c.width && !!c.height) ?? null;
        const prev = w.__tipSnapA as Record<number, number> | undefined;
        if (!prev) return {tooltipCount: tips().length, deltaPx: -1};
        const delta = (colors: Record<number, number>) => {
          let d = 0;
          const keys = new Set([...Object.keys(colors), ...Object.keys(prev)].map(Number));
          for (const k of keys) d += Math.abs((colors[k] ?? 0) - (prev[k] ?? 0));
          return d;
        };
        try {
          const deltaPx = await w.__poll(() => { const cv = tipCanvas(); return cv ? delta(w.__tipHistogram(cv)) : -1; },
            (d: number) => d > 0, 1500, 50);
          return {tooltipCount: tips().length, deltaPx};
        } catch { return {tooltipCount: tips().length, deltaPx: -1}; }
      });
      console.log(`[S3] Cell A snap=${JSON.stringify(snapA)} Cell B tooltipCount=${r.tooltipCount} deltaPx=${r.deltaPx}`);

      expect(r.tooltipCount).toBe(1);
      expect(r.deltaPx).toBeGreaterThanOrEqual(0);
      expect(r.deltaPx).toBeGreaterThan(0);
    });

    await softStep('Scenario 3 Step 6 — hover the pinned row-header name column shows a column-statistics tooltip (not an R tooltip)', async () => {
      await refreshRoot(page, geom);
      await dismissTooltip(page, geom);
      const nameX = geom.rootX + geom.pinnedW - 45;
      const nameY = geom.rootY + geom.headerH + geom.rowH * 1.5;

      await page.mouse.move(nameX, nameY + geom.rowH, {steps: 3});
      await page.mouse.move(nameX, nameY, {steps: 4});

      await page.waitForFunction(() => {
        const t = Array.from(document.querySelectorAll('.d4-tooltip'))
          .map((tip) => tip.textContent ?? '').join(' ');
        return /(min:|avg:|max:|std:|nulls:)/.test(t);
      }, undefined, {timeout: 6000});
      const r = await page.evaluate(() => {
        const text = Array.from(document.querySelectorAll('.d4-tooltip'))
          .map((tip) => tip.textContent ?? '').join(' ').replace(/\s+/g, ' ').trim();
        return {hasText: text.length > 0, text: text.slice(0, 200), hasRLine: /\bR:\s*-?\d/.test(text)};
      });
      console.log(`[S3] row-header tooltip hasText=${r.hasText} hasRLine=${r.hasRLine} text="${r.text}"`);

      expect(r.hasText).toBe(true);
      expect(r.hasRLine).toBe(false);
    });

    await softStep('Scenario 4 Step 2-4 — right-click menu has Show Pearson R, Columns > X/Y, and a Tooltip submenu with Edit.../Visible/Properties...', async () => {
      await page.mouse.move(geom.rootX + 5, geom.rootY + geom.headerH + 200);
      await v.pollValue(() => tooltipShown(page), (shown) => !shown, 400, 50);
      const c = cellCenter(geom, xiHeight, yiAge);
      await page.mouse.click(c.x, c.y, {button: 'right'});
      await waitMenuOpen(page, 700);
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

      const sub = await page.evaluate(async () => {
        const w = window as any;
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

        const laidOut = () => {
          const a = Array.from(document.querySelectorAll('.d4-menu-item-label'))
            .find((l) => (l.textContent ?? '').trim() === 'Use as Group Tooltip');
          const b = a?.closest('.d4-menu-item')?.getBoundingClientRect();
          return !!b && b.width > 0 && b.height > 0;
        };
        await w.__poll(laidOut, (ok: boolean) => ok, 800, 25);

        const anchor = Array.from(document.querySelectorAll('.d4-menu-item-label'))
          .find((l) => (l.textContent ?? '').trim() === 'Use as Group Tooltip');
        const container = anchor?.closest('.d4-menu-item')?.parentElement ?? null;
        const sibs = container
          ? Array.from(container.querySelectorAll('.d4-menu-item-label')).map((s) => (s.textContent ?? '').trim())
          : [];
        return {edit: sibs.includes('Edit...'), visible: sibs.includes('Visible'), properties: sibs.includes('Properties...')};
      });
      console.log(`[S4] tooltip submenu=${JSON.stringify(sub)}`);

      expect((sub as any).edit).toBe(true);
      expect((sub as any).visible).toBe(true);
      expect((sub as any).properties).toBe(true);
    });

    await softStep('Scenario 4 Step 4 — Tooltip > Edit... opens the "Edit Tooltip" modal; CANCEL closes it cleanly (GROK-19054)', async () => {
      const opened = await page.evaluate(async () => {
        const w = window as any;
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

        const laidOut = () => {
          const a = Array.from(document.querySelectorAll('.d4-menu-item-label'))
            .find((l) => (l.textContent ?? '').trim() === 'Use as Group Tooltip');
          const b = a?.closest('.d4-menu-item')?.getBoundingClientRect();
          return !!b && b.width > 0 && b.height > 0;
        };
        await w.__poll(laidOut, (ok: boolean) => ok, 800, 25);

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

        const item = target.item;
        await w.__poll(() => item.getBoundingClientRect().width > 0, (ok: boolean) => ok, 200, 25);
        const r2 = item.getBoundingClientRect();
        for (const t of ['mousedown', 'mouseup', 'click'])
          item.dispatchEvent(new MouseEvent(t, {bubbles: true, cancelable: true, clientX: (r2.x || 10) + 5, clientY: (r2.y || 10) + 5, view: window}));

        const dlg = await w.__poll(() => document.querySelector('.d4-dialog[name="dialog-Edit-Tooltip"]'),
          (d: Element | null) => !!d, 1000, 25);
        return {
          title: (dlg?.querySelector('.d4-dialog-header, .d4-dialog-title')?.textContent ?? '').trim(),
          present: !!dlg,
          hasCancel: !!dlg?.querySelector('[name="button-CANCEL"]'),
        };
      });
      console.log(`[S4] Edit... -> ${JSON.stringify(opened)}`);
      expect((opened as any).present).toBe(true);
      expect((opened as any).title).toBe('Edit Tooltip');

      await page.evaluate(() => {
        const dlg = document.querySelector('.d4-dialog[name="dialog-Edit-Tooltip"]');
        const cancel = dlg?.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
        if (cancel) cancel.click();
        return null;
      });
      const closed = await v.pollValue(() => page.evaluate(() =>
        document.querySelectorAll('.d4-dialog[name="dialog-Edit-Tooltip"]').length), (n) => n === 0, 500, 50);
      expect(closed).toBe(0);
    });

    await softStep('Scenario 4 Step 6-7 — Tooltip > Properties... pushes scatter props to the context panel, NO dialog', async () => {
      await page.mouse.move(geom.rootX + 5, geom.rootY + geom.headerH + 200);
      await v.pollValue(() => tooltipShown(page), (shown) => !shown, 300, 50);
      const c = cellCenter(geom, xiHeight, yiAge);
      await page.mouse.click(c.x, c.y, {button: 'right'});
      await waitMenuOpen(page, 700);
      const r = await page.evaluate(async () => {
        const w = window as any;
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

        const laidOut = () => {
          const a = Array.from(document.querySelectorAll('.d4-menu-item-label'))
            .find((l) => (l.textContent ?? '').trim() === 'Use as Group Tooltip');
          const b = a?.closest('.d4-menu-item')?.getBoundingClientRect();
          return !!b && b.width > 0 && b.height > 0;
        };
        await w.__poll(laidOut, (ok: boolean) => ok, 800, 25);
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

        const item = target.item;
        await w.__poll(() => item.getBoundingClientRect().width > 0, (ok: boolean) => ok, 200, 25);
        const r2 = item.getBoundingClientRect();
        for (const t of ['mousedown', 'mouseup', 'click'])
          item.dispatchEvent(new MouseEvent(t, {bubbles: true, cancelable: true, clientX: (r2.x || 10) + 5, clientY: (r2.y || 10) + 5, view: window}));

        const panelHasScatterProps = () => /Zoom and Filter|Marker|Row Source/.test(document.querySelector('.grok-prop-panel')?.textContent ?? '');
        const hasScatterProps = await w.__poll(panelHasScatterProps, (ok: boolean) => ok, 1000, 50);
        const dialogsAfter = document.querySelectorAll('.d4-dialog').length;
        return {dialogsBefore, dialogsAfter, hasScatterProps};
      });
      console.log(`[S4] Properties... -> ${JSON.stringify(r)}`);

      expect((r as any).dialogsAfter).toBe((r as any).dialogsBefore);
      expect((r as any).hasScatterProps).toBe(true);

      await page.keyboard.press('Escape');
      await page.evaluate(() => document.body.click());
      await waitMenuClosed(page, 300);
    });

    await softStep('Scenario 5 Step 3 — Open as table yields a table whose HEIGHT/AGE value matches the matrix (GROK-19053)', async () => {
      const expected: number = await page.evaluate(() => {
        const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
        const df = grok.shell.tv.dataFrame;
        return cp.getCorrelation(df.col('HEIGHT'), df.col('AGE'));
      });
      const beforeTables: string[] = await page.evaluate(() =>
        Array.from(grok.shell.tableViews).map((v: any) => v.table?.name));

      const c = cellCenter(geom, xiHeight, yiAge);
      await page.mouse.click(c.x, c.y, {button: 'right'});
      await waitMenuOpen(page, 700);
      await page.evaluate(async ({before}) => {
        const w = window as any;
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
        await w.__poll(() => (Array.from(grok.shell.tableViews) as any[]).some((tv: any) => !before.includes(tv.table?.name)),
          (opened: boolean) => opened, 1800, 25);
        return null;
      }, {before: beforeTables});
      const r = await page.evaluate(({before, exp, tol}) => {
        const tvs = Array.from(grok.shell.tableViews);
        const newTv = tvs.find((v: any) => !before.includes(v.table?.name)) as any;
        if (!newTv || !newTv.table) return {opened: false};
        const t = newTv.table;
        const cols: string[] = Array.from({length: t.columns.length}, (_: any, i: number) => t.columns.byIndex(i).name);

        // the opened correlation table names its label column "__name" (with "__t" alongside),
        // not "name" — resolve either so the row lookup below actually finds HEIGHT
        const nameCol = t.col('__name') ?? t.col('name');
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

      expect((r as any).opened).toBe(true);
      expect((r as any).hasHeight).toBe(true);
      expect((r as any).hasAgeCol).toBe(true);
      expect(Number.isFinite((r as any).hAge)).toBe(true);
      expect((r as any).matches).toBe(true);

      await page.evaluate(({before}) => {
        const tv = Array.from(grok.shell.tableViews).find((v: any) => !before.includes(v.table?.name)) as any;
        if (tv) tv.close();
      }, {before: beforeTables});
    });
  } finally {
    await page.evaluate(() => { const w = window as any; delete w.__corrClicks; delete w.__tipSnapA; }).catch(() => {});
    await v.closeAllAndWait(page);
  }

  v.finishSpec();
});
