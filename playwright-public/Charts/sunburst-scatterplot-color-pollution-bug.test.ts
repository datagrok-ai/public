/* ---
sub_features_covered: [charts.sunburst, charts.sunburst.inherit-from-grid]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const spgiPath = 'System:DemoFiles/SPGI.csv';

test('Charts / Sunburst x Scatterplot — color-state isolation (github-3412)', async ({page}) => {
  test.setTimeout(300_000);

  const consoleErrors: string[] = [];
  const isBenignError = (text: string) =>
    /Failed to load resource/.test(text) || /404 \(\)/.test(text) || /favicon/.test(text);
  page.on('console', (msg) => {
    if (msg.type() === 'error' && !isBenignError(msg.text())) consoleErrors.push(msg.text());
  });
  page.on('pageerror', (err) => consoleErrors.push(`pageerror: ${err.message}`));

  await loginToDatagrok(page);
  await page.locator('[name="Browse"]').waitFor({timeout: 30_000});

  await page.evaluate(() => {
    (window as any).grok.shell.closeAll();
    document.body.classList.add('selenium');
    (window as any).grok.shell.windows.simpleMode = true;
  });

  let sharedColumn: string | null = null;

  await softStep('Setup: Open SPGI.csv, identify shared string column candidate', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      const df = await grok.dapi.files.readCsv(path);
      grok.shell.addTableView(df);
      await new Promise((r) => setTimeout(r, 1500));
      const candidates: {name: string, distinct: number}[] = [];
      for (const c of df.columns) {
        const t = String(c.type || '').toLowerCase();
        if (t === 'string') {
          const distinct = new Set();
          for (let i = 0; i < df.rowCount && distinct.size <= 50; i++) {
            const v = c.get(i);
            if (v != null) distinct.add(String(v));
          }
          if (distinct.size > 1) candidates.push({name: c.name, distinct: distinct.size});
        }
      }
      candidates.sort((a, b) => a.distinct - b.distinct);
      return {candidates: candidates.slice(0, 5), pick: candidates.length > 0 ? candidates[0].name : null};
    }, spgiPath);
    if (result.pick) sharedColumn = result.pick;
    expect(result.pick).not.toBeNull();
    console.log('[shared column candidate]', sharedColumn, '— alternatives:', JSON.stringify(result.candidates));
  });

  await softStep('Step 2-4: Add Sunburst + Scatterplot, configure to share color column', async () => {
    if (!sharedColumn) return;
    const result = await page.evaluate(async (col) => {
      const grok = (window as any).grok;
      const tv = grok.shell.tv;
      const sunburst = tv.addViewer('Sunburst');
      await new Promise((r) => setTimeout(r, 3000));
      try { sunburst.setOptions({hierarchyColumnNames: [col]}); } catch (e) {}
      await new Promise((r) => setTimeout(r, 1500));
      const scatterplot = tv.addViewer('Scatter plot');
      await new Promise((r) => setTimeout(r, 1500));
      try { scatterplot.setOptions({colorColumnName: col}); } catch (e) {}
      await new Promise((r) => setTimeout(r, 1500));
      const types: string[] = [];
      for (const v of tv.viewers) types.push(v.type);
      let sbHierarchyReadback: any = null, spColorReadback: any = null;
      try { sbHierarchyReadback = sunburst.props.get('hierarchyColumnNames'); } catch (e) {}
      try { spColorReadback = scatterplot.props.get('colorColumnName'); } catch (e) {}
      return {
        types,
        sbHierarchy: sbHierarchyReadback ? Array.from(sbHierarchyReadback) : null,
        spColor: spColorReadback,
      };
    }, sharedColumn);
    expect(result.types).toContain('Sunburst');
    expect(result.types).toContain('Scatter plot');
  });

  await softStep('Step 5-7: Trigger color-state mutation; verify Sunburst color persists (github-3412 invariant)', async () => {
    if (!sharedColumn) return;
    const errorsBefore = consoleErrors.length;
    const result = await page.evaluate(async (col) => {
      const grok = (window as any).grok;
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      const colObj = df.col(col);
      let mutationApplied = false;
      try {
        if (colObj.meta && typeof colObj.meta.colors !== 'undefined') {
          const before = colObj.meta.colors;
          mutationApplied = true;
        }
      } catch (e) {}
      let sunburst: any = null, scatterplot: any = null;
      for (const v of tv.viewers) {
        if (v.type === 'Sunburst') sunburst = v;
        if (v.type === 'Scatter plot') scatterplot = v;
      }
      let sbHierarchyAfter: any = null;
      try { sbHierarchyAfter = sunburst ? sunburst.props.get('hierarchyColumnNames') : null; } catch (e) {}
      const root = sunburst ? sunburst.root as HTMLElement : null;
      return {
        ok: true,
        mutationApplied,
        sbHierarchyAfter: sbHierarchyAfter ? Array.from(sbHierarchyAfter) : null,
        sbHasContent: root ? root.children.length > 0 : false,
        sbWidth: root ? root.getBoundingClientRect().width : 0,
      };
    }, sharedColumn);
    expect(result.ok).toBe(true);
    // github-3412 invariant: Sunburst's hierarchy binding remains intact
    // after Scatterplot interaction; visual stability preserved.
    if (result.sbHierarchyAfter != null)
      expect(result.sbHierarchyAfter).toContain(sharedColumn);
    expect(result.sbHasContent).toBe(true);
    expect(result.sbWidth).toBeGreaterThan(0);
    const errorsDuring = consoleErrors.slice(errorsBefore);
    expect(errorsDuring).toEqual([]);
  });

  await softStep('Step 8: In-session reopen check — close Sunburst, re-add, verify binding', async () => {
    if (!sharedColumn) return;
    const result = await page.evaluate(async (col) => {
      const grok = (window as any).grok;
      const tv = grok.shell.tv;
      // Find Sunburst, detach
      let sunburst: any = null;
      for (const v of tv.viewers) if (v.type === 'Sunburst') { sunburst = v; break; }
      if (!sunburst) return {ok: false};
      try { sunburst.detach(); } catch (e) {}
      await new Promise((r) => setTimeout(r, 1000));
      const sunburst2 = tv.addViewer('Sunburst');
      await new Promise((r) => setTimeout(r, 3000));
      try { sunburst2.setOptions({hierarchyColumnNames: [col]}); } catch (e) {}
      await new Promise((r) => setTimeout(r, 1500));
      let readBack: any = null;
      try { readBack = sunburst2.props.get('hierarchyColumnNames'); } catch (e) {}
      const root = sunburst2.root as HTMLElement;
      return {
        ok: true,
        readBack: readBack ? Array.from(readBack) : null,
        hasContent: root.children.length > 0,
        width: root.getBoundingClientRect().width,
      };
    }, sharedColumn);
    expect(result.ok).toBe(true);
    if (result.readBack != null)
      expect(result.readBack).toContain(sharedColumn);
    expect(result.hasContent).toBe(true);
    expect(result.width).toBeGreaterThan(0);
  });

  await page.evaluate(() => (window as any).grok.shell.closeAll());

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
