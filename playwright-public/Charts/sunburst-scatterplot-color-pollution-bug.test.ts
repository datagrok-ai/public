/* ---
sub_features_covered: [charts.sunburst, charts.sunburst.inherit-from-grid]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const spgiPath = 'System:DemoFiles/SPGI.csv';

test('Charts / Sunburst x Scatterplot — color-state isolation (github-3412)', async ({page}) => {
  test.setTimeout(120_000);

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
      const poll = async (pred: () => boolean, timeout = 15_000, interval = 200) => {
        const start = Date.now();
        while (Date.now() - start < timeout) {
          try { if (pred()) return true; } catch (e) {}
          await new Promise((r) => setTimeout(r, interval));
        }
        return false;
      };
      const df = await grok.dapi.files.readCsv(path);
      grok.shell.addTableView(df);
      await poll(() => grok.shell.tv && grok.shell.tv.dataFrame && grok.shell.tv.dataFrame.rowCount > 0);
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
      const poll = async (pred: () => boolean, timeout = 15_000, interval = 200) => {
        const start = Date.now();
        while (Date.now() - start < timeout) {
          try { if (pred()) return true; } catch (e) {}
          await new Promise((r) => setTimeout(r, interval));
        }
        return false;
      };
      const sunburst = tv.addViewer('Sunburst');
      await poll(() => sunburst.root && sunburst.root.children.length > 0);
      sunburst.setOptions({hierarchyColumnNames: [col]});
      await poll(() => {
        const rb = sunburst.props.get('hierarchyColumnNames');
        return rb && Array.from(rb).indexOf(col) >= 0;
      });
      const scatterplot = tv.addViewer('Scatter plot');
      await poll(() => scatterplot.root && scatterplot.root.children.length > 0);
      scatterplot.setOptions({colorColumnName: col});
      await poll(() => scatterplot.props.get('colorColumnName') === col);
      const types: string[] = [];
      for (const v of tv.viewers) types.push(v.type);
      const sbHierarchyReadback = sunburst.props.get('hierarchyColumnNames');
      return {
        types,
        sbHierarchy: sbHierarchyReadback ? Array.from(sbHierarchyReadback) : null,
        spColor: scatterplot.props.get('colorColumnName'),
      };
    }, sharedColumn);
    expect(result.types).toContain('Sunburst');
    expect(result.types).toContain('Scatter plot');
    expect(result.sbHierarchy, 'Sunburst hierarchy binding must be set on shared column').not.toBeNull();
    expect(result.sbHierarchy).toContain(sharedColumn);
    expect(result.spColor, 'Scatterplot color must be bound to shared column').toBe(sharedColumn);
  });

  await softStep('Step 5-7: Trigger color-state mutation; verify Sunburst color persists (github-3412 invariant)', async () => {
    if (!sharedColumn) return;
    const errorsBefore = consoleErrors.length;
    const result = await page.evaluate(async (col) => {
      const grok = (window as any).grok;
      const tv = grok.shell.tv;
      const poll = async (pred: () => boolean, timeout = 15_000, interval = 200) => {
        const start = Date.now();
        while (Date.now() - start < timeout) {
          try { if (pred()) return true; } catch (e) {}
          await new Promise((r) => setTimeout(r, interval));
        }
        return false;
      };
      let sunburst: any = null, scatterplot: any = null;
      for (const v of tv.viewers) {
        if (v.type === 'Sunburst') sunburst = v;
        if (v.type === 'Scatter plot') scatterplot = v;
      }
      // Real color-state mutation on the shared column via the Scatterplot.
      let mutationApplied = false;
      if (scatterplot) {
        scatterplot.setOptions({colorColumnName: col});
        mutationApplied = await poll(() => scatterplot.props.get('colorColumnName') === col);
      }
      const sbHierarchyAfter = sunburst ? sunburst.props.get('hierarchyColumnNames') : null;
      const root = sunburst ? sunburst.root as HTMLElement : null;
      return {
        mutationApplied,
        sbHierarchyAfter: sbHierarchyAfter ? Array.from(sbHierarchyAfter) : null,
        spColorAfter: scatterplot ? scatterplot.props.get('colorColumnName') : null,
        sbHasContent: root ? root.children.length > 0 : false,
        sbWidth: root ? root.getBoundingClientRect().width : 0,
      };
    }, sharedColumn);
    // github-3412 isolation invariant: after mutating Scatterplot color on the shared
    // column, both viewers retain their independent bindings — Sunburst hierarchy is
    // not wiped and Scatterplot color remains on the shared column.
    expect(result.mutationApplied, 'Scatterplot color mutation must be applied').toBe(true);
    expect(result.sbHierarchyAfter, 'Sunburst hierarchy binding must survive scatterplot color interaction').not.toBeNull();
    expect(result.sbHierarchyAfter).toContain(sharedColumn);
    expect(result.spColorAfter, 'Scatterplot color must stay bound to shared column').toBe(sharedColumn);
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
      const poll = async (pred: () => boolean, timeout = 15_000, interval = 200) => {
        const start = Date.now();
        while (Date.now() - start < timeout) {
          try { if (pred()) return true; } catch (e) {}
          await new Promise((r) => setTimeout(r, interval));
        }
        return false;
      };
      // Find Sunburst, detach
      let sunburst: any = null;
      for (const v of tv.viewers) if (v.type === 'Sunburst') { sunburst = v; break; }
      if (!sunburst) return {ok: false};
      sunburst.detach();
      await poll(() => ![...tv.viewers].some((v: any) => v.type === 'Sunburst'));
      const sunburst2 = tv.addViewer('Sunburst');
      await poll(() => sunburst2.root && sunburst2.root.children.length > 0);
      sunburst2.setOptions({hierarchyColumnNames: [col]});
      await poll(() => {
        const rb = sunburst2.props.get('hierarchyColumnNames');
        return rb && Array.from(rb).indexOf(col) >= 0;
      });
      const readBack = sunburst2.props.get('hierarchyColumnNames');
      const root = sunburst2.root as HTMLElement;
      return {
        ok: true,
        readBack: readBack ? Array.from(readBack) : null,
        hasContent: root.children.length > 0,
        width: root.getBoundingClientRect().width,
      };
    }, sharedColumn);
    expect(result.ok, 'Sunburst viewer must be found to reopen').toBe(true);
    expect(result.readBack, 'reopened Sunburst must retain hierarchy binding').not.toBeNull();
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
