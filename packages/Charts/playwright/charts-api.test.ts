/* ---
sub_features_covered: [charts.radar, charts.sunburst, charts.tree, charts.timelines, charts.echart-base, charts.echart-base.table]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';
const aePath = 'System:AppData/Charts/ae.csv';

test('Charts / Viewer API contract', async ({page}) => {
  test.setTimeout(120_000);

  const consoleErrors: string[] = [];
  const isBenignError = (text: string) =>
    /Failed to load resource/.test(text) ||
    /404 \(\)/.test(text) ||
    /favicon/.test(text);
  page.on('console', (msg) => {
    if (msg.type() === 'error' && !isBenignError(msg.text())) consoleErrors.push(msg.text());
  });
  page.on('pageerror', (err) => consoleErrors.push(`pageerror: ${err.message}`));

  await loginToDatagrok(page);

  await page.evaluate(() => {
    document.querySelectorAll('.d4-dialog').forEach((d) => {
      const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
      if (cancel) cancel.click();
    });
    (window as any).grok.shell.closeAll();
    document.body.classList.add('selenium');
    (window as any).grok.shell.windows.simpleMode = true;
  });

  // === Scenario 1: addViewer round-trip across all 4 viewer types ===

  await softStep('Scenario 1: addViewer round-trip across Radar/Sunburst/Tree/Timelines', async () => {
    const result = await page.evaluate(async ([dPath, aPath]) => {
      const grok = (window as any).grok;
      const waitUntil = async (pred: () => any, timeoutMs = 15000, stepMs = 100) => {
        const t0 = Date.now();
        while (Date.now() - t0 < timeoutMs) {
          try { if (await pred()) return true; } catch (e) {}
          await new Promise((r) => setTimeout(r, stepMs));
        }
        return false;
      };
      const viewerReady = (v: any) => waitUntil(() => v.props.getProperties().length > 0);
      const out: any = {viewers: {}};
      const collect = async (df: any, type: string) => {
        const tv = grok.shell.addTableView(df);
        await waitUntil(() => grok.shell.tv?.dataFrame === df);
        const viewer = tv.addViewer(type);
        await viewerReady(viewer);
        const types: string[] = [];
        for (const v of tv.viewers) types.push(v.type);
        return {types, propCount: viewer.props.getProperties().length};
      };
      grok.shell.closeAll();
      await waitUntil(() => (grok.shell.tableViews?.length ?? 0) === 0);
      const dfRadar = await grok.dapi.files.readCsv(dPath);
      out.viewers.Radar = await collect(dfRadar, 'Radar');
      const dfSunburst = await grok.dapi.files.readCsv(dPath);
      out.viewers.Sunburst = await collect(dfSunburst, 'Sunburst');
      const dfTree = await grok.dapi.files.readCsv(dPath);
      out.viewers.Tree = await collect(dfTree, 'Tree');
      const dfTl = await grok.dapi.files.readCsv(aPath);
      out.viewers.Timelines = await collect(dfTl, 'Timelines');
      return out;
    }, [demogPath, aePath]);
    expect(result.viewers.Radar.types).toContain('Radar');
    expect(result.viewers.Sunburst.types).toContain('Sunburst');
    expect(result.viewers.Tree.types).toContain('Tree');
    expect(result.viewers.Timelines.types).toContain('Timelines');
    expect(result.viewers.Radar.propCount).toBeGreaterThan(0);
    expect(result.viewers.Sunburst.propCount).toBeGreaterThan(0);
    expect(result.viewers.Tree.propCount).toBeGreaterThan(0);
    expect(result.viewers.Timelines.propCount).toBeGreaterThan(0);
  });

  // === Scenario 2: getProperties + categories enumeration (Radar) ===

  await softStep('Scenario 2: Radar getProperties categories include Data/Selection/Value/Style/Legend', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      const waitUntil = async (pred: () => any, timeoutMs = 15000, stepMs = 100) => {
        const t0 = Date.now();
        while (Date.now() - t0 < timeoutMs) {
          try { if (await pred()) return true; } catch (e) {}
          await new Promise((r) => setTimeout(r, stepMs));
        }
        return false;
      };
      grok.shell.closeAll();
      await waitUntil(() => (grok.shell.tableViews?.length ?? 0) === 0);
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await waitUntil(() => grok.shell.tv?.dataFrame === df);
      const radar = tv.addViewer('Radar');
      await waitUntil(() => radar.props.getProperties().length > 0);
      const cats = new Set<string>();
      for (const p of radar.props.getProperties()) if (p.category) cats.add(p.category);
      return {categories: Array.from(cats)};
    }, demogPath);
    expect(result.categories).toEqual(expect.arrayContaining(['Data', 'Selection', 'Value', 'Style', 'Legend']));
  });

  // === Scenario 3: setOptions round-trip across all 4 viewers ===

  await softStep('Scenario 3: Radar setOptions(backgroundMinColor) round-trip', async () => {
    const result = await page.evaluate(async () => {
      const grok = (window as any).grok;
      const waitUntil = async (pred: () => any, timeoutMs = 15000, stepMs = 100) => {
        const t0 = Date.now();
        while (Date.now() - t0 < timeoutMs) {
          try { if (await pred()) return true; } catch (e) {}
          await new Promise((r) => setTimeout(r, stepMs));
        }
        return false;
      };
      const tv = grok.shell.tv;
      let radar: any = null;
      for (const v of tv.viewers) if (v.type === 'Radar') { radar = v; break; }
      if (!radar) return {ok: false};
      const newColor = 0xFF123456 | 0;
      radar.setOptions({backgroundMinColor: newColor});
      await waitUntil(() => radar.props.get('backgroundMinColor') === newColor);
      return {ok: true, newColor, echoed: radar.props.get('backgroundMinColor')};
    });
    expect(result.ok).toBe(true);
    expect(result.echoed).toBe(result.newColor);
  });

  await softStep('Scenario 3: Sunburst setOptions(hierarchyColumnNames) round-trip', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      const waitUntil = async (pred: () => any, timeoutMs = 15000, stepMs = 100) => {
        const t0 = Date.now();
        while (Date.now() - t0 < timeoutMs) {
          try { if (await pred()) return true; } catch (e) {}
          await new Promise((r) => setTimeout(r, stepMs));
        }
        return false;
      };
      const eq = (a: any, b: any[]) => a != null && Array.from(a).join(',') === b.join(',');
      grok.shell.closeAll();
      await waitUntil(() => (grok.shell.tableViews?.length ?? 0) === 0);
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await waitUntil(() => grok.shell.tv?.dataFrame === df);
      const sunburst = tv.addViewer('Sunburst');
      await waitUntil(() => sunburst.props.getProperties().length > 0);
      sunburst.setOptions({hierarchyColumnNames: ['SEX', 'RACE']});
      await waitUntil(() => eq(sunburst.props.get('hierarchyColumnNames'), ['SEX', 'RACE']));
      const cols = sunburst.props.get('hierarchyColumnNames');
      return {cols: cols ? Array.from(cols) : null};
    }, demogPath);
    expect(result.cols).toEqual(['SEX', 'RACE']);
  });

  await softStep('Scenario 3: Tree setOptions(hierarchyColumnNames) round-trip', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      const waitUntil = async (pred: () => any, timeoutMs = 15000, stepMs = 100) => {
        const t0 = Date.now();
        while (Date.now() - t0 < timeoutMs) {
          try { if (await pred()) return true; } catch (e) {}
          await new Promise((r) => setTimeout(r, stepMs));
        }
        return false;
      };
      const eq = (a: any, b: any[]) => a != null && Array.from(a).join(',') === b.join(',');
      grok.shell.closeAll();
      await waitUntil(() => (grok.shell.tableViews?.length ?? 0) === 0);
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await waitUntil(() => grok.shell.tv?.dataFrame === df);
      const tree = tv.addViewer('Tree');
      await waitUntil(() => tree.props.getProperties().length > 0);
      tree.setOptions({hierarchyColumnNames: ['CONTROL', 'SEX', 'RACE']});
      await waitUntil(() => eq(tree.props.get('hierarchyColumnNames'), ['CONTROL', 'SEX', 'RACE']));
      const hierarchy = tree.props.get('hierarchyColumnNames');
      return {hierarchy: hierarchy ? Array.from(hierarchy) : null};
    }, demogPath);
    expect(result.hierarchy).toEqual(['CONTROL', 'SEX', 'RACE']);
  });

  await softStep('Scenario 3: Timelines setOptions(colorColumnName) round-trip', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      const waitUntil = async (pred: () => any, timeoutMs = 15000, stepMs = 100) => {
        const t0 = Date.now();
        while (Date.now() - t0 < timeoutMs) {
          try { if (await pred()) return true; } catch (e) {}
          await new Promise((r) => setTimeout(r, stepMs));
        }
        return false;
      };
      grok.shell.closeAll();
      await waitUntil(() => (grok.shell.tableViews?.length ?? 0) === 0);
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await waitUntil(() => grok.shell.tv?.dataFrame === df);
      const timelines = tv.addViewer('Timelines');
      await waitUntil(() => timelines.props.getProperties().length > 0);
      timelines.setOptions({colorColumnName: 'AESOC'});
      await waitUntil(() => timelines.props.get('colorColumnName') === 'AESOC');
      return {color: timelines.props.get('colorColumnName')};
    }, aePath);
    expect(result.color).toBe('AESOC');
  });

  // === Scenario 4: getProperties surface check across all 4 viewers ===

  await softStep('Scenario 4: getProperties surface check Radar/Sunburst/Tree/Timelines', async () => {
    const result = await page.evaluate(async ([dPath, aPath]) => {
      const grok = (window as any).grok;
      const waitUntil = async (pred: () => any, timeoutMs = 15000, stepMs = 100) => {
        const t0 = Date.now();
        while (Date.now() - t0 < timeoutMs) {
          try { if (await pred()) return true; } catch (e) {}
          await new Promise((r) => setTimeout(r, stepMs));
        }
        return false;
      };
      const out: any = {};
      const addAndCollect = async (df: any, type: string) => {
        const tv = grok.shell.addTableView(df);
        await waitUntil(() => grok.shell.tv?.dataFrame === df);
        const viewer = tv.addViewer(type);
        await waitUntil(() => viewer.props.getProperties().length > 0);
        const names: string[] = [];
        for (const p of viewer.props.getProperties()) names.push(p.name);
        return names;
      };
      grok.shell.closeAll();
      await waitUntil(() => (grok.shell.tableViews?.length ?? 0) === 0);
      out.Radar = await addAndCollect(await grok.dapi.files.readCsv(dPath), 'Radar');
      out.Sunburst = await addAndCollect(await grok.dapi.files.readCsv(dPath), 'Sunburst');
      out.Tree = await addAndCollect(await grok.dapi.files.readCsv(dPath), 'Tree');
      out.Timelines = await addAndCollect(await grok.dapi.files.readCsv(aPath), 'Timelines');
      return out;
    }, [demogPath, aePath]);
    expect(result.Radar).toEqual(expect.arrayContaining(['table', 'colorColumnName', 'backgroundMinColor', 'currentRowColor', 'legendVisibility']));
    expect(result.Sunburst).toEqual(expect.arrayContaining(['hierarchyColumnNames', 'inheritFromGrid', 'includeNulls']));
    expect(result.Tree).toEqual(expect.arrayContaining(['hierarchyColumnNames', 'orient', 'onClick', 'showCounts', 'includeNulls']));
    expect(result.Timelines).toEqual(expect.arrayContaining(['splitByColumnName', 'startColumnName', 'endColumnName', 'colorColumnName', 'legendVisibility', 'marker']));
  });

  await page.evaluate(() => (window as any).grok.shell.closeAll());

  expect(consoleErrors, consoleErrors.join('\n')).toHaveLength(0);

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
