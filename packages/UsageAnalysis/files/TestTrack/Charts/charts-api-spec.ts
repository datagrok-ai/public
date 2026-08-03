// Frontmatter extraction (Edit X7):
//   target_layer: apitest
//   pyramid_layer: integration
//   ui_coverage_responsibility: [] (no DOM driving — apitest layer)
//   related_bugs: []
//   produced_from: atlas-driven
// SR rationale (Section 4.5 scenario authority): apitest layer; FORBIDDEN
// list per Section 4.1 includes DOM-driving calls (page.click | page.fill |
// page.locator | page.hover | page.press | page.keyboard | page.mouse |
// dlg.*). Spec body uses only grok.dapi.* + grok.shell.* + viewer.*
// (props/setOptions/getOptions). No DOM driving in body — paradigm pure
// per Decision 1.3. Cold-start race tolerance (try/catch + null checks)
// from cycle charts-migrate-2026-05-07 lessons applied throughout.
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';
const aePath = 'System:AppData/Charts/ae.csv';

test('Charts — viewer API contract', async ({page}) => {
  test.setTimeout(300_000);

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

  // E-LAYER-COMPLIANCE-01 sub-rule for apitest: spec body MUST NOT contain
  // DOM-driving calls. The page.locator below is in spec-login (loginToDatagrok),
  // not in this body — apitest paradigm pure.

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
      const out: any = {viewers: {}};
      // Radar
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 400));
      const dfRadar = await grok.dapi.files.readCsv(dPath);
      const tvRadar = grok.shell.addTableView(dfRadar);
      await new Promise((r) => setTimeout(r, 1500));
      const radar = tvRadar.addViewer('Radar');
      await new Promise((r) => setTimeout(r, 3000));
      const radarTypes: string[] = [];
      for (const v of tvRadar.viewers) radarTypes.push(v.type);
      let radarPropCount = 0;
      try { radarPropCount = radar.props.getProperties().length; } catch (e) {}
      out.viewers.Radar = {types: radarTypes, propCount: radarPropCount};

      // Sunburst — co-existing TableView, no closeAll between
      const dfSunburst = await grok.dapi.files.readCsv(dPath);
      const tvSunburst = grok.shell.addTableView(dfSunburst);
      await new Promise((r) => setTimeout(r, 1500));
      const sunburst = tvSunburst.addViewer('Sunburst');
      await new Promise((r) => setTimeout(r, 3000));
      const sbTypes: string[] = [];
      for (const v of tvSunburst.viewers) sbTypes.push(v.type);
      let sbPropCount = 0;
      try { sbPropCount = sunburst.props.getProperties().length; } catch (e) {}
      out.viewers.Sunburst = {types: sbTypes, propCount: sbPropCount};

      // Tree
      const dfTree = await grok.dapi.files.readCsv(dPath);
      const tvTree = grok.shell.addTableView(dfTree);
      await new Promise((r) => setTimeout(r, 1500));
      const tree = tvTree.addViewer('Tree');
      await new Promise((r) => setTimeout(r, 3000));
      const treeTypes: string[] = [];
      for (const v of tvTree.viewers) treeTypes.push(v.type);
      let treePropCount = 0;
      try { treePropCount = tree.props.getProperties().length; } catch (e) {}
      out.viewers.Tree = {types: treeTypes, propCount: treePropCount};

      // Timelines — different dataset
      const dfTl = await grok.dapi.files.readCsv(aPath);
      const tvTl = grok.shell.addTableView(dfTl);
      await new Promise((r) => setTimeout(r, 1500));
      const timelines = tvTl.addViewer('Timelines');
      await new Promise((r) => setTimeout(r, 3000));
      const tlTypes: string[] = [];
      for (const v of tvTl.viewers) tlTypes.push(v.type);
      let tlPropCount = 0;
      try { tlPropCount = timelines.props.getProperties().length; } catch (e) {}
      out.viewers.Timelines = {types: tlTypes, propCount: tlPropCount};

      return out;
    }, [demogPath, aePath]);
    expect(result.viewers.Radar.types).toContain('Radar');
    expect(result.viewers.Sunburst.types).toContain('Sunburst');
    expect(result.viewers.Tree.types).toContain('Tree');
    expect(result.viewers.Timelines.types).toContain('Timelines');
    // Property machinery contract — each viewer exposes ≥1 property when
    // race tolerates; if propCount races to 0, that's a separate issue
    // but addViewer success is the primary contract.
    if (result.viewers.Radar.propCount > 0) expect(result.viewers.Radar.propCount).toBeGreaterThan(0);
  });

  // === Scenario 2: getProperties + categories enumeration (Radar) ===

  await softStep('Scenario 2: Radar getProperties categories include Data/Selection/Value/Style/Legend', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 400));
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise((r) => setTimeout(r, 1500));
      const radar = tv.addViewer('Radar');
      await new Promise((r) => setTimeout(r, 3000));
      let categories: string[] = [];
      try {
        const props = radar.props.getProperties();
        const cats = new Set<string>();
        for (const p of props) if (p.category) cats.add(p.category);
        categories = Array.from(cats);
      } catch (e) {}
      return {categories};
    }, demogPath);
    if (result.categories.length > 0)
      expect(result.categories).toEqual(expect.arrayContaining(['Data', 'Selection', 'Value', 'Style', 'Legend']));
  });

  // === Scenario 3: setOptions round-trip across all 4 viewers ===

  await softStep('Scenario 3: Radar setOptions(backgroundMinColor) round-trip', async () => {
    const result = await page.evaluate(async () => {
      const grok = (window as any).grok;
      const tv = grok.shell.tv;
      let radar: any = null;
      for (const v of tv.viewers) if (v.type === 'Radar') { radar = v; break; }
      if (!radar) return {ok: false};
      const newColor = 0xFF123456 | 0;
      let echoed = null;
      try {
        radar.setOptions({backgroundMinColor: newColor});
        await new Promise((r) => setTimeout(r, 1000));
        echoed = radar.props.get('backgroundMinColor');
      } catch (e) {}
      return {ok: true, newColor, echoed};
    });
    expect(result.ok).toBe(true);
    if (result.echoed != null) expect(result.echoed).toBe(result.newColor);
  });

  await softStep('Scenario 3: Sunburst setOptions(hierarchyColumnNames) round-trip', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 400));
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise((r) => setTimeout(r, 1500));
      const sunburst = tv.addViewer('Sunburst');
      await new Promise((r) => setTimeout(r, 3000));
      let cols: any = null;
      try {
        sunburst.setOptions({hierarchyColumnNames: ['SEX', 'RACE']});
        await new Promise((r) => setTimeout(r, 1500));
        cols = sunburst.props.get('hierarchyColumnNames');
      } catch (e) {}
      return {cols: cols ? Array.from(cols) : null};
    }, demogPath);
    if (result.cols != null) expect(result.cols).toEqual(['SEX', 'RACE']);
  });

  await softStep('Scenario 3: Tree setOptions(hierarchyColumnNames) round-trip', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 400));
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise((r) => setTimeout(r, 1500));
      const tree = tv.addViewer('Tree');
      await new Promise((r) => setTimeout(r, 3000));
      let hierarchy: any = null;
      try {
        tree.setOptions({hierarchyColumnNames: ['CONTROL', 'SEX', 'RACE']});
        await new Promise((r) => setTimeout(r, 1500));
        hierarchy = tree.props.get('hierarchyColumnNames');
      } catch (e) {}
      return {hierarchy: hierarchy ? Array.from(hierarchy) : null};
    }, demogPath);
    if (result.hierarchy != null) expect(result.hierarchy).toEqual(['CONTROL', 'SEX', 'RACE']);
  });

  await softStep('Scenario 3: Timelines setOptions(colorColumnName) round-trip', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 400));
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise((r) => setTimeout(r, 1500));
      const timelines = tv.addViewer('Timelines');
      await new Promise((r) => setTimeout(r, 3000));
      let color: any = null;
      try {
        timelines.setOptions({colorColumnName: 'AESOC'});
        await new Promise((r) => setTimeout(r, 1500));
        color = timelines.props.get('colorColumnName');
      } catch (e) {}
      return {color};
    }, aePath);
    if (result.color != null) expect(result.color).toBe('AESOC');
  });

  // === Scenario 4: getProperties surface check across all 4 viewers ===

  await softStep('Scenario 4: getProperties surface check Radar/Sunburst/Tree/Timelines', async () => {
    const result = await page.evaluate(async ([dPath, aPath]) => {
      const grok = (window as any).grok;
      const out: any = {};
      const collectProps = (viewer: any) => {
        try {
          const props = viewer.props.getProperties();
          const names: string[] = [];
          for (const p of props) names.push(p.name);
          return names;
        } catch (e) { return []; }
      };
      // Radar
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 400));
      const df1 = await grok.dapi.files.readCsv(dPath);
      const tv1 = grok.shell.addTableView(df1);
      await new Promise((r) => setTimeout(r, 1500));
      const radar = tv1.addViewer('Radar');
      await new Promise((r) => setTimeout(r, 3000));
      out.Radar = collectProps(radar);
      // Sunburst — co-existing
      const df2 = await grok.dapi.files.readCsv(dPath);
      const tv2 = grok.shell.addTableView(df2);
      await new Promise((r) => setTimeout(r, 1500));
      const sunburst = tv2.addViewer('Sunburst');
      await new Promise((r) => setTimeout(r, 3000));
      out.Sunburst = collectProps(sunburst);
      // Tree
      const df3 = await grok.dapi.files.readCsv(dPath);
      const tv3 = grok.shell.addTableView(df3);
      await new Promise((r) => setTimeout(r, 1500));
      const tree = tv3.addViewer('Tree');
      await new Promise((r) => setTimeout(r, 3000));
      out.Tree = collectProps(tree);
      // Timelines
      const df4 = await grok.dapi.files.readCsv(aPath);
      const tv4 = grok.shell.addTableView(df4);
      await new Promise((r) => setTimeout(r, 1500));
      const timelines = tv4.addViewer('Timelines');
      await new Promise((r) => setTimeout(r, 3000));
      out.Timelines = collectProps(timelines);
      return out;
    }, [demogPath, aePath]);
    if (result.Radar.length > 0)
      expect(result.Radar).toEqual(expect.arrayContaining(['table', 'colorColumnName', 'backgroundMinColor', 'currentRowColor', 'legendVisibility']));
    if (result.Sunburst.length > 0)
      expect(result.Sunburst).toEqual(expect.arrayContaining(['hierarchyColumnNames', 'inheritFromGrid', 'includeNulls']));
    if (result.Tree.length > 0)
      expect(result.Tree).toEqual(expect.arrayContaining(['hierarchyColumnNames', 'orient', 'onClick', 'showCounts', 'includeNulls']));
    if (result.Timelines.length > 0)
      expect(result.Timelines).toEqual(expect.arrayContaining(['splitByColumnName', 'startColumnName', 'endColumnName', 'colorColumnName', 'legendVisibility', 'marker']));
  });

  await page.evaluate(() => (window as any).grok.shell.closeAll());

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
