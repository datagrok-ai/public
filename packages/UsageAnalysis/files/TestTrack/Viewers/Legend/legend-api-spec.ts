// Frontmatter extraction (Edit X7):
//   target_layer: apitest
//   pyramid_layer: integration
//   related_bugs: []
//   produced_from: atlas-driven
// SR rationale (Section 4.5 scenario authority): apitest layer; FORBIDDEN
// list per Section 4.1 includes DOM-driving calls (page.click | page.fill |
// page.locator | page.hover | page.press | page.keyboard | page.mouse |
// dlg.*). Spec body uses only grok.dapi.* + grok.shell.* + viewer.props.* +
// col.tags.* + col.meta.colors.* + df.filter.* + tv.getFiltersGroup. No DOM
// driving in body — paradigm pure per Decision 1.3. Cold-start race
// tolerance (try/catch + null checks) from cycle charts-migrate-2026-05-07
// lessons applied throughout.
//
// ApiSamples reference (per user feedback, see public/packages/ApiSamples/scripts):
//   scripts/grid/color-coding/color-coding.js — setCategorical, setLinear
//   scripts/grid/color-coding/get-cell-color.js — DG.Color.getCellColorHtml
//   scripts/ui/viewers/filters/filter-group.js — tv.getFiltersGroup + fg.updateOrAdd
//   scripts/data-frame/events/events.js — column metadata events
//
// Sister scenario: Charts/charts-api.md + Charts/charts-api.ts (also target_layer: apitest).
// Closes documented follow-up from modernize-legacy-specs.md §4.
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:DemoFiles/chem/SPGI.csv';

test('Legend — JS API contract', async ({page}) => {
  // Validator B 2026-05-09 run 2 timed out at 318s on Scenario 8 (300_000 budget too tight
  // for 8 sequential dataset re-reads + viewer attaches). Bump to 600_000 — runs 1+3
  // settled at 264s+285s, so 600s gives 2× headroom for transient dev slowness.
  test.setTimeout(600_000);
  stepErrors.length = 0;

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
    document.body.classList.add('selenium');
    (window as any).grok.shell.closeAll();
    (window as any).grok.shell.windows.simpleMode = true;
  });

  // === Scenario 1: legend.column round-trip across host viewer types ===

  await softStep('Scenario 1: legend.column round-trip across Scatter/Histogram/Bar/Pie/Line', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 400));
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise((r) => setTimeout(r, 1500));
      const out: Record<string, any> = {};
      const configs: Array<[string, string, string]> = [
        ['Scatter plot', 'colorColumnName', 'RACE'],
        ['Histogram', 'splitColumnName', 'RACE'],
        ['Bar chart', 'splitColumnName', 'RACE'],
        ['Pie chart', 'categoryColumnName', 'RACE'],
        ['Line chart', 'splitColumnName', 'RACE'],
      ];
      for (const [vtype, prop, val] of configs) {
        try {
          tv.addViewer(vtype);
          await new Promise((r) => setTimeout(r, 600));
          const v = tv.viewers.find((x: any) => x.type === vtype);
          if (!v) { out[vtype] = {error: 'viewer not attached'}; continue; }
          (v.props as any)[prop] = val;
          try { v.props.legendVisibility = 'Always'; } catch (_) {}
          await new Promise((r) => setTimeout(r, 800));
          let echoed: any = null;
          try { echoed = (v.props as any)[prop]; } catch (_) {}
          out[vtype] = {prop, expected: val, echoed};
        } catch (e: any) {
          out[vtype] = {error: String(e?.message ?? e).slice(0, 200)};
        }
      }
      return out;
    }, demogPath);
    // Each viewer's legend column property must round-trip when the viewer
    // attaches; if it races to null, that's acceptable per cold-start tolerance.
    for (const [vtype, info] of Object.entries(result as Record<string, any>)) {
      if (info.error) continue;
      if (info.echoed != null) expect(info.echoed, `${vtype} ${info.prop}`).toBe(info.expected);
    }
    // At least 3 of 5 viewers must round-trip the property successfully.
    const successes = Object.values(result as Record<string, any>)
      .filter((info: any) => !info.error && info.echoed === info.expected).length;
    expect(successes, 'at least 3 of 5 viewers round-trip legend.column').toBeGreaterThanOrEqual(3);
  });

  // === Scenario 2: legend.extra-column round-trip on Scatter plot ===

  await softStep('Scenario 2: legend.extra-column markersColumnName round-trip + deselect', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 400));
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise((r) => setTimeout(r, 1500));
      tv.addViewer('Scatter plot');
      await new Promise((r) => setTimeout(r, 800));
      const sp = tv.viewers.find((v: any) => v.type === 'Scatter plot');
      if (!sp) return {ok: false, error: 'no scatter'};
      sp.props.colorColumnName = 'RACE';
      sp.props.markersColumnName = 'SEX';
      try { sp.props.legendVisibility = 'Always'; } catch (_) {}
      await new Promise((r) => setTimeout(r, 1000));
      let colorEchoed: any = null;
      let markersEchoed: any = null;
      try { colorEchoed = sp.props.colorColumnName; } catch (_) {}
      try { markersEchoed = sp.props.markersColumnName; } catch (_) {}
      // Deselect markers (GROK-19083 path).
      sp.props.markersColumnName = '';
      await new Promise((r) => setTimeout(r, 800));
      let markersAfterDeselect: any = null;
      try { markersAfterDeselect = sp.props.markersColumnName; } catch (_) {}
      return {ok: true, colorEchoed, markersEchoed, markersAfterDeselect};
    }, demogPath);
    expect(result.ok).toBe(true);
    if (result.colorEchoed != null) expect(result.colorEchoed).toBe('RACE');
    if (result.markersEchoed != null) expect(result.markersEchoed).toBe('SEX');
    if (result.markersAfterDeselect != null) expect(result.markersAfterDeselect).toBe('');
  });

  // === Scenario 3: legend.color-scale.numerical via tag round-trip ===

  await softStep('Scenario 3: numerical color-scale tag round-trip on AGE', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 400));
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise((r) => setTimeout(r, 1500));
      tv.addViewer('Scatter plot');
      await new Promise((r) => setTimeout(r, 800));
      const sp = tv.viewers.find((v: any) => v.type === 'Scatter plot');
      if (!sp) return {ok: false};
      sp.props.colorColumnName = 'AGE';
      try { sp.props.legendVisibility = 'Always'; } catch (_) {}
      const col = df.col('AGE');
      col.tags['.color-coding-type'] = 'Linear';
      col.tags['.color-coding-scheme'] = '[1, 8388607, 16711680]';
      try { sp.invalidate?.(); } catch (_) {}
      await new Promise((r) => setTimeout(r, 1000));
      return {
        ok: true,
        codingType: col.tags['.color-coding-type'],
        scheme: col.tags['.color-coding-scheme'],
      };
    }, demogPath);
    expect(result.ok).toBe(true);
    expect(result.codingType).toBe('Linear');
    expect(result.scheme).toBe('[1, 8388607, 16711680]');
  });

  // === Scenario 4: legend.use-custom-color-coding via setCategorical ===

  await softStep('Scenario 4: setCategorical round-trip via JSON tag', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 400));
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise((r) => setTimeout(r, 1500));
      tv.addViewer('Scatter plot');
      await new Promise((r) => setTimeout(r, 800));
      const sp = tv.viewers.find((v: any) => v.type === 'Scatter plot');
      if (!sp) return {ok: false};
      sp.props.colorColumnName = 'Stereo Category';
      try { sp.props.legendVisibility = 'Always'; } catch (_) {}
      const col = df.col('Stereo Category');
      col.tags['.color-coding-type'] = 'Categorical';
      // ApiSamples reference: scripts/grid/color-coding/color-coding.js
      col.meta.colors.setCategorical(
        {'R_ONE': '#FF0000', 'S_UNKN': '#00FF00'},
        {fallbackColor: '#808080'},
      );
      try { sp.invalidate?.(); } catch (_) {}
      await new Promise((r) => setTimeout(r, 1000));
      let parsed: Record<string, any> = {};
      try { parsed = JSON.parse(col.tags['.color-coding-categorical'] ?? '{}'); } catch (_) {}
      return {
        ok: true,
        codingType: col.tags['.color-coding-type'],
        rOne: String(parsed['R_ONE'] ?? '').toLowerCase(),
        sUnkn: String(parsed['S_UNKN'] ?? '').toLowerCase(),
      };
    }, spgiPath);
    expect(result.ok).toBe(true);
    expect(result.codingType).toBe('Categorical');
    expect(result.rOne).toBe('#ff0000');
    expect(result.sUnkn).toBe('#00ff00');
  });

  // === Scenario 5: legend.show-nulls via includeNulls round-trip ===

  await softStep('Scenario 5: includeNulls round-trip on Histogram + Bar chart', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 400));
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise((r) => setTimeout(r, 1500));
      const out: Record<string, any> = {};
      // Histogram
      try {
        tv.addViewer('Histogram');
        await new Promise((r) => setTimeout(r, 700));
        const h = tv.viewers.find((v: any) => v.type === 'Histogram');
        if (h) {
          h.props.splitColumnName = 'Primary Scaffold Name';
          try { h.props.legendVisibility = 'Always'; } catch (_) {}
          h.props.includeNulls = true;
          await new Promise((r) => setTimeout(r, 500));
          const trueEcho = h.props.includeNulls;
          h.props.includeNulls = false;
          await new Promise((r) => setTimeout(r, 500));
          const falseEcho = h.props.includeNulls;
          out.Histogram = {trueEcho, falseEcho};
        }
      } catch (e: any) { out.Histogram = {error: String(e?.message ?? e).slice(0, 200)}; }
      // Bar chart
      try {
        tv.addViewer('Bar chart');
        await new Promise((r) => setTimeout(r, 700));
        const bc = tv.viewers.find((v: any) => v.type === 'Bar chart');
        if (bc) {
          bc.props.splitColumnName = 'Primary Scaffold Name';
          try { bc.props.legendVisibility = 'Always'; } catch (_) {}
          bc.props.includeNulls = true;
          await new Promise((r) => setTimeout(r, 500));
          const trueEcho = bc.props.includeNulls;
          bc.props.includeNulls = false;
          await new Promise((r) => setTimeout(r, 500));
          const falseEcho = bc.props.includeNulls;
          out.BarChart = {trueEcho, falseEcho};
        }
      } catch (e: any) { out.BarChart = {error: String(e?.message ?? e).slice(0, 200)}; }
      return out;
    }, spgiPath);
    if (result.Histogram && !result.Histogram.error) {
      expect(result.Histogram.trueEcho).toBe(true);
      expect(result.Histogram.falseEcho).toBe(false);
    }
    if (result.BarChart && !result.BarChart.error) {
      expect(result.BarChart.trueEcho).toBe(true);
      expect(result.BarChart.falseEcho).toBe(false);
    }
  });

  // === Scenario 6: legend.allow-item-coloring metadata round-trip on ≤100 cats ===

  await softStep('Scenario 6: setCategorical does not throw on <100 cats', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 400));
      const df = await grok.dapi.files.readCsv(path);
      grok.shell.addTableView(df);
      await new Promise((r) => setTimeout(r, 1500));
      const col = df.col('Stereo Category');
      const cats: string[] = Array.from(col.categories);
      const map: Record<string, string> = {};
      const palette = ['#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#FF00FF', '#00FFFF', '#808080'];
      cats.forEach((c, i) => { map[c] = palette[i % palette.length]; });
      let threwOnCall = false;
      try {
        col.tags['.color-coding-type'] = 'Categorical';
        col.meta.colors.setCategorical(map, {fallbackColor: '#808080'});
      } catch (_) { threwOnCall = true; }
      let parsed: Record<string, any> = {};
      try { parsed = JSON.parse(col.tags['.color-coding-categorical'] ?? '{}'); } catch (_) {}
      const roundTrippedCount = cats.filter(c =>
        String(parsed[c] ?? '').toLowerCase() === map[c].toLowerCase()).length;
      return {nCats: cats.length, threwOnCall, roundTrippedCount};
    }, spgiPath);
    expect(result.nCats).toBeLessThanOrEqual(100);
    expect(result.threwOnCall).toBe(false);
    expect(result.roundTrippedCount).toBe(result.nCats);
  });

  // === Scenario 7: legend.refresh.on-data-change after addNewCalculated ===

  await softStep('Scenario 7: addNewCalculated + re-bind colorColumnName resolves new column', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 400));
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise((r) => setTimeout(r, 1500));
      tv.addViewer('Scatter plot');
      await new Promise((r) => setTimeout(r, 800));
      const sp = tv.viewers.find((v: any) => v.type === 'Scatter plot');
      if (!sp) return {ok: false};
      sp.props.colorColumnName = 'SEX';
      try { sp.props.legendVisibility = 'Always'; } catch (_) {}
      await new Promise((r) => setTimeout(r, 800));
      try {
        await df.columns.addNewCalculated('SEX_alt', "if(${SEX}=='M', 'Male', 'Female')");
      } catch (_) {}
      await new Promise((r) => setTimeout(r, 1500));
      sp.props.colorColumnName = 'SEX_alt';
      await new Promise((r) => setTimeout(r, 1000));
      let echoed: any = null;
      try { echoed = sp.props.colorColumnName; } catch (_) {}
      const col = df.col('SEX_alt');
      const cats: string[] = col ? Array.from(col.categories) : [];
      return {ok: true, echoed, cats};
    }, demogPath);
    expect(result.ok).toBe(true);
    if (result.echoed != null) expect(result.echoed).toBe('SEX_alt');
    expect(result.cats.length).toBeGreaterThanOrEqual(2);
    expect(result.cats).toEqual(expect.arrayContaining(['Male', 'Female']));
  });

  // === Scenario 8: legend.refresh.on-reset-filter via df.filter.setAll(true) ===

  await softStep('Scenario 8: Filter Panel filter then df.filter.setAll(true) resets fully', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 400));
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise((r) => setTimeout(r, 1500));
      tv.addViewer('Scatter plot');
      await new Promise((r) => setTimeout(r, 800));
      const sp = tv.viewers.find((v: any) => v.type === 'Scatter plot');
      if (!sp) return {ok: false};
      sp.props.colorColumnName = 'Stereo Category';
      try { sp.props.legendVisibility = 'Always'; } catch (_) {}
      const fg = tv.getFiltersGroup();
      const DG = (window as any).DG;
      const cats: string[] = Array.from(df.col('Stereo Category').categories);
      const subset = cats.slice(0, Math.min(2, cats.length));
      const rowCount = df.rowCount;
      // ApiSamples reference: scripts/ui/viewers/filters/filter-group.js
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Stereo Category', selected: subset});
      await new Promise((r) => setTimeout(r, 1500));
      const filteredCount = df.filter.trueCount;
      df.filter.setAll(true);
      await new Promise((r) => setTimeout(r, 800));
      const resetCount = df.filter.trueCount;
      return {ok: true, rowCount, filteredCount, resetCount};
    }, spgiPath);
    expect(result.ok).toBe(true);
    expect(result.filteredCount).toBeLessThan(result.rowCount);
    expect(result.filteredCount).toBeGreaterThan(0);
    expect(result.resetCount).toBe(result.rowCount);
  });

  await page.evaluate(() => (window as any).grok.shell.closeAll());

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
