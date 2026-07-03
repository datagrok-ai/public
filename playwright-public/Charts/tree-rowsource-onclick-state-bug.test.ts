/* ---
sub_features_covered: [charts.tree, charts.tree.color-palette, charts.tree.on-click]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';

test('Charts / Tree — rowSource x onClick state machine (github-3245)', async ({page}) => {
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

  const hierarchy = ['CONTROL', 'SEX', 'RACE'];

  await softStep('Setup: Open demog.csv, add Tree, set hierarchy', async () => {
    await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      tv.addViewer('Tree');
    }, demogPath);
    await page.waitForFunction(() => {
      const grok = (window as any).grok;
      const tv = grok.shell.tv;
      if (!tv || !tv.grid) return false;
      let tree: any = null;
      for (const v of tv.viewers) if (v.type === 'Tree') { tree = v; break; }
      return !!tree && !!tree.root && tree.root.children.length > 0;
    }, null, {timeout: 30_000});

    const result = await page.evaluate((h) => {
      const grok = (window as any).grok;
      const tv = grok.shell.tv;
      let tree: any = null;
      for (const v of tv.viewers) if (v.type === 'Tree') { tree = v; break; }
      tree.setOptions({hierarchyColumnNames: h});
      const types: string[] = [];
      for (const v of tv.viewers) types.push(v.type);
      return {types};
    }, hierarchy);
    expect(result.types).toContain('Tree');

    await page.waitForFunction((h) => {
      const grok = (window as any).grok;
      const tv = grok.shell.tv;
      let tree: any = null;
      for (const v of tv.viewers) if (v.type === 'Tree') { tree = v; break; }
      if (!tree) return false;
      let hb: any = null;
      try { hb = tree.props.get('hierarchyColumnNames'); } catch (e) { return false; }
      return Array.isArray(hb) && hb.join(',') === h.join(',') && tree.root.children.length > 0;
    }, hierarchy, {timeout: 20_000});

    const readback = await page.evaluate(() => {
      const grok = (window as any).grok;
      const tv = grok.shell.tv;
      let tree: any = null;
      for (const v of tv.viewers) if (v.type === 'Tree') { tree = v; break; }
      return tree.props.get('hierarchyColumnNames');
    });
    expect(readback).toEqual(hierarchy);
  });

  // Palette: colorColumnName drives the Tree's color mapping (charts.tree.color-palette).
  await softStep('Palette: set colorColumnName=AGE, verify readback + render', async () => {
    const result = await page.evaluate(async () => {
      const grok = (window as any).grok;
      const tv = grok.shell.tv;
      let tree: any = null;
      for (const v of tv.viewers) if (v.type === 'Tree') { tree = v; break; }
      if (!tree) return {ok: false};
      let threw = false;
      try { tree.setOptions({colorColumnName: 'AGE'}); }
      catch (e) { threw = true; }
      let readback: any = null;
      for (let i = 0; i < 30; i++) {
        readback = tree.props.get('colorColumnName');
        if (readback === 'AGE' && tree.root.children.length > 0) break;
        await new Promise((r) => setTimeout(r, 100));
      }
      const root = tree.root as HTMLElement;
      return {ok: true, threw, readback, hasContent: root.children.length > 0, width: root.getBoundingClientRect().width};
    });
    expect(result.ok).toBe(true);
    expect(result.threw).toBe(false);
    expect(result.readback).toBe('AGE');
    expect(result.hasContent).toBe(true);
    expect(result.width).toBeGreaterThan(0);
  });

  const rowSources = ['All', 'Filtered', 'Selected'];
  const onClicks = ['Select', 'Filter', 'None'];
  // github-3245: onClick coerces rowSource — setting either forces rowSource = rowSourceMap[onClick].
  const rowSourceMap: Record<string, string> = {Select: 'Filtered', Filter: 'All', None: 'All'};

  const iterateCombo = async (rs: string, oc: string, suffix: string) => {
    await softStep(`Combo ${rs} x ${oc} ${suffix}`, async () => {
      const expectedRs = rowSourceMap[oc];
      const errorsBefore = consoleErrors.length;
      const result = await page.evaluate(async ([r, o, expRs]) => {
        const grok = (window as any).grok;
        const tv = grok.shell.tv;
        let tree: any = null;
        for (const v of tv.viewers) if (v.type === 'Tree') { tree = v; break; }
        if (!tree) return {ok: false};
        let setOptionsThrew = false;
        try { tree.setOptions({rowSource: r, onClick: o}); }
        catch (e) { setOptionsThrew = true; }
        let readbackRs: any = null, readbackOc: any = null;
        for (let i = 0; i < 50; i++) {
          readbackRs = tree.props.get('rowSource');
          readbackOc = tree.props.get('onClick');
          if (readbackOc === o && readbackRs === expRs && tree.root.children.length > 0) break;
          await new Promise((res) => setTimeout(res, 100));
        }
        const root = tree.root as HTMLElement;
        const rect = root.getBoundingClientRect();
        return {ok: true, setOptionsThrew, readbackRs, readbackOc,
          hasContent: root.children.length > 0, width: rect.width, height: rect.height};
      }, [rs, oc, expectedRs] as [string, string, string]);
      expect(result.ok, `${rs}/${oc} tree present`).toBe(true);
      expect(result.setOptionsThrew, `${rs}/${oc} setOptions`).toBe(false);
      expect(result.readbackOc, `${rs}/${oc} onClick readback`).toBe(oc);
      expect(result.readbackRs, `${rs}/${oc} rowSource normalized`).toBe(expectedRs);
      expect(result.hasContent, `${rs}/${oc} content`).toBe(true);
      expect(result.width, `${rs}/${oc} width`).toBeGreaterThan(0);
      const errorsDuring = consoleErrors.slice(errorsBefore);
      expect(errorsDuring).toEqual([]);
    });
  };

  // Forward iteration: 9 combinations
  for (const rs of rowSources) {
    for (const oc of onClicks) {
      await iterateCombo(rs, oc, '(forward)');
    }
  }

  // Reverse iteration: same 9 combos, reverse order
  for (const rs of [...rowSources].reverse()) {
    for (const oc of [...onClicks].reverse()) {
      await iterateCombo(rs, oc, '(reverse)');
    }
  }

  // ESC-equivalent state-clear check
  await softStep('Step 4: ESC-equivalent (df.selection.setAll(false)) — selection clears, viewer stable', async () => {
    const result = await page.evaluate(async () => {
      const grok = (window as any).grok;
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      df.selection.fireChanged();
      const tv = grok.shell.tv;
      let tree: any = null;
      for (const v of tv.viewers) if (v.type === 'Tree') { tree = v; break; }
      if (!tree) return {ok: false};
      for (let i = 0; i < 30; i++) {
        if (df.selection.trueCount === 0 && tree.root.children.length > 0) break;
        await new Promise((r) => setTimeout(r, 100));
      }
      const root = tree.root as HTMLElement;
      const rect = root.getBoundingClientRect();
      return {ok: true, trueCount: df.selection.trueCount, hasContent: root.children.length > 0, width: rect.width};
    });
    expect(result.ok).toBe(true);
    expect(result.trueCount).toBe(0);
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
