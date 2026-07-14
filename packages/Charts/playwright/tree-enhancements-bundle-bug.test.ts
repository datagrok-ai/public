/* ---
sub_features_covered: [charts.tree, charts.tree.font-size, charts.tree.include-nulls, charts.tree.layout, charts.tree.on-click, charts.tree.orient, charts.tree.show-counts]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';

test('Charts / Tree — 8-enhancement bundle regression (github-3221)', async ({page}) => {
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

  let availableProps: string[] = [];

  await softStep('Setup: Open demog.csv, add Tree, set hierarchy CONTROL/SEX/RACE', async () => {
    await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      const df = await grok.dapi.files.readCsv(path);
      grok.shell.addTableView(df);
    }, demogPath);
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 15_000});

    await page.evaluate(() => (window as any).grok.shell.tv.addViewer('Tree'));
    await page.waitForFunction(() => {
      const t = (window as any).grok.shell.tv.viewers.find((v: any) => v.type === 'Tree');
      return t && (t.root as HTMLElement).children.length > 0;
    }, null, {timeout: 15_000});

    const setThrew = await page.evaluate(() => {
      const tree = (window as any).grok.shell.tv.viewers.find((v: any) => v.type === 'Tree');
      try { tree.setOptions({hierarchyColumnNames: ['CONTROL', 'SEX', 'RACE']}); return false; }
      catch (e) { return true; }
    });
    expect(setThrew).toBe(false);
    // Let the hierarchy change re-render the tree. hierarchyColumnNames is a COLUMN_LIST
    // property that is not readable via props.get on the dev Tree build (throws
    // "Property not found"), so gate on the viewer staying rendered rather than
    // round-tripping the value.
    await page.waitForFunction(() => {
      const t = (window as any).grok.shell.tv.viewers.find((v: any) => v.type === 'Tree');
      if (!t || (t.root as HTMLElement).children.length === 0) return false;
      try { return t.props.getProperties().length > 0; } catch { return false; }
    }, null, {timeout: 30_000});

    const meta = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const tree = tv.viewers.find((v: any) => v.type === 'Tree');
      return {
        types: tv.viewers.map((v: any) => v.type) as string[],
        propNames: tree.props.getProperties().map((p: any) => p.name) as string[],
      };
    });
    expect(meta.types).toContain('Tree');
    expect(meta.propNames.length).toBeGreaterThan(0);
    availableProps = meta.propNames;
  });

  // Helper applied 8 times — exercises each enhancement
  const exerciseProperty = async (propName: string, value: any, label: string) => {
    await softStep(`Capability: ${label} (${propName})`, async () => {
      const errorsBefore = consoleErrors.length;
      const setThrew = await page.evaluate(([p, v]) => {
        const tree = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Tree');
        if (!tree) return null;
        try { const opts: any = {}; opts[p] = v; tree.setOptions(opts); return false; }
        catch (e) { return true; }
      }, [propName, value] as [string, any]);
      expect(setThrew).toBe(false);

      // github-3221 invariant: setOptions actually takes effect — poll until the property round-trips.
      await page.waitForFunction(([p, v]) => {
        const tree = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Tree');
        try { return JSON.stringify(tree.props.get(p)) === JSON.stringify(v); }
        catch { return false; }
      }, [propName, value] as [string, any], {timeout: 10_000});

      const result = await page.evaluate(([p]) => {
        const tree = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Tree');
        const root = tree.root as HTMLElement;
        return {
          readBack: tree.props.get(p),
          hasContent: root.children.length > 0,
          width: root.getBoundingClientRect().width,
        };
      }, [propName] as [string]);
      expect(result.readBack).toStrictEqual(value);
      expect(result.hasContent).toBe(true);
      expect(result.width).toBeGreaterThan(0);
      const errorsDuring = consoleErrors.slice(errorsBefore);
      expect(errorsDuring).toEqual([]);
    });
  };

  await exerciseProperty('showCounts', true, 'showCounts ON');
  await exerciseProperty('showCounts', false, 'showCounts OFF');
  await exerciseProperty('onClick', 'Filter', 'onClick=Filter');
  await exerciseProperty('onClick', 'Select', 'onClick=Select');
  await exerciseProperty('fontSize', 14, 'fontSize=14');
  await exerciseProperty('fontSize', 30, 'fontSize=30 (atlas-cited max)');
  await exerciseProperty('orient', 'LR', 'orient=LR');
  await exerciseProperty('orient', 'RL', 'orient=RL');
  await exerciseProperty('orient', 'TB', 'orient=TB (default)');
  await exerciseProperty('includeNulls', true, 'includeNulls ON');
  await exerciseProperty('includeNulls', false, 'includeNulls OFF');

  // Capability 4 (layout) — 'layout' may not be exposed on this Tree build (see riskNotes).
  if (availableProps.includes('layout'))
    await exerciseProperty('layout', 'orthogonal', 'layout=orthogonal');
  else
    test.info().annotations.push({type: 'skip', description: 'charts.tree.layout not exposed on current Tree build'});

  await exerciseProperty('showMouseOverLine', true, 'showMouseOverLine ON');

  await softStep('Final visual stability check after 8 enhancements', async () => {
    const result = await page.evaluate(() => {
      const tree = (window as any).grok.shell.tv.viewers.find((v: any) => v.type === 'Tree');
      if (!tree) return {ok: false};
      const root = tree.root as HTMLElement;
      const rect = root.getBoundingClientRect();
      return {ok: true, hasContent: root.children.length > 0, width: rect.width, height: rect.height};
    });
    expect(result.ok).toBe(true);
    expect(result.hasContent).toBe(true);
    expect(result.width).toBeGreaterThan(0);
    expect(result.height).toBeGreaterThan(0);
    expect(consoleErrors).toEqual([]);
  });

  await page.evaluate(() => (window as any).grok.shell.closeAll());

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
