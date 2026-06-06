import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';

test('Tree — 8-enhancement bundle regression (github-3221)', async ({page}) => {
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

  let availableProps: string[] = [];

  await softStep('Setup: Open demog.csv, add Tree, set hierarchy CONTROL/SEX/RACE', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise((r) => setTimeout(r, 1500));
      const tree = tv.addViewer('Tree');
      await new Promise((r) => setTimeout(r, 3000));
      try { tree.setOptions({hierarchyColumnNames: ['CONTROL', 'SEX', 'RACE']}); } catch (e) {}
      await new Promise((r) => setTimeout(r, 1500));
      const types: string[] = [];
      for (const v of tv.viewers) types.push(v.type);
      let propNames: string[] = [];
      try { propNames = tree.props.getProperties().map((p: any) => p.name); } catch (e) {}
      return {types, propNames};
    }, demogPath);
    expect(result.types).toContain('Tree');
    console.log('[Tree props]', result.propNames.join(', '));
    availableProps = result.propNames;
  });

  // Helper applied 8 times — exercises each enhancement
  const exerciseProperty = async (propName: string, value: any, label: string) => {
    await softStep(`Capability: ${label} (${propName})`, async () => {
      const errorsBefore = consoleErrors.length;
      const result = await page.evaluate(async ([p, v]) => {
        const grok = (window as any).grok;
        const tv = grok.shell.tv;
        let tree: any = null;
        for (const vw of tv.viewers) if (vw.type === 'Tree') { tree = vw; break; }
        if (!tree) return {ok: false};
        let setOptionsThrew = false;
        try {
          const opts: any = {}; opts[p] = v;
          tree.setOptions(opts);
        } catch (e) { setOptionsThrew = true; }
        await new Promise((r) => setTimeout(r, 800));
        let readBack: any = null;
        try { readBack = tree.props.get(p); } catch (e) {}
        const root = tree.root as HTMLElement;
        return {
          ok: true,
          setOptionsThrew,
          readBack,
          hasContent: root.children.length > 0,
          width: root.getBoundingClientRect().width,
        };
      }, [propName, value] as [string, any]);
      expect(result.ok).toBe(true);
      // github-3221 invariant: each capability's setOptions does not throw.
      expect(result.setOptionsThrew).toBe(false);
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

  // Capability 4 (layout) — guarded probe; 'layout' may not be exposed on dev build.
  if (availableProps.includes('layout'))
    await exerciseProperty('layout', 'orthogonal', 'layout=orthogonal');
  else
    console.warn('[SKIP] charts.tree.layout not exposed on current Tree build; capability 4 deferred per scenario step 7 conditional');

  await exerciseProperty('showMouseOverLine', true, 'showMouseOverLine ON');

  await softStep('Final visual stability check after 8 enhancements', async () => {
    const result = await page.evaluate(async () => {
      const grok = (window as any).grok;
      const tv = grok.shell.tv;
      let tree: any = null;
      for (const v of tv.viewers) if (v.type === 'Tree') { tree = v; break; }
      if (!tree) return {ok: false};
      const root = tree.root as HTMLElement;
      const rect = root.getBoundingClientRect();
      return {ok: true, hasContent: root.children.length > 0, width: rect.width, height: rect.height};
    });
    expect(result.ok).toBe(true);
    expect(result.hasContent).toBe(true);
    expect(result.width).toBeGreaterThan(0);
    expect(result.height).toBeGreaterThan(0);
  });

  await page.evaluate(() => (window as any).grok.shell.closeAll());

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
