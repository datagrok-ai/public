import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const aePath = 'System:AppData/Charts/ae.csv';

test('Sunburst — date-column hierarchy graceful handling (github-2954)', async ({page}) => {
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

  await softStep('Steps 1-2: Open ae.csv, add Sunburst', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise((r) => setTimeout(r, 1500));
      tv.addViewer('Sunburst');
      await new Promise((r) => setTimeout(r, 3000));
      const types: string[] = [];
      for (const v of tv.viewers) types.push(v.type);
      const dateColumns: string[] = [];
      for (const c of df.columns) {
        const t = String(c.type || '').toLowerCase();
        if (t.includes('date') || t.includes('time')) dateColumns.push(c.name);
      }
      return {types, dateColumns};
    }, aePath);
    expect(result.types).toContain('Sunburst');
    expect(result.dateColumns.length).toBeGreaterThan(0);
  });

  await softStep('Step 3: Configure hierarchy with date column — github-2954 invariant (no NullError)', async () => {
    const errorsBefore = consoleErrors.length;
    const result = await page.evaluate(async () => {
      const grok = (window as any).grok;
      const tv = grok.shell.tv;
      let sunburst: any = null;
      for (const v of tv.viewers) if (v.type === 'Sunburst') { sunburst = v; break; }
      if (!sunburst) return {ok: false};
      const df = tv.dataFrame;
      let dateCol: string | null = null;
      for (const c of df.columns) {
        const t = String(c.type || '').toLowerCase();
        if (t.includes('date') || t.includes('time')) { dateCol = c.name; break; }
      }
      if (!dateCol) return {ok: false, reason: 'no-date-column'};
      let setOptionsThrew = false;
      try {
        sunburst.setOptions({hierarchyColumnNames: [dateCol]});
      } catch (e) { setOptionsThrew = true; }
      await new Promise((r) => setTimeout(r, 2000));
      const root = sunburst.root as HTMLElement;
      const rect = root.getBoundingClientRect();
      let readBack: any = null;
      try { readBack = sunburst.props.get('hierarchyColumnNames'); } catch (e) {}
      return {
        ok: true,
        dateCol,
        setOptionsThrew,
        readBack: readBack ? Array.from(readBack) : null,
        hasContent: root.children.length > 0,
        width: rect.width,
        height: rect.height,
      };
    });
    expect(result.ok).toBe(true);
    // github-2954 invariant: setOptions does NOT throw with date column.
    expect(result.setOptionsThrew).toBe(false);
    expect(result.hasContent).toBe(true);
    expect(result.width).toBeGreaterThan(0);
    const errorsDuring = consoleErrors.slice(errorsBefore);
    expect(errorsDuring).toEqual([]);
  });

  await softStep('Step 6: Recovery — set hierarchy to known-good string column', async () => {
    const result = await page.evaluate(async () => {
      const grok = (window as any).grok;
      const tv = grok.shell.tv;
      let sunburst: any = null;
      for (const v of tv.viewers) if (v.type === 'Sunburst') { sunburst = v; break; }
      if (!sunburst) return {ok: false};
      try {
        sunburst.setOptions({hierarchyColumnNames: ['AETERM']});
        await new Promise((r) => setTimeout(r, 1500));
      } catch (e) { return {ok: false, err: String(e).substring(0, 100)}; }
      const root = sunburst.root as HTMLElement;
      return {ok: true, hasContent: root.children.length > 0};
    });
    expect(result.ok).toBe(true);
    expect(result.hasContent).toBe(true);
  });

  await page.evaluate(() => (window as any).grok.shell.closeAll());

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
