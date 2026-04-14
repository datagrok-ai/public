import {test, expect, chromium} from '@playwright/test';

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';
const datasetPath = 'System:DemoFiles/demog.csv';
const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Tree Viewer — Collaborative Filtering', async () => {
  const browser = await chromium.connectOverCDP('http://localhost:9222');
  const context = browser.contexts()[0];
  let page = context.pages().find(p => p.url().includes('datagrok'));
  if (!page) {
    page = await context.newPage();
    await page.goto(baseUrl, {waitUntil: 'networkidle', timeout: 60000});
    await page.waitForFunction(() => {
      try { return typeof grok !== 'undefined' && typeof grok.shell.closeAll === 'function'; }
      catch { return false; }
    }, {timeout: 45000});
  }

  // Setup: Open demog.csv, add Tree viewer, set hierarchy, open filters
  await softStep('Setup: Open demog.csv and add Tree viewer', async () => {
    const result = await page!.evaluate(async (path) => {
      document.querySelectorAll('.d4-dialog').forEach(d => {
        const cancel = d.querySelector('[name="button-CANCEL"]');
        if (cancel) cancel.click();
      });
      grok.shell.closeAll();
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      grok.shell.windows.simpleMode = false;

      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });

      const tree = tv.addViewer('Tree');
      await new Promise(r => setTimeout(r, 2000));
      tree.setOptions({hierarchyColumnNames: ['CONTROL', 'SEX', 'RACE']});
      await new Promise(r => setTimeout(r, 1000));

      tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 2000));

      return {rows: df.rowCount, viewers: tv.viewers.length};
    }, datasetPath);
    expect(result.rows).toBe(5850);
    expect(result.viewers).toBe(3); // Grid + Tree + Filters
  });

  // Step 1: Select branches false→F→Asian, false→F→Black, false→M→Asian
  await softStep('Step 1: Select three branches', async () => {
    const result = await page!.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const control = df.col('CONTROL');
      const sex = df.col('SEX');
      const race = df.col('RACE');
      df.selection.setAll(false);
      for (let i = 0; i < df.rowCount; i++) {
        const c = control.get(i);
        const s = sex.get(i);
        const r = race.get(i);
        if (c === false && s === 'F' && r === 'Asian') df.selection.set(i, true);
        if (c === false && s === 'F' && r === 'Black') df.selection.set(i, true);
        if (c === false && s === 'M' && r === 'Asian') df.selection.set(i, true);
      }
      df.selection.fireChanged();
      return {selected: df.selection.trueCount};
    });
    expect(result.selected).toBe(174);
  });

  // Step 2: Filter CONTROL=true → selected ∩ filtered = 0
  await softStep('Step 2: Filter CONTROL=true, expect overlap = 0', async () => {
    const result = await page!.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'CONTROL', selected: ['true']});
      await new Promise(r => setTimeout(r, 1000));

      const df = grok.shell.tv.dataFrame;
      let overlap = 0;
      for (let i = 0; i < df.rowCount; i++)
        if (df.filter.get(i) && df.selection.get(i)) overlap++;
      return {filtered: df.filter.trueCount, overlap};
    });
    expect(result.overlap).toBe(0);
  });

  // Step 3: Add true→F→Black to selection → selected ∩ filtered = 2
  await softStep('Step 3: Add true→F→Black, expect overlap = 2', async () => {
    const result = await page!.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const control = df.col('CONTROL');
      const sex = df.col('SEX');
      const race = df.col('RACE');
      for (let i = 0; i < df.rowCount; i++) {
        if (control.get(i) === true && sex.get(i) === 'F' && race.get(i) === 'Black')
          df.selection.set(i, true);
      }
      df.selection.fireChanged();
      await new Promise(r => setTimeout(r, 500));

      let overlap = 0;
      for (let i = 0; i < df.rowCount; i++)
        if (df.filter.get(i) && df.selection.get(i)) overlap++;
      return {selected: df.selection.trueCount, overlap};
    });
    expect(result.selected).toBe(176);
    expect(result.overlap).toBe(2);
  });

  // Step 4: Clear CONTROL filter → all rows visible, 176 selected
  await softStep('Step 4: Clear filter, expect 176 selected', async () => {
    const result = await page!.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true);
      df.filter.fireChanged();
      await new Promise(r => setTimeout(r, 500));
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'CONTROL'});
      await new Promise(r => setTimeout(r, 500));
      return {filtered: df.filter.trueCount, selected: df.selection.trueCount};
    });
    expect(result.filtered).toBe(5850);
    expect(result.selected).toBe(176);
  });

  // Cleanup
  await page!.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
