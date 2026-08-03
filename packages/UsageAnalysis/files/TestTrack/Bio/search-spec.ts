import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
test.use(specTestOptions);
test('Bio Subsequence Search filter and reset on filter_FASTA', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:AppData/Bio/tests/filter_FASTA.csv');
    grok.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(() => resolve(), 4000);
    });
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    const hasMacromolecule = cols.some((c: any) => c.semType === 'Macromolecule');
    if (hasMacromolecule) {
      for (let i = 0; i < 60; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 5000));
    }
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
  await page.evaluate(async () => {
    const probes = ['Bio:getSeqHelper', 'Bio:getMonomerLibHelper', 'Bio:getBioLib'];
    for (const fn of probes) {
      try { await (grok as any).functions.call(fn, {}); return; } catch {  }
    }
    await new Promise((r) => setTimeout(r, 3000));
  });
  await page.evaluate(async () => {
    const names = [
      'Bio:bioSubstructureFilter', 'Bio:bioSubstructureFilterPanel',
      'Bio:subsequenceSearchTopMenu', 'Bio:subsequenceSearch',
    ];
    const findAny = (): boolean => {
      for (const n of names) {
        try {
          if ((grok as any).functions.find && (grok as any).functions.find(n)) return true;
        } catch {  }
      }
      return false;
    };
    const deadline = Date.now() + 15_000;
    while (Date.now() < deadline) {
      if (findAny()) return;
      await new Promise((r) => setTimeout(r, 300));
    }
    await new Promise((r) => setTimeout(r, 1500));
  });
  let baseRows = -1;
  await softStep('Open filter_FASTA.csv (14 rows, Macromolecule fasta column)', async () => {
    const info = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macro = cols.find((c: any) => c.semType === 'Macromolecule');
      return {
        rows: df.rowCount,
        macroName: macro ? (macro as any).name : null,
        units: macro ? (macro as any).getTag('units') : null,
      };
    });
    baseRows = info.rows;
    expect(info.rows).toBeGreaterThan(0);
    expect(info.macroName).not.toBeNull();
    expect(info.units).toBe('fasta');
  });
  await softStep('Open Bio > Search > Subsequence Search', async () => {
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 400));
      document.querySelector('[name="div-Bio---Search"]')!
        .dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 400));
      (document.querySelector('[name="div-Bio---Search---Subsequence-Search-..."]') as HTMLElement).click();
    });
    await page.locator('[name="viewer-Filters"] input[placeholder="Substructure"]')
      .waitFor({timeout: 30_000});
    const info = await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      return {
        filterCount: fg.filters.length,
        column: (fg.filters[0] as any).columnName ?? null,
      };
    });
    expect(info.filterCount).toBe(1);
    expect(info.column).not.toBeNull();
  });
  const subsequence = 'RTDEVSNHTHDKPTLTWFEEIFEEYHSP';
  await softStep(`Set Sequence filter to ${subsequence} -> 1 row`, async () => {
    await page.evaluate((seq) => {
      const input = document.querySelector(
        '[name="viewer-Filters"] input[placeholder="Substructure"]',
      ) as HTMLInputElement;
      const nativeSetter = Object.getOwnPropertyDescriptor(
        window.HTMLInputElement.prototype, 'value')!.set!;
      nativeSetter.call(input, seq);
      input.dispatchEvent(new Event('input', {bubbles: true}));
      input.dispatchEvent(new Event('change', {bubbles: true}));
      input.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', code: 'Enter', bubbles: true}));
      input.dispatchEvent(new KeyboardEvent('keyup', {key: 'Enter', code: 'Enter', bubbles: true}));
    }, subsequence);
    await page.waitForFunction(() => {
      const df = grok.shell.tv.dataFrame;
      return df.filter.trueCount < df.rowCount;
    }, null, {timeout: 30_000});
    const state = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const input = document.querySelector(
        '[name="viewer-Filters"] input[placeholder="Substructure"]',
      ) as HTMLInputElement;
      return {value: input.value, filtered: df.filter.trueCount, total: df.rowCount};
    });
    expect(state.value).toBe(subsequence);
    expect(state.filtered).toBe(1);
    expect(state.total).toBe(baseRows);
  });
  await softStep('Click Reset Filter -> all rows present', async () => {
    await page.locator('[name="viewer-Filters"] .d4-filter-group-header [name="icon-arrow-rotate-left"]')
      .click();
    await page.waitForFunction(() => {
      const df = grok.shell.tv.dataFrame;
      return df.filter.trueCount === df.rowCount;
    }, null, {timeout: 15_000});
    await page.waitForFunction(() => {
      const input = document.querySelector(
        '[name="viewer-Filters"] input[placeholder="Substructure"]',
      ) as HTMLInputElement | null;
      return input !== null && input.value === '';
    }, null, {timeout: 15_000});
    const state = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const input = document.querySelector(
        '[name="viewer-Filters"] input[placeholder="Substructure"]',
      ) as HTMLInputElement | null;
      return {value: input?.value ?? null, filtered: df.filter.trueCount, total: df.rowCount};
    });
    expect(state.filtered).toBe(state.total);
    expect(state.total).toBe(baseRows);
    expect(state.value).toBe('');
  });
  finishSpec();
});
