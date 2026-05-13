import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Bio Subsequence Search on FASTA', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  // Step 1: open sample_FASTA.csv — loaded from the Bio samples folder.
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:AppData/Bio/samples/FASTA.csv');
    grok.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(() => resolve(), 3000);
    });
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    const hasBioChem = cols.some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 5000));
    }
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Open sample_FASTA.csv', async () => {
    const info = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      return {
        rows: df.rowCount,
        hasMacro: cols.some((c: any) => c.name === 'Sequence' && c.semType === 'Macromolecule'),
      };
    });
    expect(info.rows).toBe(64);
    expect(info.hasMacro).toBe(true);
  });

  // Step 2: Bio > Search > Subsequence Search... opens the filter panel with a
  // single Bio substructure filter card scoped to the Sequence column (no dialog).
  // Menu item name ends with `-...` because the label is "Subsequence Search " (trailing space).
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
      .waitFor({timeout: 15000});
    const info = await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      return {
        filterCount: fg.filters.length,
        column: (fg.filters[0] as any).columnName ?? null,
      };
    });
    expect(info.filterCount).toBe(1);
    expect(info.column).toBe('Sequence');
  });

  // Step 3: typing the sequence narrows the table to rows whose Sequence cell
  // contains the subsequence. FASTA sample has exactly one matching row.
  await softStep('Set Sequence filter to MHAILRYFIRRLFYHIFYKIYSLISKKHQSLPSDVRQF', async () => {
    await page.evaluate(() => {
      const input = document.querySelector(
        '[name="viewer-Filters"] input[placeholder="Substructure"]',
      ) as HTMLInputElement;
      const seq = 'MHAILRYFIRRLFYHIFYKIYSLISKKHQSLPSDVRQF';
      const nativeSetter = Object.getOwnPropertyDescriptor(
        window.HTMLInputElement.prototype, 'value')!.set!;
      nativeSetter.call(input, seq);
      input.dispatchEvent(new Event('input', {bubbles: true}));
      input.dispatchEvent(new Event('change', {bubbles: true}));
      input.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', code: 'Enter', bubbles: true}));
      input.dispatchEvent(new KeyboardEvent('keyup', {key: 'Enter', code: 'Enter', bubbles: true}));
    });
    await page.waitForFunction(() => {
      const df = grok.shell.tv.dataFrame;
      return df.filter.trueCount < df.rowCount;
    }, null, {timeout: 15000});
    const state = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const input = document.querySelector(
        '[name="viewer-Filters"] input[placeholder="Substructure"]',
      ) as HTMLInputElement;
      return {value: input.value, filtered: df.filter.trueCount, total: df.rowCount};
    });
    expect(state.value).toBe('MHAILRYFIRRLFYHIFYKIYSLISKKHQSLPSDVRQF');
    expect(state.filtered).toBe(1);
    expect(state.total).toBe(64);
  });

  // Step 4: Reset (↺) icon on the filter group header clears all filter values.
  // No confirmation popup is shown — substructure input becomes empty and the
  // table is back to the full row count.
  await softStep('Click Reset Filter', async () => {
    await page.locator(
      '[name="viewer-Filters"] .d4-filter-group-header [name="icon-arrow-rotate-left"]',
    ).click();
    await page.waitForFunction(() => {
      const df = grok.shell.tv.dataFrame;
      return df.filter.trueCount === df.rowCount;
    }, null, {timeout: 10000});
    const state = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const input = document.querySelector(
        '[name="viewer-Filters"] input[placeholder="Substructure"]',
      ) as HTMLInputElement | null;
      return {value: input?.value ?? null, filtered: df.filter.trueCount, total: df.rowCount};
    });
    expect(state.value).toBe('');
    expect(state.filtered).toBe(state.total);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
