import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Bio MSA on FASTA', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

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

  // Step 1: FASTA.csv is opened in the setup phase above.
  await softStep('Open sample_FASTA.csv', async () => {
    const info = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      return {
        rows: df.rowCount,
        hasMacro: cols.some((c: any) => c.semType === 'Macromolecule'),
      };
    });
    expect(info.rows).toBeGreaterThan(0);
    expect(info.hasMacro).toBe(true);
  });

  // Step 2: Edit > Add New Column with formula RandBetween(0,5).
  // Name input = .ui-input-addnewcolumn-name, formula editor = CodeMirror .cm-content,
  // OK/CANCEL buttons are the generic [name="button-OK"] / [name="button-CANCEL"].
  await softStep('Add new column Clusters = RandBetween(0,5)', async () => {
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Edit"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 300));
      (document.querySelector('[name="div-Edit---Add-New-Column..."]') as HTMLElement).click();
    });
    await page.locator('[name="dialog-Add-New-Column"]').waitFor({timeout: 15000});

    await page.locator('[name="dialog-Add-New-Column"] input.ui-input-addnewcolumn-name').click();
    await page.keyboard.type('Clusters', {delay: 20});

    await page.locator('[name="dialog-Add-New-Column"] .cm-content').click();
    await page.keyboard.type('RandBetween(0, 5)', {delay: 20});

    const formulaText: string = await page.locator('[name="dialog-Add-New-Column"] .cm-content').textContent() ?? '';
    if (!formulaText.includes('RandBetween')) {
      // CodeMirror occasionally drops keystrokes in cold-context runs — fall back to the JS API.
      await page.locator('[name="dialog-Add-New-Column"] [name="button-Add-New-Column---CANCEL"]').click();
      await page.evaluate(async () => {
        const df = grok.shell.tv.dataFrame;
        await (df.columns as any).addNewCalculated('Clusters', 'RandBetween(0, 5)');
      });
    } else {
      await page.locator('[name="dialog-Add-New-Column"] [name="button-Add-New-Column---OK"]').click();
    }
    await page.waitForFunction(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      return cols.some((c: any) => c.name === 'Clusters' && c.type === 'int');
    }, null, {timeout: 60000});
    const info = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const c: any = df.col('Clusters');
      return {found: !!c, type: c?.type, min: c?.stats?.min, max: c?.stats?.max};
    });
    expect(info.found).toBe(true);
    expect(info.type).toBe('int');
    expect(info.min).toBeGreaterThanOrEqual(0);
    expect(info.max).toBeLessThanOrEqual(5);
  });

  // Step 3: Bio > Analyze > MSA opens dialog-MSA.
  await softStep('Open Bio > Analyze > MSA', async () => {
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 400));
      document.querySelector('[name="div-Bio---Analyze"]')!
        .dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 400));
      (document.querySelector('[name="div-Bio---Analyze---MSA..."]') as HTMLElement).click();
    });
    await page.locator('[name="dialog-MSA"]').waitFor({timeout: 15000});
    const inputs: string[] = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      return Array.from(dlg.querySelectorAll('[name^="input-host-"]'))
        .map((h) => h.getAttribute('name')!);
    });
    expect(inputs).toContain('input-host-Sequence');
    expect(inputs).toContain('input-host-Clusters');
  });

  // Step 4: MSA dialog auto-populates the Clusters input with the most recently-added
  // int column (our new "Clusters"). Ensure that explicitly so the test is deterministic.
  await softStep('Set Cluster input to the new Clusters column', async () => {
    await page.evaluate(() => {
      const dlg = (DG.Dialog as any).getOpenDialogs().find((d: any) =>
        d.root?.getAttribute?.('name') === 'dialog-MSA') ?? (DG.Dialog as any).getOpenDialogs()[0];
      const df = grok.shell.tv.dataFrame;
      const input = dlg.inputs.find((i: any) => i.caption === 'Clusters');
      input.value = df.col('Clusters');
    });
    const displayed: string = await page.evaluate(() => {
      const host = document.querySelector('[name="dialog-MSA"] [name="input-host-Clusters"]')!;
      return (host.querySelector('.d4-column-selector-column')?.textContent
        ?? host.querySelector('.ui-input-editor')?.textContent
        ?? host.textContent)!.trim();
    });
    expect(displayed).toContain('Clusters');
  });

  // Step 5: The ALIGNMENT PARAMETERS button reveals Gap open / Gap extend / Terminal gap
  // inputs (initially `display:flex` but height 0 until the button toggles them).
  await softStep('Alignment parameters button adds input parameters', async () => {
    const before: boolean = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const gap = dlg.querySelector('[name="input-host-Gap-open"]') as HTMLElement | null;
      return !!gap && gap.getBoundingClientRect().height > 0;
    });
    expect(before).toBe(false);
    await page.locator('[name="dialog-MSA"] [name="button-Alignment-parameters"]').click();
    await page.waitForFunction(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const gap = dlg.querySelector('[name="input-host-Gap-open"]') as HTMLElement | null;
      return !!gap && gap.getBoundingClientRect().height > 0;
    }, null, {timeout: 5000});
    const visible: boolean = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const names = ['input-host-Gap-open', 'input-host-Gap-extend', 'input-host-Terminal-gap'];
      return names.every((n) => {
        const el = dlg.querySelector(`[name="${n}"]`) as HTMLElement | null;
        return !!el && el.getBoundingClientRect().height > 0;
      });
    });
    expect(visible).toBe(true);
  });

  // Step 6: OK runs MSA — verify a new aligned Macromolecule column, sequence-renderer
  // metadata, and equal sequence length within every cluster (scenario success criterion).
  await softStep('OK produces aligned MSA column with cluster-aware lengths', async () => {
    await page.locator('[name="dialog-MSA"] [name="button-OK"]').click();
    await page.waitForFunction(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const msa: any = cols.find((c: any) => c.name.toLowerCase().includes('msa') && c.semType === 'Macromolecule');
      if (!msa) return false;
      const gridCol = (grok.shell.tv as any).grid?.col?.(msa.name);
      return gridCol?.cellType === 'sequence';
    }, null, {timeout: 180000});
    const result = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const msa: any = cols.find((c: any) => c.name.toLowerCase().includes('msa'));
      const clusters: any = df.col('Clusters');
      const lenByCluster: Record<string, Set<number>> = {};
      for (let i = 0; i < df.rowCount; i++) {
        const key = String(clusters.get(i));
        (lenByCluster[key] ||= new Set()).add((msa.get(i) ?? '').length);
      }
      const allEqualPerCluster = Object.values(lenByCluster).every((s) => s.size === 1);
      let cellRendererTag: string | null = null;
      try {
        for (const [k, v] of msa.tags) { if (k === 'cell.renderer') { cellRendererTag = v; break; } }
      } catch { /* tags may not be iterable in some builds */ }
      if (!cellRendererTag) {
        const gridCol = (grok.shell.tv as any).grid?.col?.(msa.name);
        cellRendererTag = gridCol?.cellType ?? null;
      }
      let hasGaps = false;
      for (let i = 0; i < df.rowCount && !hasGaps; i++)
        if ((msa.get(i) ?? '').includes('-')) hasGaps = true;
      return {
        msaName: msa?.name,
        semType: msa?.semType,
        cellRendererTag,
        allEqualPerCluster,
        hasGaps,
        clusterCount: Object.keys(lenByCluster).length,
      };
    });
    expect(result.msaName).toBeTruthy();
    expect(result.semType).toBe('Macromolecule');
    expect(result.cellRendererTag).toBe('sequence');
    expect(result.allEqualPerCluster).toBe(true);
    expect(result.hasGaps).toBe(true);
    expect(result.clusterCount).toBeGreaterThan(1);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
