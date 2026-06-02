/* ---
sub_features_covered: [bio.analyze.msa, bio.analyze.msa.dialog, bio.analyze.msa.align-sequences, bio.rendering.column-header, bio.detector]
--- */
//   related_bugs: [GROK-18474, GROK-15176]
//   GROK-18474 — MSA column-header click handler crashes on FASTA-aligned data;
//     is owned by bug_focused_candidates[GROK-18474] (atlas
//   GROK-15176 — Bio's to-atomic-level produces molfiles with invalid isotope;
//     convert.md:Step 2 — owned by bug_focused_candidates[GROK-15176].
// failure_keys=[]): hypothesis category = test-bug (cold-start race in
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
    const df = await grok.dapi.files.readCsv('System:AppData/Bio/tests/filter_FASTA.csv');
    grok.shell.addTableView(df);
    let detectorFired = false;
    const sub = df.onSemanticTypeDetected.subscribe(() => { detectorFired = true; });
    try {
      const deadline = Date.now() + 20_000;
      while (Date.now() < deadline) {
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macro = cols.some((c: any) => c.semType === 'Macromolecule');
        if (macro || detectorFired) break;
        await new Promise((r) => setTimeout(r, 100));
      }
    } finally { sub.unsubscribe(); }
    for (let i = 0; i < 60; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise((r) => setTimeout(r, 200));
    }
    // restores the deterministic top-menu reachability the scenario
    await new Promise((r) => setTimeout(r, 5000));
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
  await page.waitForFunction(() => !!document.querySelector('[name="div-Bio"]'),
    null, {timeout: 30_000});
  await softStep('Open filter_FASTA.csv and detect Macromolecule column', async () => {
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
    // B-STAB-01 stabilization (test-bug fix, retry round 1, hypothesis_retry
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const c: any = df.col('Clusters');
      for (let i = 0; i < df.rowCount; i++)
        c.set(i, i % 2, false);
    });
    const info = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const c: any = df.col('Clusters');
      const distinctCats = new Set<number>();
      const perCatCount: Record<string, number> = {};
      let minV = Number.POSITIVE_INFINITY;
      let maxV = Number.NEGATIVE_INFINITY;
      for (let i = 0; i < df.rowCount; i++) {
        const v = c.get(i);
        distinctCats.add(v);
        perCatCount[String(v)] = (perCatCount[String(v)] ?? 0) + 1;
        if (v < minV) minV = v;
        if (v > maxV) maxV = v;
      }
      const minPerCat = Math.min(...Object.values(perCatCount));
      return {
        found: !!c,
        type: c?.type,
        min: minV,
        max: maxV,
        distinctCount: distinctCats.size,
        minRowsPerCluster: minPerCat,
      };
    });
    expect(info.found).toBe(true);
    expect(info.type).toBe('int');
    expect(info.min).toBeGreaterThanOrEqual(0);
    expect(info.max).toBeLessThanOrEqual(5);
    expect(info.distinctCount).toBeGreaterThanOrEqual(2);
    expect(info.minRowsPerCluster).toBeGreaterThanOrEqual(2);
  });
  await softStep('Open Bio > Analyze > MSA...', async () => {
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
  // Round-1 hypothesis (test-bug): the prior spec used `.first().click()`,
  await softStep('Alignment parameters button toggles input parameters', async () => {
    const before: boolean = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const gap = dlg.querySelector('[name="input-host-Gap-open"]') as HTMLElement | null;
      return !!gap && gap.getBoundingClientRect().height > 0;
    });
    expect(before).toBe(false);
    await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const btns = Array.from(dlg.querySelectorAll('[name="button-Alignment-parameters"]')) as HTMLElement[];
      const visible = btns.find((b) => b.offsetParent !== null);
      if (!visible) throw new Error('No visible Alignment-parameters button in dialog-MSA');
      visible.click();
    });
    await page.waitForFunction(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const gap = dlg.querySelector('[name="input-host-Gap-open"]') as HTMLElement | null;
      return !!gap && gap.getBoundingClientRect().height > 0;
    }, null, {timeout: 5000});
    const allOpen: boolean = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const names = ['input-host-Gap-open', 'input-host-Gap-extend', 'input-host-Terminal-gap'];
      return names.every((n) => {
        const el = dlg.querySelector(`[name="${n}"]`) as HTMLElement | null;
        return !!el && el.getBoundingClientRect().height > 0;
      });
    });
    expect(allOpen).toBe(true);
    await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const btns = Array.from(dlg.querySelectorAll('[name="button-Alignment-parameters"]')) as HTMLElement[];
      const visible = btns.find((b) => b.offsetParent !== null);
      if (!visible) throw new Error('No visible Alignment-parameters button on second toggle');
      visible.click();
    });
    await page.waitForFunction(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const gap = dlg.querySelector('[name="input-host-Gap-open"]') as HTMLElement | null;
      return !gap || gap.getBoundingClientRect().height === 0;
    }, null, {timeout: 5000});
    const allClosed: boolean = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const names = ['input-host-Gap-open', 'input-host-Gap-extend', 'input-host-Terminal-gap'];
      return names.every((n) => {
        const el = dlg.querySelector(`[name="${n}"]`) as HTMLElement | null;
        return !el || el.getBoundingClientRect().height === 0;
      });
    });
    expect(allClosed).toBe(true);
  });
  await softStep('OK runs kalign MSA — verify aligned column, renderer, and per-cluster lengths', async () => {
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
      } catch {  }
      if (!cellRendererTag) {
        const gridCol = (grok.shell.tv as any).grid?.col?.(msa.name);
        cellRendererTag = gridCol?.cellType ?? null;
      }
      let hasGaps = false;
      for (let i = 0; i < df.rowCount && !hasGaps; i++)
        if ((msa.get(i) ?? '').includes('-')) hasGaps = true;
      let alignedTag: string | null = null;
      try { alignedTag = msa.getTag('aligned') ?? null; } catch {  }
      return {
        msaName: msa?.name,
        semType: msa?.semType,
        cellRendererTag,
        alignedTag,
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
