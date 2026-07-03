/* ---
sub_features_covered: [bio.analyze.sequence-space.transform, bio.api.get-seq-helper, bio.detector, bio.io.fasta-handler, bio.io.save-as-fasta, bio.rendering, bio.transform.convert-notation, bio.transform.convert-notation.action]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import * as bio from '../helpers/bio';
import {
  saveAllTablesWithProvenance,
  reopenAndAssertProvenance,
  deleteProjectWithCleanup,
} from '../helpers/projects';
test.use(specTestOptions);
test('Bio macromolecule_column source-class lifecycle: detect → convert → fasta round-trip → save+reopen', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;
  const stamp = Date.now();
  const projectName = `bio-lifecycle-macromolecule-${stamp}`;
  const fastaTempPath = `System:AppData/UsageAnalysis/temp/bio-lifecycle-${stamp}.fasta`;
  let saved: {projectId: string; primaryTableInfoId: string; layoutId: string | null} | null = null;
  await loginToDatagrok(page);
  // Scenario 1 — Detect on open + convert-notation round trip
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
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
    }
  }, 'System:AppData/Bio/tests/filter_HELM.csv');
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
  await page.evaluate(async () => {
    const probes = ['Bio:getSeqHelper', 'Bio:getMonomerLibHelper', 'Bio:getBioLib'];
    for (let i = 0; i < 30; i++) {
      for (const fn of probes) {
        try { await (grok as any).functions.call(fn, {}); return; } catch { /* try next */ }
      }
      await new Promise((r) => setTimeout(r, 200));
    }
  });
  const bioCellTypes = ['sequence', 'helm', 'separator', 'biln', 'custom', 'fasta'];
  await softStep('S1.2: Macromolecule detector classifies HELM column synchronously (units=helm)', async () => {
    await page.waitForFunction((accepted: string[]) => {
      const df = grok.shell.tv?.dataFrame;
      if (!df) return false;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
      if (!macro) return false;
      const ct = (grok.shell.tv as any).grid?.col?.(macro.name)?.cellType ?? null;
      return ct !== null && accepted.indexOf(ct) >= 0;
    }, bioCellTypes, {timeout: 60_000});
    const info = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
      return {
        hasMacro: !!macro,
        units: macro?.getTag?.('units') ?? macro?.meta?.units ?? null,
        gridCellType: (grok.shell.tv as any).grid?.col?.(macro?.name)?.cellType ?? null,
      };
    });
    expect(info.hasMacro).toBe(true);
    expect(info.units).toBe('helm');
    expect(bioCellTypes, `HELM renderer must dispatch a bio sequence cellType, got ${info.gridCellType}`)
      .toContain(info.gridCellType!);
  });
  await softStep('S1.3-1.4: Convert HELM → SEPARATOR via top-menu; new Macromolecule column appears with units=separator', async () => {
    const before: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
    await page.evaluate(() => (document.querySelector('[name="div-Bio"]') as HTMLElement).click());
    await page.locator('[name="div-Bio---Transform"]').waitFor({state: 'attached', timeout: 15_000});
    await page.evaluate(() => document.querySelector('[name="div-Bio---Transform"]')!
      .dispatchEvent(new MouseEvent('mouseover', {bubbles: true})));
    await page.locator('[name="div-Bio---Transform---Convert-Sequence-Notation..."]').waitFor({state: 'attached', timeout: 15_000});
    await page.evaluate(() => (document.querySelector(
      '[name="div-Bio---Transform---Convert-Sequence-Notation..."]') as HTMLElement).click());
    const dlg = page.locator('[name="dialog-Convert-Sequence-Notation"]');
    await dlg.waitFor({timeout: 60_000});
    await dlg.locator('[name="input-host-Convert-to"] select').selectOption('separator');
    await dlg.locator('[name="input-host-Separator"] select').selectOption('-');
    await dlg.locator('[name="button-OK"]').click();
    await page.waitForFunction(
      (b) => grok.shell.tv.dataFrame.columns.length > b, before, {timeout: 60_000});
    const info: {macroCount: number, lastUnits: string | null, allUnits: string[]} =
      await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macroCols = cols.filter((c: any) => c.semType === 'Macromolecule');
        const last: any = macroCols[macroCols.length - 1];
        return {
          macroCount: macroCols.length,
          lastUnits: last?.meta?.units ?? null,
          allUnits: macroCols.map((c: any) => c.meta?.units ?? null).filter((u: any) => u !== null),
        };
      });
    expect(info.macroCount).toBeGreaterThanOrEqual(2);
    expect(info.lastUnits).toBe('separator');
    await page.waitForFunction(
      () => document.querySelectorAll('[name="dialog-Convert-Sequence-Notation"]').length === 0,
      null, {timeout: 15_000}).catch(() => {});
  });
  // Scenario 2 — Import FASTA → Export FASTA → re-import round trip
  await page.evaluate(async (path) => {
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    grok.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(() => resolve(), 4000);
    });
  }, 'System:AppData/Bio/tests/filter_FASTA.csv');
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
  await page.waitForFunction(() => {
    const df = grok.shell.tv?.dataFrame;
    if (!df) return false;
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    return cols.some((c: any) => c.semType === 'Macromolecule' &&
      (c.getTag?.('units') ?? c.meta?.units) === 'fasta');
  }, null, {timeout: 30_000});
  await softStep('S2.1: filter_FASTA.csv opens with Macromolecule semType (units=fasta)', async () => {
    const info = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
      return {
        hasMacro: !!macro,
        units: macro?.meta?.units ?? null,
        rowCount: df.rowCount,
        colName: macro?.name ?? null,
      };
    });
    expect(info.hasMacro).toBe(true);
    expect(info.units).toBe('fasta');
    expect(info.rowCount).toBeGreaterThan(0);
  });
  await softStep('S2.2-2.4: saveAsFastaDo → write → re-readCsv produces a Macromolecule round-trip', async () => {
    const result = await page.evaluate(async ({tempPath}) => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const seqCol: any = cols.find((c: any) => c.semType === 'Macromolecule');
      if (!seqCol) throw new Error('S2.2: no Macromolecule column found');
      const idCol: any = cols.find((c: any) => c.semType !== 'Macromolecule') ?? null;
      const idColList = idCol ? [idCol] : [];
      const seqHelper: any = await (grok as any).functions.call('Bio:getSeqHelper', {});
      const seqHandler: any = seqHelper.getSeqHandler(seqCol);
      const canonSeq = (handler: any, row: number): string => {
        const ss: any = handler.getSplitted(row);
        const mons: string[] = [];
        for (let p = 0; p < ss.length; p++) mons.push(ss.getOriginal(p));
        return mons.join('/');
      };
      const originalFirstSeq: string = canonSeq(seqHandler, 0);
      const fastaLines: string[] = [];
      const lineWidth = 60;
      for (let rowIdx = 0; rowIdx < seqHandler.length; rowIdx++) {
        const seqId: string = idColList.length > 0
          ? idColList.map((c: any) => c.get(rowIdx)?.toString() ?? '').join('|')
          : String(rowIdx + 1);
        const srcSS: any = seqHandler.getSplitted(rowIdx);
        const monomers: string[] = [];
        for (let p = 0; p < srcSS.length; p++)
          monomers.push(srcSS.getOriginal(p));
        const seqText: string = monomers.map((om: string) => om.length > 1 ? `[${om}]` : om).join('');
        fastaLines.push(`>${seqId}\n`);
        for (let i = 0; i < seqText.length; i += lineWidth)
          fastaLines.push(seqText.slice(i, i + lineWidth) + '\n');
      }
      const fastaText: string = fastaLines.join('');
      if (!fastaText.startsWith('>'))
        throw new Error('S2.2: exported FASTA does not start with > header');
      if (fastaText.length < 4)
        throw new Error('S2.2: exported FASTA is empty');
      let writeErr: string | null = null;
      try {
        await grok.dapi.files.writeAsText(tempPath, fastaText);
      } catch (e) {
        writeErr = String(e).slice(0, 200);
      }
      let reimportErr: string | null = null;
      let reimportedSummary: {rowCount: number, cols: number, firstColSem: string | null} | null = null;
      let reimportedFirstSeq: string | null = null;
      try {
        const dfs: any = await (grok as any).functions.call('Bio:importFasta', {fileContent: fastaText});
        const reimported: any = Array.isArray(dfs) ? dfs[0] : dfs;
        const rCols = Array.from({length: reimported.columns.length}, (_, i) => reimported.columns.byIndex(i));
        const rMacro: any = rCols.find((c: any) => c.semType === 'Macromolecule')
          ?? reimported.columns.byIndex(reimported.columns.length - 1);
        reimportedFirstSeq = canonSeq(seqHelper.getSeqHandler(rMacro), 0);
        reimportedSummary = {
          rowCount: reimported.rowCount,
          cols: reimported.columns.length,
          firstColSem: rMacro?.semType ?? null,
        };
      } catch (e) {
        reimportErr = String(e).slice(0, 200);
      }
      try { await grok.dapi.files.delete(tempPath); } catch (_) { /* best effort */ }
      return {
        fastaShape: {
          startsWithHeader: fastaText.startsWith('>'),
          lineCount: fastaText.split('\n').length,
          totalLen: fastaText.length,
        },
        writeErr,
        reimported: reimportedSummary,
        reimportErr,
        reimportedFirstSeq,
        originalFirstSeq,
        originalRowCount: df.rowCount,
      };
    }, {tempPath: fastaTempPath});
    expect(result.fastaShape.startsWithHeader).toBe(true);
    expect(result.fastaShape.totalLen).toBeGreaterThan(0);
    expect(result.writeErr, `FASTA write failed: ${result.writeErr}`).toBeNull();
    expect(result.reimportErr, `Bio:importFasta failed: ${result.reimportErr}`).toBeNull();
    expect(result.reimported).not.toBeNull();
    expect(result.reimported!.rowCount).toBe(result.originalRowCount);
    expect(result.reimported!.firstColSem).toBe('Macromolecule');
    expect(result.reimportedFirstSeq, 'FASTA round-trip must preserve first sequence content')
      .toBe(result.originalFirstSeq);
  });
  // Scenario 3 — Save project with analysis + reopen restores analysis output
  await softStep('S3.1: Open Bio | Analyze | Sequence Space with defaults — embedding columns + ScatterPlot dock', async () => {
    const baseCols: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
    await bio.openBioAnalyze(page, 'div-Bio---Analyze---Sequence-Space...');
    await page.locator('[name="dialog-Sequence-Space"] [name="button-OK"]').waitFor({timeout: 60_000});
    await page.locator('[name="dialog-Sequence-Space"] [name="button-OK"]').click();
    await page.waitForFunction(
      (base) => grok.shell.tv.dataFrame.columns.length > base &&
        Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
      baseCols, {timeout: 240_000});
    const result = await page.evaluate(() => ({
      hasScatter: Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
      cols: grok.shell.tv.dataFrame.columns.length,
      colNames: Array.from({length: grok.shell.tv.dataFrame.columns.length},
        (_, i) => grok.shell.tv.dataFrame.columns.byIndex(i).name),
    }));
    expect(result.hasScatter).toBe(true);
    expect(result.cols).toBeGreaterThan(baseCols);
    const hasEmbedX = result.colNames.some((n: string) => /^Embed_X_\d+$/.test(n));
    const hasEmbedY = result.colNames.some((n: string) => /^Embed_Y_\d+$/.test(n));
    expect(hasEmbedX).toBe(true);
    expect(hasEmbedY).toBe(true);
  });
  try {
    await softStep('S3.3: Save project with provenance (JS API path)', async () => {
      saved = await saveAllTablesWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
      expect(saved.primaryTableInfoId).toBeTruthy();
    });
    await softStep('S3.4-3.5: reopen project — embedding columns + Macromolecule semType survive', async () => {
      if (!saved) throw new Error('S3.3 did not produce a saved project');
      const result = await reopenAndAssertProvenance(page, saved.projectId);
      expect(result.tablesAfter).toBeGreaterThan(0);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
      const post = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
        const colNames = cols.map((c: any) => c.name);
        return {
          hasMacro: !!macro,
          units: macro?.meta?.units ?? null,
          hasEmbedX: colNames.some((n: string) => /^Embed_X_\d+$/.test(n)),
          hasEmbedY: colNames.some((n: string) => /^Embed_Y_\d+$/.test(n)),
          hasScatter: Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
        };
      });
      expect(post.hasMacro).toBe(true);
      expect(post.units).toBe('fasta');
      expect(post.hasEmbedX).toBe(true);
      expect(post.hasEmbedY).toBe(true);
    });
  } finally {
    // Scenario 4 — Cleanup (runs regardless of earlier failures)
    if (saved) {
      await deleteProjectWithCleanup(page, {
        projectId: saved.projectId,
        tableInfoId: saved.primaryTableInfoId,
      });
    }
    await page.evaluate(async (p) => {
      try { await grok.dapi.files.delete(p); } catch (_) { /* best effort */ }
    }, fastaTempPath).catch(() => {});
  }
  finishSpec();
});
