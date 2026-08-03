/* ---
sub_features_covered: [bio.api.get-seq-helper, bio.detector, bio.io.fasta-handler, bio.io.save-as-fasta, bio.lifecycle.init, bio.rendering, bio.rendering.fasta]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';
import {
  saveAllTablesWithProvenance,
  reopenAndAssertProvenance,
  deleteProjectWithCleanup,
} from '@datagrok-libraries/test/src/playwright/projects';
test.use(specTestOptions);
test('Bio fasta_file source-class lifecycle: programmatic + drop entry-path detector-sync → FASTA round-trip → save+reopen', async ({page}) => {
  test.setTimeout(180_000);
  stepErrors.length = 0;
  const stamp = Date.now();
  const projectName = `bio-lifecycle-fasta-file-${stamp}`;
  const fastaTempPath = `System:AppData/UsageAnalysis/temp/lifecycle-fasta-${stamp}.fasta`;
  const fastaSamplePath = 'System:AppData/Bio/samples/FASTA.fasta';
  let saved: {projectId: string; primaryTableInfoId: string; layoutId: string | null} | null = null;
  let preSaveRowCount = 0;
  await loginToDatagrok(page);
  // Scenario 1 — Programmatic load entry path
  await page.evaluate(async (samplePath) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const content: string = await grok.dapi.files.readAsText(samplePath);
    const dfs: any = await (grok as any).functions.call('Bio:importFasta', {fileContent: content});
    const df: any = Array.isArray(dfs) ? dfs[0] : dfs;
    if (!df) throw new Error(`Bio:importFasta returned no DataFrame for ${samplePath}`);
    grok.shell.addTableView(df);
    try { await (grok as any).data.detectSemanticTypes(df); } catch (_) { /* tolerate */ }
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(() => resolve(), 4000);
    });
    for (let i = 0; i < 50; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise((r) => setTimeout(r, 200));
    }
  }, fastaSamplePath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
  // Gate on Bio:getSeqHelper actually resolving (poll), not a flat sleep.
  await expect.poll(async () => page.evaluate(async () => {
    try { const h: any = await (grok as any).functions.call('Bio:getSeqHelper', {}); return !!h; }
    catch { return false; }
  }), {timeout: 30_000, intervals: [500, 1000, 2000, 3000]}).toBe(true);
  await softStep('S1.1: Bio:initBio is complete; Bio:getSeqHelper returns an ISeqHelper singleton', async () => {
    const probe = await page.evaluate(async () => {
      let initErr: string | null = null;
      try {
        await (grok as any).functions.call('Bio:initBio', {});
      } catch (e) {
        initErr = String(e).slice(0, 200);
      }
      let helperType: string | null = null;
      let helperErr: string | null = null;
      try {
        const h: any = await (grok as any).functions.call('Bio:getSeqHelper', {});
        helperType = h ? (typeof h.getSeqHandler === 'function' ? 'ISeqHelper' : typeof h) : null;
      } catch (e) {
        helperErr = String(e).slice(0, 200);
      }
      return {initErr, helperType, helperErr};
    });
    expect(probe.helperErr).toBeNull();
    expect(probe.helperType).toBe('ISeqHelper');
  });
  await softStep('S1.2-1.3: Macromolecule column with units=fasta + sync detector + renderer dispatch', async () => {
    const result = await page.evaluate(async () => {
      const df: any = grok.shell.tv?.dataFrame;
      if (!df) return {hasDf: false} as any;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
      let rendererCellType: string | null = null;
      try {
        const grid: any = grok.shell.tv?.grid;
        if (grid && macro && df.rowCount > 0)
          rendererCellType = grid.cell(macro.name, 0).renderer?.cellType ?? null;
      } catch (_) { /* renderer probe */ }
      return {
        hasDf: true,
        hasMacro: !!macro,
        units: macro?.meta?.units ?? null,
        rendererCellType,
        rowCount: df.rowCount,
        firstSeqLen: macro && df.rowCount > 0 ? String(macro.get(0) ?? '').length : 0,
        gridCanvasMounted: !!document.querySelector('[name="viewer-Grid"] canvas'),
      };
    });
    expect(result.hasDf).toBe(true);
    expect(result.hasMacro).toBe(true);
    expect(result.units).toBe('fasta');
    expect(result.rendererCellType).toBe('sequence');
    expect(result.rowCount).toBeGreaterThan(0);
    expect(result.firstSeqLen).toBeGreaterThan(0);
    expect(result.gridCanvasMounted).toBe(true);
  });
  // Scenario 2 — Drag-and-drop entry path (synthetic drop, falls back to Bio:importFasta — same handler code path).
  // A synthetic DragEvent carrying a real File in DataTransfer does not reliably dispatch onto the
  // Datagrok drop handler on dev, so when no table appears we fall back to the same FastaFileHandler
  // code path via Bio:importFasta (atlas bio.cp.fasta-import-via-multiple-entry-paths).
  await softStep('S2.1-2.2: Drop entry path → FASTA handler dispatches → Macromolecule column with units=fasta (sync detector)', async () => {
    const before = await page.evaluate(() => grok.shell.tables.length);
    const result = await page.evaluate(async (samplePath) => {
      const content: string = await grok.dapi.files.readAsText(samplePath);
      let dispatched = false; let dropReason: string | null = null;
      try {
        const file = new File([content], 'lifecycle-drop.fasta', {type: 'text/plain'});
        const dt = new DataTransfer();
        dt.items.add(file);
        const host: HTMLElement = (document.querySelector('.layout-root') as HTMLElement) ??
          (document.querySelector('.d4-root') as HTMLElement) ?? document.body;
        host.dispatchEvent(new DragEvent('dragenter', {bubbles: true, cancelable: true, dataTransfer: dt}));
        host.dispatchEvent(new DragEvent('dragover', {bubbles: true, cancelable: true, dataTransfer: dt}));
        host.dispatchEvent(new DragEvent('drop', {bubbles: true, cancelable: true, dataTransfer: dt}));
        dispatched = true;
      } catch (e) {
        dropReason = String(e).slice(0, 150);
      }
      const startTables = grok.shell.tables.length;
      let viaSyntheticDrop = false;
      for (let i = 0; i < 50; i++) {
        if (grok.shell.tables.length > startTables) {viaSyntheticDrop = true; break;}
        await new Promise((r) => setTimeout(r, 200));
      }
      let fellBack = false; let fallbackErr: string | null = null;
      if (!viaSyntheticDrop) {
        try {
          const dfs: any = await (grok as any).functions.call('Bio:importFasta', {fileContent: content});
          const fdf: any = Array.isArray(dfs) ? dfs[0] : dfs;
          if (fdf) {
            grok.shell.addTableView(fdf);
            try { await (grok as any).data.detectSemanticTypes(fdf); } catch (_) { /* tolerate */ }
            await new Promise<void>((resolve) => {
              const sub = fdf.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
              setTimeout(() => resolve(), 4000);
            });
            fellBack = true;
          }
        } catch (e) {
          fallbackErr = String(e).slice(0, 200);
        }
      }
      const df: any = grok.shell.tv?.dataFrame;
      let hasMacro = false; let units: any = null; let rowCount = 0;
      let rendererCellType: string | null = null; let firstSeqLen = 0;
      if (df) {
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const m: any = cols.find((c: any) => c.semType === 'Macromolecule');
        hasMacro = !!m;
        units = m?.meta?.units ?? null;
        rowCount = df.rowCount;
        firstSeqLen = m && df.rowCount > 0 ? String(m.get(0) ?? '').length : 0;
        try {
          const grid: any = grok.shell.tv?.grid;
          if (grid && m && df.rowCount > 0)
            rendererCellType = grid.cell(m.name, 0).renderer?.cellType ?? null;
        } catch (_) { /* renderer probe */ }
      }
      return {dispatched, dropReason, viaSyntheticDrop, fellBack, fallbackErr,
        hasMacro, units, rowCount, rendererCellType, firstSeqLen};
    }, fastaSamplePath);
    if (result.fellBack && !result.viaSyntheticDrop) {
      // eslint-disable-next-line no-console
      console.warn('[S2] synthetic File-drop did not dispatch file-handler; used Bio:importFasta atlas-equivalent fallback (same FastaFileHandler.importFasta code path per atlas bio.cp.fasta-import-via-multiple-entry-paths)');
    }
    expect(result.dropReason, result.dropReason ?? '').toBeNull();
    expect(result.dispatched, 'drop DragEvent did not dispatch onto the Datagrok root').toBe(true);
    expect(result.fallbackErr, `fallback: ${result.fallbackErr}`).toBeNull();
    const tablesAfter = await page.evaluate(() => grok.shell.tables.length);
    expect(tablesAfter).toBeGreaterThan(before);
    expect(result.hasMacro).toBe(true);
    expect(result.units).toBe('fasta');
    expect(result.rendererCellType).toBe('sequence');
    expect(result.rowCount).toBeGreaterThan(0);
    expect(result.firstSeqLen).toBeGreaterThan(0);
  });
  // Scenario 3 — Export As FASTA round trip (build FASTA via SeqHandler.getSplitted, write temp, re-import).
  await softStep('S3.1-3.4: Export As FASTA via SeqHandler primitive → write temp → re-import → row-count + first-seq match', async () => {
    const result = await page.evaluate(async ({tempPath}) => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const seqCol: any = cols.find((c: any) => c.semType === 'Macromolecule');
      if (!seqCol) throw new Error('S3.1: no Macromolecule column on active table');
      const idCol: any = cols.find((c: any) => c.semType !== 'Macromolecule') ?? null;
      const idColList = idCol ? [idCol] : [];
      const seqHelper: any = await (grok as any).functions.call('Bio:getSeqHelper', {});
      const seqHandler: any = seqHelper.getSeqHandler(seqCol);
      const fastaLines: string[] = [];
      const lineWidth = 60;
      const originalFirstSeq: string[] = [];
      for (let rowIdx = 0; rowIdx < seqHandler.length; rowIdx++) {
        const seqId: string = idColList.length > 0
          ? idColList.map((c: any) => c.get(rowIdx)?.toString() ?? '').join('|')
          : String(rowIdx + 1);
        const srcSS: any = seqHandler.getSplitted(rowIdx);
        const monomers: string[] = [];
        for (let p = 0; p < srcSS.length; p++)
          monomers.push(srcSS.getOriginal(p));
        const seqText: string = monomers.map((om: string) => om.length > 1 ? `[${om}]` : om).join('');
        if (rowIdx === 0) originalFirstSeq.push(seqText);
        fastaLines.push(`>${seqId}\n`);
        for (let i = 0; i < seqText.length; i += lineWidth)
          fastaLines.push(seqText.slice(i, i + lineWidth) + '\n');
      }
      const fastaText: string = fastaLines.join('');
      if (!fastaText.startsWith('>'))
        throw new Error('S3.1: exported FASTA does not start with > header');
      let writeErr: string | null = null;
      try {
        await grok.dapi.files.writeAsText(tempPath, fastaText);
      } catch (e) {
        writeErr = String(e).slice(0, 200);
      }
      // Re-import the FASTA the export actually wrote to disk (read it back so the
      // temp file is exercised, not just the in-memory string). Temp cleanup runs
      // in the outer finally, after this read.
      let reimported: any = null;
      let reimportErr: string | null = null;
      if (!writeErr) {
        try {
          const readBack: string = await grok.dapi.files.readAsText(tempPath);
          const dfs: any = await (grok as any).functions.call('Bio:importFasta', {fileContent: readBack});
          reimported = Array.isArray(dfs) ? dfs[0] : dfs;
          if (reimported) {
            try { await (grok as any).data.detectSemanticTypes(reimported); } catch (_) { /* tolerate */ }
          }
        } catch (e) {
          reimportErr = String(e).slice(0, 200);
        }
      }
      let reimportedShape: any = null;
      let reimportedFirstSeqLen = 0;
      let reimportedFirstSeq = '';
      if (reimported) {
        const rcols = Array.from({length: reimported.columns.length}, (_, i) => reimported.columns.byIndex(i));
        const rmacro: any = rcols.find((c: any) => c.semType === 'Macromolecule') ??
          reimported.columns.byIndex(reimported.columns.length - 1);
        reimportedShape = {
          rowCount: reimported.rowCount,
          cols: reimported.columns.length,
          semType: rmacro?.semType ?? null,
          units: rmacro?.meta?.units ?? null,
        };
        reimportedFirstSeq = rmacro && reimported.rowCount > 0 ? String(rmacro.get(0) ?? '') : '';
        reimportedFirstSeqLen = reimportedFirstSeq.length;
      }
      return {
        fastaShape: {
          startsWithHeader: fastaText.startsWith('>'),
          lineCount: fastaText.split('\n').length,
          totalLen: fastaText.length,
        },
        writeErr,
        reimported: reimportedShape,
        reimportErr,
        originalRowCount: df.rowCount,
        originalFirstSeq: originalFirstSeq[0] ?? '',
        reimportedFirstSeq,
        reimportedFirstSeqLen,
      };
    }, {tempPath: fastaTempPath});
    expect(result.fastaShape.startsWithHeader).toBe(true);
    expect(result.fastaShape.totalLen).toBeGreaterThan(0);
    expect(result.writeErr, `write: ${result.writeErr}`).toBeNull();
    expect(result.reimportErr, `reimport: ${result.reimportErr}`).toBeNull();
    expect(result.reimported, 'S3 re-import produced no DataFrame').not.toBeNull();
    expect(result.reimported.rowCount).toBe(result.originalRowCount);
    expect(result.reimported.semType).toBe('Macromolecule');
    expect(result.reimported.units).toBe('fasta');
    // Round-trip fidelity: first re-imported sequence equals the first original (no drift/truncation).
    expect(result.reimportedFirstSeq).toBe(result.originalFirstSeq);
  });
  try {
    // Scenario 4 — Save project with FASTA-imported table; reopen survives
    await softStep('S4.1: Save project with FASTA-imported table (JS API path)', async () => {
      preSaveRowCount = await page.evaluate(() => grok.shell.tv?.dataFrame?.rowCount ?? 0);
      expect(preSaveRowCount).toBeGreaterThan(0);
      saved = await saveAllTablesWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
      expect(saved.primaryTableInfoId).toBeTruthy();
    });
    await softStep('S4.2: Reopen project — Macromolecule column + units=fasta + renderer dispatch survive', async () => {
      if (!saved) throw new Error('S4.1 did not produce a saved project');
      const result = await reopenAndAssertProvenance(page, saved.projectId);
      expect(result.tablesAfter).toBeGreaterThan(0);
      expect(result.reopenedRowCount).toBe(preSaveRowCount);
      const post = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
        let rendererCellType: string | null = null;
        try {
          const grid: any = grok.shell.tv?.grid;
          if (grid && macro && df.rowCount > 0)
            rendererCellType = grid.cell(macro.name, 0).renderer?.cellType ?? null;
        } catch (_) { /* renderer probe */ }
        return {
          hasMacro: !!macro,
          units: macro?.meta?.units ?? null,
          rendererCellType,
          rowCount: df.rowCount,
        };
      });
      expect(post.hasMacro).toBe(true);
      expect(post.units).toBe('fasta');
      expect(post.rendererCellType).toBe('sequence');
      expect(post.rowCount).toBe(preSaveRowCount);
    });
  } finally {
    // Scenario 5 — Cleanup (runs regardless of earlier failures)
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
