import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import {
  saveAllTablesWithProvenance,
  reopenAndAssertProvenance,
  deleteProjectWithCleanup,
} from '../helpers/projects';
test.use(specTestOptions);
test('Bio immunum_wasm source-class lifecycle: init → IMGT numbering → save+reopen → re-run', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;
  const stamp = Date.now();
  const projectName = `bio-lifecycle-immunum-wasm-${stamp}`;
  const antibodyFixturePath = 'System:AppData/Bio/samples/antibodies.csv';
  let saved: {projectId: string; primaryTableInfoId: string; layoutId: string | null} | null = null;
  const EXPECTED_IMMUNUM_COLS = [
    'position_names', 'chain_type', 'annotations_json',
    'numbering_detail', 'numbering_map',
  ] as const;
  await loginToDatagrok(page);
  // Setup — open antibody fixture, await semType + Bio init readiness.
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
      await new Promise((r) => setTimeout(r, 5000));
    }
  }, antibodyFixturePath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
  await page.evaluate(async () => {
    const probes = ['Bio:getSeqHelper', 'Bio:getMonomerLibHelper', 'Bio:getBioLib'];
    for (const fn of probes) {
      try { await (grok as any).functions.call(fn, {}); return; } catch { /* try next */ }
    }
    await new Promise((r) => setTimeout(r, 3000));
  });
  // Scenario 1 — Trigger initBio + verify getSeqHelper resolves
  await softStep('S1.1: initBio completes; Bio:getSeqHelper resolves to a usable singleton', async () => {
    const info = await page.evaluate(async () => {
      let resolved = false;
      let helperKind: string | null = null;
      let methodNames: string[] = [];
      try {
        const helper: any = await (grok as any).functions.call('Bio:getSeqHelper', {});
        resolved = !!helper;
        helperKind = helper ? (helper.constructor?.name ?? typeof helper) : null;
        if (helper) {
          methodNames = ['getSeqHandler', 'getSeqMonomers', 'helmToAtomicLevel',
            'setUnitsToFastaColumn']
            .filter((m) => typeof helper[m] === 'function');
        }
      } catch (e) {
        return {resolved: false, helperKind: null, methodNames: [], err: String(e).slice(0, 200)};
      }
      return {resolved, helperKind, methodNames, err: null};
    });
    expect(info.resolved).toBe(true);
    expect(info.methodNames).toContain('getSeqHandler');
  });
  await softStep('S1.2: antibodies.csv opens with at least one Macromolecule column (FASTA notation)', async () => {
    const info = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macroCols = cols.filter((c: any) => c.semType === 'Macromolecule');
      const macro: any = macroCols[0];
      return {
        macroCount: macroCols.length,
        macroName: macro?.name ?? null,
        units: macro?.meta?.units ?? null,
        rowCount: df.rowCount,
      };
    });
    expect(info.macroCount).toBeGreaterThan(0);
    expect(info.rowCount).toBeGreaterThan(0);
    expect(info.units).not.toBeNull();
  });
  await softStep('S1.3-1.4: Apply Antibody Numbering (Immunum/IMGT) via top-menu; sibling column appended', async () => {
    const baseCols: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 400));
      document.querySelector('[name="div-Bio---Annotate"]')!
        .dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 400));
      (document.querySelector(
        '[name="div-Bio---Annotate---Apply-Numbering-Scheme..."]') as HTMLElement).click();
    });
    // Dialog title is "Apply Antibody Numbering" (not the menu label).
    await page.locator('[name="dialog-Apply-Antibody-Numbering"]').waitFor({timeout: 60_000});
    await page.locator('[name="dialog-Apply-Antibody-Numbering"] [name="button-OK"]').click();
    await page.waitForFunction(
      (b) => grok.shell.tv.dataFrame.columns.length > b ||
        (() => {
          const df = grok.shell.tv.dataFrame;
          const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
          return cols.some((c: any) => {
            try {
              const tag = c.getTag?.('.annotations') ?? c.tags?.get?.('.annotations');
              return typeof tag === 'string' && tag.length > 0 && tag.startsWith('[');
            } catch { return false; }
          });
        })(),
      baseCols, {timeout: 120_000});
    await page.waitForFunction(
      () => document.querySelectorAll('[name="dialog-Apply-Antibody-Numbering"]').length === 0,
      null, {timeout: 30_000}).catch(() => {});
  });
  await softStep('S1.4: Immunum engine returns 5-column result DataFrame with documented column names + non-null per-row values', async () => {
    const info = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macroCol: any = cols.find((c: any) => c.semType === 'Macromolecule');
      if (!macroCol) throw new Error('S1.4: no Macromolecule column for Immunum engine call');
      let result: any = null;
      let invokeErr: string | null = null;
      try {
        result = await (grok as any).functions.call(
          'Bio:immunumAntibodyNumbering',
          {df, seqCol: macroCol, scheme: 'imgt'},
        );
      } catch (e) {
        invokeErr = String(e).slice(0, 300);
      }
      if (!result) return {result: null, invokeErr, colNames: null, sampleRow: null};
      const colNames = Array.from({length: result.columns.length},
        (_, i) => result.columns.byIndex(i).name);
      const sampleRow: Record<string, string | null> = {};
      for (let i = 0; i < result.columns.length; i++) {
        const col = result.columns.byIndex(i);
        const v = col.get(0);
        sampleRow[col.name] = v == null ? null : String(v).slice(0, 80);
      }
      return {
        result: {rowCount: result.rowCount, colCount: result.columns.length},
        invokeErr,
        colNames,
        sampleRow,
      };
    });
    expect(info.invokeErr).toBeNull();
    expect(info.result).not.toBeNull();
    expect(info.result!.colCount).toBe(5);
    expect(info.result!.rowCount).toBeGreaterThan(0);
    for (const expectedName of EXPECTED_IMMUNUM_COLS)
      expect(info.colNames).toContain(expectedName);
    for (const expectedName of EXPECTED_IMMUNUM_COLS)
      expect(info.sampleRow![expectedName]).not.toBeNull();
  });
  try {
    // Scenario 2 — Save project with the numbering output
    await softStep('S2.1: Save project with numbering output (JS API path)', async () => {
      saved = await saveAllTablesWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
      expect(saved.primaryTableInfoId).toBeTruthy();
    });
    await softStep('S2.2: saved project is findable server-side via dapi.projects.find(id)', async () => {
      if (!saved) throw new Error('S2.1 did not produce a saved project');
      const ok = await page.evaluate(async (id) => {
        const proj = await (grok as any).dapi.projects.find(id);
        return proj != null && proj.id === id;
      }, saved.projectId);
      expect(ok).toBe(true);
    });
    // Scenario 3 — Reopen project + WASM re-load + deterministic re-run.
    await softStep('S3.1-3.2: reopen project — antibody table + Macromolecule semType survive', async () => {
      if (!saved) throw new Error('S2.1 did not produce a saved project');
      const result = await reopenAndAssertProvenance(page, saved.projectId);
      expect(result.tablesAfter).toBeGreaterThan(0);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
      const post = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
        return {
          hasMacro: !!macro,
          units: macro?.meta?.units ?? null,
          rowCount: df.rowCount,
        };
      });
      expect(post.hasMacro).toBe(true);
      expect(post.units).not.toBeNull();
      expect(post.rowCount).toBeGreaterThan(0);
    });
    await softStep('S3.3-3.4: re-run Immunum on reopened table; result shape deterministic on WASM re-load', async () => {
      const info = await page.evaluate(async () => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macroCol: any = cols.find((c: any) => c.semType === 'Macromolecule');
        if (!macroCol) throw new Error('S3.3: no Macromolecule column after reopen');
        let result: any = null;
        let invokeErr: string | null = null;
        try {
          result = await (grok as any).functions.call(
            'Bio:immunumAntibodyNumbering',
            {df, seqCol: macroCol, scheme: 'imgt'},
          );
        } catch (e) {
          invokeErr = String(e).slice(0, 300);
        }
        if (!result) return {result: null, invokeErr, colNames: null, sampleRow: null};
        const colNames = Array.from({length: result.columns.length},
          (_, i) => result.columns.byIndex(i).name);
        const sampleRow: Record<string, string | null> = {};
        for (let i = 0; i < result.columns.length; i++) {
          const col = result.columns.byIndex(i);
          const v = col.get(0);
          sampleRow[col.name] = v == null ? null : String(v).slice(0, 80);
        }
        return {
          result: {rowCount: result.rowCount, colCount: result.columns.length},
          invokeErr,
          colNames,
          sampleRow,
        };
      });
      expect(info.invokeErr).toBeNull();
      expect(info.result).not.toBeNull();
      expect(info.result!.colCount).toBe(5);
      expect(info.result!.rowCount).toBeGreaterThan(0);
      for (const expectedName of EXPECTED_IMMUNUM_COLS)
        expect(info.colNames).toContain(expectedName);
      for (const expectedName of EXPECTED_IMMUNUM_COLS)
        expect(info.sampleRow![expectedName]).not.toBeNull();
    });
  } finally {
    // Scenario 4 — Cleanup (runs regardless of earlier failures)
    if (saved) {
      await deleteProjectWithCleanup(page, {
        projectId: saved.projectId,
        tableInfoId: saved.primaryTableInfoId,
      });
    }
  }
  finishSpec();
});
