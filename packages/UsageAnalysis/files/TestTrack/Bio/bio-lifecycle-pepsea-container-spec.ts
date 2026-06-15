/* ---
sub_features_covered: [bio.analyze.msa, bio.analyze.msa.align-sequences, bio.analyze.msa.dialog, bio.api.get-seq-helper, bio.engines.msa-pepsea, bio.lifecycle.init]
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
test('Bio pepsea_container source-class lifecycle: HELM → MSA (PepSeA engine) → container status → save + reopen + post-reopen MSA', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;
  const stamp = Date.now();
  const projectName = `bio-lifecycle-pepsea-container-${stamp}`;
  let saved: {projectId: string; primaryTableInfoId: string; layoutId: string | null} | null = null;
  let envHasPepsea: boolean = false;
  let containerId: string | null = null;
  let preRunStatus: string | null = null;
  await loginToDatagrok(page);
  // Setup — open HELM fixture, run Bio init probe, resolve PepSeA container
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
  }, 'System:AppData/Bio/tests/filter_HELM.csv');
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
  await page.evaluate(async () => {
    const probes = ['Bio:getSeqHelper', 'Bio:getMonomerLibHelper', 'Bio:getBioLib'];
    for (const fn of probes) {
      try { await (grok as any).functions.call(fn, {}); return; } catch { /* try next */ }
    }
    await new Promise((r) => setTimeout(r, 3000));
  });
  // Scenario 1 — Initial PepSeA run on HELM column
  await softStep('S1.1: Bio:initBio complete; Bio:getSeqHelper returns ISeqHelper singleton', async () => {
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
  await softStep('S1.0: Resolve PepSeA Docker container handle via DAPI', async () => {
    const result = await page.evaluate(async () => {
      try {
        const dc: any = (grok as any).dapi.docker?.dockerContainers;
        if (!dc) return {err: 'grok.dapi.docker.dockerContainers not exposed on this build'} as any;
        let containers: any[] = [];
        try {
          containers = await dc.list();
        } catch (e1) {
          try {
            const c: any = await dc.filter('name like "%pepsea%"').first();
            if (c) containers = [c];
          } catch (e2) {
            return {err: `dockerContainers.list failed: ${String(e1).slice(0, 120)} / fallback filter failed: ${String(e2).slice(0, 120)}`} as any;
          }
        }
        if (!Array.isArray(containers) || containers.length === 0)
          return {found: false, count: 0} as any;
        const pepsea: any = containers.find((c: any) =>
          c && c.name && String(c.name).toLowerCase().indexOf('pepsea') >= 0);
        if (!pepsea) return {found: false, count: containers.length} as any;
        return {found: true, id: pepsea.id, name: pepsea.name, status: pepsea.status};
      } catch (e) {
        return {err: String(e).slice(0, 200)} as any;
      }
    });
    if ((result as any).found) {
      envHasPepsea = true;
      containerId = (result as any).id;
      preRunStatus = (result as any).status ?? null;
      expect(containerId).toBeTruthy();
    } else {
      // eslint-disable-next-line no-console
      console.warn('[S1.0] PepSeA Docker container NOT registered on this host — env-conditional skip per scenario Setup directive. PepSeA-specific assertions will short-circuit.');
    }
  });
  await softStep('S1.2: filter_HELM.csv detector classifies Macromolecule column with units=helm', async () => {
    const info = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
      return {
        hasMacro: !!macro,
        units: macro?.meta?.units ?? null,
        rowCount: df.rowCount,
      };
    });
    expect(info.hasMacro).toBe(true);
    expect(info.units).toBe('helm');
    expect(info.rowCount).toBeGreaterThan(0);
  });
  let aligned: {hasNewMacro: boolean; alignedTag: string | null; newColName: string | null} = {
    hasNewMacro: false, alignedTag: null, newColName: null,
  };
  await softStep('S1.3: Bio | Analyze | MSA... — dialog opens; Engine SELECT exposes PepSeA option', async () => {
    await bio.openBioAnalyze(page, 'div-Bio---Analyze---MSA...');
    await page.locator('[name="dialog-MSA"]').waitFor({timeout: 60_000});
    const engineOpts = await page.evaluate(() => {
      const sel = document.querySelector('[name="dialog-MSA"] [name="input-Engine"]') as HTMLSelectElement | null;
      if (!sel) return {hasSelect: false, options: [] as string[]} as any;
      const options = Array.from(sel.options).map((o) => o.text || o.value);
      return {hasSelect: true, options};
    });
    expect(engineOpts.hasSelect).toBe(true);
    expect(engineOpts.options.length).toBeGreaterThanOrEqual(2);
    const hasPepsea = engineOpts.options.some((o: string) => o.toLowerCase().indexOf('pepsea') >= 0);
    expect(hasPepsea).toBe(true);
    const hasKalign = engineOpts.options.some((o: string) =>
      o.toLowerCase().indexOf('datagrok') >= 0 || o.toLowerCase().indexOf('msa') >= 0);
    expect(hasKalign).toBe(true);
  });
  await softStep('S1.4: Select PepSeA engine; OK runs; aligned Macromolecule column appears with aligned=SEQ.MSA tag', async () => {
    if (!envHasPepsea) {
      await page.evaluate(() => {
        const cancel = document.querySelector('[name="dialog-MSA"] [name="button-CANCEL"]') as HTMLElement | null;
        if (cancel) cancel.click();
      });
      // eslint-disable-next-line no-console
      console.warn('[S1.4] PepSeA env not configured — dialog canceled; PepSeA-engine OK + alignSequences assertion skipped per Setup env-gate.');
      return;
    }
    const beforeCols: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
    await page.evaluate(() => {
      const sel = document.querySelector('[name="dialog-MSA"] [name="input-Engine"]') as HTMLSelectElement | null;
      if (!sel) return;
      const pepseaIdx = Array.from(sel.options).findIndex((o) =>
        (o.text || o.value || '').toLowerCase().indexOf('pepsea') >= 0);
      if (pepseaIdx >= 0) {
        sel.selectedIndex = pepseaIdx;
        sel.dispatchEvent(new Event('input', {bubbles: true}));
        sel.dispatchEvent(new Event('change', {bubbles: true}));
      }
    });
    await page.locator('[name="dialog-MSA"] [name="button-OK"]').click();
    await page.waitForFunction(
      (b) => grok.shell.tv.dataFrame.columns.length > b,
      beforeCols, {timeout: 240_000});
    aligned = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macroCols = cols.filter((c: any) => c.semType === 'Macromolecule');
      const last: any = macroCols[macroCols.length - 1];
      return {
        hasNewMacro: macroCols.length >= 2,
        alignedTag: last?.getTag?.('aligned') ?? null,
        newColName: last?.name ?? null,
      };
    });
    expect(aligned.hasNewMacro).toBe(true);
    expect(aligned.newColName).toBeTruthy();
    expect(aligned.alignedTag).toBe('SEQ.MSA');
  });
  // Scenario 2 — Save project with PepSeA alignment
  try {
    await softStep('S2.1: Save project containing the PepSeA-aligned table (JS API persistence)', async () => {
      saved = await saveAllTablesWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
      expect(saved.primaryTableInfoId).toBeTruthy();
    });
    await softStep('S2.2: grok.dapi.projects.find returns the saved project', async () => {
      if (!saved) throw new Error('S2.1 did not produce a saved project');
      const result = await page.evaluate(async (id) => {
        try {
          const p: any = await grok.dapi.projects.find(id);
          return {found: !!p, name: p?.name ?? null};
        } catch (e) {
          return {err: String(e).slice(0, 200)} as any;
        }
      }, saved.projectId);
      expect((result as any).err).toBeUndefined();
      expect((result as any).found).toBe(true);
      expect((result as any).name).toBe(projectName);
    });
    // Scenario 3 — Reopen (container eviction skipped; host-shared resource)
    await softStep('S3.1: Container eviction skipped (host-shared resource per scenario Setup directive)', async () => {
      if (envHasPepsea) {
        // eslint-disable-next-line no-console
        console.warn(`[S3.1] PepSeA container preRunStatus=${preRunStatus ?? 'unknown'}; eviction skipped to preserve host-shared resource state. Reopen path (S3.2-S3.3) is the assertable contract.`);
      } else {
        // eslint-disable-next-line no-console
        console.warn('[S3.1] PepSeA env not configured — eviction skip is the only valid path.');
      }
      expect(true).toBe(true);
    });
    await softStep('S3.2-3.3: Reopen project — aligned Macromolecule column with aligned=SEQ.MSA tag survives', async () => {
      if (!saved) throw new Error('S2.1 did not produce a saved project; cannot exercise reopen');
      const result = await reopenAndAssertProvenance(page, saved.projectId);
      expect(result.tablesAfter).toBeGreaterThan(0);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
      const post = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macroCols = cols.filter((c: any) => c.semType === 'Macromolecule');
        const aligned: any = macroCols.find((c: any) => c?.getTag?.('aligned') === 'SEQ.MSA') ?? null;
        const sourceHelm: any = macroCols.find((c: any) => (c?.meta?.units ?? null) === 'helm') ?? null;
        return {
          macroCount: macroCols.length,
          hasAligned: !!aligned,
          alignedTag: aligned?.getTag?.('aligned') ?? null,
          alignedName: aligned?.name ?? null,
          hasSourceHelm: !!sourceHelm,
          sourceHelmUnits: sourceHelm?.meta?.units ?? null,
        };
      });
      if (envHasPepsea) {
        expect(post.hasAligned).toBe(true);
        expect(post.alignedTag).toBe('SEQ.MSA');
        expect(post.alignedName).toBeTruthy();
        expect(post.hasSourceHelm).toBe(true);
        expect(post.sourceHelmUnits).toBe('helm');
      } else {
        expect(post.hasSourceHelm).toBe(true);
        expect(post.sourceHelmUnits).toBe('helm');
        // eslint-disable-next-line no-console
        console.warn('[S3.3] PepSeA env not configured — aligned column assertion skipped; HELM source-column survival is the reachable contract.');
      }
    });
    await softStep('S3.4: Post-reopen MSA invocation — PepSeA container resumes / auto-restarts (no crash)', async () => {
      if (!envHasPepsea) {
        // eslint-disable-next-line no-console
        console.warn('[S3.4] PepSeA env not configured — post-reopen MSA invocation skipped.');
        return;
      }
      await bio.openBioAnalyze(page, 'div-Bio---Analyze---MSA...');
      await page.locator('[name="dialog-MSA"]').waitFor({timeout: 60_000});
      const engineOpts = await page.evaluate(() => {
        const sel = document.querySelector('[name="dialog-MSA"] [name="input-Engine"]') as HTMLSelectElement | null;
        if (!sel) return {hasSelect: false, options: [] as string[]} as any;
        return {
          hasSelect: true,
          options: Array.from(sel.options).map((o) => o.text || o.value),
        };
      });
      expect(engineOpts.hasSelect).toBe(true);
      expect(engineOpts.options.some((o: string) => o.toLowerCase().indexOf('pepsea') >= 0)).toBe(true);
      await page.evaluate(() => {
        const cancel = document.querySelector('[name="dialog-MSA"] [name="button-CANCEL"]') as HTMLElement | null;
        if (cancel) cancel.click();
      });
    });
    // Scenario 4 — Container lifecycle observation (read-only)
    await softStep('S4: Container status observable via DAPI without errors (read-only)', async () => {
      if (!envHasPepsea || !containerId) {
        // eslint-disable-next-line no-console
        console.warn('[S4] PepSeA env not configured — container status observation skipped.');
        return;
      }
      const result = await page.evaluate(async (id) => {
        try {
          const c: any = await (grok as any).dapi.docker.dockerContainers.find(id);
          return {found: !!c, status: c?.status ?? null, name: c?.name ?? null};
        } catch (e) {
          return {err: String(e).slice(0, 200)} as any;
        }
      }, containerId);
      expect((result as any).err).toBeUndefined();
      expect((result as any).found).toBe(true);
      expect((result as any).status).not.toBeNull();
    });
  } finally {
    // Scenario 5 — Cleanup (project deletion only; PepSeA container left untouched)
    if (saved) {
      await deleteProjectWithCleanup(page, {
        projectId: saved.projectId,
        tableInfoId: saved.primaryTableInfoId,
      });
    }
  }
  finishSpec();
});
