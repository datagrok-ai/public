// GROK-17967: Mol* (Biostructure) and NGL viewers must agree on multi-ligand state. Asserted via the
// JS property contract (same wired ligandColumnName + showSelectedRowsLigands + shared selection count
// on both engines) rather than a canvas pixel-diff. Fixture: dock.csv (multi-row 'ligand' column).
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

declare const grok: any;
declare const DG: any;

const sampleDockCsv = 'System:AppData/BiostructureViewer/samples/dock.csv';

test('BiostructureViewer — GROK-17967 multi-ligand NGL/Biostructure parity regression guard', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  // Baseline environment setup.
  await page.evaluate(() => {
    document.querySelectorAll('.d4-dialog').forEach((d) => {
      const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
      if (cancel) cancel.click();
    });
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
  });
  await page.locator('[name="Browse"]').waitFor({timeout: 30_000});

  try {
    // SCENARIO 1 — Multi-ligand fixture in the Biostructure (Mol*) viewer.
    let bioMountedDiag: any = null;
    let bioPropsDiag: any = null;
    let bioSelectionCount = 0;

    await softStep('Scenario 1 step 1 — Open multi-ligand dock.csv; Biostructure viewer mounts', async () => {
      const result = await page.evaluate(async (path) => {
        grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1500));
        const df = await grok.dapi.files.readCsv(path);
        df.name = 'biostructure-bug-grok-17967-multiligand';
        const tv = grok.shell.addTableView(df);
        // Wait for semantic-type detection (best-effort).
        await new Promise((resolve) => {
          const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
          setTimeout(resolve, 5000);
        });
        await new Promise((r) => setTimeout(r, 1500));
        const vBio = tv.addViewer('Biostructure');
        await new Promise((r) => setTimeout(r, 3000));
        // Tag the ligand column Molecule3D so the picker accepts it as a wiring target.
        const ligandCol = df.col('ligand');
        if (ligandCol) ligandCol.semType = 'Molecule3D';
        await new Promise((r) => setTimeout(r, 800));
        return {
          rowCount: df.rowCount,
          hasLigandCol: !!ligandCol,
          ligandColType: ligandCol?.type,
          ligandSemType: ligandCol?.semType,
          hasBioContainer: !!document.querySelector('[name="viewer-Biostructure"]'),
          bioType: vBio?.type,
          defaultRep: vBio?.props?.get?.('representation') ?? null,
          viewerTypesAfter: tv.viewers ? Array.from(tv.viewers).map((v: any) => v.type) : [],
        };
      }, sampleDockCsv);

      await expect(page.locator('[name="viewer-Biostructure"]')).toBeVisible({timeout: 30_000});

      expect(result.hasBioContainer).toBe(true);
      expect(result.bioType).toBe('Biostructure');
      expect(result.defaultRep).toBe('cartoon');
      expect(result.hasLigandCol).toBe(true);
      expect(result.ligandSemType).toBe('Molecule3D');
      expect(result.rowCount).toBeGreaterThan(1);
      expect(result.viewerTypesAfter).toContain('Biostructure');
      bioMountedDiag = result;
    });

    await softStep('Scenario 1 steps 2-3 — Wire ligandColumnName + showSelectedRowsLigands on Biostructure', async () => {
      const result = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        if (!tv) return {ok: false, reason: 'no tv'};
        const vBio = tv.viewers ? Array.from(tv.viewers).find((v: any) => v.type === 'Biostructure') as any : null;
        if (!vBio) return {ok: false, reason: 'no Biostructure viewer'};
        vBio.setOptions({
          ligandColumnName: 'ligand',
          showSelectedRowsLigands: true,
        });
        await new Promise((r) => setTimeout(r, 1500));
        return {
          ok: true,
          ligandColAfter: vBio.props.get('ligandColumnName'),
          showSelectedAfter: vBio.props.get('showSelectedRowsLigands'),
          showCurrentDefault: vBio.props.get('showCurrentRowLigand'),
        };
      });
      expect(result.ok).toBe(true);
      // Mol* side: the wired ligandColumnName must round-trip.
      expect(result.ligandColAfter).toBe('ligand');
      expect(result.showSelectedAfter).toBe(true);
      expect(result.showCurrentDefault).toBe(true);
      bioPropsDiag = result;
    });

    await softStep('Scenario 1 step 3 — Select multiple ligand-bearing rows on shared dataframe', async () => {
      const result = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        if (!tv) return {ok: false, reason: 'no tv'};
        const df = tv.dataFrame;
        if (!df) return {ok: false, reason: 'no dataframe'};
        // Select the first 3 rows (multi-ligand selection).
        df.selection.init((i: number) => i < 3);
        await new Promise((r) => setTimeout(r, 1500));
        return {
          ok: true,
          selectedCount: df.selection.trueCount,
          rowCount: df.rowCount,
        };
      });
      expect(result.ok).toBe(true);
      expect(result.selectedCount).toBeGreaterThanOrEqual(2);
      expect(result.selectedCount).toBeLessThanOrEqual(result.rowCount);
      bioSelectionCount = result.selectedCount;
    });

    // SCENARIO 2 — Same table, add NGL viewer, mirror the property contract; assert parity.
    let nglMountedDiag: any = null;
    let nglPropsDiag: any = null;

    await softStep('Scenario 2 step 1 — Add NGL viewer to the same table; NGL mounts', async () => {
      const result = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        if (!tv) return {ok: false, reason: 'no tv'};
        const vNgl = tv.addViewer('NGL');
        await new Promise((r) => setTimeout(r, 3000));
        return {
          ok: true,
          hasNglContainer: !!document.querySelector('[name="viewer-NGL"]'),
          nglType: vNgl?.type,
          viewerTypesAfter: tv.viewers ? Array.from(tv.viewers).map((v: any) => v.type) : [],
        };
      });
      await expect(page.locator('[name="viewer-NGL"]')).toBeVisible({timeout: 30_000});
      expect(result.ok).toBe(true);
      expect(result.hasNglContainer).toBe(true);
      expect(result.nglType).toBe('NGL');
      expect(result.viewerTypesAfter).toContain('NGL');
      expect(result.viewerTypesAfter).toContain('Biostructure');
      nglMountedDiag = result;
    });

    await softStep('Scenario 2 steps 2-3 — Wire ligandColumnName + showSelectedRowsLigands on NGL', async () => {
      const result = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        if (!tv) return {ok: false, reason: 'no tv'};
        const vNgl = tv.viewers ? Array.from(tv.viewers).find((v: any) => v.type === 'NGL') as any : null;
        if (!vNgl) return {ok: false, reason: 'no NGL viewer'};
        // Wire the same property contract as Scenario 1 (parity input).
        vNgl.setOptions({
          ligandColumnName: 'ligand',
          showSelectedRowsLigands: true,
        });
        await new Promise((r) => setTimeout(r, 1500));
        return {
          ok: true,
          ligandColAfter: vNgl.props.get('ligandColumnName'),
          showSelectedAfter: vNgl.props.get('showSelectedRowsLigands'),
          showCurrentDefault: vNgl.props.get('showCurrentRowLigand'),
        };
      });
      expect(result.ok).toBe(true);
      // NGL side: the wired ligandColumnName must round-trip the same as Biostructure.
      expect(result.ligandColAfter).toBe('ligand');
      expect(result.showSelectedAfter).toBe(true);
      expect(result.showCurrentDefault).toBe(true);
      nglPropsDiag = result;
    });

    await softStep('Scenario 2 step 5 — Parity assertion (load-bearing GROK-17967 invariant)', async () => {
      // Parity on the JS-API surface: both viewers report the same ligandColumnName,
      // same showSelectedRowsLigands, and read the same shared-selection count.
      const result = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        if (!tv) return {ok: false, reason: 'no tv'};
        const df = tv.dataFrame;
        const vBio = tv.viewers ? Array.from(tv.viewers).find((v: any) => v.type === 'Biostructure') as any : null;
        const vNgl = tv.viewers ? Array.from(tv.viewers).find((v: any) => v.type === 'NGL') as any : null;
        return {
          ok: !!(vBio && vNgl && df),
          bioWired: vBio?.props?.get?.('ligandColumnName'),
          nglWired: vNgl?.props?.get?.('ligandColumnName'),
          bioShowSel: vBio?.props?.get?.('showSelectedRowsLigands'),
          nglShowSel: vNgl?.props?.get?.('showSelectedRowsLigands'),
          bioShowCur: vBio?.props?.get?.('showCurrentRowLigand'),
          nglShowCur: vNgl?.props?.get?.('showCurrentRowLigand'),
          sharedSelectionCount: df?.selection?.trueCount,
        };
      });
      expect(
        result.ok,
        `Parity precondition failed: both viewers must be present. result=${JSON.stringify(result)}`,
      ).toBe(true);

      // Mol* and NGL must agree on the wired ligandColumnName.
      expect(
        result.bioWired,
        `GROK-17967 parity violated (ligandColumnName divergence): ` +
        `Mol*=${result.bioWired}, NGL=${result.nglWired}. See ` +
        `bug-library/biostructureviewer.yaml#GROK-17967.`,
      ).toBe(result.nglWired);
      expect(result.bioWired).toBe('ligand');

      // The two engines must agree on showSelectedRowsLigands.
      expect(
        result.bioShowSel,
        `GROK-17967 parity violated (showSelectedRowsLigands divergence): ` +
        `Mol*=${result.bioShowSel}, NGL=${result.nglShowSel}.`,
      ).toBe(result.nglShowSel);
      expect(result.bioShowSel).toBe(true);

      // Both engines must report the same showCurrentRowLigand default.
      expect(result.bioShowCur).toBe(result.nglShowCur);

      // Shared selected-row count drives the per-row overlay on both engines.
      expect(result.sharedSelectionCount).toBe(bioSelectionCount);
      expect(result.sharedSelectionCount).toBeGreaterThanOrEqual(2);
    });

    await softStep('Scenario 2 step 6 — Teardown (best-effort)', async () => {
      await page.evaluate(() => {
        try { grok.shell.closeAll(); } catch (e) { /* best-effort */ }
      });
      // Reference diag values to avoid unused-var lint flags.
      if (!bioMountedDiag || !bioPropsDiag || !nglMountedDiag || !nglPropsDiag)
        return;
    });
  } finally {
    // Best-effort cleanup. No server-side state created by this scenario.
    try {
      await page.evaluate(() => { (window as any).grok?.shell?.closeAll?.(); });
    } catch (e) { /* best-effort */ }
  }

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
