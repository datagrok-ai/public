/* ---
sub_features_covered: [bio.detector, bio.rendering, bio.search.diversity, bio.search.diversity.top-menu, bio.viewers.diversity-search]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
test.use(specTestOptions);
const FASTA_DATASET_PATH = 'System:AppData/Bio/tests/filter_FASTA.csv';
const HELM_DATASET_PATH = 'System:AppData/Bio/tests/filter_HELM.csv';
test('Bio Diversity Search docks viewer + re-runs diversity selection on fresh dataset', async ({page}) => {
  test.setTimeout(240_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    grok.shell.addTableView(df);
  }, FASTA_DATASET_PATH);
  await page.waitForFunction(() => {
    const df = grok.shell.tv?.dataFrame;
    if (!df) return false;
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    return cols.some((c: any) => c.semType === 'Macromolecule');
  }, null, {timeout: 30_000});
  await page.locator('[name="viewer-Grid"] canvas').first().waitFor({state: 'attached', timeout: 30_000});
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
  await page.waitForFunction(async () => {
    const probes = ['Bio:getSeqHelper', 'Bio:getMonomerLibHelper', 'Bio:getBioLib'];
    for (const fn of probes) {
      try { await (grok as any).functions.call(fn, {}); return true; } catch {  }
    }
    return false;
  }, null, {timeout: 30_000});
  await page.evaluate(async () => {
    const candidates = ['Bio:diversitySearchTopMenu', 'Bio:diversitySearch'];
    const findAny = (names: string[]): boolean => {
      for (const n of names) {
        try {
          if ((grok as any).functions.find && (grok as any).functions.find(n)) return true;
        } catch {  }
      }
      return false;
    };
    const deadline = Date.now() + 15_000;
    while (Date.now() < deadline) {
      if (findAny(candidates)) return;
      await new Promise((r) => setTimeout(r, 300));
    }
  });
  const setupProbe = await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    const macro = cols.find((c: any) => c.semType === 'Macromolecule') as any;
    return {
      hasMacromoleculeCol: !!macro,
      macroName: macro?.name ?? null,
      macroSemType: macro?.semType ?? null,
      macroUnits: macro?.getTag?.('units') ?? null,
      rowCount: df.rowCount,
    };
  });
  expect(setupProbe.hasMacromoleculeCol,
    'atlas bio.detector contract: readCsv path MUST classify a Macromolecule column synchronously').toBe(true);
  expect(setupProbe.macroSemType).toBe('Macromolecule');
  expect(setupProbe.macroUnits,
    'atlas bio.rendering: filter_FASTA.csv MUST classify with units=fasta').toBe('fasta');
  expect(setupProbe.rowCount,
    'scenario .md Setup: ≥ 5 rows so the diversity subset is meaningful').toBeGreaterThanOrEqual(5);
  await softStep('Scenario 1.1: click Bio > Search > Diversity Search (no dialog — viewer docks directly)', async () => {
    await page.evaluate(() => (document.querySelector('[name="div-Bio"]') as HTMLElement).click());
    await page.locator('[name="div-Bio---Search"]').waitFor({state: 'attached', timeout: 15_000});
    await page.evaluate(() => {
      const group = document.querySelector('[name="div-Bio---Search"]')!;
      group.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      group.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
    });
    await page.locator('[name="div-Bio---Search---Diversity-Search"]').waitFor({state: 'attached', timeout: 15_000});
    await page.evaluate(() => (document.querySelector('[name="div-Bio---Search---Diversity-Search"]') as HTMLElement).click());
  });
  await softStep('Scenario 1.2-3: viewer Sequence Diversity Search docks; diverse subset populated', async () => {
    await page.locator('[name="viewer-Sequence-Diversity-Search"]').waitFor({timeout: 60_000});
    await page.waitForFunction(() => {
      const viewers = Array.from((grok.shell.tv as any).viewers) as any[];
      const v = viewers.find((vw) => vw.type === 'Sequence Diversity Search');
      if (!v) return false;
      const ids = (v as any).renderMolIds;
      return Array.isArray(ids) && ids.length > 0;
    }, null, {timeout: 180_000});
    const dockProbe = await page.evaluate(() => {
      const tv = grok.shell.tv as any;
      const viewers = Array.from(tv.viewers) as any[];
      const v = viewers.find((vw) => vw.type === 'Sequence Diversity Search');
      const viewerTypes = viewers.map((vw) => vw.type);
      const renderMolIds: number[] | null = (v as any)?.renderMolIds ?? null;
      const subsetSize = Array.isArray(renderMolIds) ? renderMolIds.length : null;
      const subsetIndices: number[] | null = Array.isArray(renderMolIds) ? [...renderMolIds] : null;
      const df = tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macro = cols.find((c: any) => c.semType === 'Macromolecule') as any;
      const subsetSequences: (string | null)[] | null =
        subsetIndices ? subsetIndices.map((i) => macro?.get?.(i) ?? null) : null;
      return {
        viewerPresent: !!v,
        viewerType: v?.type ?? null,
        viewerTypes,
        subsetSize,
        subsetIndices,
        subsetSequences,
        dfId: df.id ?? null,
        dfName: df.name ?? null,
        viewerDockedDOM: !!document.querySelector('[name="viewer-Sequence-Diversity-Search"]'),
      };
    });
    expect(dockProbe.viewerPresent,
      `Expected 'Sequence Diversity Search' viewer in tv.viewers; actual types: ${JSON.stringify(dockProbe.viewerTypes)}`).toBe(true);
    expect(dockProbe.viewerType).toBe('Sequence Diversity Search');
    expect(dockProbe.viewerDockedDOM,
      'viewer dock container [name="viewer-Sequence-Diversity-Search"] MUST be present in DOM').toBe(true);
    expect(dockProbe.subsetSize,
      'Viewer diversity compute MUST populate renderMolIds with >=2 entries (non-empty diverse subset).' +
      ` Observed: ${dockProbe.subsetSize}`).not.toBeNull();
    expect(dockProbe.subsetSize!).toBeGreaterThanOrEqual(2);
    expect(dockProbe.subsetSequences,
      'renderMolIds MUST resolve to macromolecule sequence values').not.toBeNull();
    const seqs1 = dockProbe.subsetSequences!.map((s) => String(s ?? ''));
    // filter_FASTA.csv legitimately contains empty macromolecule cells, and the diversity
    // selection may pick them; assert on the real (non-empty) sequences it resolves.
    const nonEmptySeqs1 = seqs1.filter((s) => s.length > 0);
    expect(nonEmptySeqs1.length,
      'diverse subset MUST resolve to at least 2 real (non-empty) macromolecule sequences').toBeGreaterThanOrEqual(2);
    expect(new Set(nonEmptySeqs1).size,
      'diverse subset MUST contain >=2 distinct sequences (not a trivial repeat or empty-content result)')
      .toBeGreaterThanOrEqual(2);
    const errorBalloonCount = await page.locator('.d4-balloon.error, .grok-balloon-error').count();
    expect(errorBalloonCount,
      'Scenario 1 Expected: "No error balloon appears" on positive path').toBe(0);
    await page.evaluate((d) => {
      (window as any).__scenario1SubsetIndices = d.indices;
      (window as any).__scenario1RowCount = grok.shell.tv.dataFrame.rowCount;
      (window as any).__scenario1DfId = d.dfId;
      (window as any).__scenario1DfName = d.dfName;
    }, {indices: dockProbe.subsetIndices, dfId: dockProbe.dfId, dfName: dockProbe.dfName});
  });
  // Scenario 2 — Reopening with a fresh dataset re-runs diversity selection.
  await softStep('Scenario 2.1: close the Diversity Search viewer', async () => {
    await page.evaluate(() => {
      const viewers = Array.from((grok.shell.tv as any).viewers) as any[];
      const v = viewers.find((vw) => vw.type === 'Sequence Diversity Search');
      if (!v) throw new Error('Scenario 2.1: expected docked Sequence Diversity Search viewer from Scenario 1');
      v.close();
    });
    await page.waitForFunction(() => {
      const viewers = Array.from((grok.shell.tv as any).viewers) as any[];
      const stillThere = viewers.some((vw) => vw.type === 'Sequence Diversity Search');
      const stillInDom = !!document.querySelector('[name="viewer-Sequence-Diversity-Search"]');
      return !stillThere && !stillInDom;
    }, null, {timeout: 30_000});
  });
  await softStep('Scenario 2.2: open filter_HELM.csv; verify HELM Macromolecule detector classifies synchronously', async () => {
    await page.evaluate(async (path) => {
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv(path);
      grok.shell.addTableView(df);
    }, HELM_DATASET_PATH);
    await page.waitForFunction(() => {
      const df = grok.shell.tv?.dataFrame;
      if (!df) return false;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      return cols.some((c: any) => c.semType === 'Macromolecule');
    }, null, {timeout: 30_000});
    await page.locator('[name="viewer-Grid"] canvas').first().waitFor({state: 'attached', timeout: 30_000});
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
    await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
    const helmProbe = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macro = cols.find((c: any) => c.semType === 'Macromolecule') as any;
      return {
        hasMacromoleculeCol: !!macro,
        macroSemType: macro?.semType ?? null,
        macroUnits: macro?.getTag?.('units') ?? null,
        rowCount: df.rowCount,
      };
    });
    expect(helmProbe.hasMacromoleculeCol,
      'Scenario 2 Expected: HELM Macromolecule detector classifies synchronously on open').toBe(true);
    expect(helmProbe.macroSemType).toBe('Macromolecule');
    expect(helmProbe.macroUnits,
      'atlas bio.detector: filter_HELM.csv MUST classify with units=helm').toBe('helm');
    expect(helmProbe.rowCount,
      'filter_HELM.csv MUST carry >=2 rows so a non-trivial diverse subset can be selected').toBeGreaterThanOrEqual(2);
  });
  await softStep('Scenario 2.3: click Bio > Search > Diversity Search against the HELM TableView', async () => {
    await page.evaluate(() => (document.querySelector('[name="div-Bio"]') as HTMLElement).click());
    await page.locator('[name="div-Bio---Search"]').waitFor({state: 'attached', timeout: 15_000});
    await page.evaluate(() => {
      const group = document.querySelector('[name="div-Bio---Search"]')!;
      group.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      group.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
    });
    await page.locator('[name="div-Bio---Search---Diversity-Search"]').waitFor({state: 'attached', timeout: 15_000});
    await page.evaluate(() => (document.querySelector('[name="div-Bio---Search---Diversity-Search"]') as HTMLElement).click());
  });
  await softStep('Scenario 2.4: viewer re-docks on HELM view; renderMolIds drawn from new dataset', async () => {
    await page.locator('[name="viewer-Sequence-Diversity-Search"]').waitFor({timeout: 60_000});
    await page.waitForFunction(() => {
      const viewers = Array.from((grok.shell.tv as any).viewers) as any[];
      const v = viewers.find((vw) => vw.type === 'Sequence Diversity Search');
      if (!v) return false;
      const ids = (v as any).renderMolIds;
      return Array.isArray(ids) && ids.length > 0;
    }, null, {timeout: 180_000});
    const scenario2Probe = await page.evaluate(() => {
      const tv = grok.shell.tv as any;
      const viewers = Array.from(tv.viewers) as any[];
      const v = viewers.find((vw) => vw.type === 'Sequence Diversity Search');
      const renderMolIds: number[] | null = (v as any)?.renderMolIds ?? null;
      const subsetIndices: number[] | null = Array.isArray(renderMolIds) ? [...renderMolIds] : null;
      const df = tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macro = cols.find((c: any) => c.semType === 'Macromolecule') as any;
      const subsetSequences: (string | null)[] | null =
        subsetIndices ? subsetIndices.map((i) => macro?.get?.(i) ?? null) : null;
      return {
        viewerPresent: !!v,
        subsetIndices,
        subsetSize: subsetIndices?.length ?? null,
        subsetSequences,
        macroUnits: macro?.getTag?.('units') ?? null,
        dfId: df.id ?? null,
        dfName: df.name ?? null,
        scenario1Indices: (window as any).__scenario1SubsetIndices ?? null,
        scenario1RowCount: (window as any).__scenario1RowCount ?? null,
        scenario1DfId: (window as any).__scenario1DfId ?? null,
        currentRowCount: df.rowCount,
        viewerDockedDOM: !!document.querySelector('[name="viewer-Sequence-Diversity-Search"]'),
      };
    });
    expect(scenario2Probe.viewerPresent,
      'Scenario 2 Expected: viewer re-docks on the new HELM TableView').toBe(true);
    expect(scenario2Probe.viewerDockedDOM,
      'viewer dock container [name="viewer-Sequence-Diversity-Search"] MUST be present in DOM (Scenario 2)').toBe(true);
    expect(scenario2Probe.subsetSize,
      `Scenario 2 viewer MUST populate renderMolIds with >=2 entries. Observed: ${scenario2Probe.subsetSize}`)
      .not.toBeNull();
    expect(scenario2Probe.subsetSize!).toBeGreaterThanOrEqual(2);
    expect(scenario2Probe.subsetIndices,
      'renderMolIds MUST be a concrete index array on the HELM compute path').not.toBeNull();
    // Re-run proof: the viewer computed over the fresh HELM Macromolecule column, not the cached FASTA one.
    expect(scenario2Probe.macroUnits,
      'Scenario 2: viewer MUST compute over the HELM column (units=helm), confirming the diversity ' +
      'selection re-ran on the fresh dataset rather than re-displaying the cached FASTA subset').toBe('helm');
    if (scenario2Probe.scenario1DfId != null && scenario2Probe.dfId != null)
      expect(scenario2Probe.dfId,
        'Scenario 2: viewer re-docked against a fresh DataFrame, not the cached Scenario-1 FASTA table')
        .not.toBe(scenario2Probe.scenario1DfId);
    expect(scenario2Probe.subsetSequences,
      'Scenario 2 renderMolIds MUST resolve to HELM sequence values').not.toBeNull();
    const seqs2 = scenario2Probe.subsetSequences!.map((s) => String(s ?? ''));
    for (const s of seqs2) {
      expect(s.length, 'each Scenario-2 subset entry MUST be a non-empty sequence').toBeGreaterThan(0);
      expect(s.includes('$') || s.includes('{'),
        `Scenario 2 subset MUST contain HELM-notation sequences drawn from the HELM dataset, ` +
        `not the cached FASTA subset. Observed: ${s}`).toBe(true);
    }
    expect(new Set(seqs2).size,
      'Scenario 2 diverse subset MUST contain >=2 distinct sequences').toBeGreaterThanOrEqual(2);
    // Secondary guard: every selected index resolves inside the HELM rowCount range.
    for (const idx of scenario2Probe.subsetIndices!) {
      expect(idx, `Scenario 2 renderMolIds index ${idx} MUST be in [0, ${scenario2Probe.currentRowCount})`)
        .toBeGreaterThanOrEqual(0);
      expect(idx).toBeLessThan(scenario2Probe.currentRowCount);
    }
    const errorBalloonCount = await page.locator('.d4-balloon.error, .grok-balloon-error').count();
    expect(errorBalloonCount,
      'Scenario 2 Expected: "No error balloon appears" on positive path').toBe(0);
  });
  finishSpec();
});
