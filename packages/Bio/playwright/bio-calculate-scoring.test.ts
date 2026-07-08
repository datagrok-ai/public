/* ---
sub_features_covered: [bio.analyze.alignment-pairwise, bio.calculate.get-region, bio.calculate.get-region.api, bio.calculate.identity, bio.calculate.seq-identity, bio.calculate.similarity]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';
test.use(specTestOptions);
const HELM_DATASET_PATH = 'System:AppData/Bio/tests/filter_HELM.csv';
const FASTA_DATASET_PATH = 'System:AppData/Bio/tests/filter_FASTA.csv';
async function openBioDataset(page: import('@playwright/test').Page, path: string): Promise<void> {
  await page.evaluate(async (datasetPath) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(datasetPath);
    grok.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(() => resolve(), 4000);
    });
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    const hasMacromolecule = cols.some((c: any) => c.semType === 'Macromolecule');
    if (hasMacromolecule) {
      const deadline = Date.now() + 20_000;
      while (Date.now() < deadline) {
        const canvas = document.querySelector('[name="viewer-Grid"] canvas') as HTMLCanvasElement | null;
        if (canvas && canvas.width > 0 && canvas.height > 0) break;
        await new Promise((r) => setTimeout(r, 200));
      }
    }
  }, path);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
  await page.evaluate(async () => {
    const probes = ['Bio:getSeqHelper', 'Bio:getMonomerLibHelper', 'Bio:getBioLib'];
    const deadline = Date.now() + 15_000;
    while (Date.now() < deadline) {
      for (const fn of probes) {
        try { await (grok as any).functions.call(fn, {}); return; } catch {  }
      }
      await new Promise((r) => setTimeout(r, 300));
    }
  });
  await page.evaluate(async () => {
    const candidates = ['Bio:Identity', 'Bio:Similarity', 'Bio:getRegion',
      'Bio:seqIdentity', 'Bio:sequenceAlignment'];
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
}
async function openBioCalculateLeaf(
  page: import('@playwright/test').Page,
  leafSelector: string,
  leafLabel: string,
): Promise<void> {
  await page.evaluate(async ({leafSel, leafName}) => {
    (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
    for (let i = 0; i < 50; i++) {
      if (document.querySelector('[name="div-Bio---Calculate"]')) break;
      await new Promise((r) => setTimeout(r, 100));
    }
    const group = document.querySelector('[name="div-Bio---Calculate"]');
    if (!group) throw new Error('[name="div-Bio---Calculate"] group anchor not found');
    group.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
    group.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
    for (let i = 0; i < 50; i++) {
      if (document.querySelector(leafSel)) break;
      await new Promise((r) => setTimeout(r, 100));
    }
    const leaf = document.querySelector(leafSel);
    if (!leaf)
      throw new Error(`${leafSel} leaf (${leafName}) not found under Bio > Calculate`);
    (leaf as HTMLElement).click();
  }, {leafSel: leafSelector, leafName: leafLabel});
}
test('Bio Calculate scoring — Identity / Similarity top-menu + seqIdentity / getRegion / sequenceAlignment API', async ({page}) => {
  test.setTimeout(180_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  // Scenario 1 — Identity scoring via top-menu (filter_HELM.csv).
  await openBioDataset(page, HELM_DATASET_PATH);
  const helmSetup = await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    const macro = cols.find((c: any) => c.semType === 'Macromolecule') as any;
    if (!macro) throw new Error('Setup: filter_HELM.csv MUST classify a Macromolecule column synchronously');
    const firstSeq = macro.get(0) as string;
    return {
      preColumnCount: df.columns.length,
      macroName: macro.name,
      macroSemType: macro.semType,
      macroUnits: macro.getTag?.('units') ?? null,
      rowCount: df.rowCount,
      firstSeq,
    };
  });
  expect(helmSetup.macroSemType,
    'atlas bio.detector contract: filter_HELM.csv MUST classify Macromolecule synchronously').toBe('Macromolecule');
  expect(helmSetup.macroUnits,
    'atlas bio.detector: filter_HELM.csv MUST classify with units=helm (HELM-monomer fingerprint paths)').toBe('helm');
  expect(helmSetup.rowCount,
    'scenario .md Setup: HELM dataset MUST have >=2 rows for cross-row scoring').toBeGreaterThanOrEqual(2);
  expect(helmSetup.firstSeq).toBeTruthy();
  await softStep('Scenario 1.1: click Bio > Calculate > Identity... (dialog opens)', async () => {
    await openBioCalculateLeaf(page,
      '[name="div-Bio---Calculate---Identity..."]', 'Identity...');
  });
  await softStep('Scenario 1.2: dialog opens; fill Reference with row 0 HELM sequence; click OK', async () => {
    await page.locator('[name="dialog-Identity"]').waitFor({timeout: 30_000});
    await page.evaluate((refValue) => {
      const host = document.querySelector('[name="dialog-Identity"] [name="input-host-Reference"]');
      if (!host) throw new Error('Identity dialog: Reference input host not found');
      const input = host.querySelector('input, textarea') as HTMLInputElement | HTMLTextAreaElement | null;
      if (!input) throw new Error('Identity dialog: Reference input element not found inside host');
      input.focus();
      input.value = refValue;
      input.dispatchEvent(new Event('input', {bubbles: true}));
      input.dispatchEvent(new Event('change', {bubbles: true}));
      input.blur();
    }, helmSetup.firstSeq);
    await page.locator('[name="dialog-Identity"] [name="button-OK"]').click();
  });
  await softStep('Scenario 1.3-4: dialog closes; new Identity score column appended; row 0 == 1.0; rows in [0,1]', async () => {
    await page.locator('[name="dialog-Identity"]').waitFor({state: 'detached', timeout: 60_000});
    await page.waitForFunction((preCount) => {
      const df = grok.shell.tv.dataFrame;
      return df.columns.length > preCount;
    }, helmSetup.preColumnCount, {timeout: 120_000});
    const identityProbe = await page.evaluate((preCount) => {
      const df = grok.shell.tv.dataFrame;
      const allCols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const newCols = allCols.slice(preCount);
      const numericCols = newCols.filter((c: any) =>
        c.type === 'double' || c.type === 'int' || c.type === 'float' || c.type === 'qnum');
      if (numericCols.length === 0) {
        return {
          addedCols: newCols.map((c: any) => ({name: c.name, type: c.type})),
          numericCount: 0,
        };
      }
      const scoreCol = numericCols[0] as any;
      const values: number[] = [];
      for (let i = 0; i < df.rowCount; i++) {
        const v = scoreCol.get(i);
        if (typeof v === 'number' && !Number.isNaN(v)) values.push(v);
      }
      const row0 = scoreCol.get(0) as number;
      return {
        addedCols: newCols.map((c: any) => ({name: c.name, type: c.type})),
        numericCount: numericCols.length,
        scoreColName: scoreCol.name,
        scoreColType: scoreCol.type,
        row0,
        min: values.length > 0 ? Math.min(...values) : null,
        max: values.length > 0 ? Math.max(...values) : null,
        nonNullCount: values.length,
      };
    }, helmSetup.preColumnCount);
    expect(identityProbe.numericCount,
      `Scenario 1 Expected: a new numeric score column MUST be appended. Added cols: ${JSON.stringify(identityProbe.addedCols)}`)
      .toBeGreaterThanOrEqual(1);
    expect(identityProbe.row0,
      `Scenario 1 Expected: row 0 self-identity MUST be ~1.0 (tol 0.01 per scoring.ts precedent); got ${identityProbe.row0}`)
      .toBeCloseTo(1, 2);
    expect(identityProbe.min!,
      `Scenario 1 Expected: min Identity score MUST be >= 0.0; got ${identityProbe.min}`).toBeGreaterThanOrEqual(0);
    expect(identityProbe.max!,
      `Scenario 1 Expected: max Identity score MUST be <= 1.0 (tolerated +0.001 rounding); got ${identityProbe.max}`)
      .toBeLessThanOrEqual(1.001);
    const errorBalloonCount = await page.locator('.d4-balloon.error, .grok-balloon-error').count();
    expect(errorBalloonCount,
      'Scenario 1 Expected: "No error balloon appears"').toBe(0);
  });
  // Scenario 2 — Similarity scoring via top-menu (same filter_HELM.csv view).
  const preSimilarityColumnCount = await page.evaluate(() =>
    grok.shell.tv.dataFrame.columns.length);
  await softStep('Scenario 2.1: click Bio > Calculate > Similarity... (dialog opens)', async () => {
    await openBioCalculateLeaf(page,
      '[name="div-Bio---Calculate---Similarity..."]', 'Similarity...');
  });
  await softStep('Scenario 2.2: Similarity dialog opens; fill Reference; click OK', async () => {
    await page.waitForFunction(() => {
      const sim = document.querySelector('[name="dialog-Similarity"]');
      const anyDialogWithRef = Array.from(document.querySelectorAll('.d4-dialog'))
        .some((d) => d.querySelector('[name="input-host-Reference"]'));
      return !!sim || anyDialogWithRef;
    }, null, {timeout: 30_000});
    await page.evaluate((refValue) => {
      let dialog = document.querySelector('[name="dialog-Similarity"]') as HTMLElement | null;
      if (!dialog) {
        const candidates = Array.from(document.querySelectorAll('.d4-dialog'))
          .filter((d) => d.querySelector('[name="input-host-Reference"]')) as HTMLElement[];
        const open = candidates.filter((d) => !d.matches('[name="dialog-Identity"]'));
        dialog = (open.length > 0 ? open[open.length - 1] : candidates[candidates.length - 1]) ?? null;
      }
      if (!dialog) throw new Error('Similarity dialog: no open dialog with [name="input-host-Reference"] found');
      const host = dialog.querySelector('[name="input-host-Reference"]');
      if (!host) throw new Error('Similarity dialog: Reference input host not found');
      const input = host.querySelector('input, textarea') as HTMLInputElement | HTMLTextAreaElement | null;
      if (!input) throw new Error('Similarity dialog: Reference input element not found inside host');
      input.focus();
      input.value = refValue;
      input.dispatchEvent(new Event('input', {bubbles: true}));
      input.dispatchEvent(new Event('change', {bubbles: true}));
      input.blur();
      const ok = dialog.querySelector('[name="button-OK"]') as HTMLElement | null;
      if (!ok) throw new Error('Similarity dialog: OK button not found');
      ok.click();
    }, helmSetup.firstSeq);
  });
  await softStep('Scenario 2.3: dialog closes; new Similarity column appended; row 0 self-similarity is the column max', async () => {
    await page.waitForFunction(() => {
      const sim = document.querySelector('[name="dialog-Similarity"]');
      const anyDialogWithRef = Array.from(document.querySelectorAll('.d4-dialog'))
        .some((d) => d.querySelector('[name="input-host-Reference"]'));
      return !sim && !anyDialogWithRef;
    }, null, {timeout: 60_000});
    await page.waitForFunction((preCount) => {
      const df = grok.shell.tv.dataFrame;
      return df.columns.length > preCount;
    }, preSimilarityColumnCount, {timeout: 120_000});
    const similarityProbe = await page.evaluate((preCount) => {
      const df = grok.shell.tv.dataFrame;
      const allCols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const newCols = allCols.slice(preCount);
      const numericCols = newCols.filter((c: any) =>
        c.type === 'double' || c.type === 'int' || c.type === 'float' || c.type === 'qnum');
      if (numericCols.length === 0) {
        return {
          addedCols: newCols.map((c: any) => ({name: c.name, type: c.type})),
          numericCount: 0,
        };
      }
      const scoreCol = numericCols[0] as any;
      const row0 = scoreCol.get(0) as number;
      const values: number[] = [];
      for (let i = 0; i < df.rowCount; i++) {
        const v = scoreCol.get(i);
        if (typeof v === 'number' && !Number.isNaN(v)) values.push(v);
      }
      return {
        addedCols: newCols.map((c: any) => ({name: c.name, type: c.type})),
        numericCount: numericCols.length,
        scoreColName: scoreCol.name,
        scoreColType: scoreCol.type,
        row0,
        nonNullCount: values.length,
        min: values.length > 0 ? Math.min(...values) : null,
        max: values.length > 0 ? Math.max(...values) : null,
        totalColumnCount: df.columns.length,
      };
    }, preSimilarityColumnCount);
    expect(similarityProbe.numericCount,
      `Scenario 2 Expected: a new numeric Similarity score column MUST be appended. Added cols: ${JSON.stringify(similarityProbe.addedCols)}`)
      .toBeGreaterThanOrEqual(1);
    expect(similarityProbe.nonNullCount,
      'Scenario 2 Expected: Similarity scores are non-null for at least 2 rows').toBeGreaterThanOrEqual(2);
    expect(similarityProbe.min!,
      `Scenario 2 Expected: Similarity scores are non-negative; got min=${similarityProbe.min}`).toBeGreaterThanOrEqual(0);
    expect(similarityProbe.row0,
      `Scenario 2 Expected: reference is row 0, so its self-similarity MUST be the column max; got row0=${similarityProbe.row0}, max=${similarityProbe.max}`)
      .toBe(similarityProbe.max);
    expect(similarityProbe.max!,
      `Scenario 2 Expected: similarity scores vary across rows (max > min); got max=${similarityProbe.max}, min=${similarityProbe.min}`)
      .toBeGreaterThan(similarityProbe.min!);
    expect(similarityProbe.totalColumnCount,
      'Scenario 2 Expected: Identity column from Scenario 1 + Similarity column from Scenario 2 coexist')
      .toBeGreaterThan(preSimilarityColumnCount);
    const errorBalloonCount = await page.locator('.d4-balloon.error, .grok-balloon-error').count();
    expect(errorBalloonCount,
      'Scenario 2 Expected: "No error balloon appears"').toBe(0);
  });
  // Scenario 3 — getRegion API on filter_FASTA.csv.
  await openBioDataset(page, FASTA_DATASET_PATH);
  const fastaSetup = await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    const macro = cols.find((c: any) => c.semType === 'Macromolecule') as any;
    if (!macro) throw new Error('Setup: filter_FASTA.csv MUST classify a Macromolecule column synchronously');
    const samples: string[] = [];
    for (let i = 0; i < Math.min(df.rowCount, 3); i++) samples.push(macro.get(i) as string);
    return {
      macroName: macro.name,
      macroSemType: macro.semType,
      macroUnits: macro.getTag?.('units') ?? null,
      rowCount: df.rowCount,
      preColumnCount: df.columns.length,
      samples,
    };
  });
  expect(fastaSetup.macroSemType).toBe('Macromolecule');
  expect(fastaSetup.macroUnits,
    'atlas bio.detector: filter_FASTA.csv MUST classify with units=fasta').toBe('fasta');
  expect(fastaSetup.rowCount).toBeGreaterThanOrEqual(2);
  await softStep('Scenario 3.1-2: invoke Bio:getRegion via grok.functions.call(sequence, start, end, name)', async () => {
    const regionResult = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macroCol = cols.find((c: any) => c.semType === 'Macromolecule') as any;
      if (!macroCol) throw new Error('Scenario 3: Macromolecule column missing on FASTA setup');
      let regionCol: any;
      let callError: string | null = null;
      try {
        // Positions are looked up by 1-indexed NAME against posList (get-region.ts:128,137-142),
        // so start must be >= '1'; '0' has no matching position name and is correctly rejected.
        regionCol = await (grok as any).functions.call('Bio:getRegion', {
          sequence: macroCol,
          start: '1',
          end: '4',
          name: 'region_0_4',
        });
      } catch (e: any) {
        callError = e?.message ?? String(e);
      }
      if (callError || !regionCol)
        return {ok: false, error: callError ?? 'getRegion returned null'};
      const sampleCells: string[] = [];
      const len = regionCol.length ?? regionCol.rowCount ?? df.rowCount;
      for (let i = 0; i < Math.min(len, 3); i++) sampleCells.push(regionCol.get(i) as string);
      const result = {
        ok: true,
        name: regionCol.name as string,
        type: regionCol.type as string,
        length: len as number,
        semType: regionCol.semType as string | undefined,
        units: regionCol.getTag?.('units') as string | undefined,
        sampleCells,
      };
      df.columns.add(regionCol);
      return result;
    });
    expect(regionResult.ok,
      `Scenario 3 Expected: Bio:getRegion call MUST succeed; got: ${regionResult.error}`).toBe(true);
    expect(regionResult.name,
      'Scenario 3 Expected: returned column name MUST be the `name` parameter "region_0_4"').toBe('region_0_4');
    expect(regionResult.type,
      'Scenario 3 Expected: getRegion returns DG.Column<string>').toBe('string');
    expect(regionResult.length,
      'Scenario 3 Expected: returned column row count matches source row count').toBe(fastaSetup.rowCount);
    expect(regionResult.sampleCells.length).toBeGreaterThan(0);
    for (let i = 0; i < regionResult.sampleCells.length; i++) {
      const cell = regionResult.sampleCells[i];
      const source = fastaSetup.samples[i];
      expect(typeof cell === 'string' && cell.length > 0,
        `Scenario 3 Expected: region cell for row ${i} is a non-empty positional slice; got ${JSON.stringify(cell)}`)
        .toBe(true);
      expect(source.startsWith(cell),
        `Scenario 3 Expected: region_0_4 cell for row ${i} MUST be the leading slice (first monomers) of the source sequence "${source}"; got "${cell}"`)
        .toBe(true);
      expect(cell.length,
        `Scenario 3 Expected: leading region for row ${i} is shorter than the full source "${source}"; got "${cell}"`)
        .toBeLessThan(source.length);
    }
    expect(regionResult.semType,
      'Scenario 3 Expected: returned column MUST preserve Macromolecule semType').toBe('Macromolecule');
    expect(regionResult.units,
      'Scenario 3 Expected: returned column MUST preserve units=fasta').toBe('fasta');
    const errorBalloonCount = await page.locator('.d4-balloon.error, .grok-balloon-error').count();
    expect(errorBalloonCount,
      'Scenario 3 Expected: "No error balloon appears"').toBe(0);
  });
  // Scenario 4 — seqIdentity API (single-pair + empty-input contract).
  await openBioDataset(page, HELM_DATASET_PATH);
  const helmSeq4Setup = await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    const macro = cols.find((c: any) => c.semType === 'Macromolecule') as any;
    if (!macro) throw new Error('Scenario 4 setup: Macromolecule column missing on HELM re-open');
    return {
      row0: macro.get(0) as string,
      row1: macro.get(1) as string,
      rowCount: df.rowCount,
    };
  });
  expect(helmSeq4Setup.rowCount).toBeGreaterThanOrEqual(2);
  expect(helmSeq4Setup.row0).toBeTruthy();
  expect(helmSeq4Setup.row1).toBeTruthy();
  await softStep('Scenario 4.1: Bio:sequenceIdentityScoring(table, col, ref=row0Helm) returns row 0 ~1.0 (self-identity)', async () => {
    const r = await page.evaluate(async (ref) => {
      try {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macroCol = cols.find((c: any) => c.semType === 'Macromolecule') as any;
        if (!macroCol)
          return {ok: false, error: 'Scenario 4.1 setup: Macromolecule column missing on HELM re-open'};
        const scoresCol = await (grok as any).functions.call('Bio:sequenceIdentityScoring',
          {table: df, macromolecule: macroCol, reference: ref});
        if (!scoresCol)
          return {ok: false, error: 'Bio:sequenceIdentityScoring returned null'};
        return {ok: true, value: scoresCol.get(0) as number | null};
      } catch (e: any) {
        return {ok: false, error: e?.message ?? String(e)};
      }
    }, helmSeq4Setup.row0);
    expect(r.ok,
      `Scenario 4.1: sequenceIdentityScoring self-call MUST succeed; got: ${r.error}`).toBe(true);
    expect(r.value,
      `Scenario 4.1 Expected: row 0 self-identity ~1.0 (tol 0.01 per scoring.ts precedent); got ${r.value}`)
      .toBeCloseTo(1, 2);
  });
  await softStep('Scenario 4.2: Bio:seqIdentity(seq="", ref=row0) returns null (empty-input contract)', async () => {
    const r = await page.evaluate(async ([ref]) => {
      try {
        const v = await (grok as any).functions.call('Bio:seqIdentity', {seq: '', ref});
        return {ok: true, value: v as number | null};
      } catch (e: any) {
        return {ok: false, error: e?.message ?? String(e)};
      }
    }, [helmSeq4Setup.row0]);
    expect(r.ok,
      `Scenario 4.2: seqIdentity empty-seq call MUST succeed (not throw); got: ${r.error}`).toBe(true);
    const emptyInputNoValue = r.value === null || r.value === undefined;
    expect(emptyInputNoValue,
      `Scenario 4.2 Expected: seqIdentity(seq='', ref) MUST return null/undefined ` +
      `(empty-value branch in calculateScoresWithEmptyValues); got ${JSON.stringify(r.value)}`)
      .toBe(true);
  });
  await softStep('Scenario 4.3: Bio:sequenceIdentityScoring(table, col, ref=row1Helm) → row 0 score in [0,1]', async () => {
    const r = await page.evaluate(async (ref) => {
      try {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macroCol = cols.find((c: any) => c.semType === 'Macromolecule') as any;
        if (!macroCol)
          return {ok: false, error: 'Scenario 4.3 setup: Macromolecule column missing'};
        const scoresCol = await (grok as any).functions.call('Bio:sequenceIdentityScoring',
          {table: df, macromolecule: macroCol, reference: ref});
        if (!scoresCol)
          return {ok: false, error: 'Bio:sequenceIdentityScoring returned null'};
        return {ok: true, value: scoresCol.get(0) as number | null};
      } catch (e: any) {
        return {ok: false, error: e?.message ?? String(e)};
      }
    }, helmSeq4Setup.row1);
    expect(r.ok,
      `Scenario 4.3: cross-row sequenceIdentityScoring MUST succeed; got: ${r.error}`).toBe(true);
    expect(r.value,
      'Scenario 4.3 Expected: cross-row score is a non-null number').not.toBeNull();
    expect(typeof r.value).toBe('number');
    expect(r.value!,
      `Scenario 4.3 Expected: cross-row Identity score >= 0; got ${r.value}`).toBeGreaterThanOrEqual(0);
    expect(r.value!,
      `Scenario 4.3 Expected: cross-row Identity score <= 1.0 (tolerated +0.001 rounding); got ${r.value}`)
      .toBeLessThanOrEqual(1.001);
    if (helmSeq4Setup.row0 !== helmSeq4Setup.row1)
      expect(r.value!,
        `Scenario 4.3 Expected: distinct sequences give a non-trivial score strictly below the self-identity 1.0; got ${r.value}`)
        .toBeLessThan(1);
    const errorBalloonCount = await page.locator('.d4-balloon.error, .grok-balloon-error').count();
    expect(errorBalloonCount,
      'Scenario 4 Expected: "No error balloon appears across the three invocations"').toBe(0);
  });
  // Scenario 5 — sequenceAlignment API (global + local).
  await openBioDataset(page, FASTA_DATASET_PATH);
  const fastaSeq5Setup = await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    const macro = cols.find((c: any) => c.semType === 'Macromolecule') as any;
    if (!macro) throw new Error('Scenario 5 setup: Macromolecule column missing on FASTA re-open');
    return {
      row0: macro.get(0) as string,
      row1: macro.get(1) as string,
    };
  });
  expect(fastaSeq5Setup.row0).toBeTruthy();
  expect(fastaSeq5Setup.row1).toBeTruthy();
  await softStep('Scenario 5.1: Bio:sequenceAlignment global (Needleman-Wunsch / BLOSUM62) returns non-null result', async () => {
    const r = await page.evaluate(async ([seq1, seq2]) => {
      try {
        const result = await (grok as any).functions.call('Bio:sequenceAlignment', {
          alignType: 'Global alignment',
          alignTable: 'BLOSUM62',
          gap: -10,
          seq1,
          seq2,
        });
        const sampleKeys = result && typeof result === 'object' ? Object.keys(result).slice(0, 10) : [];
        return {
          ok: true,
          isNull: result == null,
          isObject: result != null && typeof result === 'object',
          stringFields: result && typeof result === 'object'
            ? Object.entries(result)
                .filter(([, v]) => typeof v === 'string')
                .map(([k, v]) => ({k, len: (v as string).length}))
            : [],
          sampleKeys,
        };
      } catch (e: any) {
        return {ok: false, error: e?.message ?? String(e)};
      }
    }, [fastaSeq5Setup.row0, fastaSeq5Setup.row1]);
    expect(r.ok,
      `Scenario 5.1: Bio:sequenceAlignment global call MUST succeed; got: ${r.error}`).toBe(true);
    expect(r.isNull,
      'Scenario 5 Expected: global alignment returns a non-null result').toBe(false);
    expect(r.isObject,
      'Scenario 5 Expected: global alignment returns an object').toBe(true);
    const maxInputLen = Math.max(fastaSeq5Setup.row0.length, fastaSeq5Setup.row1.length);
    const longEnoughString = r.stringFields!.some((f: any) => f.len >= maxInputLen);
    expect(longEnoughString,
      `Scenario 5 Expected: result contains a string field whose length >= max input length (${maxInputLen}); ` +
      `got: ${JSON.stringify(r.stringFields)}`).toBe(true);
  });
  await softStep('Scenario 5.3: Bio:sequenceAlignment local (Smith-Waterman / BLOSUM45) returns non-null result', async () => {
    const r = await page.evaluate(async ([seq1, seq2]) => {
      try {
        const result = await (grok as any).functions.call('Bio:sequenceAlignment', {
          alignType: 'Local alignment',
          alignTable: 'BLOSUM45',
          gap: -10,
          seq1,
          seq2,
        });
        const stringFields = result && typeof result === 'object'
          ? Object.entries(result)
              .filter(([, v]) => typeof v === 'string')
              .map(([k, v]) => ({k, len: (v as string).length}))
          : [];
        return {
          ok: true,
          isNull: result == null,
          isObject: result != null && typeof result === 'object',
          stringFields,
        };
      } catch (e: any) {
        return {ok: false, error: e?.message ?? String(e)};
      }
    }, [fastaSeq5Setup.row0, fastaSeq5Setup.row1]);
    expect(r.ok,
      `Scenario 5.3: Bio:sequenceAlignment local call MUST succeed; got: ${r.error}`).toBe(true);
    expect(r.isNull,
      'Scenario 5 Expected: local alignment (Smith-Waterman) returns a non-null result').toBe(false);
    expect(r.isObject).toBe(true);
    const hasNonEmptyString = r.stringFields!.some((f: any) => f.len > 0);
    expect(hasNonEmptyString,
      `Scenario 5 Expected: local alignment result contains at least one non-empty string field; ` +
      `got: ${JSON.stringify(r.stringFields)}`).toBe(true);
    const errorBalloonCount = await page.locator('.d4-balloon.error, .grok-balloon-error').count();
    expect(errorBalloonCount,
      'Scenario 5 Expected: "No error balloon appears"').toBe(0);
  });
  finishSpec();
});
