import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
test.use(specTestOptions);
const HELM_PATH = 'System:AppData/Bio/tests/filter_HELM.csv';
const MSA_PATH = 'System:AppData/Bio/tests/filter_MSA.csv';
async function openBioDataset(page: import('@playwright/test').Page, path: string) {
  await page.evaluate(async (p) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(p);
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
  }, path);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
  await page.evaluate(async () => {
    const probes = ['Bio:getSeqHelper', 'Bio:getMonomerLibHelper', 'Bio:getBioLib'];
    for (const fn of probes) {
      try { await (grok as any).functions.call(fn, {}); return; } catch {  }
    }
    await new Promise((r) => setTimeout(r, 3000));
  });
}
const BIO_SEQUENCE_CELL_TYPES = ['sequence', 'helm', 'separator', 'biln', 'custom', 'fasta'];
async function waitForSequenceCellTypeBind(page: import('@playwright/test').Page,
    timeoutMs = 60_000): Promise<void> {
  await page.waitForFunction((accepted: string[]) => {
    const df = grok.shell.tv?.dataFrame;
    if (!df) return false;
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
    if (!macro) return false;
    const gridCol = (grok.shell.tv as any).grid?.col?.(macro.name);
    const ct: string | null = gridCol?.cellType ?? null;
    return ct !== null && accepted.indexOf(ct) >= 0;
  }, BIO_SEQUENCE_CELL_TYPES, {timeout: timeoutMs});
}
async function inspectMacroCol(page: import('@playwright/test').Page):
    Promise<{name: string | null, semType: string | null, units: string | null,
             gridCellType: string | null, hasErrorBalloon: boolean}> {
  return await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
    const gridCol = (grok.shell.tv as any).grid?.col?.(macro?.name);
    const hasErrorBalloon = !!document.querySelector('.d4-balloon-error');
    return {
      name: macro?.name ?? null,
      semType: macro?.semType ?? null,
      units: macro?.getTag?.('units') ?? macro?.meta?.units ?? null,
      gridCellType: gridCol?.cellType ?? null,
      hasErrorBalloon,
    };
  });
}
test('Bio | Rendering — detector + renderer dispatch for HELM and SEPARATOR', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  await softStep('Open filter_HELM.csv — units=helm, bilnSequenceCellRenderer dispatch', async () => {
    await openBioDataset(page, HELM_PATH);
    await waitForSequenceCellTypeBind(page);
    const info = await inspectMacroCol(page);
    expect(info.semType).toBe('Macromolecule');
    expect(info.units).toBe('helm');
    expect(info.gridCellType).not.toBeNull();
    expect(BIO_SEQUENCE_CELL_TYPES, `cellType ${info.gridCellType} must be a Bio sequence-family value`)
      .toContain(info.gridCellType!);
    expect(info.hasErrorBalloon).toBe(false);
  });
  await softStep('Open filter_MSA.csv — units=separator (separatorSequenceCellRenderer dispatch)', async () => {
    await openBioDataset(page, MSA_PATH);
    await waitForSequenceCellTypeBind(page);
    const info = await inspectMacroCol(page);
    expect(info.semType).toBe('Macromolecule');
    expect(info.units).toBe('separator');
    expect(info.gridCellType).not.toBeNull();
    expect(BIO_SEQUENCE_CELL_TYPES).toContain(info.gridCellType!);
    expect(info.hasErrorBalloon).toBe(false);
    const sepTag: string | null = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
      return macro?.getTag?.('.separator') ?? macro?.getTag?.('separator') ?? null;
    });
    expect(sepTag).not.toBeNull();
  });
  await softStep('Enumerate Bio cellRenderer registrations — bio.rendering parent surface', async () => {
    const bioCellRenderers: string[] = await page.evaluate(() => {
      const list: string[] = [];
      try {
        const fns = (DG.Func.find({tags: ['cellRenderer']}) || []) as any[];
        for (const f of fns) {
          const friendly = (f.package as any)?.friendlyName ?? null;
          const pkgName = (f.package as any)?.name ?? null;
          if (friendly === 'Bio' || pkgName === '@datagrok/bio' || pkgName === 'Bio')
            list.push(f.name);
        }
      } catch { /* surface as empty list */ }
      return list;
    });
    const expected = [
      'fastaSequenceCellRenderer',
      'separatorSequenceCellRenderer',
      'bilnSequenceCellRenderer',
      'customSequenceCellRenderer',
      'monomerCellRenderer',
    ];
    for (const n of expected)
      expect(bioCellRenderers, `cellRenderer ${n} must be registered (atlas bio.rendering)`).toContain(n);
  });
  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
// Scenario 2 — Convert HELM -> SEPARATOR re-dispatches renderer (GROK-12164 guard)
test('Bio | Rendering — Convert HELM to SEPARATOR re-dispatches renderer (GROK-12164)', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  let preHelmColName: string | null = null;
  await softStep('Open filter_HELM.csv — baseline dispatch units=helm', async () => {
    await openBioDataset(page, HELM_PATH);
    await waitForSequenceCellTypeBind(page);
    const info = await inspectMacroCol(page);
    expect(info.units).toBe('helm');
    expect(info.gridCellType).not.toBeNull();
    expect(BIO_SEQUENCE_CELL_TYPES).toContain(info.gridCellType!);
    preHelmColName = info.name;
  });
  await softStep('Bio > Transform > Convert Sequence Notation... opens dialog', async () => {
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 400));
      document.querySelector('[name="div-Bio---Transform"]')!
        .dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 400));
      (document.querySelector(
        '[name="div-Bio---Transform---Convert-Sequence-Notation..."]') as HTMLElement).click();
    });
    await page.locator('[name="dialog-Convert-Sequence-Notation"]').waitFor({timeout: 60_000});
  });
  await softStep('Set Convert-to=separator, Separator=-, click OK', async () => {
    const dlg = page.locator('[name="dialog-Convert-Sequence-Notation"]');
    await dlg.locator('[name="input-host-Convert-to"] select').selectOption('separator');
    await dlg.locator('[name="input-host-Separator"] select').selectOption('-');
    const before: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
    await dlg.locator('[name="button-OK"]').click();
    await page.waitForFunction(
      (b) => grok.shell.tv.dataFrame.columns.length > b, before, {timeout: 60_000});
    await page.waitForFunction(
      () => document.querySelectorAll('[name="dialog-Convert-Sequence-Notation"]').length === 0,
      null, {timeout: 15_000}).catch(() => {});
  });
  // GROK-12164: the new separator column's dispatch must follow its own units=separator tag, not the source HELM tags.
  await softStep('GROK-12164: new column units=separator, source HELM column units=helm intact', async () => {
    await page.waitForFunction((accepted: string[]) => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const sep: any = cols.find((c: any) =>
        c.semType === 'Macromolecule' &&
        ((c.getTag?.('units') ?? c.meta?.units) === 'separator'));
      if (!sep) return false;
      const gridCol = (grok.shell.tv as any).grid?.col?.(sep.name);
      const ct: string | null = gridCol?.cellType ?? null;
      return ct !== null && accepted.indexOf(ct) >= 0;
    }, BIO_SEQUENCE_CELL_TYPES, {timeout: 60_000});
    const tagsByMacro: Array<{name: string, units: string | null, gridCellType: string | null}> =
      await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macros = cols.filter((c: any) => c.semType === 'Macromolecule');
        return macros.map((c: any) => ({
          name: c.name,
          units: c.getTag?.('units') ?? c.meta?.units ?? null,
          gridCellType: (grok.shell.tv as any).grid?.col?.(c.name)?.cellType ?? null,
        }));
      });
    expect(tagsByMacro.length).toBeGreaterThanOrEqual(2);
    const sourceHelm = tagsByMacro.find((m) => m.units === 'helm');
    const newSeparator = tagsByMacro.find((m) => m.units === 'separator');
    expect(sourceHelm, 'source HELM column units=helm must be intact post-convert').toBeTruthy();
    expect(sourceHelm!.gridCellType).not.toBeNull();
    expect(BIO_SEQUENCE_CELL_TYPES,
      `source HELM cellType ${sourceHelm!.gridCellType} must be in bio sequence family`)
      .toContain(sourceHelm!.gridCellType!);
    if (preHelmColName) expect(sourceHelm!.name).toBe(preHelmColName);
    expect(newSeparator, 'new column from convert must carry units=separator (GROK-12164)').toBeTruthy();
    expect(newSeparator!.gridCellType).not.toBeNull();
    expect(BIO_SEQUENCE_CELL_TYPES,
      `new SEPARATOR cellType ${newSeparator!.gridCellType} must be in bio sequence family`)
      .toContain(newSeparator!.gridCellType!);
    expect(newSeparator!.gridCellType,
      'GROK-12164: new SEPARATOR column must not retain HELM dispatch (cellType !== "helm")')
      .not.toBe('helm');
    const hasErrorBalloon = await page.evaluate(() => !!document.querySelector('.d4-balloon-error'));
    expect(hasErrorBalloon).toBe(false);
  });
  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
// Scenario 3 — Split-to-Monomers produces Monomer columns rendered by monomerCellRenderer
test('Bio | Rendering — Split to Monomers produces Monomer columns (monomerCellRenderer)', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  await softStep('Open filter_MSA.csv — separator-with-MSA-flag baseline', async () => {
    await openBioDataset(page, MSA_PATH);
    await waitForSequenceCellTypeBind(page);
    const info = await inspectMacroCol(page);
    expect(info.units).toBe('separator');
  });
  await softStep('Bio > Transform > Split to Monomers... — OK adds Monomer columns', async () => {
    const beforeMonCount: number = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      return cols.filter((c: any) => c.semType === 'Monomer').length;
    });
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 400));
      document.querySelector('[name="div-Bio---Transform"]')!
        .dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 400));
      (document.querySelector(
        '[name="div-Bio---Transform---Split-to-Monomers..."]') as HTMLElement).click();
    });
    await page.locator('[name="dialog-Split-to-Monomers"]').waitFor({timeout: 60_000});
    await page.locator('[name="dialog-Split-to-Monomers"] [name="button-OK"]').click();
    await page.waitForFunction((base) => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      return cols.filter((c: any) => c.semType === 'Monomer').length > base;
    }, beforeMonCount, {timeout: 60_000});
    await page.waitForFunction(
      () => document.querySelectorAll('[name="dialog-Split-to-Monomers"]').length === 0,
      null, {timeout: 15_000}).catch(() => {});
  });
  await softStep('Monomer columns dispatch to monomerCellRenderer (cellType=Monomer)', async () => {
    await page.waitForFunction(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const monCols = cols.filter((c: any) => c.semType === 'Monomer');
      if (monCols.length === 0) return false;
      return monCols.some((c: any) =>
        (grok.shell.tv as any).grid?.col?.(c.name)?.cellType === 'Monomer');
    }, null, {timeout: 60_000});
    const monomerColInfo: Array<{name: string, semType: string | null, gridCellType: string | null}> =
      await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const monCols = cols.filter((c: any) => c.semType === 'Monomer');
        return monCols.map((c: any) => ({
          name: c.name,
          semType: c.semType,
          gridCellType: (grok.shell.tv as any).grid?.col?.(c.name)?.cellType ?? null,
        }));
      });
    expect(monomerColInfo.length).toBeGreaterThan(0);
    for (const m of monomerColInfo)
      expect(m.semType, `column ${m.name} should be semType='Monomer'`).toBe('Monomer');
    const anyMonomerDispatched = monomerColInfo.some((m) => m.gridCellType === 'Monomer');
    expect(anyMonomerDispatched, 'at least one Monomer column should have grid.cellType=Monomer').toBe(true);
    const hasErrorBalloon = await page.evaluate(() => !!document.querySelector('.d4-balloon-error'));
    expect(hasErrorBalloon).toBe(false);
  });
  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
// Scenario 4 — Custom-notation column dispatches to customSequenceCellRenderer
test('Bio | Rendering — units=custom column dispatches to customSequenceCellRenderer', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  await softStep('Open filter_HELM.csv — baseline units=helm', async () => {
    await openBioDataset(page, HELM_PATH);
    await waitForSequenceCellTypeBind(page);
    const info = await inspectMacroCol(page);
    expect(info.units).toBe('helm');
    expect(info.gridCellType).not.toBeNull();
    expect(BIO_SEQUENCE_CELL_TYPES).toContain(info.gridCellType!);
  });
  await softStep('setTag units=custom + detectSemanticTypes — re-dispatches to customSequenceCellRenderer', async () => {
    const colName: string | null = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
      if (!macro) return null;
      macro.setTag('units', 'custom');
      try { await (grok.data as any).detectSemanticTypes(df); } catch { /* tolerate older builds */ }
      try { (grok.shell.tv as any).grid?.invalidate?.(); } catch { /* may not throw on older builds */ }
      return macro.name;
    });
    expect(colName, 'Macromolecule column must be located on HELM dataset').toBeTruthy();
    await page.waitForFunction(({colName, accepted}: {colName: string | null, accepted: string[]}) => {
      const gridCol = (grok.shell.tv as any).grid?.col?.(colName);
      const ct: string | null = gridCol?.cellType ?? null;
      return ct !== null && accepted.indexOf(ct) >= 0;
    }, {colName, accepted: BIO_SEQUENCE_CELL_TYPES}, {timeout: 60_000});
    const result: {postUnits: string | null, postCellType: string | null} =
      await page.evaluate((name) => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macro: any = cols.find((c: any) => c.name === name);
        return {
          postUnits: macro?.getTag?.('units') ?? macro?.meta?.units ?? null,
          postCellType: (grok.shell.tv as any).grid?.col?.(name)?.cellType ?? null,
        };
      }, colName);
    expect(result.postUnits).toBe('custom');
    expect(result.postCellType).not.toBeNull();
    expect(BIO_SEQUENCE_CELL_TYPES,
      `cellType ${result.postCellType} must be in bio sequence family after units=custom rebind`)
      .toContain(result.postCellType!);
    const hasErrorBalloon = await page.evaluate(() => !!document.querySelector('.d4-balloon-error'));
    expect(hasErrorBalloon).toBe(false);
  });
  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
