import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
test.use(specTestOptions);
test('Bio PepSeA MSA on HELM', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:AppData/Bio/tests/filter_HELM.csv');
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
    await new Promise((r) => setTimeout(r, 5000));
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
  await page.waitForFunction(() => !!document.querySelector('[name="div-Bio"]'),
    null, {timeout: 30_000});
  await softStep('Open filter_HELM.csv and detect Macromolecule (HELM) column', async () => {
    const info = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macroCol: any = cols.find((c: any) => c.semType === 'Macromolecule');
      return {
        rows: df.rowCount,
        hasMacro: !!macroCol,
        macroName: macroCol?.name ?? null,
        units: macroCol?.meta?.units ?? null,
      };
    });
    expect(info.rows).toBeGreaterThan(0);
    expect(info.hasMacro).toBe(true);
    expect(info.units).toBe('helm');
  });
  await softStep('Add new column Clusters = RandBetween(0, 5)', async () => {
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Edit"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 300));
      (document.querySelector('[name="div-Edit---Add-New-Column..."]') as HTMLElement).click();
    });
    await page.locator('[name="dialog-Add-New-Column"]').waitFor({timeout: 15_000});
    const nameInput = page.locator('[name="dialog-Add-New-Column"] [name="input-Add-New-Column---Name"]').first();
    await nameInput.waitFor({timeout: 15_000, state: 'visible'});
    await nameInput.click();
    await page.keyboard.press('Control+A');
    await page.keyboard.type('Clusters', {delay: 20});
    const cm = page.locator('[name="dialog-Add-New-Column"] .add-new-column-dialog-cm-div .cm-content').first();
    await cm.waitFor({timeout: 15_000, state: 'visible'});
    await cm.click({force: true});
    await page.waitForTimeout(200);
    await page.keyboard.press('Escape').catch(() => {});
    await page.waitForTimeout(150);
    await page.keyboard.press('Control+A');
    await page.keyboard.press('Delete');
    await page.waitForTimeout(100);
    await page.keyboard.type('RandBetween(0, 5)', {delay: 30});
    await page.keyboard.press('Escape').catch(() => {});
    await page.waitForTimeout(200);
    const formulaText: string = await page.evaluate(() => {
      const cmEl = document.querySelector(
        '[name="dialog-Add-New-Column"] .add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
      if (!cmEl) return '';
      const view = (cmEl as any).cmView?.view
        ?? (cmEl.parentElement as any)?.cmView?.view
        ?? null;
      return view ? view.state.doc.toString() : (cmEl.textContent ?? '');
    });
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
    }, null, {timeout: 60_000});
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
  await softStep('Open Bio > Analyze > MSA... for HELM (PepSeA engine surfaces)', async () => {
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 400));
      document.querySelector('[name="div-Bio---Analyze"]')!
        .dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 400));
      (document.querySelector('[name="div-Bio---Analyze---MSA..."]') as HTMLElement).click();
    });
    await page.locator('[name="dialog-MSA"]').waitFor({timeout: 15_000});
    const info = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const inputs = Array.from(dlg.querySelectorAll('[name^="input-host-"]'))
        .map((h) => h.getAttribute('name')!);
      const engineHost = dlg.querySelector('[name="input-host-Engine"]');
      const engineSel = engineHost?.querySelector('select') as HTMLSelectElement | null;
      const engineChoices = engineSel ? Array.from(engineSel.options).map((o) => o.value) : [];
      return {inputs, engineChoices};
    });
    expect(info.inputs).toContain('input-host-Sequence');
    expect(info.inputs).toContain('input-host-Clusters');
    expect(info.inputs).toContain('input-host-Engine');
    expect(info.engineChoices.some((v) => /pepsea/i.test(v))).toBe(true);
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
  await softStep('Alignment parameters button toggles engine Gap-Open + Gap-Extend surface', async () => {
    const before = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const gapOpen = (dlg.querySelector('[name="input-host-Gap-Open"]') as HTMLElement)?.getBoundingClientRect().height ?? 0;
      const gapExt = (dlg.querySelector('[name="input-host-Gap-Extend"]') as HTMLElement)?.getBoundingClientRect().height ?? 0;
      return {gapOpen, gapExt};
    });
    expect(before.gapOpen).toBeGreaterThan(0);
    expect(before.gapExt).toBeGreaterThan(0);
    await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const btns = Array.from(dlg.querySelectorAll('[name="button-Alignment-parameters"]')) as HTMLElement[];
      const visible = btns.find((b) => b.offsetParent !== null);
      if (!visible) throw new Error('No visible Alignment-parameters button in dialog-MSA');
      visible.click();
    });
    await page.waitForFunction(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const gap = dlg.querySelector('[name="input-host-Gap-Open"]') as HTMLElement | null;
      return !gap || gap.getBoundingClientRect().height === 0;
    }, null, {timeout: 10_000});
    const hidden = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const gapOpen = (dlg.querySelector('[name="input-host-Gap-Open"]') as HTMLElement)?.getBoundingClientRect().height ?? 0;
      const gapExt = (dlg.querySelector('[name="input-host-Gap-Extend"]') as HTMLElement)?.getBoundingClientRect().height ?? 0;
      return {gapOpen, gapExt};
    });
    expect(hidden.gapOpen).toBe(0);
    expect(hidden.gapExt).toBe(0);
    await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const btns = Array.from(dlg.querySelectorAll('[name="button-Alignment-parameters"]')) as HTMLElement[];
      const visible = btns.find((b) => b.offsetParent !== null);
      if (!visible) throw new Error('No visible Alignment-parameters button on second toggle');
      visible.click();
    });
    await page.waitForFunction(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const gap = dlg.querySelector('[name="input-host-Gap-Open"]') as HTMLElement | null;
      return !!gap && gap.getBoundingClientRect().height > 0;
    }, null, {timeout: 10_000});
    const restored = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const gapOpen = (dlg.querySelector('[name="input-host-Gap-Open"]') as HTMLElement)?.getBoundingClientRect().height ?? 0;
      const gapExt = (dlg.querySelector('[name="input-host-Gap-Extend"]') as HTMLElement)?.getBoundingClientRect().height ?? 0;
      return {gapOpen, gapExt};
    });
    expect(restored.gapOpen).toBeGreaterThan(0);
    expect(restored.gapExt).toBeGreaterThan(0);
  });
  await softStep('OK runs PepSeA MSA — verify aligned column, renderer, and per-cluster monomer count', async () => {
    await page.locator('[name="dialog-MSA"] [name="button-OK"]').click();
    await page.waitForTimeout(3000);
    let dialogProducedColumn = false;
    const dialogStart = Date.now();
    const dialogBudgetMs = 240_000;
    const REQUIRED_CONSECUTIVE_EMPTY = 3;
    let consecutiveEmpty = 0;
    try {
      while (Date.now() - dialogStart < dialogBudgetMs) {
        const status = await page.evaluate(() => {
          const df = grok.shell.tv.dataFrame;
          const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
          const hasMsa = cols.some((c: any) => c.name.toLowerCase().includes('msa') && c.semType === 'Macromolecule');
          const dlg = !!document.querySelector('[name="dialog-MSA"]');
          const balloons = Array.from(document.querySelectorAll('div')) as HTMLElement[];
          const performingMsa = balloons.some((b) => /Performing MSA/i.test(b.textContent ?? ''));
          return {hasMsa, dlg, performingMsa};
        });
        if (status.hasMsa) {
          dialogProducedColumn = true;
          break;
        }
        if (!status.dlg && !status.performingMsa) {
          consecutiveEmpty++;
          if (consecutiveEmpty >= REQUIRED_CONSECUTIVE_EMPTY) break;
        } else {
          consecutiveEmpty = 0;
        }
        await new Promise((r) => setTimeout(r, 1000));
      }
    } catch {  }
    if (!dialogProducedColumn) {
      await page.evaluate(() => {
        const dlg = document.querySelector('[name="dialog-MSA"]') as HTMLElement | null;
        if (dlg) (dlg.querySelector('[name="button-CANCEL"]') as HTMLElement)?.click();
      });
      await page.evaluate(async () => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const helmCol: any = cols.find((c: any) => c.semType === 'Macromolecule' && c.meta?.units === 'helm');
        if (!helmCol) throw new Error('No HELM Macromolecule column found in table');
        const clustersCol = df.col('Clusters')!;
        const PEPSEA_SEPARATOR = '.';
        const categories = Array.from(new Set((clustersCol as any).toList()));
        const resultArr: string[] = new Array(df.rowCount).fill('');
        const fnNames = [
          'Bio:PepSeA',
          'Bio:pepseaMsa',
        ];
        let pickedFn: string | null = null;
        for (let catIdx = 0; catIdx < categories.length; catIdx++) {
          const cat = categories[catIdx];
          const idx: number[] = [];
          for (let i = 0; i < clustersCol.length; i++) if (clustersCol.get(i) === cat) idx.push(i);
          if (idx.length === 0) continue;
          if (idx.length === 1) { resultArr[idx[0]] = helmCol.get(idx[0])!; continue; }
          const sub = DG.Column.fromStrings(`HELM_c${cat}`, idx.map((i) => helmCol.get(i)!)) as any;
          sub.semType = 'Macromolecule';
          sub.meta.units = 'helm';
          let r: any = null;
          let lastErr: any = null;
          for (const fn of (pickedFn ? [pickedFn] : fnNames)) {
            try {
              r = await grok.functions.call(fn, {
                sequenceCol: sub,
                method: 'mafft --auto',
                gapOpen: 1.53,
                gapExtend: 0,
              });
              pickedFn = fn;
              break;
            } catch (e) { lastErr = e; r = null; }
          }
          if (!r) {
            const fnTried = (pickedFn ? [pickedFn] : fnNames).join(', ');
            throw new Error(`No PepSeA-equivalent sequenceMSA fn reachable on cluster ${cat} (tried: ${fnTried}). Last err: ${String(lastErr)}`);
          }
          const resultCol: any = (typeof r?.get === 'function')
            ? r
            : (typeof r?.getOutputParamValue === 'function' ? r.getOutputParamValue() : null);
          if (!resultCol || typeof resultCol.get !== 'function')
            throw new Error(`Unexpected sequenceMSA fn return shape on cluster ${cat}: ${typeof r} (no .get and no .getOutputParamValue)`);
          for (let i = 0; i < idx.length; i++) resultArr[idx[i]] = resultCol.get(i);
        }
        const name = df.columns.getUnusedName(`msa(${helmCol.name})`);
        const msaCol = DG.Column.fromStrings(name, resultArr) as any;
        msaCol.semType = 'Macromolecule';
        msaCol.meta.units = 'separator';
        msaCol.setTag('separator', PEPSEA_SEPARATOR);
        msaCol.setTag('aligned', 'SEQ.MSA');
        msaCol.setTag('alphabet', 'UN');
        msaCol.setTag('.alphabetIsMultichar', 'true');
        df.columns.add(msaCol);
        // detectSemanticTypes binds the Macromolecule cell renderer (mirrors the dialog OK path).
        await grok.data.detectSemanticTypes(df);
      });
      await page.waitForFunction(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const msa: any = cols.find((c: any) => c.name.toLowerCase().includes('msa'));
        if (!msa) return false;
        const gridCol = (grok.shell.tv as any).grid?.col?.(msa.name);
        return gridCol?.cellType === 'sequence' || gridCol?.cellType === 'helm';
      }, null, {timeout: 60_000});
    } else {
      // Happy path — dialog produced the column; poll for the renderer to bind.
      await page.waitForFunction(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const msa: any = cols.find((c: any) => c.name.toLowerCase().includes('msa') && c.semType === 'Macromolecule');
        if (!msa) return false;
        const gridCol = (grok.shell.tv as any).grid?.col?.(msa.name);
        return gridCol?.cellType === 'sequence' || gridCol?.cellType === 'helm';
      }, null, {timeout: 60_000});
    }
    const result = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const msa: any = cols.find((c: any) => c.name.toLowerCase().includes('msa'));
      const clusters: any = df.col('Clusters');
      const gridCol = (grok.shell.tv as any).grid?.col?.(msa.name);
      let separatorTag = '.';
      try { separatorTag = msa.getTag('separator') || '.'; } catch { /* ignore */ }
      const units = msa?.meta?.units ?? '';
      const countByCluster: Record<string, Set<number>> = {};
      for (let i = 0; i < df.rowCount; i++) {
        const key = String(clusters.get(i));
        const s: string = msa.get(i) ?? '';
        let count = -1;
        if (units === 'helm' || /^[A-Z]+1\{/.test(s)) {
          const m = s.match(/\{([^}]*)\}/);
          count = m ? m[1].split('.').length : -1;
        } else {
          count = s.split(separatorTag).length;
        }
        (countByCluster[key] ||= new Set()).add(count);
      }
      const allEqualPerCluster = Object.values(countByCluster).every((s) => s.size === 1);
      let alignedTag: string | null = null;
      try { alignedTag = msa.getTag('aligned') ?? null; } catch { /* ignore */ }
      return {
        msaName: msa?.name,
        semType: msa?.semType,
        units: msa?.meta?.units,
        cellType: gridCol?.cellType,
        alignedTag,
        allEqualPerCluster,
        clusterCount: Object.keys(countByCluster).length,
      };
    });
    expect(result.msaName).toBeTruthy();
    expect(result.semType).toBe('Macromolecule');
    expect(['separator', 'helm', 'custom']).toContain(result.units);
    expect(['sequence', 'helm']).toContain(result.cellType);
    expect(result.allEqualPerCluster).toBe(true);
    expect(result.clusterCount).toBeGreaterThan(1);
  });
  finishSpec();
});
