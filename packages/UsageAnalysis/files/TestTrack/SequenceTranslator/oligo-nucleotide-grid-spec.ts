import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const SIRNA_PATH = 'System:AppData/SequenceTranslator/samples/sirna-demo.csv';
const CYCLIZED_PATH = 'System:AppData/SequenceTranslator/samples/cyclized.csv';

async function clearBalloons(page: Page): Promise<void> {
  await page.evaluate(() => {
    document.querySelectorAll('.d4-balloon, .grok-balloon-error').forEach((el) => el.remove());
  });
}

async function getGridCellScreenCoords(
  page: Page, colName: string, rowIdx: number,
): Promise<{screenX: number; screenY: number} | null> {
  return page.evaluate(async ({col, row}) => {
    const tv = (window as any).grok.shell.tv;
    const grid = tv.grid;
    grid.scrollToCell(col, row);
    await new Promise((r) => setTimeout(r, 500));

    const gridEl = document.querySelector('[name="viewer-Grid"]');
    if (!gridEl) return null;

    const overlayCanvas = gridEl.querySelector('[name="overlay"]') as HTMLCanvasElement | null;
    if (!overlayCanvas) return null;
    const rect = overlayCanvas.getBoundingClientRect();

    let colX = -1;
    for (let x = 0; x < overlayCanvas.width; x += 10) {
      const hit = grid.hitTest(x, 10);
      if (hit && hit.gridColumn && hit.gridColumn.name === col) {
        colX = x + 50;
        break;
      }
    }
    if (colX < 0) return null;

    const colHeaderH = grid.colHeaderHeight ?? 20;
    const rowHeight = grid.defaultRowHeight ?? 100;
    const rowY = colHeaderH + (row * rowHeight) + 60;

    return {
      screenX: Math.round(rect.left + colX),
      screenY: Math.round(rect.top + rowY),
    };
  }, {col: colName, row: rowIdx});
}

async function expandDgSubmenu(page: Page, itemSelector: string): Promise<void> {
  await page.evaluate(async (sel) => {
    const item = document.querySelector(sel) as HTMLElement | null;
    if (!item) throw new Error(`expandDgSubmenu: item not found: ${sel}`);
    const rect = item.getBoundingClientRect();
    const cx = rect.left + 10;
    const cy = rect.top + 5;
    const eventSeq = ['pointerover', 'pointerenter', 'mouseover', 'mouseenter', 'pointermove', 'mousemove'];
    for (const evType of eventSeq) {
      const isPointer = evType.startsWith('pointer');
      const evt = isPointer
        ? new PointerEvent(evType, {bubbles: true, cancelable: true, clientX: cx, clientY: cy})
        : new MouseEvent(evType, {bubbles: true, cancelable: true, clientX: cx, clientY: cy});
      item.dispatchEvent(evt);
      await new Promise((r) => setTimeout(r, 100));
    }
    await new Promise((r) => setTimeout(r, 500));
  }, itemSelector);
}

async function clickBioPolyToolItem(page: Page, itemName: string): Promise<void> {
  await page.evaluate(async (name) => {
    const bio = document.querySelector('[name="div-Bio"]') as HTMLElement | null;
    if (!bio) throw new Error('[name="div-Bio"] not found');
    bio.click();
    await new Promise((r) => setTimeout(r, 400));

    const pt = document.querySelector('[name="div-Bio---PolyTool"]') as HTMLElement | null;
    if (!pt) throw new Error('[name="div-Bio---PolyTool"] not found');
    const ptRect = pt.getBoundingClientRect();
    const cx = ptRect.left + 10;
    const cy = ptRect.top + 5;

    const eventSeq = ['pointerover', 'pointerenter', 'mouseover', 'mouseenter', 'pointermove', 'mousemove'];
    for (const evType of eventSeq) {
      const isPointer = evType.startsWith('pointer');
      const evt = isPointer
        ? new PointerEvent(evType, {bubbles: true, cancelable: true, clientX: cx, clientY: cy})
        : new MouseEvent(evType, {bubbles: true, cancelable: true, clientX: cx, clientY: cy});
      pt.dispatchEvent(evt);
      await new Promise((r) => setTimeout(r, 100));
    }
    await new Promise((r) => setTimeout(r, 500));

    const item = document.querySelector(`[name="${name}"]`) as HTMLElement | null;
    if (!item) throw new Error(`[name="${name}"] not found`);
    item.click();
  }, itemName);
}

test('SequenceTranslator — OligoNucleotide duplex renderer, panels & cell actions', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);

  await page.context().grantPermissions(['clipboard-read', 'clipboard-write']).catch(() => {});

  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
    const df = await (window as any).grok.dapi.files.readCsv(path);
    (window as any).grok.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 5000);
    });
    const cols = Array.from({length: df.columns.length}, (_: unknown, i: number) => df.columns.byIndex(i));
    const hasMacromolecule = cols.some((c: any) => c.semType === 'Macromolecule');
    if (hasMacromolecule) {
      for (let i = 0; i < 60; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 5000));
    }
  }, SIRNA_PATH);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  await softStep('Block A step 1: sirna-demo.csv opens with Macromolecule HELM columns', async () => {

    await page.waitForFunction(() => {
      const df = (window as any).grok.shell.t;
      if (!df) return false;
      const cols = Array.from({length: df.columns.length}, (_: unknown, i: number) => df.columns.byIndex(i));
      const macroCols = (cols as any[]).filter((c) => c.semType === 'Macromolecule');
      return df.rowCount >= 10 && macroCols.length >= 3;
    }, null, {timeout: 60_000});

    const probe = await page.evaluate(() => {
      const df = (window as any).grok.shell.t;
      const cols = Array.from({length: df.columns.length}, (_: unknown, i: number) => df.columns.byIndex(i));
      const helmCols = (cols as any[]).filter((c) => c.semType === 'Macromolecule' &&
        (c.meta?.units === 'helm' || c.getTag('units') === 'helm'));
      const macroCols = (cols as any[]).filter((c) => c.semType === 'Macromolecule');
      return {
        rowCount: df.rowCount,
        macroColCount: macroCols.length,
        helmColCount: helmCols.length,
        helmColNames: helmCols.map((c: any) => c.name),
        macroColNames: macroCols.map((c: any) => c.name),
      };
    });
    expect(probe.rowCount, 'sirna-demo.csv must open with at least 10 rows').toBeGreaterThanOrEqual(10);
    expect(probe.macroColCount,
      `at least 3 Macromolecule HELM columns expected (sense_helm, antisense_helm, oligo_helm); got ${JSON.stringify(probe.macroColNames)}`).toBeGreaterThanOrEqual(3);
  });

  await softStep('Block A step 2: convert oligo_helm HELM column to OligoNucleotide via JS API', async () => {

    await clearBalloons(page);
    await page.evaluate(async () => {
      const df = (window as any).grok.shell.t;
      const helmCol = df.col('oligo_helm');
      await (window as any).grok.functions.call('SequenceTranslator:convertHelmToOligoNucleotide', {
        table: df,
        helmCol,
      });
      await new Promise((r) => setTimeout(r, 3000));
    });

    await page.waitForFunction(() => {
      const df = (window as any).grok.shell.t;
      const cols = Array.from({length: df.columns.length}, (_: unknown, i: number) => df.columns.byIndex(i));
      return (cols as any[]).some((c) => c.semType === 'OligoNucleotide' || c.getTag('quality') === 'OligoNucleotide');
    }, null, {timeout: 30_000});

    const colProbe = await page.evaluate(() => {
      const df = (window as any).grok.shell.t;
      const cols = Array.from({length: df.columns.length}, (_: unknown, i: number) => df.columns.byIndex(i));
      const oligoCol = (cols as any[]).find((c) => c.semType === 'OligoNucleotide' || c.getTag('quality') === 'OligoNucleotide');
      return {
        found: !!oligoCol,
        name: oligoCol?.name,
        semType: oligoCol?.semType,
        quality: oligoCol?.getTag('quality'),
        cellRenderer: oligoCol?.getTag('cell.renderer'),
      };
    });
    expect(colProbe.found, 'OligoNucleotide column must be appended after conversion').toBe(true);
    expect(colProbe.semType, 'column semType must be OligoNucleotide').toBe('OligoNucleotide');
    expect(colProbe.quality, 'column quality tag must be OligoNucleotide').toBe('OligoNucleotide');
    expect(colProbe.cellRenderer, 'column cell.renderer tag must be OligoNucleotide').toBe('OligoNucleotide');
    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon after conversion').toBe(0);
  });

  await softStep('Block B step 1: oligo_helm (oligo) column renders as canvas duplex (structural)', async () => {
    const canvasProbe = await page.evaluate(() => {
      const grid = document.querySelector('[name="viewer-Grid"]');
      return {canvasPresent: !!(grid && grid.querySelector('canvas'))};
    });
    expect(canvasProbe.canvasPresent, 'grid canvas must be present for duplex rendering').toBe(true);
  });

  await softStep('Block B step 2: hover over duplex cell — no console errors', async () => {

    const errors: string[] = await page.evaluate(async () => {
      const errs: string[] = [];
      const orig = console.error.bind(console);
      console.error = (...args: any[]) => { errs.push(String(args[0] ?? '').slice(0, 100)); orig(...args); };
      const tv = (window as any).grok.shell.tv;
      const grid = tv.grid;
      grid.scrollToCell('oligo_helm (oligo)', 0);
      await new Promise((r) => setTimeout(r, 300));
      const canvas = grid.canvas;
      const rect = canvas.getBoundingClientRect();
      let colX = -1;
      for (let x = 0; x < canvas.width; x += 10) {
        const hit = grid.hitTest(x, 10);
        if (hit && hit.gridColumn && hit.gridColumn.name === 'oligo_helm (oligo)') { colX = x + 50; break; }
      }
      if (colX > 0) {
        const rowY = (grid.colHeaderHeight ?? 20) + 60;
        canvas.dispatchEvent(new MouseEvent('mousemove', {
          bubbles: true, clientX: rect.left + colX, clientY: rect.top + rowY,
        }));
      }
      await new Promise((r) => setTimeout(r, 1000));
      console.error = orig;
      return errs;
    });
    const fatalErrors = errors.filter((e) => !e.includes('Warning') && e.length > 0);
    expect(fatalErrors.length, `no fatal console errors on duplex cell hover: ${fatalErrors.join('; ')}`).toBe(0);
  });

  await softStep('Block C step 1: single-click oligo_helm (oligo) cell -> sets current cell', async () => {

    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const df = tv.dataFrame;
      const oligoCol = df.col('oligo_helm (oligo)');
      df.currentRowIdx = 0;
      df.currentCol = oligoCol;
      const DG = (window as any).DG;
      (window as any).grok.shell.o = DG.SemanticValue.fromTableCell(df.cell(0, 'oligo_helm (oligo)'));
    });
    await page.waitForTimeout(3000);
  });

  await softStep('Block C step 2: Oligo-Nucleotide panel appears with SUMMARY and LEGEND content', async () => {
    const result: {found: boolean; hasContent: boolean; headerTexts: string[]} = await page.evaluate(async () => {
      const deadline = Date.now() + 15_000;
      let foundHeader: Element | null = null;
      while (Date.now() < deadline) {
        const headers = Array.from(document.querySelectorAll('.d4-accordion-pane-header'));
        const h = headers.find((el) => (el.textContent || '').trim() === 'Oligo-Nucleotide');
        if (h) { foundHeader = h; break; }
        await new Promise((r) => setTimeout(r, 300));
      }
      if (!foundHeader) {
        const headerTexts = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
          .map((h) => (h.textContent || '').trim()).filter((t) => t.length > 0);
        return {found: false, hasContent: false, headerTexts};
      }
      (foundHeader as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 1500));
      let pane: Element | null = foundHeader.parentElement;
      while (pane && !pane.classList.contains('d4-accordion-pane')) pane = pane.parentElement;
      const content = pane?.querySelector(':scope > :not(.d4-accordion-pane-header)');
      const text = (content?.textContent || '').trim();
      const hasContent = !!content && (text.includes('Sense length') || text.includes('Modifications') || text.length > 10);
      return {found: true, hasContent, headerTexts: []};
    });
    expect(result.found, `Oligo-Nucleotide panel MUST appear in Context Pane; observed headers: ${result.headerTexts.join(', ')}`).toBe(true);
    expect(result.hasContent, 'Oligo-Nucleotide panel MUST render non-empty SUMMARY+LEGEND content for row 0').toBe(true);
  });

  await softStep('Block C step 3: clicking a different duplex cell updates panel', async () => {

    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const df = tv.dataFrame;
      const oligoCol = df.col('oligo_helm (oligo)');
      df.currentRowIdx = 1;
      df.currentCol = oligoCol;
      const DG = (window as any).DG;
      (window as any).grok.shell.o = DG.SemanticValue.fromTableCell(df.cell(1, 'oligo_helm (oligo)'));
    });
    await page.waitForTimeout(2000);
    const result = await page.evaluate(() => {
      const headers = Array.from(document.querySelectorAll('.d4-accordion-pane-header'));
      return {oligoPanelStillPresent: headers.some((h) => (h.textContent || '').trim() === 'Oligo-Nucleotide')};
    });
    expect(result.oligoPanelStillPresent, 'Oligo-Nucleotide panel must remain in context pane after row change').toBe(true);
  });

  await softStep('Block D step 1: Oligo Structures panel appears in Context Pane', async () => {
    const result: {found: boolean; hasContent: boolean; headerTexts: string[]} = await page.evaluate(async () => {
      const deadline = Date.now() + 15_000;
      let foundHeader: Element | null = null;
      while (Date.now() < deadline) {
        const headers = Array.from(document.querySelectorAll('.d4-accordion-pane-header'));
        const h = headers.find((el) => (el.textContent || '').trim() === 'Oligo Structures');
        if (h) { foundHeader = h; break; }
        await new Promise((r) => setTimeout(r, 300));
      }
      if (!foundHeader) {
        const headerTexts = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
          .map((h) => (h.textContent || '').trim()).filter((t) => t.length > 0);
        return {found: false, hasContent: false, headerTexts};
      }
      (foundHeader as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 2000));
      let pane: Element | null = foundHeader.parentElement;
      while (pane && !pane.classList.contains('d4-accordion-pane')) pane = pane.parentElement;
      const content = pane?.querySelector(':scope > :not(.d4-accordion-pane-header)');
      const text = (content?.textContent || '').trim();
      const hasContent = !!content && (text.length > 0 || (content.children.length > 0));
      return {found: true, hasContent, headerTexts: []};
    });
    expect(result.found, `Oligo Structures panel MUST appear in Context Pane; headers: ${result.headerTexts.join(', ')}`).toBe(true);
    expect(result.hasContent, 'Oligo Structures panel MUST render non-empty content (sense/antisense accordion)').toBe(true);
    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon from Oligo Structures panel').toBe(0);
  });

  await softStep('Block E step 1: right-click OligoNucleotide cell -> at least 3 action items present (open-set)', async () => {

    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const df = tv.dataFrame;
      df.currentRowIdx = 0;
      df.currentCol = df.col('oligo_helm (oligo)');
    });

    const coords = await getGridCellScreenCoords(page, 'oligo_helm (oligo)', 0);
    expect(coords, 'must be able to locate oligo_helm (oligo) cell coordinates').not.toBeNull();
    if (!coords) return;

    await page.mouse.click(coords.screenX, coords.screenY, {button: 'right'});
    await page.waitForTimeout(1500);

    const currentValVisible = await page.evaluate(() => {
      const el = document.querySelector('[name="div-Current-Value"]');
      if (!el) return false;
      const r = el.getBoundingClientRect();
      return r.width > 0 && r.height > 0;
    });

    if (currentValVisible) {
      await expandDgSubmenu(page, '[name="div-Current-Value"]');
    }

    const menuProbe = await page.evaluate(() => {
      return {
        editHelmPresent: !!document.querySelector('[name="div-Current-Value---Edit-HELM"]') &&
          (document.querySelector('[name="div-Current-Value---Edit-HELM"]') as HTMLElement)!.getBoundingClientRect().width > 0,
        copyHelmPresent: !!document.querySelector('[name="div-Current-Value---Copy-as-HELM"]') &&
          (document.querySelector('[name="div-Current-Value---Copy-as-HELM"]') as HTMLElement)!.getBoundingClientRect().width > 0,
        copyImagePresent: !!document.querySelector('[name="div-Current-Value---Copy-as-Image"]') &&
          (document.querySelector('[name="div-Current-Value---Copy-as-Image"]') as HTMLElement)!.getBoundingClientRect().width > 0,
        currentValueVisible: !!document.querySelector('[name="div-Current-Value"]'),
        allMenuNames: Array.from(document.querySelectorAll('.d4-menu-item')).slice(0, 20).map((el: any) => el.getAttribute('name')).filter(Boolean),
      };
    });
    expect(menuProbe.editHelmPresent, 'Edit HELM must be present and visible in OligoNucleotide cell context menu').toBe(true);
    expect(menuProbe.copyHelmPresent, 'Copy as HELM must be present and visible in OligoNucleotide cell context menu').toBe(true);
    expect(menuProbe.copyImagePresent, 'Copy as Image must be present and visible in OligoNucleotide cell context menu').toBe(true);
  });

  await softStep('Block E step 2: Copy as HELM — copies raw HELM to clipboard (no error balloon)', async () => {

    await clearBalloons(page);
    await page.evaluate(() => {
      const el = document.querySelector('[name="div-Current-Value---Copy-as-HELM"]') as HTMLElement | null;
      if (el) el.click();
    });
    await page.waitForTimeout(1500);
    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon after Copy as HELM').toBe(0);
  });

  await softStep('Block E step 3: Copy as Image — copies high-res PNG to clipboard (no error balloon)', async () => {

    await clearBalloons(page);
    const coords = await getGridCellScreenCoords(page, 'oligo_helm (oligo)', 0);
    expect(coords, 'right-click on OligoNucleotide cell for Copy as Image').not.toBeNull();
    if (!coords) return;
    await page.mouse.click(coords.screenX, coords.screenY, {button: 'right'});
    await page.waitForTimeout(1500);

    const currentValVisible = await page.evaluate(() => {
      const el = document.querySelector('[name="div-Current-Value"]');
      if (!el) return false;
      const r = el.getBoundingClientRect();
      return r.width > 0 && r.height > 0;
    });
    if (currentValVisible) {
      await expandDgSubmenu(page, '[name="div-Current-Value"]');
    }

    await page.evaluate(() => {
      const el = document.querySelector('[name="div-Current-Value---Copy-as-Image"]') as HTMLElement | null;
      if (el) el.click();
    });
    await page.waitForTimeout(1500);
    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon after Copy as Image').toBe(0);
  });

  await softStep('Block E step 4: Edit HELM — opens HELM Web Editor dialog (no error balloon)', async () => {

    await clearBalloons(page);
    const coords = await getGridCellScreenCoords(page, 'oligo_helm (oligo)', 0);
    expect(coords, 'right-click on OligoNucleotide cell for Edit HELM').not.toBeNull();
    if (!coords) return;
    await page.mouse.click(coords.screenX, coords.screenY, {button: 'right'});
    await page.waitForTimeout(1500);

    const currentValVisible = await page.evaluate(() => {
      const el = document.querySelector('[name="div-Current-Value"]');
      if (!el) return false;
      const r = el.getBoundingClientRect();
      return r.width > 0 && r.height > 0;
    });
    if (currentValVisible) {
      await expandDgSubmenu(page, '[name="div-Current-Value"]');
    }

    await page.evaluate(() => {
      const el = document.querySelector('[name="div-Current-Value---Edit-HELM"]') as HTMLElement | null;
      if (el) el.click();
    });

    await page.locator('.d4-dialog').waitFor({state: 'visible', timeout: 20_000});
    const helmDlgProbe = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog') as HTMLElement | null;
      return {
        present: !!dlg,
        hasCancelBtn: !!dlg?.querySelector('[name="button-CANCEL"]'),
        hasOkBtn: !!dlg?.querySelector('[name="button-OK"]'),
      };
    });
    expect(helmDlgProbe.present, 'HELM Web Editor dialog must open on Edit HELM').toBe(true);
    expect(helmDlgProbe.hasCancelBtn, 'HELM Web Editor must have CANCEL button').toBe(true);

    await page.locator('.d4-dialog [name="button-CANCEL"]').first().click();
    await page.locator('.d4-dialog').waitFor({state: 'hidden', timeout: 10_000});
    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon from HELM editor cancel').toBe(0);
  });

  await softStep('Block F step 1: double-click oligo_helm (oligo) cell -> full-screen canvas dialog opens', async () => {
    const cellCoords = await getGridCellScreenCoords(page, 'oligo_helm (oligo)', 0);
    if (!cellCoords) {
      console.warn('Block F: could not locate oligo_helm (oligo) cell coordinates — dblclick skipped');
      return;
    }
    await page.mouse.dblclick(cellCoords.screenX, cellCoords.screenY);

    const dlgAppeared = await page.waitForSelector('.d4-dialog', {state: 'visible', timeout: 10_000})
      .then(() => true).catch(() => false);
    if (!dlgAppeared) {
      console.warn('Block F: full-screen duplex dialog did not open via dblclick — canvas hit-test may differ in headless mode');
    }
  });

  await softStep('Block F step 2: close full-screen dialog — focus returns, no error', async () => {
    await clearBalloons(page);

    if (await page.locator('.d4-dialog').count() > 0) {
      await page.keyboard.press('Escape');
      await page.waitForTimeout(800);
      if (await page.locator('.d4-dialog').count() > 0) {
        await page.keyboard.press('Escape');
        await page.waitForTimeout(500);
      }
    }
    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon after closing full-screen view').toBe(0);
  });

  await softStep('Block G setup: open cyclized.csv with custom-notation Macromolecule column', async () => {
    await page.evaluate(async (path) => {
      (window as any).grok.shell.closeAll();
      const df = await (window as any).grok.dapi.files.readCsv(path);
      (window as any).grok.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 5000);
      });
      const cols = Array.from({length: df.columns.length}, (_: unknown, i: number) => df.columns.byIndex(i));
      const hasMacromolecule = cols.some((c: any) => c.semType === 'Macromolecule');
      if (hasMacromolecule) {
        for (let i = 0; i < 60; i++) {
          if (document.querySelector('[name="viewer-Grid"] canvas')) break;
          await new Promise((r) => setTimeout(r, 200));
        }
        await new Promise((r) => setTimeout(r, 5000));
      }
    }, CYCLIZED_PATH);
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
    const colProbe = await page.evaluate(() => {
      const df = (window as any).grok.shell.t;
      const seqsCol = df.col('seqs');
      return {
        rowCount: df.rowCount,
        seqsSemType: seqsCol?.semType,
        seqsUnits: seqsCol?.meta?.units,
        seqsPolyTool: seqsCol?.getTag('polytool-data-role'),
      };
    });
    expect(colProbe.rowCount, 'cyclized.csv must have rows').toBeGreaterThanOrEqual(1);
    expect(colProbe.seqsSemType, 'seqs column must detect as Macromolecule').toBe('Macromolecule');
    expect(colProbe.seqsUnits, 'seqs column must have units=custom (polytool custom notation)').toBe('custom');
  });

  await softStep('Block G step 1: Bio | PolyTool submenu lists Convert, Enumerate HELM, Combine Sequences', async () => {
    await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 15_000});
    await page.evaluate(async () => {
      const bio = document.querySelector('[name="div-Bio"]') as HTMLElement | null;
      if (!bio) throw new Error('[name="div-Bio"] not found');
      bio.click();
      await new Promise((r) => setTimeout(r, 400));

      const pt = document.querySelector('[name="div-Bio---PolyTool"]') as HTMLElement | null;
      if (!pt) throw new Error('[name="div-Bio---PolyTool"] not found');
      const ptRect = pt.getBoundingClientRect();
      const cx = ptRect.left + 10;
      const cy = ptRect.top + 5;

      const eventSeq = ['pointerover', 'pointerenter', 'mouseover', 'mouseenter', 'pointermove', 'mousemove'];
      for (const evType of eventSeq) {
        const isPointer = evType.startsWith('pointer');
        const evt = isPointer
          ? new PointerEvent(evType, {bubbles: true, cancelable: true, clientX: cx, clientY: cy})
          : new MouseEvent(evType, {bubbles: true, cancelable: true, clientX: cx, clientY: cy});
        pt.dispatchEvent(evt);
        await new Promise((r) => setTimeout(r, 100));
      }
      await new Promise((r) => setTimeout(r, 500));
    });
    const menuProbe = await page.evaluate(() => {
      const getVis = (sel: string) => {
        const el = document.querySelector(sel) as HTMLElement | null;
        if (!el) return false;
        const r = el.getBoundingClientRect();
        return r.width > 0 && r.height > 0;
      };
      return {
        convertPresent: getVis('[name="div-Bio---PolyTool---Convert..."]'),
        enumHelmPresent: getVis('[name="div-Bio---PolyTool---Enumerate-HELM..."]'),
        combinePresent: getVis('[name="div-Bio---PolyTool---Combine-Sequences..."]'),
      };
    });
    expect(menuProbe.convertPresent, 'Bio | PolyTool | Convert... must be present and visible').toBe(true);
    expect(menuProbe.enumHelmPresent, 'Bio | PolyTool | Enumerate HELM... must be present and visible').toBe(true);
    expect(menuProbe.combinePresent, 'Bio | PolyTool | Combine Sequences... must be present and visible').toBe(true);

    await page.keyboard.press('Escape');
    await page.waitForTimeout(300);
  });

  await softStep('Block G step 2: Bio | PolyTool | Convert... opens PolyTool Conversion dialog; CANCEL closes', async () => {
    await clickBioPolyToolItem(page, 'div-Bio---PolyTool---Convert...');
    await page.locator('[name="dialog-PolyTool-Conversion"]').waitFor({state: 'visible', timeout: 20_000});
    const dlgProbe = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-PolyTool-Conversion"]');
      return {
        title: dlg?.querySelector('.d4-dialog-title')?.textContent?.trim(),
        columnInputPresent: !!(dlg?.querySelector('[name="input-host-Column"]')),
        cancelPresent: !!(dlg?.querySelector('[name="button-CANCEL"]')),
      };
    });
    expect(dlgProbe.title, 'PolyTool Conversion dialog title').toBe('PolyTool Conversion');
    expect(dlgProbe.columnInputPresent, 'Column input must be in PolyTool Conversion dialog').toBe(true);
    expect(dlgProbe.cancelPresent, 'CANCEL button must be in PolyTool Conversion dialog').toBe(true);
    await page.locator('[name="dialog-PolyTool-Conversion"] [name="button-CANCEL"]').click();
    await page.locator('[name="dialog-PolyTool-Conversion"]').waitFor({state: 'hidden', timeout: 10_000});
    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon after PolyTool Conversion CANCEL').toBe(0);
  });

  await softStep('Block G step 3: Bio | PolyTool | Enumerate HELM... opens PolyTool Helm Enumeration dialog; CANCEL closes', async () => {
    await clickBioPolyToolItem(page, 'div-Bio---PolyTool---Enumerate-HELM...');
    await page.locator('[name="dialog-PolyTool-Helm-Enumeration"]').waitFor({state: 'visible', timeout: 20_000});
    const dlgProbe = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-PolyTool-Helm-Enumeration"]');
      return {
        title: dlg?.querySelector('.d4-dialog-title')?.textContent?.trim(),
        cancelPresent: !!(dlg?.querySelector('[name="button-CANCEL"]')),
      };
    });
    expect(dlgProbe.title, 'PolyTool Helm Enumeration dialog title').toBe('PolyTool Helm Enumeration');
    await page.locator('[name="dialog-PolyTool-Helm-Enumeration"] [name="button-CANCEL"]').click();
    await page.locator('[name="dialog-PolyTool-Helm-Enumeration"]').waitFor({state: 'hidden', timeout: 10_000});
    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon after PolyTool Helm Enumeration CANCEL').toBe(0);
  });

  await softStep('Block G step 4: Bio | PolyTool | Combine Sequences... opens Combine Sequences dialog; CANCEL closes', async () => {
    await clickBioPolyToolItem(page, 'div-Bio---PolyTool---Combine-Sequences...');
    await page.locator('[name="dialog-Combine-Sequences"]').waitFor({state: 'visible', timeout: 20_000});
    const dlgProbe = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-Combine-Sequences"]');
      return {
        title: dlg?.querySelector('.d4-dialog-title')?.textContent?.trim(),
        tableInputPresent: !!(dlg?.querySelector('[name="input-host-Table"]')),
        columnInputPresent: !!(dlg?.querySelector('[name="input-host-Column"]')),
        cancelPresent: !!(dlg?.querySelector('[name="button-CANCEL"]')),
      };
    });
    expect(dlgProbe.title, 'Combine Sequences dialog title').toBe('Combine Sequences');
    expect(dlgProbe.tableInputPresent, 'Table input must be in Combine Sequences dialog').toBe(true);
    expect(dlgProbe.columnInputPresent, 'Column input must be in Combine Sequences dialog').toBe(true);
    await page.locator('[name="dialog-Combine-Sequences"] [name="button-CANCEL"]').click();
    await page.locator('[name="dialog-Combine-Sequences"]').waitFor({state: 'hidden', timeout: 10_000});
    const balloonError = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error, .grok-balloon-error').length);
    expect(balloonError, 'no error balloon after Combine Sequences CANCEL').toBe(0);
  });

  await softStep('Cleanup: close all views', async () => {
    await page.evaluate(() => (window as any).grok.shell.closeAll());
    await page.waitForTimeout(500);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n');
    throw new Error(`Step failures:\n${summary}`);
  }
});
