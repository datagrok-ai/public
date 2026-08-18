/* ---
realizes: [trellisplot.cp.split-and-pick-inner, trellisplot.int.split-columns-drive-inner-viewer-grid, trellisplot.int.viewer-type-change-control-panel-axes]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import {saveProjectViaUI, deleteProjectWithCleanup} from '../../helpers/projects';
import * as v from '../../helpers/viewers';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

const isBenignError = (text: string) =>
  /Failed to load resource/.test(text) || /404 \(\)/.test(text) || /favicon/.test(text) ||
  /Unable to find element in cloned iframe/.test(text) ||
  /NullError: method not found: '\w+' on null/.test(text) ||
  /ProjectMeta\.publish/.test(text) || /project_meta\.dart/.test(text);

async function buildDemogTrellis(page: Page): Promise<void> {
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch {}
    try { grok.shell.windows.simpleMode = true; } catch {}
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    grok.shell.addTableView(df);

    await new Promise((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
      setTimeout(resolve, 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});
  await v.addViewerByIcon(page, 'trellis-plot', 'Trellis-plot', 15000);

  await v.waitForViewerRendered(page, 'Trellis plot', 900);
}

async function reopenAndReadTrellis(page: Page, projectId: string): Promise<{
  hasTrellis: boolean; rowSource: string | null; onClick: string | null;
  viewerType: string | null; x: string[]; y: string[]; viewerTypes: string[];
  xCats: number; yCats: number; cells: number; selected: number;
  rootW: number; rootH: number;
}> {
  const reading = await page.evaluate(async (id) => {
    grok.shell.closeAll();

    await new Promise((r) => setTimeout(r, 1500));
    const proj = await grok.dapi.projects.find(id);
    await proj.open();
    let types: string[] = [];
    let tp: any = null;
    let tpView: any = null;
    for (let t = 0; t < 20; t++) {
      await new Promise((r) => setTimeout(r, 1000));
      types = [];
      tp = null;
      tpView = null;
      for (const view of grok.shell.tableViews)
        for (const vw of view.viewers) {
          types.push(vw.type);
          if (vw.type === 'Trellis plot') { tp = vw; tpView = view; }
        }
      if (tp) break;
    }

    let cells = 0;
    for (let t = 0; t < 10 && tp; t++) {
      try { cells = tp.root ? tp.root.querySelectorAll('.d4-trellis-plot-cell').length : 0; }
      catch (_) { cells = 0; }
      if (cells === 16) break;
      await new Promise((r) => setTimeout(r, 1000));
    }

    let selected = -1;
    try { selected = tpView ? tpView.dataFrame.selection.trueCount : -1; } catch (_) { selected = -1; }
    let xCats = -1;
    let yCats = -1;
    if (tp) {
      try { xCats = tp.xCategoriesCount; } catch (_) {  }
      try { yCats = tp.yCategoriesCount; } catch (_) {  }
    }

    let rootW = -1;
    let rootH = -1;
    if (tp && tp.root) {
      const b = tp.root.getBoundingClientRect();
      rootW = b.width;
      rootH = b.height;
    }
    return {
      hasTrellis: !!tp,
      rowSource: tp ? tp.props.rowSource : null,
      onClick: tp ? tp.props.onClick : null,
      viewerType: tp ? tp.props.viewerType : null,
      x: tp ? [...tp.props.xColumnNames] as string[] : [],
      y: tp ? [...tp.props.yColumnNames] as string[] : [],
      viewerTypes: types, xCats, yCats, cells, selected, rootW, rootH,
    };
  }, projectId);

  console.log(`[trellis reopen] rowSource=${reading.rowSource} onClick=${reading.onClick} ` +
    `viewerType=${reading.viewerType} x=${JSON.stringify(reading.x)} y=${JSON.stringify(reading.y)} ` +
    `cats=${reading.xCats}x${reading.yCats} selected=${reading.selected} cells=${reading.cells} ` +
    `root=${reading.rootW}x${reading.rootH}`);
  return reading;
}

async function expectRestoredEmptyGridWithLiveness(page: Page,
  r: {cells: number; rowSource: string | null; selected: number;
    x: string[]; rootW: number; rootH: number}): Promise<void> {

  expect(r.rowSource,
    `restored trellis reads Row Source = ${r.rowSource}, but the empty-grid rule graded here only holds for 'Selected'`).toBe('Selected');
  expect(r.selected,
    `the reopened view came back with ${r.selected} rows selected — every reopen path in this spec was measured to restore an EMPTY selection [DOM 2026-08-12], so a non-zero count means the recorded fact changed and the grading rule has to be re-derived before this grid can be graded`).toBe(0);
  expect(r.cells,
    `restored trellis painted ${r.cells} cells under Row Source = Selected with an EMPTY selection — that state has no rows to plot, so the grid must be empty; root measured ${r.rootW}x${r.rootH}`).toBe(0);

  const revived = await page.evaluate(async () => {
    let tp: any = null;
    let tpView: any = null;
    for (const view of grok.shell.tableViews)
      for (const vw of view.viewers) if (vw.type === 'Trellis plot') { tp = vw; tpView = view; }
    if (!tp) return {cells: -1, selected: -1, x: [] as string[]};

    const settled = new Promise<void>((res) => {
      let sub: any = null;
      try { sub = tp.onViewerRendered.subscribe(() => { sub.unsubscribe(); res(); }); }
      catch (_) {  }
      setTimeout(() => { try { sub?.unsubscribe(); } catch (_) {} res(); }, 2500);
    });
    tpView.dataFrame.selection.setAll(true);
    await settled;
    let cells = 0;
    for (let t = 0; t < 10; t++) {
      try { cells = tp.root ? tp.root.querySelectorAll('.d4-trellis-plot-cell').length : 0; }
      catch (_) { cells = 0; }
      if (cells === 16) break;
      await new Promise((res) => setTimeout(res, 1000));
    }
    return {cells, selected: tpView.dataFrame.selection.trueCount as number,
      x: [...tp.props.xColumnNames] as string[]};
  });
  expect(revived.selected,
    'selecting every row left the restored view with an empty selection — the liveness witness never ran').toBeGreaterThan(0);
  expect(revived.cells,
    `restored trellis painted ${revived.cells} cells after all ${revived.selected} rows were selected — the grid owes the full restored 4 x 4 product here, so 0 means the empty grid was NOT the empty selection (the restored viewer does not render at all) and an intermediate count means the split columns were lost on the round-trip`).toBe(16);

  expect(revived.x).toEqual(r.x);
}

test('Trellis plot: split columns, inner-type switching, persistence', async ({page}) => {
  test.setTimeout(1_500_000);
  page.setDefaultTimeout(120_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error' && !isBenignError(m.text())) consoleErrors.push(m.text()); });

  await loginToDatagrok(page);

  await page.waitForTimeout(5000);

  const addPageErrBefore = pageErrors.length;
  const addConsoleErrBefore = consoleErrors.length;

  await buildDemogTrellis(page);

  await softStep('Scenario 1 Step 1', async () => {
    await expect(page.locator('[name="viewer-Trellis-plot"]')).toHaveCount(1);

    expect(pageErrors.slice(addPageErrBefore),
      'adding the Trellis Plot viewer raised an uncaught page error (github-964 smoke guard)').toEqual([]);
    expect(consoleErrors.slice(addConsoleErrBefore),
      'adding the Trellis Plot viewer raised a non-benign console error (github-964 smoke guard)').toEqual([]);
  });

  const cardinalities = await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    return {
      sex: df.col('SEX').categories.length,
      race: df.col('RACE').categories.length,
      control: df.col('CONTROL').categories.length,
    };
  });
  expect(cardinalities).toEqual({sex: 2, race: 4, control: 2});

  const cellLocator = page.locator('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell');

  await softStep('Scenario 1 Step 4', async () => {
    const counts = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;

      const settled = new Promise<void>((res) => {
        let sub: any = null;
        try { sub = tp.onViewerRendered.subscribe(() => { sub.unsubscribe(); res(); }); }
        catch (_) {  }
        setTimeout(() => { try { sub?.unsubscribe(); } catch (_) {} res(); }, 1500);
      });
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      await settled;
      return {xCat: tp.xCategoriesCount, yCat: tp.yCategoriesCount};
    });

    await expect(cellLocator).toHaveCount(counts.xCat * counts.yCat);
    await expect(cellLocator).toHaveCount(8);
    expect(counts).toEqual({xCat: 2, yCat: 4});
  });

  await softStep('Scenario 1 Step 6', async () => {

    const backdrop = page.locator('.d4-column-selector-backdrop');
    const openAddXPopup = async (): Promise<boolean> => {
      const plus = page.locator('[name="viewer-Trellis-plot"] [name="add-x-column"]').first();

      const box = await plus.boundingBox({timeout: 8000}).catch(() => null);
      if (!box) return false;
      await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
      await page.mouse.down();
      await page.mouse.up();
      return backdrop.waitFor({timeout: 6000}).then(() => true).catch(() => false);
    };

    const findRowByHover = async (column: RegExp) => {
      const r = await page.evaluate(() => {
        const b = document.querySelector('.d4-column-selector-backdrop')!.getBoundingClientRect();
        return {left: b.left, top: b.top, height: b.height};
      });
      let rowY: number | null = null;
      const maxDy = Math.min(Math.max(r.height - 6, 54), 240);
      for (let dy = 6; dy <= maxDy; dy += 6) {
        await page.mouse.move(r.left + 60, r.top + dy);

        await page.waitForTimeout(220);
        const tip = await page.evaluate(() => {
          const el = document.querySelector('.d4-tooltip') as HTMLElement | null;
          return el && el.getBoundingClientRect().width > 0 ? el.innerText : '';
        });

        if (column.test(tip)) { rowY = r.top + dy; break; }
      }
      return {hx: r.left + 60, rowY};
    };
    const readXState = () => page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      return {x: [...tp.props.xColumnNames] as string[],
        cells: document.querySelectorAll('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell').length};
    });

    try {
      expect(await openAddXPopup(),
        'the (+) add-X-column control did not open the column selector popup').toBe(true);
      const first = await findRowByHover(/\bCONTROL\b/);
      expect(first.rowY,
        'CONTROL row was not found by tooltip scan of the popup — no blind-position fallback is taken').not.toBeNull();

      const dwell: {x: string[]; cells: number}[] = [];
      for (let i = 0; i < 4; i++) {
        await page.mouse.move(first.hx, first.rowY!);

        await page.waitForTimeout(220);
        dwell.push(await readXState());
      }

      expect(dwell.some((s) => s.cells !== 8 || s.x.length !== 1),
        'hover preview never rebuilt the grid during the dwell — preview-regression signal (GROK-19673 class)').toBe(true);

      await page.keyboard.press('Escape');
      const closedOnFirstEsc = await backdrop
        .waitFor({state: 'detached', timeout: 2000}).then(() => true).catch(() => false);
      if (!closedOnFirstEsc) {
        await page.keyboard.press('Escape');
        await backdrop.waitFor({state: 'detached', timeout: 5000});
      }
      await expect(cellLocator).toHaveCount(8);
      expect((await readXState()).x).toEqual(['SEX']);

      expect(await openAddXPopup(),
        'the (+) add-X-column control did not reopen the column selector popup').toBe(true);
      const second = await findRowByHover(/\bCONTROL\b/);
      expect(second.rowY).not.toBeNull();
      await page.mouse.click(second.hx, second.rowY!);

      await v.waitForViewerRendered(page, 'Trellis plot', 900);
      const committed = await page.evaluate(() => {
        const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
        return {x: [...tp.props.xColumnNames], xCat: tp.xCategoriesCount, yCat: tp.yCategoriesCount};
      });
      expect(committed.x).toEqual(['SEX', 'CONTROL']);
      expect({xCat: committed.xCat, yCat: committed.yCat}).toEqual({xCat: 4, yCat: 4});
      await expect(cellLocator).toHaveCount(16);
    } finally {

      await page.keyboard.press('Escape').catch(() => {});
      await page.evaluate(() =>
        document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}))).catch(() => {});
    }
  });

  await softStep('Scenario 1 Step 7', async () => {

    const BORDER = 1;
    const HEADER_H = 20;
    const ROW_H = 16;
    const AUTOSIZE_PAD = 10;

    const backdrop = page.locator('.d4-column-selector-backdrop');
    const readXState = () => page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      return {x: [...tp.props.xColumnNames] as string[],
        xCat: tp.xCategoriesCount as number, yCat: tp.yCategoriesCount as number};
    });

    const slot = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]');
      const plus = root ? root.querySelector('[name="add-x-column"]') : null;
      const host = plus ? plus.parentElement : null;
      if (!host) return null;
      const slots = Array.from(host.children).filter((el) =>
        el !== plus && el.classList.contains('d4-column-selector')) as HTMLElement[];
      const labels = slots.map((el) =>
        (el.querySelector('.d4-column-selector-column') as HTMLElement | null)?.innerText?.trim() ?? '');
      const hits = labels.filter((l) => l === 'CONTROL').length;
      if (hits !== 1) return {labels, hits, cx: -1, cy: -1};
      const b = slots[labels.indexOf('CONTROL')].getBoundingClientRect();
      return {labels, hits, cx: b.left + b.width / 2, cy: b.top + b.height / 2};
    });
    expect(slot,
      'the X selectors host ([name="add-x-column"] parent) was not found — Step 6 must have committed CONTROL first').not.toBeNull();
    expect(slot!.hits,
      `expected exactly one CONTROL slot among the X selectors, found labels ${JSON.stringify(slot?.labels)} — nothing is clicked`).toBe(1);

    await page.mouse.move(slot!.cx, slot!.cy);
    await page.mouse.down();
    await page.mouse.up();
    const opened = await backdrop.waitFor({timeout: 6000}).then(() => true).catch(() => false);
    expect(opened, 'the CONTROL X-slot did not open the column selector popup').toBe(true);

    try {
      const rect = await page.evaluate(() => {
        const b = document.querySelector('.d4-column-selector-backdrop')!.getBoundingClientRect();
        return {left: b.left, top: b.top, width: b.width, height: b.height};
      });

      const modelRows = (rect.height - (2 * BORDER + HEADER_H + AUTOSIZE_PAD)) / ROW_H;
      expect(Number.isInteger(modelRows) && modelRows >= 2,
        `column picker measured ${rect.width}x${rect.height}: its height is not 32 + 16*rows, so the layout model that locates the blank first row no longer holds — re-derive the offsets from these numbers`).toBe(true);

      const hx = rect.left + Math.min(60, rect.width / 2);
      const hy = rect.top + BORDER + HEADER_H + ROW_H / 2;

      await page.mouse.move(hx, hy);

      await page.waitForTimeout(900);
      const previewed = await readXState();

      expect(previewed.x,
        `hovering the blank first row at (${Math.round(hx)}, ${Math.round(hy)}) of a ${rect.width}x${rect.height} popup left the X split columns as ${JSON.stringify(previewed.x)} — the pointer did not land on the blank row, so nothing is committed`).toEqual(['SEX']);

      await page.mouse.click(hx, hy);
      await backdrop.waitFor({state: 'detached', timeout: 6000});

      await v.waitForViewerRendered(page, 'Trellis plot', 900);

      const after = await readXState();

      expect(after.x).toEqual(['SEX']);
      expect({xCat: after.xCat, yCat: after.yCat}).toEqual({xCat: 2, yCat: 4});
      await expect(cellLocator).toHaveCount(8);
    } finally {

      await page.keyboard.press('Escape').catch(() => {});
      await page.evaluate(() =>
        document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}))).catch(() => {});
    }
  });

  await softStep('Scenario 2 Step 3', async () => {

    const perType = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;

      function cellSig(cellIdx: number): number | null {
        const cell = root.querySelectorAll('.d4-trellis-plot-cell')[cellIdx];
        const cv = cell?.querySelector('canvas') as HTMLCanvasElement | null;
        if (!cv) return null;
        try {
          const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
          let sig = 0;
          for (let i = 0; i < img.length; i += 4)
            sig = (sig * 31 + ((img[i] << 16) | (img[i + 1] << 8) | img[i + 2])) % 2147483647;
          return sig;
        } catch { return null; }
      }

      const typeIcons: Record<string, string> = {
        'Scatter plot': 'icon-scatter-plot', 'Bar chart': 'icon-bar-chart',
        'Histogram': 'icon-histogram', 'Pie chart': 'icon-pie-chart',
      };
      const out: Record<string, {eventType: string | null; distinct: boolean}> = {};
      for (const t of ['Scatter plot', 'Bar chart', 'Histogram', 'Pie chart']) {
        let eventType: string | null = null;
        const sub = tp.onEvent('d4-trellis-plot-viewer-type-changed').subscribe((arg: any) => {
          eventType = (typeof arg === 'string' ? arg : (arg?.args?.viewerType ?? null));
        });
        const vs = root.querySelector('[name="viewer selector"]') as HTMLElement;
        vs.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0}));

        await new Promise((r) => setTimeout(r, 600));
        const item = document.querySelector(`.d4-combo-drop-down [name="${typeIcons[t]}"]`);

        const settled = new Promise<void>((res) => {
          let rsub: any = null;
          try { rsub = tp.onViewerRendered.subscribe(() => { rsub.unsubscribe(); res(); }); }
          catch (_) {  }
          setTimeout(() => { try { rsub?.unsubscribe(); } catch (_) {} res(); }, 1500);
        });
        (item?.closest('.d4-list-item') as HTMLElement | null)?.click();
        await settled;
        sub?.unsubscribe?.();

        let distinct = false;
        const deadline = Date.now() + 3000;
        do {
          const populated: number[] = [];
          const cells = root.querySelectorAll('.d4-trellis-plot-cell');
          for (let i = 0; i < cells.length && populated.length < 2; i++)
            if (cells[i].querySelector('canvas')) populated.push(i);
          if (populated.length === 2) {
            const a = cellSig(populated[0]);
            const b = cellSig(populated[1]);
            distinct = a !== null && b !== null && a !== b;
          }
          if (!distinct) await new Promise((r) => setTimeout(r, 150));
        } while (!distinct && Date.now() < deadline);
        out[t] = {eventType, distinct};
      }
      return out;
    });
    for (const t of ['Scatter plot', 'Bar chart', 'Histogram', 'Pie chart']) {
      expect(perType[t].eventType).toBe(t);
      expect(perType[t].distinct).toBe(true);
    }
  });

  await softStep('Scenario 2 Step 5', async () => {

    const result = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;

      const settle = (capMs: number, act: () => void) => new Promise<void>((res) => {
        let rsub: any = null;
        try { rsub = tp.onViewerRendered.subscribe(() => { rsub.unsubscribe(); res(); }); }
        catch (_) {  }
        setTimeout(() => { try { rsub?.unsubscribe(); } catch (_) {} res(); }, capMs);
        act();
      });

      function populatedCellIdxs(limit: number): number[] {
        const idxs: number[] = [];
        const cells = root.querySelectorAll('.d4-trellis-plot-cell');
        for (let i = 0; i < cells.length && idxs.length < limit; i++)
          if (cells[i].querySelector('canvas')) idxs.push(i);
        return idxs;
      }
      function cellHash(cellIdx: number): number | null {
        const cell = root.querySelectorAll('.d4-trellis-plot-cell')[cellIdx];
        const cv = cell?.querySelector('canvas') as HTMLCanvasElement | null;
        if (!cv) return null;
        try {
          const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
          let h = 0;
          for (let i = 0; i < img.length; i += 4)
            h = (h * 31 + ((img[i] << 16) | (img[i + 1] << 8) | img[i + 2])) % 2147483647;
          return h;
        } catch { return null; }
      }

      await settle(1500, () => {
        tp.props.viewerType = 'Pie chart';
        tp.props.xColumnNames = ['SEX'];
        tp.props.yColumnNames = ['RACE'];
      });

      await settle(1200, () => {
        tp.props.yColumnNames = [];
        tp.props.viewerType = 'Bar chart';
      });

      await settle(1500, () => { tp.props.yColumnNames = ['RACE']; });
      const idxs = populatedCellIdxs(2);
      const barA = cellHash(idxs[0]);
      const barB = cellHash(idxs[1]);

      await settle(1800, () => { tp.props.viewerType = 'Pie chart'; });

      let pieA = cellHash(idxs[0]);
      let pieB = cellHash(idxs[1]);
      const deadline = Date.now() + 3000;
      while ((pieA === barA || pieB === barB || pieA === pieB) && Date.now() < deadline) {
        await new Promise((r) => setTimeout(r, 150));
        pieA = cellHash(idxs[0]);
        pieB = cellHash(idxs[1]);
      }

      return {
        hashesRead: barA !== null && barB !== null && pieA !== null && pieB !== null,

        distinctAfter: pieA !== null && pieB !== null && pieA !== pieB,
        reRendered: barA !== pieA && barB !== pieB,
      };
    });
    expect(result.hashesRead).toBe(true);
    expect(result.distinctAfter).toBe(true);
    expect(result.reRendered).toBe(true);
  });

  await softStep('Scenario 2 Step 7', async () => {
    const result = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.viewerType = 'Pie chart';
      tp.props.showControlPanel = true;

      await new Promise((r) => setTimeout(r, 800));

      function selectorVisible(): boolean {
        const sel = root.querySelector('[name="viewer selector"]') as HTMLElement | null;
        if (!sel) return false;
        const r = sel.getBoundingClientRect();
        return r.width > 0 && r.height > 0;
      }
      const visibleBefore = selectorVisible();

      tp.props.showControlPanel = false;

      await new Promise((r) => setTimeout(r, 1000));
      const visibleAfterHide = selectorVisible();
      const typeAfterHide = tp.props.viewerType;

      tp.props.showControlPanel = true;

      await new Promise((r) => setTimeout(r, 800));
      return {visibleBefore, visibleAfterHide, typeAfterHide};
    });
    expect(result.visibleBefore).toBe(true);

    expect(result.visibleAfterHide).toBe(false);

    expect(result.typeAfterHide).toBe('Pie chart');
  });

  await softStep('Scenario 3 Step 6', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const trellises = () => Array.from(tv.viewers).filter((x: any) => x.type === 'Trellis plot') as any[];

      const readPair = () => trellises().map((v: any) => ({
        x: [...v.props.xColumnNames] as string[],
        y: [...v.props.yColumnNames] as string[],
        type: v.props.viewerType as string,
        xCat: v.xCategoriesCount as number,
        yCat: v.yCategoriesCount as number,
      })).sort((a, b) => a.type.localeCompare(b.type));

      const onEvent = (stream: any, capMs: number, act: () => void) => new Promise<void>((res) => {
        let esub: any = null;
        try { esub = stream.subscribe(() => { esub.unsubscribe(); res(); }); }
        catch (_) {  }
        setTimeout(() => { try { esub?.unsubscribe(); } catch (_) {} res(); }, capMs);
        act();
      });

      const settleViewer = (viewer: any, capMs: number, act: () => void) => new Promise<void>((res) => {
        let rsub: any = null;
        try { rsub = viewer.onViewerRendered.subscribe(() => { rsub.unsubscribe(); res(); }); }
        catch (_) {  }
        setTimeout(() => { try { rsub?.unsubscribe(); } catch (_) {} res(); }, capMs);
        act();
      });

      const tp = trellises()[0];
      await settleViewer(tp, 1500, () => {
        tp.props.xColumnNames = ['SEX', 'CONTROL'];
        tp.props.yColumnNames = ['RACE'];
        tp.props.viewerType = 'Pie chart';
      });

      let other: any = null;
      await onEvent(grok.events.onViewerAdded, 2000, () => { other = tv.addViewer('Trellis plot'); });
      await settleViewer(other, 2000, () => {
        other.props.xColumnNames = ['RACE'];
        other.props.yColumnNames = ['SEX'];
        other.props.viewerType = 'Bar chart';
      });
      const savedPair = readPair();

      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;

      await new Promise((r) => setTimeout(r, 1200));

      try {

        await onEvent(grok.events.onViewerAdded, 1000, () => { tv.addViewer('Scatter plot'); });
        const viewersBefore = Array.from(tv.viewers).map((x: any) => x.type);

        const saved = await grok.dapi.layouts.find(layoutId);

        await onEvent(grok.events.onViewLayoutApplied, 4000, () => { tv.loadLayout(saved); });
        const viewersAfterLoad = Array.from(tv.viewers).map((x: any) => x.type);
        const restoredPair = readPair();

        const pie = trellises().find((v: any) => v.props.viewerType === 'Pie chart');
        const bar = trellises().find((v: any) => v.props.viewerType === 'Bar chart');
        let crossTalk: string | null = null;
        if (pie && bar) {
          await settleViewer(pie, 1500, () => { pie.props.viewerType = 'Histogram'; });
          crossTalk = bar.props.viewerType;
          await settleViewer(pie, 1500, () => { pie.props.viewerType = 'Pie chart'; });
        }

        bar?.close?.();
        await new Promise((r) => setTimeout(r, 2500));
        const survivors = readPair();
        return {viewersBefore, viewersAfterLoad, savedPair, restoredPair, crossTalk, survivors};
      } finally {

        await grok.dapi.layouts.find(layoutId)
          .then((l: any) => l && grok.dapi.layouts.delete(l)).catch(() => {});
      }
    });
    expect(result.viewersBefore).toContain('Scatter plot');
    expect(result.viewersAfterLoad).not.toContain('Scatter plot');
    expect(result.viewersAfterLoad).toContain('Trellis plot');
    expect(result.savedPair).toEqual([
      {x: ['RACE'], y: ['SEX'], type: 'Bar chart', xCat: 4, yCat: 2},
      {x: ['SEX', 'CONTROL'], y: ['RACE'], type: 'Pie chart', xCat: 4, yCat: 4},
    ]);
    expect(result.restoredPair).toEqual(result.savedPair);
    expect(result.crossTalk,
      'changing the inner type of one restored trellis changed the other — the GROK-15494 shared-state leak').toBe('Bar chart');

    expect(result.survivors).toEqual([
      {x: ['SEX', 'CONTROL'], y: ['RACE'], type: 'Pie chart', xCat: 4, yCat: 4},
    ]);
    await expect(cellLocator).toHaveCount(16);
  });

  const nameSelected = `zz-trellis-p0-selected-${Date.now()}`;
  const nameEmpty = `zz-trellis-p0-empty-${Date.now()}`;
  let idSelected: string | null = null;
  let idEmpty: string | null = null;

  let savedSelectionCount = -1;
  try {

    await softStep('Scenario 3 Step 7', async () => {
      const pre = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const tp = Array.from(tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
        tp.props.onClick = 'Select';

        await new Promise((r) => setTimeout(r, 800));
        tv.dataFrame.selection.setAll(false);

        await new Promise((r) => setTimeout(r, 300));

        const captured: {mc: Record<string, string> | null} = {mc: null};
        const sub = tp.onEvent('d4-trellis-plot-current-cell-changed').subscribe((arg: any) => {
          const mc = arg?.args?.matchCondition ?? arg?.matchCondition ?? null;
          if (mc) captured.mc = Object.fromEntries(Object.entries(mc).map(([k, v]) => [k, String(v)]));
        });

        const cell = document.querySelector('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell') as HTMLElement;
        const rect = cell.getBoundingClientRect();
        const o = {bubbles: true, cancelable: true, view: window, button: 0,
          clientX: rect.left + 8, clientY: rect.top + 8};

        const selSettled = new Promise<void>((res) => {
          let ssub: any = null;
          try { ssub = tv.dataFrame.onSelectionChanged.subscribe(() => { ssub.unsubscribe(); res(); }); }
          catch (_) {  }
          setTimeout(() => { try { ssub?.unsubscribe(); } catch (_) {} res(); }, 1200);
        });
        cell.dispatchEvent(new MouseEvent('mousedown', o));
        cell.dispatchEvent(new MouseEvent('mouseup', o));
        cell.dispatchEvent(new MouseEvent('click', o));
        await selSettled;
        sub?.unsubscribe?.();
        const selected = tv.dataFrame.selection.trueCount;

        let expectedRows = -1;
        if (captured.mc) {
          const df = tv.dataFrame;
          const pairs = Object.entries(captured.mc).map(([c, v]) => [df.col(c), v] as [any, string]);
          expectedRows = 0;
          for (let i = 0; i < df.rowCount; i++)
            if (pairs.every(([col, v]) => String(col.get(i)) === v)) expectedRows++;
        }

        const rsSettled = new Promise<void>((res) => {
          let rsub: any = null;
          try { rsub = tp.onViewerRendered.subscribe(() => { rsub.unsubscribe(); res(); }); }
          catch (_) {  }
          setTimeout(() => { try { rsub?.unsubscribe(); } catch (_) {} res(); }, 1000);
        });
        tp.props.rowSource = 'Selected';
        await rsSettled;
        return {selected, expectedRows, matchCondition: captured.mc,
          rowSource: tp.props.rowSource, onClick: tp.props.onClick};
      });

      expect(pre.matchCondition,
        'the corner click did not fire d4-trellis-plot-current-cell-changed — the click never reached the trellis cell handler').not.toBeNull();

      expect(pre.expectedRows,
        `clicked combination ${JSON.stringify(pre.matchCondition)} holds no rows — pick a populated probe cell`).toBeGreaterThan(0);
      expect(pre.selected).toBe(pre.expectedRows);
      expect(pre.onClick).toBe('Select');
      expect(pre.rowSource).toBe('Selected');

      savedSelectionCount = pre.selected;

      const pageErrBefore = pageErrors.length;
      const saved = await saveProjectViaUI(page, nameSelected);
      idSelected = saved.projectId;
      expect(idSelected).toBeTruthy();

      expect(pageErrors.slice(pageErrBefore)).toEqual([]);
    });

    await softStep('Scenario 3 Step 11', async () => {
      const pageErrBefore = pageErrors.length;
      const errorsBefore = consoleErrors.length;
      const r = await reopenAndReadTrellis(page, idSelected!);

      expect(r.hasTrellis).toBe(true);
      expect(r.viewerTypes).toContain('Trellis plot');
      expect(pageErrors.slice(pageErrBefore)).toEqual([]);
      expect(consoleErrors.slice(errorsBefore)).toEqual([]);
      expect(r.x).toEqual(['SEX', 'CONTROL']);
      expect(r.y).toEqual(['RACE']);
      expect(r.viewerType).toBe('Pie chart');
      expect(r.onClick).toBe('Select');
      expect(r.rowSource).toBe('Selected');
      expect({xCats: r.xCats, yCats: r.yCats}).toEqual({xCats: 4, yCats: 4});

      expect(savedSelectionCount,
        'Step 7 did not record a live selection at save time — without it Step 11 cannot say anything about whether a selection survives the project round-trip').toBeGreaterThan(0);
      expect(r.selected,
        `the project was saved with ${savedSelectionCount} rows selected and reopened with ${r.selected} — the recorded behaviour is that a selection does NOT survive the project round-trip [DOM 2026-08-12], so a non-zero count means the product now persists it and the grading of both this step and Step 12 has to be re-derived`).toBe(0);

      await expectRestoredEmptyGridWithLiveness(page, r);

      const allSource = await page.evaluate(async () => {
        let tp: any = null;
        let tpView: any = null;
        for (const view of grok.shell.tableViews)
          for (const vw of view.viewers) if (vw.type === 'Trellis plot') { tp = vw; tpView = view; }
        if (!tp) return {clearedCells: -1, cells: -1, selected: -1,
          rowSource: null as string | null, x: [] as string[]};

        const clearSettled = new Promise<void>((res) => {
          let csub: any = null;
          try { csub = tp.onViewerRendered.subscribe(() => { csub.unsubscribe(); res(); }); }
          catch (_) {  }
          setTimeout(() => { try { csub?.unsubscribe(); } catch (_) {} res(); }, 2000);
        });
        tpView.dataFrame.selection.setAll(false);
        await clearSettled;
        let clearedCells = 0;
        try { clearedCells = tp.root ? tp.root.querySelectorAll('.d4-trellis-plot-cell').length : 0; }
        catch (_) { clearedCells = -1; }
        tp.props.rowSource = 'All';
        let cells = 0;
        for (let t = 0; t < 10; t++) {
          await new Promise((res) => setTimeout(res, 1000));
          try { cells = tp.root ? tp.root.querySelectorAll('.d4-trellis-plot-cell').length : 0; }
          catch (_) { cells = 0; }
          if (cells === 16) break;
        }
        return {clearedCells, cells, selected: tpView.dataFrame.selection.trueCount as number,
          rowSource: tp.props.rowSource as string, x: [...tp.props.xColumnNames] as string[]};
      });
      expect(allSource.clearedCells,
        `clearing the selection left ${allSource.clearedCells} cells under Row Source = Selected — the grid was expected to empty again, and without that the 16 asserted below could simply be the previous witness still on screen`).toBe(0);
      expect(allSource.selected,
        'the selection was not cleared before the row-source switch — the All probe would then be graded on a still-selected dataset').toBe(0);
      expect(allSource.rowSource).toBe('All');
      expect(allSource.cells,
        `restored trellis painted ${allSource.cells} cells under Row Source = All with NOTHING selected — that state plots the whole dataset, so 0 means the restored viewer only renders through a selection and an intermediate count means the split columns were lost on the round-trip (8 = CONTROL dropped)`).toBe(16);
      expect(allSource.x,
        'the grid filled under Row Source = All but with different split columns — it was rebuilt from defaults, not from the restored configuration').toEqual(['SEX', 'CONTROL']);
    });

    await softStep('Scenario 3 Step 12', async () => {
      await buildDemogTrellis(page);
      const pre = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const tp = Array.from(tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
        tv.dataFrame.selection.setAll(false);
        tp.props.xColumnNames = ['SEX', 'CONTROL'];
        tp.props.yColumnNames = ['RACE'];
        tp.props.viewerType = 'Pie chart';

        tp.props.onClick = 'Select';

        await new Promise((r) => setTimeout(r, 500));

        const rsSettled = new Promise<void>((res) => {
          let rsub: any = null;
          try { rsub = tp.onViewerRendered.subscribe(() => { rsub.unsubscribe(); res(); }); }
          catch (_) {  }
          setTimeout(() => { try { rsub?.unsubscribe(); } catch (_) {} res(); }, 1500);
        });
        tp.props.rowSource = 'Selected';
        await rsSettled;
        return {selected: tv.dataFrame.selection.trueCount, rowSource: tp.props.rowSource};
      });

      expect(pre.selected).toBe(0);
      expect(pre.rowSource).toBe('Selected');

      const saved = await saveProjectViaUI(page, nameEmpty);
      idEmpty = saved.projectId;
      expect(idEmpty).toBeTruthy();

      const pageErrBefore = pageErrors.length;
      const errorsBefore = consoleErrors.length;
      const r = await reopenAndReadTrellis(page, idEmpty!);
      expect(r.hasTrellis).toBe(true);
      expect(pageErrors.slice(pageErrBefore)).toEqual([]);
      expect(consoleErrors.slice(errorsBefore)).toEqual([]);
      expect(r.rowSource).toBe('Selected');
      expect(r.x).toEqual(['SEX', 'CONTROL']);
      expect(r.y).toEqual(['RACE']);
      expect({xCats: r.xCats, yCats: r.yCats}).toEqual({xCats: 4, yCats: 4});

      expect(r.selected,
        `the reopened project came back with ${r.selected} rows selected — it was saved with none, so this step no longer covers the empty-selection restore shape`).toBe(0);
      await expectRestoredEmptyGridWithLiveness(page, r);
    });
  } finally {

    if (idSelected) await deleteProjectWithCleanup(page, {projectId: idSelected});
    if (idEmpty) await deleteProjectWithCleanup(page, {projectId: idEmpty});
    await page.evaluate(() => grok.shell.closeAll()).catch(() => {});
  }

  v.finishSpec();
});
