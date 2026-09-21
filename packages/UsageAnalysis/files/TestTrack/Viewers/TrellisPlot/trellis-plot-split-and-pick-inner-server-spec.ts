/* ---
realizes: [trellisplot.cp.split-and-pick-inner, trellisplot.int.split-columns-drive-inner-viewer-grid]
--- */
import {expect, Page} from '@playwright/test';
import {test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import {saveProjectViaApi, deleteProjectWithCleanup} from '../../helpers/projects';
import * as v from '../../helpers/viewers';

declare const grok: any;

// Scenario 3 of the split-and-pick scenario: a layout and two projects round-tripped through the
// server. Scenarios 1-2 are trellis-plot-split-and-pick-inner-spec.ts on the local lane.
test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

const isBenignError = (text: string) =>
  /Failed to load resource/.test(text) || /404 \(\)/.test(text) || /favicon/.test(text) ||
  /Unable to find element in cloned iframe/.test(text) ||
  /NullError: method not found: '\w+' on null/.test(text) ||
  /ProjectMeta\.publish/.test(text) || /project_meta\.dart/.test(text);

async function buildDemogTrellis(page: Page): Promise<void> {
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
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
    const w = window as any;
    grok.shell.closeAll();
    await w.__poll(() => Array.from(grok.shell.tableViews).length, (n: number) => n === 0, 1500, 50);
    const proj = await grok.dapi.projects.find(id);
    await proj.open();
    const find = () => {
      let tp: any = null;
      let tpView: any = null;
      const types: string[] = [];
      for (const view of grok.shell.tableViews)
        for (const vw of view.viewers) {
          types.push(vw.type);
          if (vw.type === 'Trellis plot') { tp = vw; tpView = view; }
        }
      return {tp, tpView, types};
    };
    const found = await w.__poll(find, (f: any) => !!f.tp, 20000, 100);
    const {tp, tpView, types} = found;

    const cellCount = () => { try { return tp?.root ? tp.root.querySelectorAll('.d4-trellis-plot-cell').length : 0; } catch (_) { return 0; } };
    // the restored grid is either empty (Row Source = Selected, nothing selected) or 16 cells; a
    // quiet gap on the count covers both without waiting the full cap for the empty case
    const cells = tp ? await w.__settledFor(cellCount, 500, 4000, 50) : 0;

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
    const w = window as any;
    let tp: any = null;
    let tpView: any = null;
    for (const view of grok.shell.tableViews)
      for (const vw of view.viewers) if (vw.type === 'Trellis plot') { tp = vw; tpView = view; }
    if (!tp) return {cells: -1, selected: -1, x: [] as string[]};

    const cellCount = () => { try { return tp.root ? tp.root.querySelectorAll('.d4-trellis-plot-cell').length : 0; } catch (_) { return 0; } };
    tpView.dataFrame.selection.setAll(true);
    const cells = await w.__poll(cellCount, (n: number) => n === 16, 5000, 100);
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
  test.setTimeout(300_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  const onPageError = (e: Error) => { pageErrors.push(String(e)); };
  const onConsole = (m: any) => { if (m.type() === 'error' && !isBenignError(m.text())) consoleErrors.push(m.text()); };
  page.on('pageerror', onPageError);
  page.on('console', onConsole);

  await openDatagrok(page);
  await buildDemogTrellis(page);

  const cellLocator = page.locator('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell');

  await softStep('Scenario 3 Step 6', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const tv = grok.shell.tv;
      const trellises = () => Array.from(tv.viewers).filter((x: any) => x.type === 'Trellis plot') as any[];

      const readPair = () => trellises().map((vw: any) => ({
        x: [...vw.props.xColumnNames] as string[],
        y: [...vw.props.yColumnNames] as string[],
        type: vw.props.viewerType as string,
        xCat: vw.xCategoriesCount as number,
        yCat: vw.yCategoriesCount as number,
      })).sort((a, b) => a.type.localeCompare(b.type));

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
      await w.__settled('grok.events.onViewerAdded', () => { other = tv.addViewer('Trellis plot'); }, 2000);
      await settleViewer(other, 2000, () => {
        other.props.xColumnNames = ['RACE'];
        other.props.yColumnNames = ['SEX'];
        other.props.viewerType = 'Bar chart';
      });
      const savedPair = readPair();

      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;

      try {
        await w.__settled('grok.events.onViewerAdded', () => { tv.addViewer('Scatter plot'); }, 1000);
        const viewersBefore = Array.from(tv.viewers).map((x: any) => x.type);

        const saved = await w.__findSaved(() => grok.dapi.layouts.find(layoutId), 2000);

        await w.__settled('grok.events.onViewLayoutApplied', () => { tv.loadLayout(saved); }, 4000);
        const viewersAfterLoad = Array.from(tv.viewers).map((x: any) => x.type);
        const restoredPair = readPair();

        const pie = trellises().find((vw: any) => vw.props.viewerType === 'Pie chart');
        const bar = trellises().find((vw: any) => vw.props.viewerType === 'Bar chart');
        let crossTalk: string | null = null;
        if (pie && bar) {
          await settleViewer(pie, 1500, () => { pie.props.viewerType = 'Histogram'; });
          crossTalk = bar.props.viewerType;
          await settleViewer(pie, 1500, () => { pie.props.viewerType = 'Pie chart'; });
        }

        bar?.close?.();
        await w.__poll(() => trellises().length, (n: number) => n === 1, 2500, 50);
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
        const w = window as any;
        const tv = grok.shell.tv;
        const tp = Array.from(tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
        tp.props.onClick = 'Select';
        await w.__settled('df.onSelectionChanged', () => tv.dataFrame.selection.setAll(false), 300);

        const captured: {mc: Record<string, string> | null} = {mc: null};
        const sub = tp.onEvent('d4-trellis-plot-current-cell-changed').subscribe((arg: any) => {
          const mc = arg?.args?.matchCondition ?? arg?.matchCondition ?? null;
          if (mc) captured.mc = Object.fromEntries(Object.entries(mc).map(([k, val]) => [k, String(val)]));
        });

        const cell = document.querySelector('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell') as HTMLElement;
        const rect = cell.getBoundingClientRect();
        const o = {bubbles: true, cancelable: true, view: window, button: 0,
          clientX: rect.left + 8, clientY: rect.top + 8};

        await w.__settled('df.onSelectionChanged', () => {
          cell.dispatchEvent(new MouseEvent('mousedown', o));
          cell.dispatchEvent(new MouseEvent('mouseup', o));
          cell.dispatchEvent(new MouseEvent('click', o));
        }, 1200);
        sub?.unsubscribe?.();
        const selected = tv.dataFrame.selection.trueCount;

        let expectedRows = -1;
        if (captured.mc) {
          const df = tv.dataFrame;
          const pairs = Object.entries(captured.mc).map(([c, val]) => [df.col(c), val] as [any, string]);
          expectedRows = 0;
          for (let i = 0; i < df.rowCount; i++)
            if (pairs.every(([col, val]) => String(col.get(i)) === val)) expectedRows++;
        }

        await w.__settled('viewer:Trellis plot.onViewerRendered', () => { tp.props.rowSource = 'Selected'; }, 1000);
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
      const saved = await saveProjectViaApi(page, nameSelected);
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
        const w = window as any;
        let tp: any = null;
        let tpView: any = null;
        for (const view of grok.shell.tableViews)
          for (const vw of view.viewers) if (vw.type === 'Trellis plot') { tp = vw; tpView = view; }
        if (!tp) return {clearedCells: -1, cells: -1, selected: -1,
          rowSource: null as string | null, x: [] as string[]};

        const cellCount = () => { try { return tp.root ? tp.root.querySelectorAll('.d4-trellis-plot-cell').length : 0; } catch (_) { return -1; } };
        tpView.dataFrame.selection.setAll(false);
        const clearedCells = await w.__poll(cellCount, (n: number) => n === 0, 2000, 50);
        tp.props.rowSource = 'All';
        const cells = await w.__poll(cellCount, (n: number) => n === 16, 5000, 100);
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
        const w = window as any;
        const tv = grok.shell.tv;
        const tp = Array.from(tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
        tv.dataFrame.selection.setAll(false);
        await w.__settled('viewer:Trellis plot.onViewerRendered', () => {
          tp.props.xColumnNames = ['SEX', 'CONTROL'];
          tp.props.yColumnNames = ['RACE'];
          tp.props.viewerType = 'Pie chart';
          tp.props.onClick = 'Select';
        }, 1500);
        await w.__settled('viewer:Trellis plot.onViewerRendered', () => { tp.props.rowSource = 'Selected'; }, 1500);
        return {selected: tv.dataFrame.selection.trueCount, rowSource: tp.props.rowSource};
      });

      expect(pre.selected).toBe(0);
      expect(pre.rowSource).toBe('Selected');

      const saved = await saveProjectViaApi(page, nameEmpty);
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
    page.off('pageerror', onPageError);
    page.off('console', onConsole);
    await v.closeAllAndWait(page);
  }

  v.finishSpec();
});
