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

// Benign dev-build noise, verified live against a Row Source = All control: the NullError
// balloon on project reopen and the caught ProjectMeta.publish trace on save. GROK-19902 is an
// UNCAUGHT page error on the Selected REOPEN path, so its guard is the pageerror floor.
const isBenignError = (text: string) =>
  /Failed to load resource/.test(text) || /404 \(\)/.test(text) || /favicon/.test(text) ||
  /Unable to find element in cloned iframe/.test(text) ||
  /NullError: method not found: '\w+' on null/.test(text) ||
  /ProjectMeta\.publish/.test(text) || /project_meta\.dart/.test(text);

// Re-called mid-test by Scenario 3 Step 12, which needs a PROJECT-LESS workspace: otherwise the
// ribbon Save re-saves the project Step 11 reopened instead of creating that step's own project.
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
  await page.waitForTimeout(1500);
}

// GROK-19902 crashed on the Row Source = Selected reopen path, so the grading surface is
// "restored with its persisted config and no uncaught page error", never a same-property echo.
// xCategoriesCount / yCategoriesCount hold in both the empty and the filled state.
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
    // Settle window: waits for the full 4 x 4 product, not the first non-zero count, so a
    // half-built grid is never read as a lost split column. Both reopen paths run it out (the
    // restored selection is always empty [DOM 2026-08-12]) — that IS the "stayed empty" evidence.
    let cells = 0;
    for (let t = 0; t < 10 && tp; t++) {
      try { cells = tp.root ? tp.root.querySelectorAll('.d4-trellis-plot-cell').length : 0; }
      catch (_) { cells = 0; }
      if (cells === 16) break;
      await new Promise((r) => setTimeout(r, 1000));
    }
    // Sampled next to the cell count, not earlier: under Row Source = Selected the selection
    // decides how many cells are owed, so the two numbers must come from the same moment.
    let selected = -1;
    try { selected = tpView ? tpView.dataFrame.selection.trueCount : -1; } catch (_) { selected = -1; }
    let xCats = -1;
    let yCats = -1;
    if (tp) {
      try { xCats = tp.xCategoriesCount; } catch (_) { /* Dart getter may throw pre-render */ }
      try { yCats = tp.yCategoriesCount; } catch (_) { /* Dart getter may throw pre-render */ }
    }
    // Diagnostics only — the box is NOT a grading input: 573x1000 was measured next to a
    // legitimately empty grid. It rides into the failure messages so a future regression
    // arrives with numbers.
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
  // Printed, never asserted: on a passing run the three numbers that drive the grading
  // (rowSource, selection count, cell count) appear in no assertion message, and their
  // absence once forced the reading to be reconstructed from the pw:api call log.
  console.log(`[trellis reopen] rowSource=${reading.rowSource} onClick=${reading.onClick} ` +
    `viewerType=${reading.viewerType} x=${JSON.stringify(reading.x)} y=${JSON.stringify(reading.y)} ` +
    `cats=${reading.xCats}x${reading.yCats} selected=${reading.selected} cells=${reading.cells} ` +
    `root=${reading.rootW}x${reading.rootH}`);
  return reading;
}

// Graded against the DATA STATE, never the viewer's box size: Row Source = Selected with an
// empty selection has no rows to plot, so a zero-cell grid is CORRECT — intended behaviour,
// operator ruling 2026-08-12 (no ticket). Any other state owes 4 x 4 = 16 (8 = CONTROL dropped).
//
// Only the empty case is carried: both round-trips restore an EMPTY selection [DOM 2026-08-12],
// so a "restored selection is non-empty" branch would be dead code. The state is re-checked
// below instead — a caller landing anywhere else fails loudly.
//
// Measured live, fresh viewer, X = [SEX, CONTROL] / Y = [RACE] / Pie chart [DOM 2026-08-12]:
//   Filtered, 0 selected  -> 16 cells      Selected, 0 selected -> 0 cells
//   Selected, all 5850    -> 16 cells      All, 0 selected      -> 16 cells
//   grid container 505x837 in every state; restored root 573x1000 while painting nothing
//
// WITHDRAWN as false by those numbers: "narrow strip below the render threshold" and
// "0 is admissible only under a zero-area box".
async function expectRestoredEmptyGridWithLiveness(page: Page,
  r: {cells: number; rowSource: string | null; selected: number;
    x: string[]; rootW: number; rootH: number}): Promise<void> {
  // Precondition of the rule below: an empty grid is only correct in exactly this state.
  expect(r.rowSource,
    `restored trellis reads Row Source = ${r.rowSource}, but the empty-grid rule graded here only holds for 'Selected'`).toBe('Selected');
  expect(r.selected,
    `the reopened view came back with ${r.selected} rows selected — every reopen path in this spec was measured to restore an EMPTY selection [DOM 2026-08-12], so a non-zero count means the recorded fact changed and the grading rule has to be re-derived before this grid can be graded`).toBe(0);
  expect(r.cells,
    `restored trellis painted ${r.cells} cells under Row Source = Selected with an EMPTY selection — that state has no rows to plot, so the grid must be empty; root measured ${r.rootW}x${r.rootH}`).toBe(0);

  // Liveness witness — separates "correctly empty" from "restored dead". Driven through the
  // SELECTION rather than by rewriting rowSource, so the restored Selected wiring is exercised
  // instead of bypassed.
  const revived = await page.evaluate(async () => {
    let tp: any = null;
    let tpView: any = null;
    for (const view of grok.shell.tableViews)
      for (const vw of view.viewers) if (vw.type === 'Trellis plot') { tp = vw; tpView = view; }
    if (!tp) return {cells: -1, selected: -1, x: [] as string[]};
    tpView.dataFrame.selection.setAll(true);
    await new Promise((res) => setTimeout(res, 2500));
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
  // The filled grid must come from the RESTORED split columns, not from a viewer that rebuilt
  // itself with defaults on the way.
  expect(revived.x).toEqual(r.x);
}

test('Trellis plot: split columns, inner-type switching, persistence', async ({page}) => {
  test.setTimeout(1_500_000);
  page.setDefaultTimeout(120_000);

  // Attached before any viewer add or project action: errors raised during a project reopen
  // fire after the open call resolves, so a later attach would miss them.
  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error' && !isBenignError(m.text())) consoleErrors.push(m.text()); });

  await loginToDatagrok(page);
  await page.waitForTimeout(5000);

  // Window marks for the github-964 floor below, taken AFTER login so the bootstrap noise
  // preceding the viewer add falls outside the graded window.
  const addPageErrBefore = pageErrors.length;
  const addConsoleErrBefore = consoleErrors.length;

  await buildDemogTrellis(page);

  // #### Scenario 1 Step 1 (github-964): the viewer is addable at all AND the add raises no
  // error. Deferring the error half to the Step 7 / 11 / 12 reopen floors would leave an error
  // raised while ADDING the viewer — what github-964 broke — permanently unread.
  await softStep('Scenario 1 Step 1', async () => {
    await expect(page.locator('[name="viewer-Trellis-plot"]')).toHaveCount(1);
    // Gated on the window opened just above, never on the whole collector: login and
    // platform-bootstrap noise already sits in both arrays, and the console side is filtered a
    // second time by isBenignError at push time.
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

  // #### Scenario 1 Step 4: cell count equals the category product (SEX 2 x RACE 4 = 8)
  await softStep('Scenario 1 Step 4', async () => {
    const counts = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      await new Promise((r) => setTimeout(r, 1500));
      return {xCat: tp.xCategoriesCount, yCat: tp.yCategoriesCount};
    });
    // The DOM count may be asserted against the Dart-computed product only while both axis
    // counts stay under the 7.5-category viewport clamp (SEX 2, CONTROL 2, RACE 4 all do);
    // past the clamp the grid renders fewer cells than the product.
    await expect(cellLocator).toHaveCount(counts.xCat * counts.yCat);
    await expect(cellLocator).toHaveCount(8);
    expect(counts).toEqual({xCat: 2, yCat: 4});
  });

  // #### Scenario 1 Step 6: CONTROL added as the second X column through the (+) control, with
  // the GROK-19673 guard on the SAME actuation — hover previews, Escape reverts, only the click
  // commits, and that commit is the 16 cells [DOM 2026-08-12]. A props write proves nothing here.
  await softStep('Scenario 1 Step 6', async () => {
    // The (+) `[name="add-x-column"]` (trellis_plot.md:19,43) is a ColumnComboBox in its empty
    // `d4-column-selector-cross` form: it opens `.d4-column-selector-backdrop` on MOUSEDOWN, not
    // on click (column_combo_box.md pitfall 1) — hence the real mouse down/up below.
    const backdrop = page.locator('.d4-column-selector-backdrop');
    const openAddXPopup = async (): Promise<boolean> => {
      const plus = page.locator('[name="viewer-Trellis-plot"] [name="add-x-column"]').first();
      // Bounded explicitly: a missing (+) must fail in seconds, not on the 120s page default.
      const box = await plus.boundingBox({timeout: 8000}).catch(() => null);
      if (!box) return false;
      await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
      await page.mouse.down();
      await page.mouse.up();
      return backdrop.waitFor({timeout: 6000}).then(() => true).catch(() => false);
    };
    // Canvas ColumnGrid — rows are NOT DOM elements, so the row is found by real-mouse scan of
    // the UNFILTERED list reading `.d4-tooltip` [SRC column_grid.dart initColumnTooltips]; never
    // the search box, which corrupts the Escape semantics under test. No arithmetic fallback.
    //
    // The backdrop rect is re-read on every call: after a dock relayout the popup can open
    // elsewhere and a stale rect scans empty space.
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
        // First matching y only — the singleton tooltip lingers after the pointer leaves the
        // row, so later scan steps would false-match below it.
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

      // HOVER phase. The preview WRITES THROUGH to the look while the popup is open
      // (props.xColumnNames transiently reads the hovered column), so non-commit CANNOT be
      // asserted during the dwell — the Escape rollback in (ii) is.
      const dwell: {x: string[]; cells: number}[] = [];
      for (let i = 0; i < 4; i++) {
        await page.mouse.move(first.hx, first.rowY!);
        await page.waitForTimeout(220);
        dwell.push(await readXState());
      }
      // (i) Preview is alive: no dwell sample moving off the committed 8 cells or binding the
      // hovered column as a second split is the preview-regression signal.
      expect(dwell.some((s) => s.cells !== 8 || s.x.length !== 1),
        'hover preview never rebuilt the grid during the dwell — preview-regression signal (GROK-19673 class)').toBe(true);

      // (ii) Escape WITH the pointer still on the row — the DECISIVE non-commit signal: the
      // preview rolls back while the row is still hovered. The second Escape is a fallback for
      // the first being swallowed by the popup's focused search input.
      await page.keyboard.press('Escape');
      const closedOnFirstEsc = await backdrop
        .waitFor({state: 'detached', timeout: 2000}).then(() => true).catch(() => false);
      if (!closedOnFirstEsc) {
        await page.keyboard.press('Escape');
        await backdrop.waitFor({state: 'detached', timeout: 5000});
      }
      await expect(cellLocator).toHaveCount(8);
      expect((await readXState()).x).toEqual(['SEX']);

      // (iii) COMMIT is the CLICK: the row has to be re-found after the reopen (the fresh rect
      // rule above), and the new binding must STAY after the popup closes — a preview that
      // survives only while hovered is not a commit.
      expect(await openAddXPopup(),
        'the (+) add-X-column control did not reopen the column selector popup').toBe(true);
      const second = await findRowByHover(/\bCONTROL\b/);
      expect(second.rowY).not.toBeNull();
      await page.mouse.click(second.hx, second.rowY!);
      await page.waitForTimeout(1200);
      const committed = await page.evaluate(() => {
        const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
        return {x: [...tp.props.xColumnNames], xCat: tp.xCategoriesCount, yCat: tp.yCategoriesCount};
      });
      expect(committed.x).toEqual(['SEX', 'CONTROL']);
      expect({xCat: committed.xCat, yCat: committed.yCat}).toEqual({xCat: 4, yCat: 4});
      await expect(cellLocator).toHaveCount(16);
    } finally {
      // Runs even on a guard failure: a surviving popup is dismissed so one broken guard cannot
      // cascade into every following scenario. The column binding is deliberately NOT forced
      // back — writing it would be the very actuation under test.
      await page.keyboard.press('Escape').catch(() => {});
      await page.evaluate(() =>
        document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}))).catch(() => {});
    }
  });

  // #### Scenario 1 Step 7: remove the second X column through the BLANK first row of the picker
  // opened from the CONTROL slot — it drops exactly that one split column (refdoc "Removing a
  // Split Column"). With Step 6's (+) add this closes one UI round-trip, no props write.
  //
  // NOT the whole-axis route (right-click a split selector -> `Reset X columns`): that clears
  // every X column at once, a different operation than dropping one.
  await softStep('Scenario 1 Step 7', async () => {
    // Layout model of the canvas picker — every constant is a product-source value:
    //   BORDER 1 (column_grid.dart:379)                     ROW 16 (column_grid.dart:361-363)
    //   HEADER 20 horzColLabelsHeight (grid_look.dart:89)   PAD 10 autoSize (grid_core.dart:1718)
    //
    // => height is 2 + 20 + 16 * rows + 10, blank FIRST body row at [top + 21, top + 37]. The
    //    live 129x160 rect [DOM 2026-08-11] is exactly this model at 8 body rows (blank +
    //    demog's 7 eligible categorical columns).
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

    // The CONTROL slot is found STRUCTURALLY, never by position: `_lookToSelectors` puts the X
    // selectors and the (+) in one host div, so the slot is the (+)'s sibling labelled CONTROL —
    // which excludes the legend's `div-column-combobox-category` mirror and the Y selectors.
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
      // Gate the model BEFORE the pointer moves: the blank-row offset holds only while the popup
      // measures 2 + 20 + 16 * rows + 10. Otherwise fail with the measured numbers — the input
      // for re-deriving the offsets, NOT a licence to nudge the click until it lands.
      const modelRows = (rect.height - (2 * BORDER + HEADER_H + AUTOSIZE_PAD)) / ROW_H;
      expect(Number.isInteger(modelRows) && modelRows >= 2,
        `column picker measured ${rect.width}x${rect.height}: its height is not 32 + 16*rows, so the layout model that locates the blank first row no longer holds — re-derive the offsets from these numbers`).toBe(true);

      const hx = rect.left + Math.min(60, rect.width / 2);
      const hy = rect.top + BORDER + HEADER_H + ROW_H / 2;

      // ONE move (page.mouse.move jumps in a single mousemove), straight onto the blank row: the
      // removal fires on HOVER (allowPreview assigns a null column — trellis_plot_core.dart:785,
      // 809-812), so a live eight-row walk emptied the axis instead of dropping one [DOM 2026-08-11].
      await page.mouse.move(hx, hy);
      await page.waitForTimeout(900);
      const previewed = await readXState();
      // Targeting check, not the evidence: a miss (the header above the row, or USUBJID below
      // it) is dismissed with Escape by the finally, never clicked.
      expect(previewed.x,
        `hovering the blank first row at (${Math.round(hx)}, ${Math.round(hy)}) of a ${rect.width}x${rect.height} popup left the X split columns as ${JSON.stringify(previewed.x)} — the pointer did not land on the blank row, so nothing is committed`).toEqual(['SEX']);

      // The removal has to SURVIVE the popup close — a preview that only holds while the row is
      // hovered is not a commit.
      await page.mouse.click(hx, hy);
      await backdrop.waitFor({state: 'detached', timeout: 6000});
      await page.waitForTimeout(1500);

      const after = await readXState();
      // Both witnesses are mandatory: a DOM count back at 8 says nothing about WHICH column
      // was dropped.
      expect(after.x).toEqual(['SEX']);
      expect({xCat: after.xCat, yCat: after.yCat}).toEqual({xCat: 2, yCat: 4});
      await expect(cellLocator).toHaveCount(8);
    } finally {
      // Dismisses a surviving popup (Escape restores the column the popup opened with, so a
      // mis-targeted hover is rolled back). The split columns are deliberately NOT written back
      // — that write is the actuation under test.
      await page.keyboard.press('Escape').catch(() => {});
      await page.evaluate(() =>
        document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}))).catch(() => {});
    }
  });

  // #### Scenario 2 Step 3: each inner type renders distinct per-cell data
  await softStep('Scenario 2 Step 3', async () => {
    // Two independent signals per type: (a) 'd4-trellis-plot-viewer-type-changed' fires with the
    // new type name — a prop set->get echo does NOT prove the switch was processed; (b) two cells
    // with DIFFERENT split categories render DIFFERENT pixels (GROK-19633 / GROK-19890 class).
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

      // The event is a Dart event-bus event and does NOT fire on a JS-API viewerType assignment —
      // only the UI-driven switch emits it. Hence `[name="viewer selector"]`: a combo expanding on
      // synthetic mousedown, items `.d4-combo-drop-down .d4-list-item` [DOM 2026-08-05].
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
        (item?.closest('.d4-list-item') as HTMLElement | null)?.click();
        await new Promise((r) => setTimeout(r, 1500));
        sub?.unsubscribe?.();

        const populated: number[] = [];
        const cells = root.querySelectorAll('.d4-trellis-plot-cell');
        for (let i = 0; i < cells.length && populated.length < 2; i++)
          if (cells[i].querySelector('canvas')) populated.push(i);
        const a = cellSig(populated[0]);
        const b = cellSig(populated[1]);
        out[t] = {eventType, distinct: a !== null && b !== null && a !== b};
      }
      return out;
    });
    for (const t of ['Scatter plot', 'Bar chart', 'Histogram', 'Pie chart']) {
      expect(perType[t].eventType).toBe(t);
      expect(perType[t].distinct).toBe(true);
    }
  });

  // #### Scenario 2 Step 5: switch through a no-data config and back re-renders (github-3015)
  await softStep('Scenario 2 Step 5', async () => {
    // The stale-render class (github-3015) is the returned Pie-chart cells still showing the
    // INTERMEDIATE Bar-chart frame. A global non-white pixel count cannot catch it, so the
    // evidence is each cell's own canvas moving off the Bar-chart frame it was hashed at.
    const result = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;

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

      tp.props.viewerType = 'Pie chart';
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      await new Promise((r) => setTimeout(r, 1500));

      // An empty Y split is the no-data configuration the type switch must survive.
      tp.props.yColumnNames = [];
      tp.props.viewerType = 'Bar chart';
      await new Promise((r) => setTimeout(r, 1200));

      // Repopulated but still on Bar chart: the exact frame the stale-render class would leave
      // behind if the return to Pie chart did not re-render.
      tp.props.yColumnNames = ['RACE'];
      await new Promise((r) => setTimeout(r, 1500));
      const idxs = populatedCellIdxs(2);
      const barA = cellHash(idxs[0]);
      const barB = cellHash(idxs[1]);

      tp.props.viewerType = 'Pie chart';
      await new Promise((r) => setTimeout(r, 1800));
      const pieA = cellHash(idxs[0]);
      const pieB = cellHash(idxs[1]);

      return {
        hashesRead: barA !== null && barB !== null && pieA !== null && pieB !== null,
        // Anti-vacuity companion: two identically blank cells would satisfy reRendered just as
        // well, so per-cell binding is checked alongside it.
        distinctAfter: pieA !== null && pieB !== null && pieA !== pieB,
        reRendered: barA !== pieA && barB !== pieB,
      };
    });
    expect(result.hashesRead).toBe(true);
    expect(result.distinctAfter).toBe(true);
    expect(result.reRendered).toBe(true);
  });

  // #### Scenario 2 Step 7: hide control panel; type selector goes away, viewerType preserved
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
    // Hiding the panel does NOT remove the selector node — it collapses to a zero-size box under
    // a display:none parent, so a presence check would pass and only size tells.
    expect(result.visibleAfterHide).toBe(false);
    // Read through a different channel than the one hidden: the inner type must survive the
    // panel going away, not merely be unreadable.
    expect(result.typeAfterHide).toBe('Pie chart');
  });

  // #### Scenario 3 Step 6: layout save/re-apply restores BOTH trellis plots' exact
  // configuration and they stay INDEPENDENT afterwards (GROK-15494). A single-viewer probe
  // cannot see that bug — a shared look satisfies it — hence the second trellis.
  await softStep('Scenario 3 Step 6', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const trellises = () => Array.from(tv.viewers).filter((x: any) => x.type === 'Trellis plot') as any[];
      // Sorted by inner type: loadLayout need not restore the two viewers in creation order.
      const readPair = () => trellises().map((v: any) => ({
        x: [...v.props.xColumnNames] as string[],
        y: [...v.props.yColumnNames] as string[],
        type: v.props.viewerType as string,
        xCat: v.xCategoriesCount as number,
        yCat: v.yCategoriesCount as number,
      })).sort((a, b) => a.type.localeCompare(b.type));

      const tp = trellises()[0];
      tp.props.xColumnNames = ['SEX', 'CONTROL'];
      tp.props.yColumnNames = ['RACE'];
      tp.props.viewerType = 'Pie chart';
      await new Promise((r) => setTimeout(r, 1500));

      // Second trellis: transposed split columns and a different inner type, both under the
      // viewport clamp (RACE 4 x SEX 2), so its restored product is checkable too.
      const other = tv.addViewer('Trellis plot');
      await new Promise((r) => setTimeout(r, 2000));
      other.props.xColumnNames = ['RACE'];
      other.props.yColumnNames = ['SEX'];
      other.props.viewerType = 'Bar chart';
      await new Promise((r) => setTimeout(r, 2000));
      const savedPair = readPair();

      // The LAYOUT half persists through the JS API by rule, not omission: the grok-browser skill
      // ("Saving and restoring a layout") prescribes grok.dapi.layouts over the SAVE ribbon, which
      // does not guarantee the layout can be found and re-applied. The PROJECT half is the opposite.
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise((r) => setTimeout(r, 1200));

      try {
        // The standalone Scatter plot is the marker: absent from the saved layout, so
        // re-applying that layout must remove it again.
        tv.addViewer('Scatter plot');
        await new Promise((r) => setTimeout(r, 1000));
        const viewersBefore = Array.from(tv.viewers).map((x: any) => x.type);

        const saved = await grok.dapi.layouts.find(layoutId);
        tv.loadLayout(saved);
        await new Promise((r) => setTimeout(r, 4000));
        const viewersAfterLoad = Array.from(tv.viewers).map((x: any) => x.type);
        const restoredPair = readPair();

        // GROK-15494's actual contract: mutate ONE restored trellis and read the OTHER — a
        // shared-look leak shows up exactly here.
        const pie = trellises().find((v: any) => v.props.viewerType === 'Pie chart');
        const bar = trellises().find((v: any) => v.props.viewerType === 'Bar chart');
        let crossTalk: string | null = null;
        if (pie && bar) {
          pie.props.viewerType = 'Histogram';
          await new Promise((r) => setTimeout(r, 1500));
          crossTalk = bar.props.viewerType;
          pie.props.viewerType = 'Pie chart';
          await new Promise((r) => setTimeout(r, 1500));
        }

        // Closed before the project half: the reopen path reads "the" restored trellis, and two
        // of them make that read ambiguous.
        bar?.close?.();
        await new Promise((r) => setTimeout(r, 2500));
        const survivors = readPair();
        return {viewersBefore, viewersAfterLoad, savedPair, restoredPair, crossTalk, survivors};
      } finally {
        // The layout probe delete survives any mid-step throw (find/loadLayout).
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
    // The surviving single trellis keeps the configuration the project half depends on.
    expect(result.survivors).toEqual([
      {x: ['SEX', 'CONTROL'], y: ['RACE'], type: 'Pie chart', xCat: 4, yCat: 4},
    ]);
    await expect(cellLocator).toHaveCount(16);
  });

  // Whole-project persistence (Steps 7-12). The probe project ids are declared out here so the
  // trailing finally can delete them whichever step fails — a mid-step throw would otherwise
  // orphan real projects on the server.
  const nameSelected = `zz-trellis-p0-selected-${Date.now()}`;
  const nameEmpty = `zz-trellis-p0-empty-${Date.now()}`;
  let idSelected: string | null = null;
  let idEmpty: string | null = null;
  // Carried from Step 7 into Step 11: the rows selected when the first project was saved.
  // Grading the reopened selection against this number is the only thing that keeps Step 11's
  // reading from degenerating into a repeat of Step 12's empty state.
  let savedSelectionCount = -1;
  try {
    // #### Scenario 3 Step 7: wire the Selected mode (On Click = Select, a real
    // corner-click selection, Row Source = Selected), then drive the real ribbon Save.
    await softStep('Scenario 3 Step 7', async () => {
      const pre = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const tp = Array.from(tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
        tp.props.onClick = 'Select';
        await new Promise((r) => setTimeout(r, 800));
        tv.dataFrame.selection.setAll(false);
        await new Promise((r) => setTimeout(r, 300));
        // The click fires d4-trellis-plot-current-cell-changed with args.matchCondition, a
        // column -> category map. Reading it makes the expected selection count EXACT without
        // re-deriving the cell's position arithmetically.
        const captured: {mc: Record<string, string> | null} = {mc: null};
        const sub = tp.onEvent('d4-trellis-plot-current-cell-changed').subscribe((arg: any) => {
          const mc = arg?.args?.matchCondition ?? arg?.matchCondition ?? null;
          if (mc) captured.mc = Object.fromEntries(Object.entries(mc).map(([k, v]) => [k, String(v)]));
        });
        // Corner-click (~8px inset): a centered click is intercepted by the inner Pie chart
        // canvas and never reaches the trellis cell handler. Synthetic MouseEvents do drive this
        // path (live-verified: trueCount 0 -> >0).
        const cell = document.querySelector('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell') as HTMLElement;
        const rect = cell.getBoundingClientRect();
        const o = {bubbles: true, cancelable: true, view: window, button: 0,
          clientX: rect.left + 8, clientY: rect.top + 8};
        cell.dispatchEvent(new MouseEvent('mousedown', o));
        cell.dispatchEvent(new MouseEvent('mouseup', o));
        cell.dispatchEvent(new MouseEvent('click', o));
        await new Promise((r) => setTimeout(r, 1200));
        sub?.unsubscribe?.();
        const selected = tv.dataFrame.selection.trueCount;
        // Structural row count of the clicked combination (no filter is active here).
        let expectedRows = -1;
        if (captured.mc) {
          const df = tv.dataFrame;
          const pairs = Object.entries(captured.mc).map(([c, v]) => [df.col(c), v] as [any, string]);
          expectedRows = 0;
          for (let i = 0; i < df.rowCount; i++)
            if (pairs.every(([col, v]) => String(col.get(i)) === v)) expectedRows++;
        }
        tp.props.rowSource = 'Selected';
        await new Promise((r) => setTimeout(r, 1000));
        return {selected, expectedRows, matchCondition: captured.mc,
          rowSource: tp.props.rowSource, onClick: tp.props.onClick};
      });
      // The click has to happen BEFORE rowSource is switched to Selected: from that point the
      // grid plots the selection, so with nothing selected there is no cell left to click
      // [DOM 2026-08-12]. Hence select first, switch after.
      expect(pre.matchCondition,
        'the corner click did not fire d4-trellis-plot-current-cell-changed — the click never reached the trellis cell handler').not.toBeNull();
      // Anti-vacuity: a structurally empty combination would let "selected == expected" pass as
      // 0 == 0, the silent no-op this step exists to catch.
      expect(pre.expectedRows,
        `clicked combination ${JSON.stringify(pre.matchCondition)} holds no rows — pick a populated probe cell`).toBeGreaterThan(0);
      expect(pre.selected).toBe(pre.expectedRows);
      expect(pre.onClick).toBe('Select');
      expect(pre.rowSource).toBe('Selected');
      // Recorded after the assertions above so Step 11's baseline can only be a verified number.
      savedSelectionCount = pre.selected;

      const pageErrBefore = pageErrors.length;
      const saved = await saveProjectViaUI(page, nameSelected);
      idSelected = saved.projectId;
      expect(idSelected).toBeTruthy();
      // Only the UNCAUGHT floor applies to the save: the publish stack trace this flow logs is
      // benign (see isBenignError). The GROK-19902 crash guard is at Step 11.
      expect(pageErrors.slice(pageErrBefore)).toEqual([]);
    });

    // #### Scenario 3 Steps 8-11: Close All, reopen the saved project, read the restored
    // Trellis configuration back off the restored viewer.
    //
    // OWNS: the project was saved with a LIVE selection, so this is the section's only reading of
    // a selection across a project round-trip. It does not survive [DOM 2026-08-12], and that is
    // INTENDED — operator ruling 2026-08-12, no ticket. Step 12 saves nothing selected.
    await softStep('Scenario 3 Step 11', async () => {
      const pageErrBefore = pageErrors.length;
      const errorsBefore = consoleErrors.length;
      const r = await reopenAndReadTrellis(page, idSelected!);
      // GROK-19902 guard: the crash lived on the Row Source = Selected REOPEN path, so this is
      // where the uncaught-page-error floor earns its keep.
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
      // The save-time count is sanity-checked first so the comparison below can never degenerate
      // into 0 == 0 if the Step 7 corner click ever stops selecting anything.
      expect(savedSelectionCount,
        'Step 7 did not record a live selection at save time — without it Step 11 cannot say anything about whether a selection survives the project round-trip').toBeGreaterThan(0);
      expect(r.selected,
        `the project was saved with ${savedSelectionCount} rows selected and reopened with ${r.selected} — the recorded behaviour is that a selection does NOT survive the project round-trip [DOM 2026-08-12], so a non-zero count means the product now persists it and the grading of both this step and Step 12 has to be re-derived`).toBe(0);
      // The empty grid is the CONSEQUENCE of that loss: Row Source = Selected came back intact,
      // and with no selection there is nothing left to plot.
      await expectRestoredEmptyGridWithLiveness(page, r);

      // Second probe, NOT a repeat of the witness above: that one proves the viewer draws WITH a
      // selection, this one that it draws WITHOUT one — the empty grid is ambiguous until both
      // have run. Runs AFTER the loss is read, so rewriting rowSource cannot mask it.
      const allSource = await page.evaluate(async () => {
        let tp: any = null;
        let tpView: any = null;
        for (const view of grok.shell.tableViews)
          for (const vw of view.viewers) if (vw.type === 'Trellis plot') { tp = vw; tpView = view; }
        if (!tp) return {clearedCells: -1, cells: -1, selected: -1,
          rowSource: null as string | null, x: [] as string[]};
        // Clearing the selection first is the anti-vacuity half: the previous witness left every
        // row selected, so without it the grid is already full and the 16 below says nothing
        // about the source switch.
        tpView.dataFrame.selection.setAll(false);
        await new Promise((res) => setTimeout(res, 2000));
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

    // #### Scenario 3 Step 12: the round-trip whose SAVE-TIME state is an EMPTY selection under
    // Row Source = Selected — the viewport-restore path GROK-19902 crashed on. A project-less
    // workspace is rebuilt first so the ribbon Save creates a second project.
    //
    // Distinct from Step 11 by the state under TEST, not the state after the reopen: both land on
    // "Selected + nothing selected", but there the emptiness is a LOSS, here a faithful CARRY —
    // and only a carry holds the "empty at save time" precondition the crash needs.
    await softStep('Scenario 3 Step 12', async () => {
      await buildDemogTrellis(page);
      const pre = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const tp = Array.from(tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
        tv.dataFrame.selection.setAll(false);
        tp.props.xColumnNames = ['SEX', 'CONTROL'];
        tp.props.yColumnNames = ['RACE'];
        tp.props.viewerType = 'Pie chart';
        // onClick = Select and rowSource = Selected coexist — only the Filter + Filtered
        // combination self-corrects (live-verified).
        tp.props.onClick = 'Select';
        await new Promise((r) => setTimeout(r, 500));
        tp.props.rowSource = 'Selected';
        await new Promise((r) => setTimeout(r, 1500));
        return {selected: tv.dataFrame.selection.trueCount, rowSource: tp.props.rowSource};
      });
      // Nothing selected at save time is the exact shape the crash needed — an accidental
      // non-empty selection here silently retires the guard.
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
      // The emptiness read back here is the saved state being CARRIED, not lost (that is Step
      // 11's reading). A selection that reappeared would mean this project no longer restores
      // the shape the crash lived on — hence the assert.
      expect(r.selected,
        `the reopened project came back with ${r.selected} rows selected — it was saved with none, so this step no longer covers the empty-selection restore shape`).toBe(0);
      await expectRestoredEmptyGridWithLiveness(page, r);
    });
  } finally {
    // #### Scenario 3 Step 13: both probe projects are removed even on failure.
    if (idSelected) await deleteProjectWithCleanup(page, {projectId: idSelected});
    if (idEmpty) await deleteProjectWithCleanup(page, {projectId: idEmpty});
    await page.evaluate(() => grok.shell.closeAll()).catch(() => {});
  }

  v.finishSpec();
});
