/**
 * F0 — platform verification gate for the per-subject results grid
 * (`nca-results-grid-columns` AC-0.1 .. AC-0.5).
 *
 * **Why this file exists.** The build plan's F2 (pin the Subject column) rests on
 * `GridColumn.pin()`, an API that EXISTS in the typings but has a documented,
 * version-dependent defect surface (Datagrok Community, Jul 2024 – Jan 2025):
 * background-colour loss on the pinned column, pin state not surviving a reload,
 * and sort desync. The research cycle's own R1 review refuted the "the typings
 * say it exists, therefore it works" inference. So nothing in F2 may be built
 * until these run green on the REAL host.
 *
 * **The host under test is the EMBEDDED one, and it is the only live one.**
 * `swapOrDockResults` (nca-mapping-pane.ts) builds `DGNS.Viewer.grid(results)`
 * and docks it — it never opens a `TableView` on the results frame. The sibling
 * path `openResultsView` has ZERO callers anywhere in `package/` (blueprint
 * probe, exhaustive search), so the docked-TableView mode is unreachable in the
 * shipped app. These tests therefore reproduce the embedded construction exactly
 * — `Viewer.grid` + `dockManager.dock(FILL, mainGridNode)` — rather than the
 * cheaper `addTableView(results).grid`, which would verify a host no user meets.
 *
 * **Measured, not inferred (rule 16).** Every assertion below reads platform
 * state back through a real API (`meta.colors.getColors()`, `GridColumn.idx`,
 * `DataFrame.toCsv`) on a frame produced by the REAL pipeline (theoph through
 * `runPipelineOrThrow`), not a synthetic stub.
 *
 * **Text-colour trial update (2026-09-02).** The results grid moved from background
 * FILL colour-coding (`col.meta.colors`) to per-cell TEXT colour, applied at render
 * time in `configureResultsGrid`'s `onCellPrepare` from the SAME `buildColourRules`
 * bands. Two consequences the colour assertions below now pin: the documented DG
 * `pin()` background-formatting-loss defect no longer applies (there is no stored fill
 * for pin() to drop — text colour is re-applied every render), and
 * `configureResultsDataFrame` must register NO `col.meta.colors` (so the banded columns
 * carry ≤ 1 distinct background colour — the check that the fill was actually removed).
 */

import {category, test, expect} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {
  configureResultsDataFrame, configureResultsGrid,
  colourBandProjections, applyResultsColumnProfile, resultsGridVisibleOrder, RESULTS_EXTENDED_TAG,
  ACCEPT_REASONS_COL,
} from '../views/results-view';
import {DEFAULT_RULES} from '../adapter/defaults';
import {nca} from '@datagrok-libraries/sci-comp';
import {detectPowerGridSupport, powerGridDrawsPills} from '../views/powergrid-support';
import {
  loadFixture, withDoseMg, theophContract, theophScalars, runPipelineOrThrow,
  indomethContract, indomethScalars, setupMultirouteRoute, ratMultirouteContract, ratScalars,
} from './test-utils';

const TIMEOUT = 60000;

/** The columns this grid paints, DERIVED from the one declaration that paints
 *  them — never a restated list. Three independent readers have now restated the
 *  banded set and produced three different wrong counts (GAP-23 said 6,
 *  nca-package-architecture.md said 4, the truth was 7); AD-2 makes deriving the
 *  house rule, and this guard follows it. `crossColumn` bands (MRT / VSS) are
 *  EXCLUDED here: they carry no DG rule, so `meta.colors` holds nothing for them
 *  — they are painted per-cell in `onCellPrepare` and cannot be read back this
 *  way. Their survival is covered by the DOM-level checks instead. */
function dgPaintedColumns(): string[] {
  const out = new Set<string>();
  for (const b of colourBandProjections())
    if (b.grid === 'results' && b.ruleKey !== null) out.add(b.column);
  return [...out];
}

/** Build the results grid the way production builds it: an embedded
 *  `Viewer.grid` FILL-docked onto the host view's main grid node — a TAB of the
 *  source-data grid, per the 2026-08-13 layout. The `swapOrDockResults`
 *  construction, reproduced rather than approximated. */
async function embeddedResultsGrid(): Promise<{
  view: DG.TableView; grid: DG.Grid; results: DG.DataFrame; close: () => void;
}> {
  const src = withDoseMg(await loadFixture('theoph.csv'));
  const {results} = runPipelineOrThrow(
    src, theophContract(src), theophScalars(), DEFAULT_RULES);
  configureResultsDataFrame(results);

  const view = grok.shell.addTableView(src);
  const grid = DG.Viewer.grid(results);
  configureResultsGrid(grid);
  const mainGridNode = view.dockManager.findNode(view.grid.root) ?? null;
  view.dockManager.dock(grid.root, DG.DOCK_TYPE.FILL, mainGridNode, 'NCA results');

  return {
    view, grid, results,
    close: () => {
      try {
        grid.close();
      } catch {
        // viewer already detached — teardown is best-effort
      }
      try {
        view.close();
      } catch {
        // view already closed — teardown is best-effort
      }
    },
  };
}

/** The COUNT of distinct per-row background colours DG paints for a column. With the results
 *  grid on TEXT colour (not background fill), a banded column carries NO background coding, so
 *  this is ≤ 1 (uniform default / empty) — the check that the fill was actually removed. Before
 *  the text-colour switch it was 2–3 (one ARGB per tier). Guarded: a column with no coding set
 *  may return an empty array. */
function distinctBackgroundColours(df: DG.DataFrame, name: string): number {
  const col = df.col(name);
  if (col === null) return 0;
  try {
    return new Set(Array.from(col.meta.colors.getColors())).size;
  } catch {
    return 0; // no coding registered — nothing to read back
  }
}

category('Results grid platform gate (F0)', () => {
  test('AC-0.1: pin() holds SUBJECT on the EMBEDDED grid (the only live host)', async () => {
    const h = await embeddedResultsGrid();
    try {
      const gc = h.grid.col('SUBJECT');
      expect(gc !== null, true, 'SUBJECT grid column resolves on the embedded grid');

      const idxBefore = gc!.idx;
      const canvasBefore = h.grid.root.querySelectorAll('canvas').length;

      // The call under test. A throw here is the fallback trigger: F2 degrades to
      // column ORDERING only (never a custom pinning renderer — Datagrok-native
      // viewers only).
      gc!.pin();

      const gcAfter = h.grid.col('SUBJECT');
      expect(gcAfter !== null, true,
        'SUBJECT still resolves through grid.col() after pin() — a pinned column ' +
        'that falls out of the column list would break every later gc.visible write');

      const idxAfter = gcAfter!.idx;
      const canvasAfter = h.grid.root.querySelectorAll('canvas').length;

      // Pinning must not disturb the frame. The whole compact-profile design
      // assumes the DataFrame is untouched by grid-side presentation.
      expect(h.results.columns.contains('SUBJECT'), true,
        'the underlying DataFrame still carries SUBJECT after pin()');

      gcAfter!.unpin();
      const idxUnpinned = h.grid.col('SUBJECT')?.idx ?? -1;

      // DECISIVE CHECK — "it did not throw" is barely stronger than reading the
      // typings, which is precisely the inference the research cycle refuted.
      // SUBJECT already sits at the left edge, so pinning it can move nothing and
      // its idx proves nothing either way. Pin a MIDDLE column instead: if the
      // platform genuinely pins, the column must relocate into the pinned zone at
      // the left; if pin() is a silent no-op on an embedded host, the index will
      // not move. This is the one measurement that distinguishes the two.
      const mid = h.grid.col('HALF_LIFE');
      expect(mid !== null, true, 'HALF_LIFE resolves (the relocation probe column)');
      const midBefore = mid!.idx;
      mid!.pin();
      const midAfterIdx = h.grid.col('HALF_LIFE')?.idx ?? -1;
      const movedLeft = midAfterIdx < midBefore;
      h.grid.col('HALF_LIFE')?.unpin();
      const midRestored = h.grid.col('HALF_LIFE')?.idx ?? -1;

      expect(movedLeft, true,
        `pin() did NOT relocate a middle column (HALF_LIFE idx ${midBefore} -> ` +
        `${midAfterIdx}) on the embedded Viewer.grid host. That means pin() is a ` +
        'silent no-op here and requirement (1) cannot be met by the native API — ' +
        'ship F2 as column ORDERING only and record the version evidence. Do NOT ' +
        'build a custom pinning renderer (Datagrok-native viewers only).');

      return `pin OK | SUBJECT idx before=${idxBefore} after=${idxAfter} ` +
        `unpinned=${idxUnpinned} | canvases ${canvasBefore}->${canvasAfter} | ` +
        `RELOCATION PROBE HALF_LIFE idx ${midBefore}->${midAfterIdx} ` +
        `(moved left=${movedLeft}), after unpin=${midRestored}`;
    } finally {
      h.close();
    }
  }, {timeout: TIMEOUT});

  test('AC-0.3: banded columns carry NO background fill (text-colour trial) + pin is safe', async () => {
    // The results grid moved from background colour-coding to per-cell TEXT colour (applied at
    // render time in configureResultsGrid's onCellPrepare from the SAME buildColourRules bands).
    // This test pins two consequences:
    //  (1) `configureResultsDataFrame` registers NO `col.meta.colors` — so each banded column shows
    //      ≤ 1 distinct background colour. A future accidental re-add of the fill (which would then
    //      DOUBLE-encode with the text colour) is caught here.
    //  (2) the documented DG pin() background-formatting-loss defect no longer applies: there is no
    //      stored fill for pin() to drop (text colour is re-applied every render), so pin/unpin is
    //      safe with nothing to lose.
    const h = await embeddedResultsGrid();
    try {
      const banded = dgPaintedColumns();
      expect(banded.length > 0, true,
        'the results grid declares at least one banded column (from buildColourRules)');

      for (const n of banded) {
        expect(distinctBackgroundColours(h.results, n) <= 1, true,
          `${n}: >1 distinct BACKGROUND colour — a fill is registered on this column, but the grid ` +
          'now colours the TEXT (configureResultsDataFrame must register no col.meta.colors here).');
      }

      // pin/unpin remains safe with no stored formatting to lose.
      const gc = h.grid.col('SUBJECT');
      expect(gc !== null, true, 'SUBJECT resolves');
      gc!.pin();
      gc!.unpin();
      expect(h.grid.col('SUBJECT') !== null, true, 'SUBJECT still resolves after pin/unpin');
      return `${banded.length} banded columns carry no background fill (text-colour); pin/unpin safe: ` +
        banded.join(', ');
    } finally {
      h.close();
    }
  }, {timeout: TIMEOUT});

  test('AC-0.2: the pin survives the re-Compute dataFrame swap', async () => {
    // `swapOrDockResults` re-binds the SAME grid instance to a fresh frame on every
    // Compute (`existingGrid.dataFrame = results`), then re-runs
    // configureResultsGrid. Grid-side state is documented as resettable across that
    // swap, which is exactly why configureResultsGrid is re-applied — so the real
    // question is whether re-applying the pin after the swap works, not whether the
    // pin survives on its own.
    const h = await embeddedResultsGrid();
    try {
      h.grid.col('SUBJECT')?.pin();

      // A second, independently computed frame — the re-Compute shape.
      const src2 = withDoseMg(await loadFixture('theoph.csv'));
      const {results: results2} = runPipelineOrThrow(
        src2, theophContract(src2), theophScalars(), DEFAULT_RULES);
      configureResultsDataFrame(results2);
      h.grid.dataFrame = results2;
      configureResultsGrid(h.grid);

      const gc = h.grid.col('SUBJECT');
      expect(gc !== null, true, 'SUBJECT resolves on the swapped-in frame');
      // Re-pinning after a swap must be safe whether or not the pin persisted:
      // `configureResultsGrid` will call it unconditionally on every Compute.
      gc!.pin();
      gc!.pin(); // idempotence — the swap path calls configure* on every Compute
      expect(h.grid.col('SUBJECT') !== null, true,
        'SUBJECT still resolves after a repeated pin() on the swapped frame');

      // No background colour-coding to re-register on the swapped frame — the grid colours the TEXT
      // in configureResultsGrid's onCellPrepare (re-run on every Compute), so the banded columns
      // carry ≤ 1 distinct background colour here too (the text-colour trial invariant).
      const banded = dgPaintedColumns();
      for (const n of banded) {
        expect(distinctBackgroundColours(results2, n) <= 1, true,
          `${n}: a background fill re-appeared on the swapped-in frame (should be text-colour only)`);
      }
      return `swap + re-pin OK (idempotent); ${banded.length} banded columns text-colour (no fill)`;
    } finally {
      h.close();
    }
  }, {timeout: TIMEOUT});

  test('AC-0.5: does copy / export carry columns the grid hides?', async () => {
    // The F6 branch point, and a thesis-#4 question rather than a plumbing one: a
    // scientist who compacts the grid to scan and then exports must not silently
    // lose %AUCextrap / %AUMCextrap / adj-R² / span ratio — the fragility columns
    // the rest of this plan works to keep on screen. The typings say `toCsv`'s
    // `visibleColumnsOnly` defaults to false and that passing a grid "takes visible
    // columns from the grid"; both halves are MEASURED here, because the same
    // typings assert `Column.get()` returns null for a missing cell and that is
    // false (guardrails: three read APIs, two null behaviours).
    const h = await embeddedResultsGrid();
    try {
      const probe = 'PCT_EXTRAP';
      const gc = h.grid.col(probe);
      expect(gc !== null, true, `${probe} resolves`);
      gc!.visible = false;

      const header = (csv: string): string[] =>
        (csv.split('\n')[0] ?? '').split(',').map((s) => s.trim().replace(/^"|"$/g, ''));

      const plain = header(h.results.toCsv());
      const withGrid = header(h.results.toCsv({}, h.grid));
      const visibleOnly = header(h.results.toCsv({visibleColumnsOnly: true}, h.grid));

      // The one hard assertion: the FRAME-level export (no grid) must be complete.
      // F6's decoupling depends on there being an export path that ignores grid
      // visibility — if even the plain export dropped hidden columns, hiding a
      // column would be data loss and F1 itself would need re-examination.
      expect(plain.includes(probe), true,
        `${probe} is missing from a frame-level toCsv() even though it was only ` +
        'hidden on the GRID — grid visibility is leaking into the frame export, ' +
        'which would make the compact profile lossy (F6 blocking).');

      gc!.visible = true;
      return `hidden=${probe} | plain toCsv carries it=${plain.includes(probe)} ` +
        `| toCsv(_, grid) carries it=${withGrid.includes(probe)} ` +
        `| toCsv({visibleColumnsOnly}, grid) carries it=${visibleOnly.includes(probe)} ` +
        `| cols: plain=${plain.length} grid=${withGrid.length} visibleOnly=${visibleOnly.length}`;
    } finally {
      h.close();
    }
  }, {timeout: TIMEOUT});

  test('AC-0.4: record the installed platform / PowerGrid versions against the known fix versions', async () => {
    // Known pin() fix versions from the Datagrok Community threads: v1.22.0 /
    // PowerGrid v1.4.6 (pin state not surviving reload), v1.23 / PowerGrid v1.5.0
    // (formatting loss), PowerGrid v1.4.11 (row-height duplication). Recorded so a
    // future regression can be attributed to a version move rather than re-litigated.
    let powerGrid = 'not-installed';
    try {
      const pkgs = await grok.dapi.packages.filter('PowerGrid').list();
      const pg = pkgs.find((p) => p.name === 'PowerGrid') ?? pkgs[0];
      if (pg !== undefined) powerGrid = `${pg.name}@${pg.version}`;
    } catch (e) {
      powerGrid = `lookup-failed: ${e instanceof Error ? e.message : String(e)}`;
    }
    // Any successful lookup is informational; the gate is the behavioural tests
    // above. Asserting a version number here would make the suite fail on an
    // upgrade that changed nothing observable.
    expect(typeof powerGrid, 'string', 'version probe returns a recordable string');
    return `PowerGrid=${powerGrid} | fix versions of record: 1.22.0/PG-1.4.6 ` +
      '(reload persistence), 1.23/PG-1.5.0 (formatting loss), PG-1.4.11 (row height)';
  }, {timeout: TIMEOUT});

  test('AC-1.2: hiding a column is GRID-side only — the frame keeps it, and it still sorts', async () => {
    // The profile hides by setting `gc.visible`, never by removing a column from
    // the frame. That distinction is what keeps a hidden-then-revealed column
    // sortable and a frame-level export complete; if the profile ever "simplified"
    // into dropping columns, both properties would fail silently and this catches it.
    const h = await embeddedResultsGrid();
    try {
      const probe = 'LAMBDAZ_INTERCEPT'; // hidden by the compact profile
      applyResultsColumnProfile(h.grid, false);
      const gc = h.grid.col(probe);
      expect(gc !== null, true, `${probe} still has a GRID column while hidden`);
      expect(gc!.visible, false, `${probe} is hidden by the compact profile`);
      expect(h.results.col(probe) !== null, true,
        `${probe} must remain in the FRAME while hidden — the profile is grid-side`);

      // Sorting BY a hidden column must work, because the data is still there.
      h.grid.sort([probe]);
      applyResultsColumnProfile(h.grid, true); // reveal
      expect(h.grid.col(probe)!.visible, true, `${probe} comes back on the Full profile`);
      expect(h.results.col(probe)!.length, h.results.rowCount,
        'the revealed column still carries a full column of values');
    } finally {
      h.close();
    }
  }, {timeout: TIMEOUT});

  test('a carried Smiles detected as Molecule raises the rows only while the structure is visible', async () => {
    // User-reported 2026-11-03: the platform re-detects a carried `Smiles` as `Molecule` and raises
    // every results row to the molecule height (109 px) although Compact hides the column.
    const h = await embeddedResultsGrid();
    try {
      const smiles = h.results.columns.addNewString('Smiles');
      smiles.init(() => 'COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1');
      smiles.semType = 'Molecule';
      h.grid.props.rowHeight = 109; // what the platform set on the reported stand
      applyResultsColumnProfile(h.grid, false);
      expect(h.grid.col('Smiles')!.visible, false, 'Compact hides the carried Smiles');
      expect(h.grid.props.rowHeight, 28, 'Compact: no visible structure → the default row height');
      applyResultsColumnProfile(h.grid, true);
      expect(h.grid.col('Smiles')!.visible, true, 'Full shows the carried Smiles');
      expect(h.grid.props.rowHeight, 109, 'Full: the structure is visible → the platform molecule height');
      applyResultsColumnProfile(h.grid, false);
      expect(h.grid.props.rowHeight, 28, 'back to Compact → default again');
    } finally {
      h.close();
    }
  }, {timeout: TIMEOUT});

  // The Compact control (on the results grid's banner strip) reads `Compact` and the
  // stored flag is `extended` (`extended === !compact`); this test pins the STORAGE,
  // which is what survives the swap, so it is written in the stored convention.
  test('AC-1.9: the profile choice rides on view.temp, so a re-Compute swap cannot reset it', async () => {
    // The defect this design avoids: `swapOrDockResults` replaces the DataFrame on
    // EVERY Compute, so the study grid's `df.temp` mechanism — correct there,
    // because that frame persists — would silently drop the scientist's choice
    // once per run. Both halves are asserted: the view-side tag survives the swap,
    // and a frame-side tag demonstrably does not.
    const h = await embeddedResultsGrid();
    try {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      const bag = h.view.temp as Record<string, any>;
      bag[RESULTS_EXTENDED_TAG] = true;
      h.results.temp['probeFrameTag'] = true;

      // Re-Compute: a NEW results frame is swapped into the SAME grid.
      const src2 = withDoseMg(await loadFixture('theoph.csv'));
      const {results: results2} = runPipelineOrThrow(
        src2, theophContract(src2), theophScalars(), DEFAULT_RULES);
      configureResultsDataFrame(results2);
      h.grid.dataFrame = results2;

      const extended = (bag[RESULTS_EXTENDED_TAG] as boolean | null) === true;
      expect(extended, true, 'the view-held profile choice survives the frame swap');
      expect(results2.temp['probeFrameTag'] === true, false,
        'a FRAME-held tag does NOT survive — which is exactly why this state lives ' +
        'on view.temp rather than df.temp');

      configureResultsGrid(h.grid, undefined, extended);
      const visible = h.grid.dataFrame.columns.names()
        .filter((n) => h.grid.col(n)?.visible === true);
      expect(visible.includes('LAMBDAZ_INTERCEPT'), true,
        'after the swap the grid is still on the Full profile the scientist chose');
    } finally {
      h.close();
    }
  }, {timeout: TIMEOUT});

  test('the compact profile RENDERS its columns through the full production sequence', async () => {
    // Written to chase a report of the shipped grid showing only Subject + Status,
    // and kept as a guard because it reproduces what the unit tests cannot: the
    // real `runCompute` order — configure the grid, pin the subject, mutate the
    // frame with carried metadata AFTER the grid exists, re-apply the profile —
    // on the real indometh (all-IV-bolus) frame, docked the way the app docks it.
    //
    // It asserts RENDERED WIDTHS, not `gc.visible`. An earlier revision measured
    // the flag, reported a clean 17, and proved nothing about what is drawn: the
    // two are independent, and only the width half is what the scientist sees.
    // INDOMETH, not theoph: the reported screenshot is the all-IV-bolus indometh
    // analysis, and route drives `routeSuppressed`. Reproduce the reported frame,
    // not a convenient one.
    const src = await loadFixture('indometh.csv');
    const {results} = runPipelineOrThrow(
      src, indomethContract(src), indomethScalars(), DEFAULT_RULES);
    configureResultsDataFrame(results);
    const view = grok.shell.addTableView(src);
    const grid = DG.Viewer.grid(results);
    configureResultsGrid(grid);
    const mainGridNode = view.dockManager.findNode(view.grid.root) ?? null;
    view.dockManager.dock(grid.root, DG.DOCK_TYPE.FILL, mainGridNode, 'NCA results');
    const h = {
      view, grid, results,
      close: () => {
        try {
          grid.close();
        } catch {/* best-effort */}
        try {
          view.close();
        } catch {/* best-effort */}
      },
    };
    try {
      const order = resultsGridVisibleOrder(h.results, false);
      const afterConfigure = h.results.columns.names()
        .filter((n) => h.grid.col(n)?.visible === true);

      // `appendGroupingColumns` adds carried metadata to the frame IN PLACE, after
      // the grid was built and configured.
      h.results.columns.add(
        DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Species',
          new Array<string>(h.results.rowCount).fill('Rat')));
      h.results.columns.add(
        DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Strain',
          new Array<string>(h.results.rowCount).fill('Sprague-Dawley')));
      const afterAppend = h.results.columns.names()
        .filter((n) => h.grid.col(n)?.visible === true);

      applyResultsColumnProfile(h.grid, false);
      const afterReapply = h.results.columns.names()
        .filter((n) => h.grid.col(n)?.visible === true);

      // `grid.col(n)` resolving to null is the other candidate: a column added to
      // the frame after grid construction may have no GridColumn yet.
      const nullGridCols = h.results.columns.names()
        .filter((n) => h.grid.col(n) === null);

      // RENDER, not flag. The earlier revision of this probe measured
      // `gc.visible === true` and reported a clean 17 while the shipped grid drew
      // two columns — the two are independent, and only this half matches what the
      // scientist sees.
      await new Promise((r) => window.setTimeout(r, 400));
      const widths = afterReapply.map((n) => `${n}:${h.grid.col(n)?.width ?? 'na'}`);
      const zeroWidth = afterReapply.filter((n) => (h.grid.col(n)?.width ?? 0) <= 0);
      const totalWidth = afterReapply.reduce((s, n) => s + (h.grid.col(n)?.width ?? 0), 0);

      expect(nullGridCols.length, 0,
        `every frame column resolves to a GridColumn; unresolved: ${nullGridCols.join(',')}`);
      expect(afterReapply.length, order.length,
        `the grid shows exactly the compact partition (${order.length} columns), not ` +
        `${afterReapply.length}. Subject+Status alone is the reported defect shape.`);
      expect(afterReapply.length > 2, true,
        'more than Subject + Status survive the sequence');
      expect(zeroWidth.length, 0,
        `no visible column renders at zero width; collapsed: ${zeroWidth.join(',')}`);
      expect(afterReapply[0], 'SUBJECT', 'identity leads the rendered grid');
      // The two intermediate readings are asserted too, so a regression is located
      // at the step that caused it rather than only at the end.
      expect(afterConfigure.length, order.length,
        `configureResultsGrid alone already yields the compact partition; got ` +
        `${afterConfigure.length}`);
      expect(afterAppend.length, order.length + 2,
        `columns added to the frame after the grid exists default to VISIBLE ` +
        `(${order.length} + 2 carried), which is exactly why runCompute re-applies ` +
        `the profile afterwards; got ${afterAppend.length}`);

      // The repaint nudge changes a real width and reverts it on the next
      // macrotask. Re-applying the profile is something `runCompute` does more
      // than once per run, so the revert must leave NO drift — otherwise the
      // Status column creeps a pixel wider on every Compute.
      const statusW = grid.col('STATUS')?.width ?? -1;
      applyResultsColumnProfile(h.grid, false);
      applyResultsColumnProfile(h.grid, false);
      applyResultsColumnProfile(h.grid, false);
      await new Promise((r) => window.setTimeout(r, 200));
      expect(grid.col('STATUS')?.width ?? -1, statusW,
        'three more profile applications leave the nudged column exactly where it ' +
        'was — the deferred revert must not drift');

      return `rendered=${afterReapply.length} afterConfigure=${afterConfigure.length} ` +
        `afterAppend=${afterAppend.length} statusW=${statusW} ` +
        `rootW=${Math.round(h.grid.root.getBoundingClientRect().width)} ` +
        `totalColW=${totalWidth} widths=[${widths.join(' ')}]`;
    } finally {
      h.close();
    }
  }, {timeout: TIMEOUT});
});

category('Results grid: λz reason pills, live (acceptance-reason-pills)', () => {
  test('the Violated column resolves to the platform Tags renderer, one fill per criterion', async () => {
    const h = await embeddedResultsGrid();
    try {
      const col = h.results.col(ACCEPT_REASONS_COL)!;
      let row = -1;
      for (let i = 0; i < col.length && row < 0; i++) if (col.getString(i) !== '') row = i;
      expect(row >= 0, true, 'theoph carries at least one Flagged profile');
      // The `Tags` renderer lives in PowerGrid, which the platform loads on first use: the very
      // first read can still name the default `string` renderer. Measure both reads.
      const first = h.grid.cell(ACCEPT_REASONS_COL, row).renderer?.name;
      const tStart = Date.now();
      let renderer = h.grid.cell(ACCEPT_REASONS_COL, row).renderer;
      while (renderer?.name !== 'Tags' && Date.now() - tStart < 5000) {
        await DG.delay(100);
        renderer = h.grid.cell(ACCEPT_REASONS_COL, row).renderer;
      }
      const loadMs = Date.now() - tStart;
      expect(renderer?.name, 'Tags', `the cell renders through Tags (PowerGrid), got ${renderer?.name}`);
      const colours = h.grid.col(ACCEPT_REASONS_COL)!.categoryColors;
      const labels = col.getString(row).split(',');
      for (const l of labels)
        expect(colours[l] !== undefined, true, `pill '${l}' has a fill`);
      // acceptance-reason-pill-colors: the four criterion fills of the user's palette (AC-1), each
      // keyed to its live label, and the renderer's contrast text resolves dark on every one (AC-3).
      const fills = Object.values(colours).map((c) => typeof c === 'string' ? DG.Color.fromHtml(c) : c);
      expect(new Set(fills.map((c) => DG.Color.toHtml(c).toUpperCase())).size, 4, 'four distinct fills');
      for (const hex of ['#BFE4F7', '#F5DB98', '#EDC6DD', '#F3C6A6'])
        expect(fills.some((c) => DG.Color.toHtml(c).toUpperCase() === hex), true, `fill ${hex} present`);
      for (const c of fills) {
        const text = DG.Color.getContrastColor(c);
        const lum = 0.299 * DG.Color.r(text) + 0.587 * DG.Color.g(text) + 0.114 * DG.Color.b(text);
        expect(lum < 128, true, `pill text on ${DG.Color.toHtml(c)} is dark`);
      }
      return `row ${row}: '${col.getString(row)}' → first=${first}, ${renderer?.name} after ${loadMs} ms; ` +
        `${Object.keys(colours).length} keyed labels: ` +
        Object.entries(colours).map(([l, c]) => `${l}=${DG.Color.toHtml(c as number)}`).join(', ');
    } finally {
      h.close();
    }
  }, {timeout: TIMEOUT});

  test('real paint: does the cell outside the pills keep its background? (PowerGrid GROK-20888 fix)', async () => {
    // PowerGrid ≤ 1.8.4 filled each pill without `beginPath()`, so the grid's own cell rectangle
    // rode along in the path and the WHOLE cell took the last pill's colour, hiding every pill but
    // the last (measured 2026-10-19 on acceptance_check; fixed upstream in GROK-20888, v.next).
    // This reads the painted pixels, so it fails on a server that still has the leaking build.
    const h = await embeddedResultsGrid();
    try {
      const col = h.results.col(ACCEPT_REASONS_COL)!;
      let row = -1;
      for (let i = 0; i < col.length && row < 0; i++) if (col.getString(i) !== '') row = i;
      expect(row >= 0, true, 'theoph carries at least one Flagged profile');
      h.grid.col(ACCEPT_REASONS_COL)!.width = 400; // room to spare right of the last pill
      const t0 = Date.now();
      while (h.grid.cell(ACCEPT_REASONS_COL, row).renderer?.name !== 'Tags' && Date.now() - t0 < 5000)
        await DG.delay(100);
      h.grid.invalidate();
      await DG.delay(500);
      const b = h.grid.cell(ACCEPT_REASONS_COL, row).bounds;
      const dpr = window.devicePixelRatio || 1;
      const px = h.grid.canvas.getContext('2d')!
        .getImageData(Math.round((b.right - 6) * dpr), Math.round((b.top + b.height / 2) * dpr), 1, 1).data;
      const pill = h.grid.col(ACCEPT_REASONS_COL)!.categoryColors[col.getString(row).split(',').pop()!];
      const isPill = px[0] === DG.Color.r(pill) && px[1] === DG.Color.g(pill) && px[2] === DG.Color.b(pill);
      // REPORT-ONLY until a released PowerGrid carries the fix: the
      // released 1.8.4 still leaks, so an assertion here would be red for a reason outside this
      // package. Flip to `expect(isPill, false, …)` once that release is on the target servers.
      return isPill ?
        `LEAKING PowerGrid build: row ${row} right-edge pixel is the pill colour — pills unreadable on this server` :
        `row ${row} '${col.getString(row)}': right-edge pixel rgb(${px[0]},${px[1]},${px[2]}) — pills paint cleanly`;
    } finally {
      h.close();
    }
  }, {timeout: TIMEOUT});

  test('app-start PowerGrid check reads the loaded version (the Violated default)', async () => {
    const s = detectPowerGridSupport();
    expect(s.version !== null, true, 'PowerGrid is installed on the test server — its version must be found');
    expect(s.pills, powerGridDrawsPills(s.version!), 'pills default follows the version gate');
    return `PowerGrid ${s.version} → Violated column ${s.pills ? 'ON' : 'OFF'} by default`;
  }, {timeout: TIMEOUT});

  test('rule 16: reasons per Flagged profile on the fixtures (measured)', async () => {
    const theoph = withDoseMg(await loadFixture('theoph.csv'));
    const indo = await loadFixture('indometh.csv');
    const rat = await loadFixture('03_rat_multiroute.csv');
    const ratOpts = setupMultirouteRoute(rat);
    const frames: Array<[string, DG.DataFrame]> = [
      ['theoph', runPipelineOrThrow(theoph, theophContract(theoph), theophScalars(), DEFAULT_RULES).results],
      ['indometh', runPipelineOrThrow(indo, indomethContract(indo), indomethScalars(), DEFAULT_RULES).results],
      ['rat_multiroute', runPipelineOrThrow(
        rat, ratMultirouteContract(rat), ratScalars(nca.ROUTE_PO), DEFAULT_RULES, ratOpts).results],
    ];
    const out: string[] = [];
    for (const [name, df] of frames) {
      configureResultsDataFrame(df);
      const col = df.col(ACCEPT_REASONS_COL)!;
      const hist = [0, 0, 0, 0, 0];
      const labels = new Map<string, number>();
      for (let i = 0; i < col.length; i++) {
        const v = col.getString(i);
        const parts = v === '' ? [] : v.split(',');
        hist[parts.length]++;
        for (const p of parts) labels.set(p, (labels.get(p) ?? 0) + 1);
      }
      expect(hist.reduce((a, b) => a + b, 0), df.rowCount, `${name}: every row counted`);
      out.push(`${name} n=${df.rowCount} reasons/row 0..4=[${hist.join(',')}] ` +
        `{${[...labels].map(([k, v]) => `${k}:${v}`).join('; ')}}`);
    }
    return out.join(' | ');
  }, {timeout: TIMEOUT});
});
