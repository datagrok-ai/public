import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {category, test, expect, before, awaitCheck, delay} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';
import {createTableView} from './utils';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {ChemSimilarityViewer} from '../analysis/chem-similarity-viewer';
import {ChemDiversityViewer} from '../analysis/chem-diversity-viewer';
import {Subscription} from 'rxjs';

// Tests for the similarity viewer's select-all-similar / collaborative filter controls, plus edge cases.
category('similarity filter/select', () => {
  before(async () => {
    grok.shell.closeAll();
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  /** The filter-panel label this viewer pushes (cutoff rounded, in sync with _filterSummary). */
  function simLabel(v: ChemSimilarityViewer): string {
    return `Most similar structures: ${v.distanceMetric} ≥ ${Number(v.cutoff.toFixed(3))}`;
  }

  /** True when `bs` is non-empty and every set bit lies within the selection [0, half). */
  function confinedToSelection(bs: DG.BitSet | null, half: number, rowCount: number): boolean {
    if (!bs || bs.trueCount === 0)
      return false;
    for (let i = half; i < rowCount; i++) {
      if (bs.get(i))
        return false;
    }
    return true;
  }

  /** Opens the similarity viewer on sar-small (default cutoff 0.01 → a non-empty similar set). */
  async function openViewer(): Promise<{tv: DG.TableView, viewer: ChemSimilarityViewer, df: DG.DataFrame}> {
    const tv = await createTableView('tests/sar-small_test.csv');
    const viewer = (await tv.dataFrame.plot.fromType('Chem Similarity Search')) as ChemSimilarityViewer;
    // Gate on the computed set (not the hot renderCompleted Subject) to avoid a fast-machine flake.
    await awaitCheck(() => (viewer.similarSetBitset?.trueCount ?? 0) > 0, 'similar set was not computed', 10000);
    return {tv, viewer, df: tv.dataFrame};
  }

  test('filter toggle on then off restores the table', async () => {
    // Toggling the filter on then off must filter to the similar set, then fully restore the table.
    const {tv, viewer, df} = await openViewer();
    try {
      const similarCount = viewer.similarSetBitset!.trueCount;
      const label = simLabel(viewer);
      const gridBefore = viewer.root.querySelector('.chem-viewer-grid'); // toggling must not rebuild the panel
      viewer.filterBtn!.click();
      await awaitCheck(() => df.filter.trueCount === similarCount, 'table was not filtered to the similar set', 5000);
      expect(viewer.filterActive, true);
      expect(viewer.filterBtn!.getAttribute('aria-pressed'), 'true');
      expect(df.rows.filters.includes(label), true);
      viewer.filterBtn!.click();
      await awaitCheck(() => df.filter.trueCount === df.rowCount, 'table was not restored after toggling off', 5000);
      expect(viewer.filterActive, false);
      expect(viewer.filterBtn!.getAttribute('aria-pressed'), 'false');
      expect(df.rows.filters.includes(label), false);
      await delay(300); // let any debounced re-render settle
      // Cards must not have been cleared + rebuilt by either toggle.
      expect(viewer.root.querySelector('.chem-viewer-grid') === gridBefore, true);
    } finally {tv.close();}
  });

  test('empty set disables buttons and guards', async () => {
    // An empty similar set must disable the buttons and guard select / keyboard activation.
    const {tv, viewer, df} = await openViewer();
    try {
      // cutoff 1 → no exact duplicate of row 0 → empty similar set.
      viewer.setOptions({cutoff: 1});
      await awaitCheck(() => (viewer.similarSetBitset?.trueCount ?? 0) === 0, 'similar set did not become empty', 5000);
      expect(viewer.selectBtn!.classList.contains('chem-similarity-action-disabled'), true);
      expect(viewer.filterBtn!.classList.contains('chem-similarity-action-disabled'), true);

      // Select must not wipe an existing selection when there is nothing similar.
      df.selection.set(0, true, false);
      df.selection.set(1, true);
      const selBefore = df.selection.trueCount;
      viewer.selectBtn!.click();
      expect(df.selection.trueCount, selBefore);

      // Keyboard activation of the disabled filter button must not toggle it on with an empty set.
      viewer.filterBtn!.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', bubbles: true}));
      expect(viewer.filterActive, false);
      expect(df.filter.trueCount, df.rowCount);
    } finally {tv.close();}
  });

  test('empty set while filtered: stays usable and clears the label', async () => {
    // Emptying the set while filtered must un-filter, drop the label, and keep the toggle turn-off-able.
    const {tv, viewer, df} = await openViewer();
    try {
      viewer.filterBtn!.click();
      const label = simLabel(viewer); // captured at the original cutoff, before it empties the set
      await awaitCheck(() => viewer.filterActive && df.rows.filters.includes(label), 'filter did not apply', 5000);
      viewer.setOptions({cutoff: 1}); // empties the similar set while the toggle is still on
      await awaitCheck(() => (viewer.similarSetBitset?.trueCount ?? 0) === 0, 'set did not empty', 5000);
      // Table un-filters and the now-meaningless label is dropped.
      await awaitCheck(() => df.filter.trueCount === df.rowCount && !df.rows.filters.includes(label),
        'stale label/filter lingered after the set emptied', 5000);
      // Filter button must stay clickable (disabled class absent) so a real mouse can turn it off; select disabled.
      expect(viewer.filterBtn!.classList.contains('chem-similarity-action-disabled'), false);
      expect(viewer.selectBtn!.classList.contains('chem-similarity-action-disabled'), true);
      // Deactivation must still work on an empty set (the guard blocks only activation).
      viewer.filterBtn!.click();
      await awaitCheck(() => !viewer.filterActive && df.filter.trueCount === df.rowCount,
        'filter could not be turned off when the set is empty', 5000);
    } finally {tv.close();}
  });

  test('detach restores table', async () => {
    // Detaching the viewer while filtered must restore the table and drop the label.
    const {tv, viewer, df} = await openViewer();
    try {
      const label = simLabel(viewer);
      viewer.filterBtn!.click();
      await awaitCheck(() => df.filter.trueCount < df.rowCount, 'table was not filtered', 5000);
      viewer.detach();
      await awaitCheck(() => df.filter.trueCount === df.rowCount, 'detach did not restore the table', 5000);
      expect(df.rows.filters.includes(label), false);
    } finally {tv.close();}
  });

  // NOTE: "Reset all filters" clearing the similarity filter (via grok.events.onResetFilterRequest) is not
  // covered — there is no JS-callable way to fire that platform event from a test.

  test('enable filter when another filter is already active', async () => {
    // Turning on the similarity filter when an external filter is already active must AND with it.
    const {tv, viewer, df} = await openViewer();
    let sub: Subscription | null = null;
    try {
      const half = Math.floor(df.rowCount / 2);
      sub = df.onRowsFiltering.subscribe(() => {
        for (let i = half; i < df.rowCount; i++)
          df.filter.set(i, false, false);
      });
      df.rows.requestFilter();
      const simSet = viewer.similarSetBitset!;
      let expected = 0;
      for (let i = 0; i < df.rowCount; i++) {
        if (simSet.get(i) && i < half)
          expected++;
      }
      viewer.filterBtn!.click();
      await awaitCheck(() => viewer.filterActive && df.filter.trueCount === expected,
        'similarity filter did not AND with the pre-existing filter', 5000);
    } finally {
      sub?.unsubscribe();
      tv.close();
    }
  });

  test('diversity viewer keeps its limit cap', async () => {
    // setOptions bypasses the UI slider max; the O(n^2) diversity limit must stay clamped to 50.
    const tv = await createTableView('tests/sar-small_test.csv');
    try {
      const viewer = (await tv.dataFrame.plot.fromType('Chem Diversity Search')) as ChemDiversityViewer;
      await awaitCheck(() => viewer.renderMolIds.length > 0, 'diversity viewer did not render', 10000);
      viewer.setOptions({limit: 10000});
      expect(viewer.limit, 50); // clamped synchronously in onPropertyChanged
    } finally {tv.close();}
  });

  test('filter is collaborative', async () => {
    // A second filter applied after the similarity filter must AND with it (similar set ∩ upper half).
    const {tv, viewer, df} = await openViewer();
    let sub: Subscription | null = null;
    try {
      const simSet = viewer.similarSetBitset!;
      viewer.filterBtn!.click();
      await awaitCheck(() => viewer.filterActive && df.filter.trueCount === simSet.trueCount,
        'similarity filter was not applied', 5000);

      const half = Math.floor(df.rowCount / 2);
      let expected = 0;
      for (let i = 0; i < df.rowCount; i++) {
        if (simSet.get(i) && i < half)
          expected++;
      }
      sub = df.onRowsFiltering.subscribe(() => {
        for (let i = half; i < df.rowCount; i++)
          df.filter.set(i, false, false);
      });
      df.rows.requestFilter();
      expect(df.filter.trueCount, expected);
    } finally {
      sub?.unsubscribe();
      tv.close();
    }
  });

  test('clearing an external filter re-expands the similar set', async () => {
    // Hiding similar rows shrinks the set; clearing that filter must re-expand it back to the full set.
    const {tv, viewer, df} = await openViewer();
    let sub: Subscription | null = null;
    try {
      viewer.filterBtn!.click();
      await awaitCheck(() => viewer.filterActive, 'filter did not activate', 5000);
      // Snapshot the full similar set (no other filter active yet).
      const simIdx: number[] = [];
      for (let i = 0; i < df.rowCount; i++) {
        if (viewer.similarSetBitset!.get(i))
          simIdx.push(i);
      }
      const fullCount = simIdx.length;
      expect(fullCount > 1, true); // need at least 2 similar rows to observe a shrink

      // Hide the latter half of the similar rows → the set re-searches within the rest and shrinks.
      const hidden = new Set(simIdx.slice(Math.floor(fullCount / 2)));
      sub = df.onRowsFiltering.subscribe(() => {
        for (const i of hidden)
          df.filter.set(i, false, false);
      });
      df.rows.requestFilter();
      await awaitCheck(() => (viewer.similarSetBitset?.trueCount ?? 0) < fullCount,
        'similar set did not shrink under the external filter', 5000);

      // Clear the external filter. The net df.filter is unchanged, so the viewer must re-expand via
      // its onRowsFiltering handler (onFilterChanged stays silent).
      sub.unsubscribe();
      sub = null;
      df.rows.requestFilter();
      await awaitCheck(() => viewer.similarSetBitset!.trueCount === fullCount,
        'similar set did not re-expand after clearing the external filter', 5000);
      await awaitCheck(() => df.filter.trueCount === fullCount,
        'table did not widen back to the full similar set', 5000);
    } finally {
      sub?.unsubscribe();
      tv.close();
    }
  });

  test('select acts on the full similar set, not just the displayed cards', async () => {
    // Cap displayed cards below the full set; select must act on every similar row, not just the shown 3.
    const {tv, viewer, df} = await openViewer();
    try {
      viewer.setOptions({limit: 3});
      await awaitCheck(() => (viewer.molCol?.length ?? 99) <= 3 && (viewer.similarSetBitset?.trueCount ?? 0) > 3,
        'display did not cap below the full similar set', 5000);
      const total = viewer.similarSetBitset!.trueCount;
      viewer.selectBtn!.click();
      expect(df.selection.trueCount, total); // acts on all `total`, not the 3 shown cards
    } finally {tv.close();}
  });

  test('filter recheck recomputes even when followCurrentRow is off', async () => {
    // Regression: with followCurrentRow off, an external filter change must still rebuild the set (not AND a stale mask).
    const {tv, viewer, df} = await openViewer();
    let sub: Subscription | null = null;
    try {
      viewer.setOptions({followCurrentRow: false});
      viewer.filterBtn!.click();
      await awaitCheck(() => viewer.filterActive && df.filter.trueCount > 0, 'filter did not activate', 5000);
      await awaitCheck(() => !viewer.isComputing, 'initial search did not settle', 8000);
      const before = viewer.similarSetBitset!.trueCount;
      // Narrow the population with an external filter (keep only the first third of rows).
      const keep = Math.max(1, Math.floor(df.rowCount / 3));
      sub = df.onRowsFiltering.subscribe(() => {
        for (let i = keep; i < df.rowCount; i++)
          df.filter.set(i, false, false);
      });
      df.rows.requestFilter();
      // The set must shrink to the kept rows, proving the recheck recomputed.
      await awaitCheck(() => (viewer.similarSetBitset?.trueCount ?? before) < before,
        'similar set was not rebuilt against the narrowed population (followCurrentRow=off recompute gap)', 8000);
      // The reference must stay pinned (followCurrentRow off must not jump it to the current row).
      expect(viewer.targetMoleculeIdx, 0);
    } finally {
      sub?.unsubscribe();
      tv.close();
    }
  });

  test('an in-progress substructure search gates the filter recheck', async () => {
    // Anti-storm gate: while the watchdog _externalSearchTimer is live, _scheduleFilterRecheck must schedule
    // nothing; once cleared it schedules normally. Asserted synchronously (the event wiring is Puppeteer-flaky).
    const {tv, viewer} = await openViewer();
    const g = viewer as unknown as {
      _externalSearchTimer: ReturnType<typeof setTimeout> | null;
      _filterRecheckTimer: ReturnType<typeof setTimeout> | null;
      _scheduleFilterRecheck(): void;
    };
    try {
      if (g._filterRecheckTimer) {
        clearTimeout(g._filterRecheckTimer);
        g._filterRecheckTimer = null;
      }
      if (g._externalSearchTimer) {
        clearTimeout(g._externalSearchTimer);
        g._externalSearchTimer = null;
      }
      // Watchdog live → recheck suppressed.
      g._externalSearchTimer = setTimeout(() => {}, 30000);
      g._scheduleFilterRecheck();
      expect(g._filterRecheckTimer, null);
      // Watchdog cleared → recheck schedules normally.
      clearTimeout(g._externalSearchTimer);
      g._externalSearchTimer = null;
      g._scheduleFilterRecheck();
      expect(g._filterRecheckTimer !== null, true);
    } finally {
      if (g._filterRecheckTimer)
        clearTimeout(g._filterRecheckTimer);
      if (g._externalSearchTimer)
        clearTimeout(g._externalSearchTimer);
      tv.close();
    }
  });

  test('filter works in rowSource = Selected and FilteredSelected', async () => {
    // Both modes must confine the similar set and the applied filter to the selected rows.
    for (const mode of ['Selected', 'FilteredSelected']) {
      const {tv, viewer, df} = await openViewer();
      try {
        df.selection.setAll(false, false);
        const half = Math.max(1, Math.floor(df.rowCount / 2));
        for (let i = 0; i < half; i++)
          df.selection.set(i, true, false);
        viewer.setOptions({rowSource: mode});
        // Poll until the recompute confines the set (don't race the rowSource-change render).
        await awaitCheck(() => confinedToSelection(viewer.similarSetBitset, half, df.rowCount),
          `${mode} similar set did not confine to the selection`, 12000);
        viewer.filterBtn!.click();
        await awaitCheck(() => viewer.filterActive && confinedToSelection(df.filter, half, df.rowCount),
          `filter did not confine to the selection in ${mode} mode`, 5000);
      } finally {tv.close();}
    }
  });

  test('selecting refreshes card highlights in place (no panel rebuild)', async () => {
    // A selection change must restyle the existing cards in place, not clear + rebuild the panel (flicker).
    const {tv, viewer, df} = await openViewer();
    try {
      const gridBefore = viewer.root.querySelector('.chem-viewer-grid'); // the live card container
      expect(gridBefore !== null, true);
      viewer.selectBtn!.click(); // selects the similar set → the cards must show as selected
      await awaitCheck(() => df.selection.trueCount > 0, 'selection was not made', 5000);
      await delay(300); // let the debounced selection.onChanged render run
      // Same DOM node ⇒ refreshed in place, not rebuilt.
      expect(viewer.root.querySelector('.chem-viewer-grid') === gridBefore, true);
      // ...and a card picked up the new selection highlight.
      expect(viewer.root.querySelector('.chem-viewer-grid .d4-selected') !== null, true);
    } finally {tv.close();}
  });
});
