import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {category, test, expect, before, awaitCheck, delay} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';
import {createTableView} from './utils';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {ChemSimilarityViewer} from '../analysis/chem-similarity-viewer';
import {ChemDiversityViewer} from '../analysis/chem-diversity-viewer';
import {Subscription} from 'rxjs';

// Tests for the similarity viewer's "select all similar" and collaborative "filter to similar set"
// controls (the two action icons next to the metrics link), plus the empty-set and detach edge cases.
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

  /** True when `bs` is non-empty and every set bit lies within the selection [0, half) — i.e. the
   * Selected / FilteredSelected modes confined the set to the selected rows. */
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
    // Gate on the computed set, not the hot renderCompleted Subject (a fast machine can fire it before
    // the test subscribes, causing a timeout) — awaitCheck polls and avoids that flake.
    await awaitCheck(() => (viewer.similarSetBitset?.trueCount ?? 0) > 0, 'similar set was not computed', 10000);
    return {tv, viewer, df: tv.dataFrame};
  }

  test('filter toggle on then off restores the table', async () => {
    // (Basic "select all similar" is covered by the stronger "select acts on the full similar set" test.)
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
      // The cards (same similar set throughout) must not have been cleared + rebuilt by either toggle.
      expect(viewer.root.querySelector('.chem-viewer-grid') === gridBefore, true);
    } finally {tv.close();}
  });

  test('empty set disables buttons and guards', async () => {
    const {tv, viewer, df} = await openViewer();
    try {
      // cutoff 1 with a reference that has no exact duplicate → the similar set is empty.
      // (Relies on sar-small_test.csv row 0 having no Tanimoto-1.0 twin — keep that data file as-is.)
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

      // Keyboard activation of the (disabled) filter button must not toggle it on with an empty set.
      viewer.filterBtn!.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', bubbles: true}));
      expect(viewer.filterActive, false);
      expect(df.filter.trueCount, df.rowCount);
    } finally {tv.close();}
  });

  test('empty set while filtered: stays usable and clears the label', async () => {
    const {tv, viewer, df} = await openViewer();
    try {
      viewer.filterBtn!.click();
      const label = simLabel(viewer); // captured at the original cutoff, before it empties the set
      await awaitCheck(() => viewer.filterActive && df.rows.filters.includes(label), 'filter did not apply', 5000);
      viewer.setOptions({cutoff: 1}); // empties the similar set while the toggle is still on
      await awaitCheck(() => (viewer.similarSetBitset?.trueCount ?? 0) === 0, 'set did not empty', 5000);
      // Table un-filters (handler skips the AND) and the now-meaningless label is dropped.
      await awaitCheck(() => df.filter.trueCount === df.rowCount && !df.rows.filters.includes(label),
        'stale label/filter lingered after the set emptied', 5000);
      // The filter button must stay clickable (not pointer-events:none) so a REAL mouse can still turn it off
      // (a synthetic .click() ignores pointer-events, so assert the disabled class is absent); select disabled.
      expect(viewer.filterBtn!.classList.contains('chem-similarity-action-disabled'), false);
      expect(viewer.selectBtn!.classList.contains('chem-similarity-action-disabled'), true);
      // Deactivation must still work on an empty set (the guard blocks only activation).
      viewer.filterBtn!.click();
      await awaitCheck(() => !viewer.filterActive && df.filter.trueCount === df.rowCount,
        'filter could not be turned off when the set is empty', 5000);
    } finally {tv.close();}
  });

  test('detach restores table', async () => {
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

  // NOTE: "Reset all filters" clearing the similarity filter is implemented via the same
  // grok.events.onResetFilterRequest subscription used by chem-substructure-filter / scaffold-tree
  // filters. It isn't covered by an automated test because there is no JS-callable way to FIRE that
  // platform event from a test — df.resetFilter() exists but resets the filter without raising the event.

  test('stale filter label is stripped on attach', async () => {
    const tv = await createTableView('tests/sar-small_test.csv');
    try {
      const df = tv.dataFrame;
      // A label left behind by a restored project (toggle state itself is session-only) must not linger.
      df.rows.filters.push('Most similar structures: Tanimoto ≥ 0.5');
      const viewer = (await df.plot.fromType('Chem Similarity Search')) as ChemSimilarityViewer;
      await awaitCheck(() => (viewer.similarSetBitset?.trueCount ?? 0) > 0, 'viewer did not compute', 10000);
      expect(df.rows.filters.includes('Most similar structures: Tanimoto ≥ 0.5'), false);
    } finally {tv.close();}
  });

  test('enable filter when another filter is already active', async () => {
    const {tv, viewer, df} = await openViewer();
    let sub: Subscription | null = null;
    try {
      // An external filter is active BEFORE the similarity filter is turned on (reverse order —
      // the case the other-filters probe exists for).
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
    const tv = await createTableView('tests/sar-small_test.csv');
    try {
      const viewer = (await tv.dataFrame.plot.fromType('Chem Diversity Search')) as ChemDiversityViewer;
      await awaitCheck(() => viewer.renderMolIds.length > 0, 'diversity viewer did not render', 10000);
      // setOptions bypasses the UI slider's max; the O(n^2) diversity selection must stay bounded by 50
      // (whereas the similarity viewer raises its own cap to 10000).
      viewer.setOptions({limit: 10000});
      expect(viewer.limit, 50); // clamped synchronously in onPropertyChanged
    } finally {tv.close();}
  });

  test('filter is collaborative', async () => {
    const {tv, viewer, df} = await openViewer();
    let sub: Subscription | null = null;
    try {
      const simSet = viewer.similarSetBitset!;
      viewer.filterBtn!.click();
      await awaitCheck(() => viewer.filterActive && df.filter.trueCount === simSet.trueCount,
        'similarity filter was not applied', 5000);

      // A second filter that hides the table's lower half should AND with the similarity mask:
      // the result is (similar set ∩ upper half), proving the viewer adapts to other filters.
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

      // Hide the latter half of the similar rows with an external filter → the set is re-searched
      // within the remaining rows and shrinks; our mask is now a subset of that filter.
      const hidden = new Set(simIdx.slice(Math.floor(fullCount / 2)));
      sub = df.onRowsFiltering.subscribe(() => {
        for (const i of hidden)
          df.filter.set(i, false, false);
      });
      df.rows.requestFilter();
      await awaitCheck(() => (viewer.similarSetBitset?.trueCount ?? 0) < fullCount,
        'similar set did not shrink under the external filter', 5000);

      // Clear the external filter. Because our mask is a subset of it, the NET df.filter is unchanged,
      // so onFilterChanged stays silent — the viewer must re-expand via its onRowsFiltering handler.
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
    const {tv, viewer, df} = await openViewer();
    try {
      // Cap the DISPLAYED cards well below the full above-cutoff set (the core "see/act on more than
      // what's shown" promise): the search df (molCol) is limited, similarSetBitset is not.
      viewer.setOptions({limit: 3});
      await awaitCheck(() => (viewer.molCol?.length ?? 99) <= 3 && (viewer.similarSetBitset?.trueCount ?? 0) > 3,
        'display did not cap below the full similar set', 5000);
      const total = viewer.similarSetBitset!.trueCount;
      viewer.selectBtn!.click();
      expect(df.selection.trueCount, total); // acts on all `total`, not the 3 shown cards
    } finally {tv.close();}
  });

  test('rows added while filtering rebuild the mask', async () => {
    const {tv, viewer, df} = await openViewer();
    try {
      viewer.filterBtn!.click();
      await awaitCheck(() => viewer.filterActive && df.filter.trueCount < df.rowCount, 'filter was not applied', 5000);
      const beforeLen = viewer.similarSetBitset!.length;
      // Add a row → rowCount grows; the mask sized to the old count must be rebuilt (otherwise the
      // length-mismatch guard would silently turn the filter into a no-op).
      df.rows.addNew([df.col('smiles')!.get(0)]);
      await awaitCheck(() => viewer.similarSetBitset!.length === df.rowCount && df.rowCount > beforeLen,
        'mask was not rebuilt to the new row count', 8000);
      // The rebuilt mask must actually be APPLIED, not just rebuilt — table == the similar set.
      await awaitCheck(() => df.filter.trueCount === viewer.similarSetBitset!.trueCount,
        'rebuilt mask was not applied to the table', 5000);
      expect(viewer.filterActive, true);
    } finally {tv.close();}
  });

  test('filter works in rowSource = All', async () => {
    const {tv, viewer, df} = await openViewer();
    try {
      // In All mode the similar set is filter-independent, but the toggle must still filter the table to it.
      viewer.setOptions({rowSource: 'All'});
      await awaitCheck(() => (viewer.similarSetBitset?.trueCount ?? 0) > 0, 'set not computed in All mode', 5000);
      const simCount = viewer.similarSetBitset!.trueCount;
      viewer.filterBtn!.click();
      await awaitCheck(() => viewer.filterActive && df.filter.trueCount === simCount,
        'filter did not apply in All mode', 5000);
      viewer.filterBtn!.click();
      await awaitCheck(() => !viewer.filterActive && df.filter.trueCount === df.rowCount,
        'filter did not clear in All mode', 5000);
    } finally {tv.close();}
  });

  test('reference row change does not progressively shrink the filtered set', async () => {
    // The other-filters probe exists so moving the reference row re-searches against the SAME population
    // (the full table here — no other filter), not the already-narrowed df.filter. Without it, each move
    // would shrink the set further; returning to a row would then give a smaller count than the first visit.
    const {tv, viewer, df} = await openViewer();
    try {
      viewer.filterBtn!.click();
      await awaitCheck(() => viewer.filterActive && df.filter.trueCount > 0, 'filter did not activate', 5000);
      const countForRow = async (row: number): Promise<number> => {
        df.currentRowIdx = row;
        await awaitCheck(() => viewer.targetMoleculeIdx === row && !viewer.isComputing,
          `reference did not settle on row ${row}`, 8000);
        await delay(400); // let the debounced filter contribution apply
        return df.filter.trueCount;
      };
      const first = await countForRow(0);
      await countForRow(1);
      await countForRow(2);
      const back = await countForRow(0);
      // Same reference, same population → identical count. A shrinking probe would make `back` < `first`.
      expect(back, first);
    } finally {tv.close();}
  });

  test('filter recheck recomputes even when followCurrentRow is off', async () => {
    // Regression: with followCurrentRow off + Filtered rowSource + filter on, an external filter change must
    // still rebuild the similar set against the new population. The recompute gate used to skip it, leaving
    // the contribution handler ANDing a stale mask.
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
      // The set must shrink to reflect only the kept rows — proves the recheck recomputed despite
      // followCurrentRow being off.
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
    // The anti-storm gate (see _subscribeSearchProgress): a substructure search's progress events keep the
    // watchdog _externalSearchTimer live while results stream, and _scheduleFilterRecheck early-returns
    // (schedules nothing) while it is — so the per-batch requestFilter() storm launches NO similarity search
    // until the search finishes. This asserts that gate synchronously, the deterministic core of the throttle
    // (the event→timer wiring is verified separately; driving it through the global bus is Puppeteer-flaky).
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
      // Search streaming (watchdog live) → recheck suppressed: nothing scheduled.
      g._externalSearchTimer = setTimeout(() => {}, 30000);
      g._scheduleFilterRecheck();
      expect(g._filterRecheckTimer, null);
      // Search finished (watchdog cleared) → the recheck schedules normally.
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
    // Selected falls through to the base (raw selection); FilteredSelected is in getRowSourceIndexes's probe
    // branch AND intersects with the selection. Both must confine the similar set + applied filter to it.
    for (const mode of ['Selected', 'FilteredSelected']) {
      const {tv, viewer, df} = await openViewer();
      try {
        df.selection.setAll(false, false);
        const half = Math.max(1, Math.floor(df.rowCount / 2));
        for (let i = 0; i < half; i++)
          df.selection.set(i, true, false);
        viewer.setOptions({rowSource: mode});
        // Poll until the recompute confines the set (don't race the rowSource-change render, which would
        // otherwise read the stale Filtered-mode set built over all rows).
        await awaitCheck(() => confinedToSelection(viewer.similarSetBitset, half, df.rowCount),
          `${mode} similar set did not confine to the selection`, 12000);
        viewer.filterBtn!.click();
        await awaitCheck(() => viewer.filterActive && confinedToSelection(df.filter, half, df.rowCount),
          `filter did not confine to the selection in ${mode} mode`, 5000);
      } finally {tv.close();}
    }
  });

  test('two similarity viewers on one dataframe converge without a cascade', async () => {
    // The _probingDataframes guard stops a probe from looking like an external change to a sibling viewer.
    // Two collaborative filters on the same table must converge to a stable df.filter (bounded by the
    // _sameMask dedupe), not cascade forever — a runaway cascade would never let trueCount settle.
    const tv = await createTableView('tests/sar-small_test.csv');
    try {
      const df = tv.dataFrame;
      const v1 = (await df.plot.fromType('Chem Similarity Search')) as ChemSimilarityViewer;
      const v2 = (await df.plot.fromType('Chem Similarity Search')) as ChemSimilarityViewer;
      await awaitCheck(() => (v1.similarSetBitset?.trueCount ?? 0) > 0 && (v2.similarSetBitset?.trueCount ?? 0) > 0,
        'both viewers did not compute', 15000);
      v1.filterBtn!.click();
      await awaitCheck(() => v1.filterActive, 'v1 filter did not activate', 5000);
      v2.filterBtn!.click();
      await awaitCheck(() => v2.filterActive, 'v2 filter did not activate', 5000);
      // Poll until df.filter.trueCount holds steady across consecutive checks — i.e. the two viewers
      // reached the mutually-consistent fixed point. Times out (test fails) if it never stops changing.
      let last = -1;
      let stable = 0;
      await awaitCheck(() => {
        const c = df.filter.trueCount;
        stable = c === last ? stable + 1 : 0;
        last = c;
        return stable >= 5;
      }, 'two collaborative filters did not converge (possible recheck cascade)', 15000);
      expect(df.filter.trueCount <= df.rowCount, true);
    } finally {tv.close();}
  });

  test('selecting refreshes card highlights in place (no panel rebuild)', async () => {
    // The cards reflect selection (d4-selected); a selection change must restyle the EXISTING cards rather
    // than clear + rebuild the whole panel (which flashes every card and re-pops its molecule — the flicker).
    const {tv, viewer, df} = await openViewer();
    try {
      const gridBefore = viewer.root.querySelector('.chem-viewer-grid'); // the live card container
      expect(gridBefore !== null, true);
      viewer.selectBtn!.click(); // selects the similar set → the cards must show as selected
      await awaitCheck(() => df.selection.trueCount > 0, 'selection was not made', 5000);
      await delay(300); // let the debounced selection.onChanged render run
      // Same DOM node ⇒ the panel was refreshed in place, not rebuilt.
      expect(viewer.root.querySelector('.chem-viewer-grid') === gridBefore, true);
      // ...and a card actually picked up the new selection highlight.
      expect(viewer.root.querySelector('.chem-viewer-grid .d4-selected') !== null, true);
    } finally {tv.close();}
  });

  test('filter-panel label tracks the built cutoff, not a mid-search live value', async () => {
    // Directly exercises the _bitsetCutoff fix: _filterSummary() must describe the cutoff the mask was BUILT
    // with, not a live this.cutoff a slider drag moved before the next search settled. (The "cutoff change
    // updates the filter label" test reconstructs the label from live props, so it can't catch this.)
    const {tv, viewer} = await openViewer();
    const g = viewer as unknown as {cutoff: number; _bitsetCutoff: number; _filterSummary(): string};
    try {
      const built = Number(g._bitsetCutoff.toFixed(3)); // what the current similar set was built with
      // Move the LIVE cutoff to a clearly different value; read the label synchronously, before any new
      // search can settle and update _bitsetCutoff.
      g.cutoff = Number((built + 0.5).toFixed(3));
      const label = g._filterSummary();
      expect(label.includes(String(built)), true); // shows the BUILT cutoff
      expect(label.includes(String(g.cutoff)), false); // not the live drag value
    } finally {tv.close();}
  });
});
