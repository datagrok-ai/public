import {Subscription} from 'rxjs';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// A freshly-(re)mounted Filters viewer can asynchronously reapply a stale per-column categorical
// selection over a DataFrame's .filter BitSet shortly after construction (a Datagrok platform
// behavior, not something this app triggers) — this is how long deferredFilterReset/
// withPreservedFilters wait before re-asserting the intended filter state, past that clobber window.
const FILTER_REMOUNT_SETTLE_MS = 200;

/**
 * Owns every mounted Grid/Filters viewer for one Reaction Enumerator view instance, plus every
 * deferred filter-reset timer scheduled against them. `view.subs` does not unsubscribe for this
 * app-hosted view type — cleanup runs off `grok.events.onViewRemoved` instead (filtered by
 * `view.id`, so multiple open Reaction Enumerator tabs never cancel each other's timers/viewers).
 */
export class MountedViewerRegistry {
  private readonly mountedViewers = new Map<HTMLElement, DG.Viewer[]>();
  private readonly pendingTimers = new Set<ReturnType<typeof setTimeout>>();
  private closeSub: Subscription;

  constructor(private readonly view: DG.View) {
    this.closeSub = grok.events.onViewRemoved.subscribe((closedView) => {
      if (closedView.id !== this.view.id) return;
      for (const host of this.mountedViewers.keys()) this.close(host);
      for (const id of this.pendingTimers) clearTimeout(id);
      this.pendingTimers.clear();
      this.closeSub.unsubscribe();
    });
  }

  /** Plain viewers (grid, filters) leak (host.innerHTML='' only drops the DOM node) unless closed
   * explicitly — close the ones tracked for `host` before mounting a replacement. */
  close(host: HTMLElement): void {
    const prev = this.mountedViewers.get(host);
    if (!prev) return;
    this.mountedViewers.delete(host);
    for (const v of prev) {
      // close() no-ops for standalone viewers (see /ui skill) — detach() + root.remove() instead.
      try {v.detach(); v.root.remove();} catch (e) {
        if (!(e instanceof TypeError)) console.warn('Could not close previous viewer:', e);
      }
    }
  }

  private applyGridColumnSizing(grid: DG.Grid, extendLast = true): void {
    try {
      grid.setColumnsWidthType(DG.ColumnWidthType.Optimal);
      if (extendLast) grid.props.extendLastColumn = true;
    } catch {/* setColumnsWidthType not available on older Dart builds */}
  }

  // A single RAF after appending a resizable ui.splitH isn't reliably enough for a real clientWidth
  // yet — ui.onSizeChanged (ResizeObserver-backed) fires once real layout completes; size once,
  // then unsubscribe. Tracked in view.subs in case the view closes before that first size arrives.
  private sizeSplitOnceLaidOut(a: HTMLElement, b: HTMLElement, computeAWidth: (total: number) => number): void {
    const sub = ui.onSizeChanged(a).subscribe(() => {
      const total = a.clientWidth + b.clientWidth;
      if (total === 0) return;
      const aWidth = computeAWidth(total);
      a.style.width = aWidth + 'px';
      a.style.flexGrow = String(aWidth / total);
      b.style.width = (total - aWidth) + 'px';
      b.style.flexGrow = String((total - aWidth) / total);
      sub.unsubscribe();
    });
    this.view.subs.push(sub);
  }

  // ChemicalReaction has no meaningful substructure-filter semantics for a whole reaction template,
  // so it's simply left out of the filters below (not a workaround for anything — DG.Viewer.filters
  // only builds filters for the columns it's explicitly given).
  mountDf(host: HTMLElement, df: DG.DataFrame, withFilters: boolean,
    opts?: {rowHeight?: number; extendLastColumn?: boolean}): void {
    this.close(host);
    host.innerHTML = '';
    const grid = DG.Viewer.grid(df);
    grid.props.rowHeight = opts?.rowHeight ?? 75;
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    if (!withFilters) {
      this.mountedViewers.set(host, [grid]);
      host.appendChild(grid.root);
      this.applyGridColumnSizing(grid, opts?.extendLastColumn ?? true);
      return;
    }
    const filterStates = df.columns.toList()
      // `~`-prefixed columns are internal/technical (e.g. `~SMILES.Pattern`, added by the
      // substructure filter's own fingerprint cache) — never surface a filter for those.
      .filter((col) => col.semType !== 'ChemicalReaction' && !col.name.startsWith('~'))
      .map((col) => ({
        // A numeric column with zero variance (every value equal — common right after subsetting to
        // the exact value it was filtered on) makes the histogram filter widget clobber the whole
        // dataframe's .filter to all-false, both on construction and again if the user later toggles
        // it histogram<->categorical — a separate bug from, and not fixed by, the subset-clone filter
        // reset above. A categorical filter has no degenerate-range code path, so route zero-variance
        // numeric columns there instead of ever constructing the buggy histogram widget.
        type: col.isNumerical ?
          (col.stats.min === col.stats.max ? DG.FILTER_TYPE.CATEGORICAL : DG.FILTER_TYPE.HISTOGRAM) :
          col.semType === 'Molecule' ? DG.FILTER_TYPE.SUBSTRUCTURE : DG.FILTER_TYPE.CATEGORICAL,
        column: col.name,
      }));
    const filtersViewer = DG.Viewer.filters(df, {filters: filterStates});
    this.mountedViewers.set(host, [grid, filtersViewer]);
    filtersViewer.root.style.width = '100%';
    filtersViewer.root.style.height = '100%';
    filtersViewer.root.style.overflow = 'auto';
    const gridBox = ui.div([grid.root], {style: {flex: '1 1 0', minWidth: '0', height: '100%', overflow: 'hidden'}});
    // width must track the wrapper sizeSplitOnceLaidOut resizes below (split.children[0]), not a fixed
    // px value — a fixed width here left a dead gap between this box and the divider whenever the
    // wrapper's computed width (total * 0.25, capped 260) exceeded it.
    const filtersBox = ui.div([filtersViewer.root], {style: {flex: '0 0 auto', width: '100%', height: '100%',
      overflow: 'auto', borderRight: '1px solid var(--grey-2)'}});
    const split = ui.splitH([filtersBox, gridBox], {style: {width: '100%', height: '100%', minHeight: '0'}}, true);
    host.appendChild(split);
    this.applyGridColumnSizing(grid, opts?.extendLastColumn ?? true);
    // ui.splitH ignores child width/flex style on first layout; its resize handler reads flexGrow
    // off the wrapper boxes it creates (children[0]/[2], skipping the divider at [1]) — set size there.
    const filtersWrap = split.children[0] as HTMLElement;
    const gridWrap = split.children[2] as HTMLElement;
    if (gridWrap && filtersWrap)
      this.sizeSplitOnceLaidOut(filtersWrap, gridWrap, (total) => Math.min(260, Math.round(total * 0.25)));
  }

  // Shared tracked-timer lifecycle behind both deferredFilterReset and withPreservedFilters — a
  // closing view cancels every pending id via `pendingTimers` (see the constructor), so both callers
  // route through this one bookkeeping path instead of repeating it.
  private scheduleTimer(fn: () => void): void {
    const id = setTimeout(() => {
      this.pendingTimers.delete(id);
      fn();
    }, FILTER_REMOUNT_SETTLE_MS);
    this.pendingTimers.add(id);
  }

  // Shared by every "this table should show every row" reset after a mount: doSubsetBySelection,
  // doRestoreFullTable, subsetStepBySelection, useAllForStep.
  deferredFilterReset(df: DG.DataFrame): void {
    this.scheduleTimer(() => df.filter.setAll(true, true));
  }

  // Rebuilding the quick inputs from a loaded YAML config (see syncConfigToQuickInputs) cascades
  // into a full remount of the Reactions/BBs/Reagents grids — and each freshly-mounted Filters
  // viewer can silently reapply a stale, much-earlier categorical selection instead of the table's
  // actual current filter, since the reapply happens in the filter widget's own async
  // post-construction step (same root cause as deferredFilterReset above). Snapshot each table's
  // filter before the cascade and restore it after, once that async reapply has had time to fire.
  withPreservedFilters(dfs: DG.DataFrame[], fn: () => void): void {
    const saved = dfs.map((d) => ({df: d, mask: d.filter.clone()}));
    fn();
    this.scheduleTimer(() => {
      for (const {df, mask} of saved) df.filter.copyFrom(mask, true);
    });
  }
}
