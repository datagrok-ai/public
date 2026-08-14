import {Subscription} from 'rxjs';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// A freshly-mounted Filters viewer asynchronously reapplies a stale per-column selection over the
// DataFrame's .filter shortly after construction; this is how long to wait past that window.
const FILTER_REMOUNT_SETTLE_MS = 200;

/**
 * Owns every mounted Grid/Filters viewer for one Reaction Enumerator view, plus the deferred
 * filter-reset timers scheduled against them, and drains `view.subs` itself.
 *
 * `View.detach()` never calls `ViewBase.detach()`, so `view.subs` would never unsubscribe on its
 * own for this view type — cleanup runs off `grok.events.onViewRemoved` instead, filtered by
 * `view.id` so several open enumerator tabs don't cancel each other's viewers and timers.
 */
export class MountedViewerRegistry {
  private readonly mountedViewers = new Map<HTMLElement, DG.Viewer[]>();
  private readonly pendingTimers = new Set<ReturnType<typeof setTimeout>>();
  private readonly closeCallbacks: (() => void)[] = [];
  private closeSub: Subscription;

  constructor(private readonly view: DG.View) {
    this.closeSub = grok.events.onViewRemoved.subscribe((closedView) => {
      if (closedView.id !== this.view.id) return;
      for (const host of this.mountedViewers.keys()) this.close(host);
      for (const id of this.pendingTimers) clearTimeout(id);
      this.pendingTimers.clear();
      for (const cb of this.closeCallbacks) cb();
      this.view.subs.forEach((s) => s.unsubscribe());
      this.closeSub.unsubscribe();
    });
  }

  /** Lets other classes react to this view closing — e.g. an in-flight computation cancelling
   * itself instead of mounting a viewer into an already-torn-down view. */
  onClose(callback: () => void): void {
    this.closeCallbacks.push(callback);
  }

  /** Viewers leak unless closed explicitly — `host.innerHTML = ''` only drops the DOM node. */
  close(host: HTMLElement): void {
    const prev = this.mountedViewers.get(host);
    if (!prev) return;
    this.mountedViewers.delete(host);
    for (const v of prev) {
      // Viewer.close() no-ops for a viewer never docked in a TableView; only detach() (which
      // unsubscribes from the DataFrame) plus removing the node actually releases it.
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

  // A single RAF after appending a resizable ui.splitH isn't reliably enough for a real clientWidth;
  // ui.onSizeChanged fires once real layout completes. Size once, then unsubscribe.
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
      // ChemicalReaction has no meaningful substructure-filter semantics for a whole template, and
      // `~`-prefixed columns are internal (e.g. the substructure filter's own fingerprint cache).
      .filter((col) => col.semType !== 'ChemicalReaction' && !col.name.startsWith('~'))
      .map((col) => ({
        // A zero-variance numeric column (common right after subsetting) makes the histogram filter
        // widget clobber the whole frame's .filter to all-false. The categorical filter has no
        // degenerate-range path, so route those there instead of ever building the histogram.
        type: col.isNumerical ?
          (col.stats.min === col.stats.max ? DG.FILTER_TYPE.CATEGORICAL : DG.FILTER_TYPE.HISTOGRAM) :
          col.semType === DG.SEMTYPE.MOLECULE ? DG.FILTER_TYPE.SUBSTRUCTURE : DG.FILTER_TYPE.CATEGORICAL,
        column: col.name,
      }));
    const filtersViewer = DG.Viewer.filters(df, {filters: filterStates});
    this.mountedViewers.set(host, [grid, filtersViewer]);
    filtersViewer.root.style.width = '100%';
    filtersViewer.root.style.height = '100%';
    filtersViewer.root.style.overflow = 'auto';
    const gridBox = ui.div([grid.root], {style: {flex: '1 1 0', minWidth: '0', height: '100%', overflow: 'hidden'}});
    // Width must stay relative: the wrapper sizeSplitOnceLaidOut resizes below varies, and a fixed
    // px value here would drift out of sync with it and open a gap against the divider.
    const filtersBox = ui.div([filtersViewer.root], {style: {flex: '0 0 auto', width: '100%', height: '100%',
      overflow: 'auto', borderRight: '1px solid var(--grey-2)'}});
    const split = ui.splitH([filtersBox, gridBox], {style: {width: '100%', height: '100%', minHeight: '0'}}, true);
    host.appendChild(split);
    this.applyGridColumnSizing(grid, opts?.extendLastColumn ?? true);
    // ui.splitH ignores child width/flex on first layout; its resize handler reads flexGrow off the
    // wrapper boxes it creates (children[0]/[2], skipping the divider at [1]).
    const filtersWrap = split.children[0] as HTMLElement;
    const gridWrap = split.children[2] as HTMLElement;
    if (gridWrap && filtersWrap)
      this.sizeSplitOnceLaidOut(filtersWrap, gridWrap, (total) => Math.min(260, Math.round(total * 0.25)));
  }

  private scheduleTimer(fn: () => void): void {
    const id = setTimeout(() => {
      this.pendingTimers.delete(id);
      fn();
    }, FILTER_REMOUNT_SETTLE_MS);
    this.pendingTimers.add(id);
  }

  /** Re-asserts "show every row" after a mount, past the Filters clobber window. */
  deferredFilterReset(df: DG.DataFrame): void {
    this.scheduleTimer(() => df.filter.setAll(true, true));
  }

  /** Snapshots each table's filter across a cascade that remounts the grids (a YAML load), then
   * restores it once the same clobber window has passed. */
  withPreservedFilters(dfs: DG.DataFrame[], fn: () => void): void {
    const saved = dfs.map((d) => ({df: d, mask: d.filter.clone()}));
    fn();
    this.scheduleTimer(() => {
      for (const {df, mask} of saved) df.filter.copyFrom(mask, true);
    });
  }
}
