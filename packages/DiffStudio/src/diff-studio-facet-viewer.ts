import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {buildFacetGrid} from './facet-grid';
import {MAX_FACET_GRAPHS_COUNT} from './ui-constants';

import '../css/app-styles.css';

/** Faceted line-chart grid for Diff Studio multi-variable solutions: one chart per output
 *  variable, sharing the argument (X) column. Hosts platform line charts so rendering matches
 *  the standalone Diff Studio app 1:1. Registered as a viewer (via PackageFunctions) so RFV can
 *  render it by name; `meta.showInGallery: false` keeps it out of the Add-viewer gallery. */
export class DiffStudioFacetViewer extends DG.JsViewer {
  xColumnName: string;
  segmentColumnName: string;
  maxGraphs: number;

  /** Line charts of the currently mounted grid, reused across runs to avoid flicker. */
  private plots: DG.Viewer[] = [];
  /** Signature of the mounted grid's layout; a change forces a full rebuild. */
  private layoutKey: string | null = null;

  constructor() {
    super();
    // The JsViewer root is a `ui-box` hosted inside a `display:block` <dg-viewer>; without a
    // definite height it collapses to 0 in RFV (its split-based content is itself a flex tree
    // that only fills once the root has a height). Pin it to fill the host.
    this.root.classList.add('diff-studio-facet-viewer');
    this.xColumnName = this.string('xColumnName', '');
    this.segmentColumnName = this.string('segmentColumnName', '');
    this.maxGraphs = this.int('maxGraphs', MAX_FACET_GRAPHS_COUNT);
  }

  /** Re-render on attach and whenever the viewer is resized. */
  onTableAttached(): void {
    this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe(() => this.render()));
    this.render();
  }

  /** Rebuild the grid when a configurable property changes. */
  onPropertyChanged(property: DG.Property | null): void {
    super.onPropertyChanged(property);
    this.render();
  }

  detach(): void {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  /** Render the grid from the current data frame and properties. When the layout is unchanged
   *  (same columns and config), the existing line charts are kept and only their data frame is
   *  swapped — rebuilding from scratch on every run causes a visible flicker. */
  render(): void {
    if (!this.dataFrame) {
      ui.empty(this.root);
      this.plots = [];
      this.layoutKey = null;
      return;
    }

    const key = [
      this.dataFrame.columns.names().join('\t'),
      this.xColumnName,
      this.segmentColumnName,
      this.maxGraphs,
    ].join('|');

    // Same layout as the mounted grid: reuse the charts, just point them at the new data.
    if (this.layoutKey === key && this.plots.length > 0) {
      this.plots.forEach((plot) => plot.dataFrame = this.dataFrame!);
      return;
    }

    ui.empty(this.root);

    const {root, plots} = buildFacetGrid(this.dataFrame, {
      xColumnName: this.xColumnName || undefined,
      segmentColumnName: this.segmentColumnName || undefined,
      maxGraphs: this.maxGraphs,
    });

    root.classList.add('diff-studio-facet-grid');
    this.root.appendChild(root);

    this.plots = plots;
    this.layoutKey = key;
  } // render
}
