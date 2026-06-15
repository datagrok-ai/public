import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {buildFacetGrid} from './facet-grid';
import {MAX_FACET_GRAPHS_COUNT} from './ui-constants';

import '../css/app-styles.css';

/** Faceted line-chart grid for Diff Studio multi-variable solutions: one chart per output
 *  variable, sharing the argument (X) column. Hosts platform line charts so rendering matches
 *  the standalone Diff Studio app 1:1. Exposed as a viewer so RFV can render it by name. */
@grok.decorators.viewer({
  name: 'DiffStudio Facet',
  description: 'Faceted grid of line charts, one per output variable, for Diff Studio solutions',
})
export class DiffStudioFacetViewer extends DG.JsViewer {
  xColumnName: string;
  segmentColumnName: string;
  maxGraphs: number;

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

  /** Replace the grid contents from the current data frame and properties. */
  render(): void {
    ui.empty(this.root);

    if (!this.dataFrame)
      return;

    const {root} = buildFacetGrid(this.dataFrame, {
      xColumnName: this.xColumnName || undefined,
      segmentColumnName: this.segmentColumnName || undefined,
      maxGraphs: this.maxGraphs,
    });

    root.classList.add('diff-studio-facet-grid');
    this.root.appendChild(root);
  }
}
