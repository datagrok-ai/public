/** Singleton docked output preview panel.
 *
 * After a flow run completes, this module creates a docked panel at the bottom
 * of the Datagrok shell showing output previews as tabs.  DataFrames get a
 * full grid viewer; other supported types get lightweight representations.
 *
 * The panel is a singleton: if already docked, content is updated in-place.
 * Users can close individual tabs or the whole panel — re-running the flow
 * recreates everything. */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

/** Output types we know how to preview */
type PreviewableOutput =
  | {type: 'dataframe'; name: string; value: DG.DataFrame}
  | {type: 'viewer'; name: string; value: DG.Viewer}
  | {type: 'graphics'; name: string; value: string} // image data as string
  | {type: 'widget'; name: string; value: DG.Widget}
  | {type: 'primitive'; name: string; value: any};

export class OutputPreviewPanel {
  private rootNode: DG.DockNode | null = null;
  /** Track ALL docked elements so close() can remove every tab */
  private dockedElements: HTMLElement[] = [];
  /** The view root to dock relative to */
  private viewRoot: HTMLElement | null = null;

  /** Set the Flow view root so panels dock relative to it */
  setViewRoot(root: HTMLElement): void {
    this.viewRoot = root;
  }

  /** Clear previous outputs and show new ones.
   *  @param outputs — named output values from fc.outputs
   *  @param typeHints — optional map of output name → declared DG type (e.g. 'graphics', 'dataframe') */
  showOutputs(outputs: Record<string, any>, typeHints?: Record<string, string>): void {
    this.close();

    const entries = Object.entries(outputs);
    if (entries.length === 0) return;

    const previews = entries
      .map(([name, value]) => this.classifyOutput(name, value, typeHints?.[name]))
      .filter((p): p is PreviewableOutput => p !== null);

    if (previews.length === 0) {
      for (const [name, value] of entries)
        console.log(`Flow output [${name}]:`, value);
      return;
    }

    // Find the flow view's dock node to dock relative to
    const refNode = this.viewRoot ?
      grok.shell.dockManager.findNode(this.viewRoot) ?? null :
      null;

    // Dock the first output at the bottom of the flow view
    const first = previews[0];
    const firstEl = this.buildPreview(first);
    this.dockedElements.push(firstEl);

    try {
      this.rootNode = grok.shell.dockManager.dock(
        firstEl, DG.DOCK_TYPE.DOWN, refNode, first.name, 0.4,
      );
    } catch (e) {
      console.warn('OutputPreview: failed to dock first output', e);
      return;
    }

    // Dock remaining outputs as tabs (FILL) onto the first panel
    for (let i = 1; i < previews.length; i++) {
      const p = previews[i];
      const el = this.buildPreview(p);
      this.dockedElements.push(el);
      try {
        grok.shell.dockManager.dock(
          el, DG.DOCK_TYPE.FILL, this.rootNode, p.name,
        );
      } catch (e) {
        console.warn(`OutputPreview: failed to dock output "${p.name}"`, e);
      }
    }
  }

  /** Close all docked output panels */
  close(): void {
    for (const el of this.dockedElements) {
      try {
        const node = grok.shell.dockManager.findNode(el);
        if (node)
          grok.shell.dockManager.close(node);
      } catch {/* already closed by user */}
    }
    this.rootNode = null;
    this.dockedElements = [];
  }

  private classifyOutput(name: string, value: any, typeHint?: string): PreviewableOutput | null {
    if (value == null) return null;

    // DataFrame
    if (value.rowCount !== undefined && value.columns !== undefined)
      return {type: 'dataframe', name, value: value as DG.DataFrame};

    // Viewer
    if (value instanceof DG.Viewer)
      return {type: 'viewer', name, value};

    // Widget
    if (value instanceof DG.Widget)
      return {type: 'widget', name, value};

    // Graphics — detected by type hint or by content inspection
    if (typeHint === 'graphics' && typeof value === 'string')
      return {type: 'graphics', name, value};
    if (typeof value === 'string' && (value.startsWith('data:image') || value.startsWith('<svg')))
      return {type: 'graphics', name, value};

    // Primitives (string, number, bool) — only if simple enough to display
    if (typeof value === 'string' || typeof value === 'number' || typeof value === 'boolean')
      return {type: 'primitive', name, value};

    // Unknown complex types — just log
    console.log(`Flow output [${name}]:`, value);
    return null;
  }

  private buildPreview(preview: PreviewableOutput): HTMLElement {
    switch (preview.type) {
    case 'dataframe':
      return this.buildDataframePreview(preview.name, preview.value);
    case 'viewer':
      return preview.value.root;
    case 'widget':
      return preview.value.root;
    case 'graphics':
      return this.buildGraphicsPreview(preview.value);
    case 'primitive':
      return this.buildPrimitivePreview(preview.name, preview.value);
    }
  }

  private buildDataframePreview(name: string, df: DG.DataFrame): HTMLElement {
    const container = ui.div([], {style: {width: '100%', height: '100%'}});
    const grid = DG.Viewer.grid(df);
    container.appendChild(grid.root);
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    return container;
  }

  private buildGraphicsPreview(imageData: string): HTMLElement {
    const container = ui.div([], {style: {width: '100%', height: '100%', overflow: 'auto',
      position: 'relative',
      backgroundPosition: 'left',
      backgroundRepeat: 'no-repeat',
      backgroundSize: 'contain',
    }});
    if (imageData.startsWith('<svg'))
      container.innerHTML = imageData;
    else
      container.style.backgroundImage = `url('data:image/png;base64,${imageData}')`;

    return container;
  }

  private buildPrimitivePreview(name: string, value: any): HTMLElement {
    const container = ui.div([]);
    container.appendChild(ui.divText(`${name} = ${JSON.stringify(value)}`));
    container.style.padding = '12px';
    container.style.fontSize = '14px';
    return container;
  }
}
