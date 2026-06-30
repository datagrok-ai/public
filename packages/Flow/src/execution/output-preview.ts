/** Lazy single docked output panel.
 *
 *  Behavior contract (post-rework):
 *  - Nothing is auto-docked when a run completes.
 *  - The panel is created only when the user clicks a node that has captured
 *    output values.
 *  - Subsequent clicks on completed nodes reuse the same dock node and just
 *    swap its content.
 *  - The panel is closed when the graph changes or a new run starts (values
 *    become stale anyway).
 *
 *  Rendering is delegated to {@link buildValuePanel} from `value-inspector.ts`
 *  — keeps the DataFrame / Column / graphics / primitive layout consistent
 *  whether values are read here or (previously) in the property panel. */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {NodeExecState} from './execution-state';
import {buildValuePreviews, hasRenderablePreview} from './value-inspector';
import {setTid} from '../utils/test-ids';

export class OutputPreviewPanel {
  private rootNode: DG.DockNode | null = null;
  private hostEl: HTMLElement | null = null;
  private viewRoot: HTMLElement | null = null;

  /** Set the Flow view root so the panel docks relative to it (bottom edge). */
  setViewRoot(root: HTMLElement): void {
    this.viewRoot = root;
  }

  /** Show the runtime values for one node. If the node has no values
   *  captured yet, this is a no-op (the panel is not opened). On the first
   *  call with a value-bearing state, the panel is lazily docked at the
   *  bottom of the view; subsequent calls update its content in place. */
  showForNode(node: {id: string; label: string}, state: NodeExecState | undefined): void {
    if (!state) return;
    // Status, duration, error, primitives — those go in the property panel
    // (see `buildExecutionMeta`). The docked panel only opens when there's
    // something *renderable*: DataFrame grid, column sample, graphics image.
    if (!hasRenderablePreview(state)) return;

    const inner = buildValuePreviews(state);
    inner.style.padding = '8px 12px';

    // Reuse the open dock — but only if it's still actually docked. The user may
    // have manually closed the panel (the dock manager won't tell us), leaving
    // our refs stale; if so, drop them and re-dock below so the click reopens it.
    if (this.hostEl && this.rootNode) {
      const stillDocked = (() => {
        try {
          return !!grok.shell.dockManager.findNode(this.hostEl);
        } catch {
          return false;
        }
      })();
      if (stillDocked) {
        this.hostEl.innerHTML = '';
        this.hostEl.appendChild(inner);
        return;
      }
      this.hostEl = null;
      this.rootNode = null;
    }

    this.hostEl = setTid(ui.div([inner], {style: {
      width: '100%', height: '100%', overflow: 'auto',
    }}), 'output-panel');
    const refNode = this.viewRoot ?
      grok.shell.dockManager.findNode(this.viewRoot) ?? null :
      null;
    try {
      this.rootNode = grok.shell.dockManager.dock(
        this.hostEl, DG.DOCK_TYPE.DOWN, refNode, 'Node Output', 0.4,
      );
    } catch (e) {
      console.warn('OutputPreview: failed to dock', e);
      this.hostEl = null;
    }
  }

  /** Close the panel if open. Used on graph change / new run. */
  close(): void {
    if (this.hostEl) {
      try {
        const node = grok.shell.dockManager.findNode(this.hostEl);
        if (node) grok.shell.dockManager.close(node);
      } catch { /* already closed by user */ }
    }
    this.hostEl = null;
    this.rootNode = null;
  }
}
