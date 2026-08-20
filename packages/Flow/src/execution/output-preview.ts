/** In-view bottom output panel — a pane of the Flow view's vertical splitter,
 *  never a dock-manager dock (a dock outlives and overlaps views). */

import * as ui from 'datagrok-api/ui';

import {NodeExecState} from './execution-state';
import {buildValuePreviews, hasRenderablePreview} from './value-inspector';
import {setTid} from '../utils/test-ids';

export type OutputPanelState = 'hidden' | 'minimized' | 'expanded';

/** Height of the header strip — the whole panel when minimized. */
const HEADER_HEIGHT = 30;
const DEFAULT_EXPANDED_HEIGHT = 260;

export class OutputPreviewPanel {
  readonly root: HTMLElement;
  private readonly contentEl: HTMLElement;
  private readonly nodeLabelEl: HTMLElement;
  private readonly caretEl: HTMLElement;
  private readonly pinEl: HTMLElement;

  private state: OutputPanelState = 'hidden';
  /** Remembered for the view's lifetime — new content never pops the panel back up. */
  private userMinimized = false;
  private expandedHeight = DEFAULT_EXPANDED_HEIGHT;
  private enabled: boolean;

  /** `ExecutionState.setNodeStatus` always builds a fresh state object, so
   *  reference identity of the state IS value identity — same pair means the
   *  preview must not rebuild (grids re-mount, scroll resets). */
  private lastNodeId: string | null = null;
  private lastState: NodeExecState | null = null;

  /** Pinned node: other clicks don't switch the preview; survives `clear()`. */
  private pinnedId: string | null = null;

  private updating = false;

  /** Called when the user clicks the gear on a viewer preview. */
  onEditViewer?: (node: {id: string; label: string}, viewer: unknown) => void;

  /** Fired on every hidden/minimized/expanded transition. */
  onStateChanged?: (state: OutputPanelState) => void;

  /** Fired on a user unpin via the header icon (not on programmatic `unpin()`). */
  onUnpinned?: () => void;

  constructor(options: {enabled?: boolean} = {}) {
    this.enabled = options.enabled !== false;

    // Only the caret toggles — a fully clickable header would sit right under
    // the splitter divider and swallow near-miss resize clicks.
    this.caretEl = setTid(ui.div([], 'ff-output-panel-caret'), 'output-panel-caret');
    this.caretEl.addEventListener('click', () => this.toggle());
    ui.tooltip.bind(this.caretEl, () => this.state === 'minimized' ? 'Expand outputs' : 'Minimize outputs');
    this.nodeLabelEl = setTid(ui.div([], 'ff-output-panel-node'), 'output-panel-node');
    this.pinEl = setTid(ui.iconFA('thumbtack', () => this.togglePin()), 'output-panel-pin');
    this.pinEl.classList.add('ff-output-panel-pin');
    ui.tooltip.bind(this.pinEl, () => this.pinnedId != null ?
      'Unpin — clicking a node switches the preview again' :
      'Pin this preview — keep it while you adjust other nodes');
    const title = ui.div([], 'ff-output-panel-title');
    title.textContent = 'Preview';
    const header = setTid(
      ui.div([title, this.nodeLabelEl, this.pinEl, this.caretEl], 'ff-output-panel-header'), 'output-panel-header');

    this.contentEl = setTid(ui.div([], 'ff-output-panel-content'), 'output-panel-content');

    this.root = setTid(ui.box(), 'output-panel');
    this.root.classList.add('ff-output-panel');
    // The canvas (flex: 1 1 0) absorbs all remaining space; the pane's height IS its size.
    this.root.style.flex = '0 0 auto';
    this.root.appendChild(header);
    this.root.appendChild(this.contentEl);
    this.applyState();
    this.updatePinVisual();
  }

  get pinnedNodeId(): string | null {
    return this.pinnedId;
  }

  togglePin(): void {
    if (this.pinnedId != null) {
      this.unpin();
      this.onUnpinned?.();
      return;
    }
    if (this.lastNodeId != null) {
      this.pinnedId = this.lastNodeId;
      this.updatePinVisual();
    }
  }

  unpin(): void {
    this.pinnedId = null;
    this.updatePinVisual();
  }

  /** Overlay a spinner on the kept (stale) content while the pinned node recomputes. */
  markUpdating(message = 'Recalculating...'): void {
    if (!this.enabled || this.state === 'hidden' || this.lastNodeId == null) return;
    this.updating = true;
    ui.setUpdateIndicator(this.contentEl, true, message);
  }

  clearUpdating(): void {
    if (!this.updating) return;
    this.updating = false;
    ui.setUpdateIndicator(this.contentEl, false);
  }

  private updatePinVisual(): void {
    this.pinEl.dataset.pinned = this.pinnedId != null ? 'true' : 'false';
    // A pinned-but-cleared panel keeps the icon so the user can still unpin.
    this.pinEl.style.display = this.lastNodeId != null || this.pinnedId != null ? '' : 'none';
  }

  get panelState(): OutputPanelState {
    return this.state;
  }

  get isEnabled(): boolean {
    return this.enabled;
  }

  get currentNodeId(): string | null {
    return this.lastNodeId;
  }

  /** Off = never shows (embedded hosts). */
  setEnabled(enabled: boolean): void {
    this.enabled = enabled;
    if (!enabled) this.clear();
  }

  /** Show the runtime values for one node; no-op when disabled or nothing renderable. */
  showForNode(node: {id: string; label: string}, state: NodeExecState | undefined): void {
    if (!this.enabled || !state) return;
    if (this.pinnedId != null && node.id !== this.pinnedId) return;
    if (!hasRenderablePreview(state)) return;

    // Same node, same captured state → rebuilding would re-mount the grids and
    // reset their scroll.
    if (this.state !== 'hidden' && node.id === this.lastNodeId && state === this.lastState) {
      this.clearUpdating();
      return;
    }

    this.clearUpdating();
    const inner = buildValuePreviews(state, (viewer) => this.onEditViewer?.(node, viewer));
    inner.style.padding = '8px 12px';
    this.contentEl.innerHTML = '';
    this.contentEl.appendChild(inner);
    this.nodeLabelEl.textContent = node.label;
    this.lastNodeId = node.id;
    this.lastState = state;
    this.updatePinVisual();

    if (this.state === 'hidden')
      this.setState(this.userMinimized ? 'minimized' : 'expanded');
  }

  /** Rebuild the shown node's preview from its kept state — for when what a
   *  block renders depends on outside state (the in-node preview toggled). */
  refresh(): void {
    const nodeId = this.lastNodeId;
    const state = this.lastState;
    if (nodeId == null || state == null) return;
    this.lastState = null; // defeat the identity gate — same state must rebuild
    this.showForNode({id: nodeId, label: this.nodeLabelEl.textContent ?? ''}, state);
  }

  minimize(): void {
    this.userMinimized = true;
    if (this.state === 'expanded') {
      const h = this.root.offsetHeight;
      if (h > HEADER_HEIGHT + 10) this.expandedHeight = h;
      this.setState('minimized');
    }
  }

  expand(): void {
    this.userMinimized = false;
    if (this.state === 'minimized') this.setState('expanded');
  }

  toggle(): void {
    if (this.state === 'minimized') this.expand();
    else if (this.state === 'expanded') this.minimize();
  }

  /** Empty and hide; the user's minimized preference AND the pin survive. */
  clear(): void {
    this.clearUpdating();
    this.contentEl.innerHTML = '';
    this.nodeLabelEl.textContent = '';
    this.lastNodeId = null;
    this.lastState = null;
    this.setState('hidden');
    this.updatePinVisual();
  }

  private setState(state: OutputPanelState): void {
    if (state === this.state) return;
    this.state = state;
    this.applyState();
    this.onStateChanged?.(state);
  }

  private applyState(): void {
    const s = this.root.style;
    if (this.state === 'hidden') {
      s.display = 'none';
      return;
    }
    s.display = '';
    // `ui.splitV` keeps rewriting `style.height` on every pane; min/max clamp
    // the rendered size so its resize handling can't drift a minimized strip.
    if (this.state === 'minimized') {
      s.height = `${HEADER_HEIGHT}px`;
      s.maxHeight = `${HEADER_HEIGHT}px`;
      this.contentEl.style.display = 'none';
      this.caretEl.textContent = '▴';
      s.minHeight = `${HEADER_HEIGHT}px`;
    } else {
      s.height = `${this.expandedHeight}px`;
      s.maxHeight = '';
      this.contentEl.style.display = '';
      this.caretEl.textContent = '▾';
      s.minHeight = `${180}px`;
    }
  }
}
