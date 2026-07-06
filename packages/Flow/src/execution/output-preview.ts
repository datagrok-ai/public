/** In-view bottom output panel.
 *
 *  Behavior contract (splitter rework — replaces the old dock-manager panel):
 *  - The panel is a permanent pane of the Flow view's vertical splitter
 *    (`ui.splitV`), owned by the view — never docked to the platform dock
 *    manager, so it can never leak over other views or outlive its flow.
 *  - It is not closable. It starts fully hidden; the first renderable output
 *    (clicking a completed node, or a run's focus node) expands it. The header
 *    caret minimizes it to a slim strip; the choice is remembered — while
 *    minimized, later clicks update the content but never pop the panel back
 *    up. Restoring is always an explicit header click.
 *  - `clear()` empties and hides it (graph change / new run — values are
 *    stale), preserving the user's minimized preference.
 *  - Embedded hosts (the creation-script dialog) construct it disabled: it
 *    never shows there, only in the real editor view.
 *
 *  Rendering is delegated to {@link buildValuePreviews} from
 *  `value-inspector.ts` — keeps the DataFrame / Column / graphics / primitive
 *  layout consistent with the context panel. */

import * as ui from 'datagrok-api/ui';

import {NodeExecState} from './execution-state';
import {buildValuePreviews, hasRenderablePreview} from './value-inspector';
import {setTid} from '../utils/test-ids';

export type OutputPanelState = 'hidden' | 'minimized' | 'expanded';

/** Height of the header strip — the whole panel when minimized. */
const HEADER_HEIGHT = 30;
/** Default expanded height; replaced by the real height once the user resizes
 *  (captured on minimize, and the splitter divider writes `style.height`). */
const DEFAULT_EXPANDED_HEIGHT = 260;

export class OutputPreviewPanel {
  /** The splitter pane (`ui.box`) the view mounts as the bottom `splitV` item. */
  readonly root: HTMLElement;
  private readonly contentEl: HTMLElement;
  private readonly nodeLabelEl: HTMLElement;
  private readonly caretEl: HTMLElement;

  private state: OutputPanelState = 'hidden';
  /** The user minimized the panel — remembered for the view's lifetime so
   *  showing new content never pops the panel back up uninvited. */
  private userMinimized = false;
  private expandedHeight = DEFAULT_EXPANDED_HEIGHT;
  /** Disabled panels never show — embedded hosts (dialogs) construct it so. */
  private enabled: boolean;

  /** What the panel currently renders. `ExecutionState.setNodeStatus` always
   *  builds a fresh state object, so reference identity of the state IS value
   *  identity — re-clicking the same node with unchanged results must not
   *  rebuild the preview (grids re-mount, scroll resets, the panel jumps). */
  private lastNodeId: string | null = null;
  private lastState: NodeExecState | null = null;

  /** Called when the user clicks "Edit settings" on a viewer preview — the host
   *  shows the viewer in the context panel and captures its option changes. */
  onEditViewer?: (node: {id: string; label: string}, viewer: unknown) => void;

  /** Fired on every hidden/minimized/expanded transition — the view syncs the
   *  splitter divider (resizing a hidden or minimized pane makes no sense). */
  onStateChanged?: (state: OutputPanelState) => void;

  constructor(options: {enabled?: boolean} = {}) {
    this.enabled = options.enabled !== false;

    // Only the caret toggles — a fully clickable header would sit right under
    // the splitter divider and swallow near-miss resize clicks.
    this.caretEl = setTid(ui.div([], 'ff-output-panel-caret'), 'output-panel-caret');
    this.caretEl.addEventListener('click', () => this.toggle());
    ui.tooltip.bind(this.caretEl, () => this.state === 'minimized' ? 'Expand outputs' : 'Minimize outputs');
    this.nodeLabelEl = setTid(ui.div([], 'ff-output-panel-node'), 'output-panel-node');
    const title = ui.div([], 'ff-output-panel-title');
    title.textContent = 'Outputs';
    const header = setTid(
      ui.div([title, this.nodeLabelEl, this.caretEl], 'ff-output-panel-header'), 'output-panel-header');

    this.contentEl = setTid(ui.div([], 'ff-output-panel-content'), 'output-panel-content');

    this.root = setTid(ui.box(), 'output-panel');
    this.root.classList.add('ff-output-panel');
    // Pin the pane so the canvas (flex: 1 1 0) absorbs all remaining space:
    // the pane's height IS its size, in every state.
    this.root.style.flex = '0 0 auto';
    this.root.appendChild(header);
    this.root.appendChild(this.contentEl);
    this.applyState();
  }

  get panelState(): OutputPanelState {
    return this.state;
  }

  get isEnabled(): boolean {
    return this.enabled;
  }

  /** Turn the panel on/off. Off = never shows (embedded hosts). Turning it on
   *  does not show anything by itself — the next renderable output does. */
  setEnabled(enabled: boolean): void {
    this.enabled = enabled;
    if (!enabled) this.clear();
  }

  /** Show the runtime values for one node. No-op when disabled or when the
   *  node has nothing renderable (DataFrame grid, column sample, graphics,
   *  widget/viewer root). First content ever → the panel expands (or appears
   *  minimized if the user minimized it earlier); afterwards content updates
   *  in place and the current state is respected. */
  showForNode(node: {id: string; label: string}, state: NodeExecState | undefined): void {
    if (!this.enabled || !state) return;
    // Status, duration, error, primitives — those go in the property panel
    // (see `buildExecutionMeta`). This panel only shows rich values.
    if (!hasRenderablePreview(state)) return;

    // Same node, same captured state, panel visible → the content is already
    // right; rebuilding would re-mount the grids and reset their scroll.
    if (this.state !== 'hidden' && node.id === this.lastNodeId && state === this.lastState) return;

    const inner = buildValuePreviews(state, (viewer) => this.onEditViewer?.(node, viewer));
    inner.style.padding = '8px 12px';
    this.contentEl.innerHTML = '';
    this.contentEl.appendChild(inner);
    this.nodeLabelEl.textContent = node.label;
    this.lastNodeId = node.id;
    this.lastState = state;

    if (this.state === 'hidden')
      this.setState(this.userMinimized ? 'minimized' : 'expanded');
  }

  /** Collapse to the header strip. Remembered: subsequent content never
   *  auto-expands — only an explicit `expand()` (header click) does. */
  minimize(): void {
    this.userMinimized = true;
    if (this.state === 'expanded') {
      // Keep the height the user (or the splitter divider) last gave it.
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

  /** Empty and hide the panel (graph change / new run — values are stale).
   *  The user's minimized preference survives. */
  clear(): void {
    this.contentEl.innerHTML = '';
    this.nodeLabelEl.textContent = '';
    this.lastNodeId = null;
    this.lastState = null;
    this.setState('hidden');
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
    // `ui.splitV`'s container-resize handler keeps rewriting `style.height` on
    // every pane; min/max clamp the rendered size so a minimized strip stays a
    // strip and an expanded pane can't be squeezed below its header.
    s.minHeight = `${HEADER_HEIGHT}px`;
    if (this.state === 'minimized') {
      s.height = `${HEADER_HEIGHT}px`;
      s.maxHeight = `${HEADER_HEIGHT}px`;
      this.contentEl.style.display = 'none';
      this.caretEl.textContent = '▴';
    } else {
      s.height = `${this.expandedHeight}px`;
      s.maxHeight = '';
      this.contentEl.style.display = '';
      this.caretEl.textContent = '▾';
    }
  }
}
