/** Toolbox Suggestions pane — the bottom ~30% of the function browser.
 *
 *  Renders the ranked next steps from the suggestion engine
 *  ([suggestion-engine.ts](../suggest/suggestion-engine.ts)) and refreshes,
 *  debounced, whenever the host reports a context change (selection, graph
 *  edit, run completion, a preview-cell click). The header caret minimizes it
 *  to a slim strip; the choice persists in localStorage.
 *
 *  Items behave like toolbox nodes: DOUBLE-click hands the suggestion back to
 *  the host (`onAccept` — creates the node, prefills its inputs, wires the
 *  suggested connections), and they are HTML5-draggable onto the canvas — the
 *  full suggestion travels as JSON under {@link FF_SUGGEST_MIME} so a drop
 *  still gets the wiring and prefill, just at the drop point. */

import * as ui from 'datagrok-api/ui';

import {Suggestion} from '../suggest/suggestion-engine';
import {setTid} from '../utils/test-ids';

/** DataTransfer type for dragging a suggestion onto the canvas. The payload is
 *  the whole `Suggestion` as JSON (a plain-data object) — unlike the browser's
 *  `FF_DRAG_MIME`, which carries only a typeName. */
export const FF_SUGGEST_MIME = 'application/x-funcflow-suggestion';

const LS_KEY = 'funcflow-suggestions-collapsed';
const REFRESH_DEBOUNCE_MS = 250;

export class SuggestionPane {
  readonly root: HTMLElement;
  private readonly listEl: HTMLElement;
  private readonly caretEl: HTMLElement;
  private readonly countEl: HTMLElement;
  private collapsed = false;
  private timer: ReturnType<typeof setTimeout> | null = null;
  private refreshSeq = 0;
  /** Last rendered set — exposed for tests. */
  suggestions: Suggestion[] = [];

  constructor(
    private readonly provide: () => Promise<Suggestion[]>,
    private readonly onAccept: (s: Suggestion) => void,
  ) {
    try {
      this.collapsed = localStorage.getItem(LS_KEY) === 'true';
    } catch {/* storage blocked */}

    this.caretEl = setTid(ui.div([], 'ff-suggest-pane-caret'), 'suggest-pane-caret');
    this.caretEl.addEventListener('click', () => this.setCollapsed(!this.collapsed));
    ui.tooltip.bind(this.caretEl, () => this.collapsed ? 'Expand suggestions' : 'Minimize suggestions');

    const title = ui.div([], 'ff-suggest-pane-title');
    title.textContent = 'Suggestions';
    this.countEl = setTid(ui.div([], 'ff-suggest-pane-count'), 'suggest-pane-count');
    const header = setTid(
      ui.div([title, this.countEl, this.caretEl], 'ff-suggest-pane-header'), 'suggest-pane-header');
    header.addEventListener('click', (e) => {
      if (e.target === header || e.target === title) this.setCollapsed(!this.collapsed);
    });

    this.listEl = setTid(ui.div([], 'ff-suggest-pane-list'), 'suggest-pane-list');

    this.root = setTid(ui.div([header, this.listEl], 'ff-suggest-pane'), 'suggest-pane');
    this.applyCollapsed();
  }

  private setCollapsed(v: boolean): void {
    this.collapsed = v;
    try {
      localStorage.setItem(LS_KEY, String(v));
    } catch {/* storage blocked */}
    this.applyCollapsed();
    if (!v) this.refresh();
  }

  private applyCollapsed(): void {
    this.root.dataset.collapsed = String(this.collapsed);
    this.caretEl.textContent = this.collapsed ? '▸' : '▾';
  }

  /** Debounced recompute + render. Safe to call from any context signal. */
  refresh(): void {
    if (this.timer !== null) clearTimeout(this.timer);
    this.timer = setTimeout(() => {
      this.timer = null;
      void this.refreshNow();
    }, REFRESH_DEBOUNCE_MS);
  }

  /** Immediate recompute (tests / expand). Stale async results are dropped. */
  async refreshNow(): Promise<void> {
    if (this.collapsed) return;
    const seq = ++this.refreshSeq;
    let items: Suggestion[] = [];
    try {
      items = await this.provide();
    } catch {/* collection failed (mid-teardown) — render empty */}
    if (seq !== this.refreshSeq) return;
    this.suggestions = items;
    this.render();
  }

  private render(): void {
    this.listEl.innerHTML = '';
    this.countEl.textContent = this.suggestions.length ? String(this.suggestions.length) : '';
    if (this.suggestions.length === 0) {
      const empty = ui.div([], 'ff-suggest-pane-empty');
      empty.textContent = 'Select a node — or run the flow — to see what fits next.';
      this.listEl.appendChild(empty);
      return;
    }
    for (const s of this.suggestions) {
      const label = ui.div([], 'ff-suggest-pane-item-label');
      label.textContent = s.label;
      const reason = ui.div([], 'ff-suggest-pane-item-reason');
      reason.textContent = s.reason;
      const item = setTid(ui.div([label, reason], 'ff-suggest-pane-item'), 'suggest-pane-item');
      item.addEventListener('dblclick', () => this.onAccept(s));
      item.draggable = true;
      item.addEventListener('dragstart', (ev) => {
        if (!ev.dataTransfer) return;
        ev.dataTransfer.setData(FF_SUGGEST_MIME, JSON.stringify(s));
        ev.dataTransfer.setData('text/plain', s.label);
        ev.dataTransfer.effectAllowed = 'copy';
        item.classList.add('funcflow-func-item-dragging');
      });
      item.addEventListener('dragend', () => item.classList.remove('funcflow-func-item-dragging'));
      ui.tooltip.bind(item, () =>
        `${s.label} — ${s.reason}. Double-click to add${s.wire.length ? ' and connect' : ''}, or drag onto the canvas.`);
      this.listEl.appendChild(item);
    }
  }

  destroy(): void {
    if (this.timer !== null) clearTimeout(this.timer);
  }
}
