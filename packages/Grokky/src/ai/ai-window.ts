import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {AIPanel} from './panel';
import {getAIPanelToggleSubscription, isToggleKey} from '../utils';

type AnyPanel = AIPanel<any, any>;

export class AIWindowManager {
  private static _instance: AIWindowManager | null = null;
  static get instance(): AIWindowManager {
    if (!AIWindowManager._instance)
      AIWindowManager._instance = new AIWindowManager();
    return AIWindowManager._instance;
  }

  private readonly byView = new WeakMap<DG.ViewBase, AnyPanel>();
  /** Panels owned by a single view (e.g. DB query assistants) — disposed when that view closes. */
  private readonly owned = new WeakSet<AnyPanel>();
  private shell: AnyPanel | null = null;
  private current: AnyPanel | null = null;
  private initialized = false;

  init(shell: AnyPanel): void {
    this.shell = shell;
    if (this.initialized)
      return;
    this.initialized = true;

    // The mounted panel deliberately survives view switches: one panel, one session.
    // It learns the current view at prompt time (buildWorkspaceContext), not by remounting.

    grok.events.onViewRemoved.subscribe((v) => this.unregister(v));

    // Keep the empty-state suggestion cards in sync with the current view.
    grok.events.onCurrentViewChanged.subscribe(() => {
      if (grok.shell.windows.showAI)
        this.current?.activate(false);
    });

    getAIPanelToggleSubscription().subscribe((v) => this.toggle((v as DG.ViewBase) ?? grok.shell.v));

    document.addEventListener('keydown', (e: KeyboardEvent) => {
      if (isToggleKey(e)) {
        e.preventDefault();
        this.toggle(grok.shell.v);
      }
    });
  }

  register(view: DG.ViewBase, panel: AnyPanel, opts: {owned?: boolean} = {}): void {
    this.byView.set(view, panel);
    if (opts.owned)
      this.owned.add(panel);
  }

  show(view: DG.ViewBase | null = grok.shell.v): void {
    this.open(this.panelFor(view), true);
  }

  showPanel(panel: AnyPanel, focus = true): void {
    this.open(panel, focus);
  }

  toggle(view: DG.ViewBase | null): void {
    const panel = this.panelFor(view);
    const showing = grok.shell.windows.showAI && grok.shell.windows.ai.contains(panel.root);
    showing ? this.close() : this.open(panel, true);
  }

  private panelFor(view: DG.ViewBase | null): AnyPanel {
    return (view && this.byView.get(view)) || this.current || this.shell!;
  }

  private open(panel: AnyPanel, focus: boolean): void {
    const previouslyFocused = document.activeElement as HTMLElement | null;
    this.mount(panel);
    panel.activate(focus);
    if (!focus)
      previouslyFocused?.focus();
  }

  private close(): void {
    grok.shell.windows.showAI = false;
  }

  private mount(panel: AnyPanel): void {
    const win = grok.shell.windows.ai;
    if (this.current && this.current !== panel && win.contains(this.current.root))
      win.removeChild(this.current.root);
    if (!win.contains(panel.root))
      win.appendChild(panel.root);
    grok.shell.windows.showAI = true;
    this.current = panel;
  }

  private unregister(view: DG.ViewBase): void {
    const panel = this.byView.get(view);
    if (!panel)
      return;
    this.byView.delete(view);
    if (!this.owned.has(panel))
      return; // shared panel (shell/scripting singleton) — keep it alive
    if (this.current === panel) {
      this.current = null;
      if (grok.shell.windows.showAI && this.shell)
        this.open(this.shell, false);
    }
    panel.dispose();
  }
}
