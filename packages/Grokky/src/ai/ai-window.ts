import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {AIPanel} from './panel';
import {getAIPanelToggleSubscription, isToggleKey} from '../utils';

type AnyPanel = AIPanel<any, any>;

const near = (a: number, b: number, tolerance = 2): boolean => Math.abs(a - b) <= tolerance;

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
  private props: HTMLElement | null = null;
  private propsShown = false;
  private stacked = false;
  private moved: DG.DockContainer | null = null;

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

    const follow = () => Promise.resolve().then(() => this.followProperties());
    follow();
    grok.shell.windows.onPanelVisibilityChanged.subscribe(follow);
    grok.shell.dockManager.onClosed.subscribe(follow);

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

  private followProperties(): void {
    const {windows, dockManager: dm} = grok.shell;
    const ai = dm.findNode(windows.ai);
    if (!windows.showAI) {
      if (ai)
        dm.close(ai);
      this.moved?.destroy();
      this.moved = null;
    }
    const opened = windows.showContextPanel && !this.propsShown;
    this.propsShown = windows.showContextPanel;
    const propsNode = this.propsShown ? dm.findNode(windows.context.root!) : undefined;
    const previous = this.props;
    this.props = propsNode?.container.containerElement ?? null;
    if (!ai || !windows.showAI)
      return;
    const a = ai.container.containerElement.getBoundingClientRect();
    if (!a.width)
      return;
    if (this.props) {
      const p = this.props.getBoundingClientRect();
      const beside = opened && near(p.top, a.top) && near(p.bottom, a.bottom) &&
        (near(p.left, a.right, 12) || near(a.left, p.right, 12));
      if (beside)
        this.moved = dm.dock(windows.ai, DG.DOCK_TYPE.DOWN, propsNode!, 'AI', 0.3).container;
      this.stacked = beside || (near(p.left, a.left) && near(p.right, a.right));
    } else if (previous && this.stacked) {
      this.stacked = false;
      const width = parseFloat(previous.style.width);
      const area = dm.element.getBoundingClientRect();
      if (width && !near(a.width, width) && near(a.right, area.right) && near(a.height, area.height))
        this.moved = dm.dock(windows.ai, DG.DOCK_TYPE.RIGHT, null, 'AI', width / area.width).container;
    }
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
