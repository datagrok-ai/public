/** Spotfire-style output tabs: one lazily-created detached TableView per table output, switched from
 *  the status bar. Runtime identity is the output node id; layouts persist by `paramName`. */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {setTid} from '../utils/test-ids';

export interface OutputTabInfo {
  nodeId: string;
  paramName: string;
}

export interface OutputTab {
  nodeId: string;
  paramName: string;
  /** Updated in place, never rebuilt — rebuilding the pressed chip would swallow its click. */
  chip: HTMLElement;
  chipLabel: HTMLElement;
  /** Display-toggled; never unmounted while the tab exists — the TableView DOM must persist. */
  pane: HTMLElement;
  emptyEl: HTMLElement;
  tv: DG.TableView | null;
  df: DG.DataFrame | null;
  stale: boolean;
  /** Applied once after the tv exists and the pane is visible (a hidden view defers forever). */
  pendingLayout: string | null;
  galleryConfigured: boolean;
  /** Captured ONCE, unwrapped — `setRibbonPanels` moves elements, so re-reading after a swap returns empty husks. */
  ribbonPanels?: HTMLElement[][];
  creating?: boolean;
}

export interface OutputViewsCallbacks {
  /** Fired on tab switches and when the active tab's tv appears; `null` = the Canvas tab. */
  onActiveTabChanged: (tab: OutputTab | null) => void;
  runFlow: () => void;
}

export const CANVAS_TAB = 'canvas';

export class OutputViewsManager {
  private readonly tabs = new Map<string, OutputTab>();
  private active: string = CANVAS_TAB;
  private readonly canvasChip: HTMLElement;
  private pendingByParam: Record<string, string> = {};
  private destroyed = false;

  constructor(
    private readonly contentHost: HTMLElement,
    private readonly canvasPane: HTMLElement,
    private readonly tabStripHost: HTMLElement,
    private readonly cb: OutputViewsCallbacks,
  ) {
    this.canvasChip = this.buildChipEl('Canvas', 'canvas');
    this.canvasChip.dataset.active = 'true';
    this.canvasChip.onclick = (): void => this.activate(CANVAS_TAB);
    ui.tooltip.bind(this.canvasChip, 'The flow canvas');
    this.tabStripHost.appendChild(this.canvasChip);
  }

  get activeKey(): string {return this.active;}

  getTab(nodeId: string): OutputTab | undefined {return this.tabs.get(nodeId);}
  getTabs(): OutputTab[] {return [...this.tabs.values()];}

  syncTabs(infos: OutputTabInfo[]): void {
    if (this.destroyed) return;
    const seen = new Set<string>();
    for (const info of infos) {
      seen.add(info.nodeId);
      let tab = this.tabs.get(info.nodeId);
      if (!tab) {
        tab = this.buildTab(info);
        this.tabs.set(info.nodeId, tab);
        this.tabStripHost.appendChild(tab.chip);
        this.contentHost.appendChild(tab.pane);
      }
      else if (tab.paramName !== info.paramName) {
        tab.paramName = info.paramName;
        if (tab.tv == null && tab.pendingLayout == null && this.pendingByParam[info.paramName])
          tab.pendingLayout = this.pendingByParam[info.paramName];
        this.refreshChip(tab);
      }
    }
    for (const [id, tab] of [...this.tabs]) {
      if (!seen.has(id)) this.destroyTab(id, tab);
    }
    // Only touch the DOM when the order differs — moving the pressed chip mid-click would swallow its click.
    const want = infos.map((i) => this.tabs.get(i.nodeId)?.chip).filter((c): c is HTMLElement => c != null);
    const have = Array.from(this.tabStripHost.children).filter((c) => c !== this.canvasChip);
    if (want.some((c, i) => c !== have[i]))
      for (const chip of want) this.tabStripHost.appendChild(chip);
  }

  setValue(nodeId: string, df: DG.DataFrame): void {
    const tab = this.tabs.get(nodeId);
    if (!tab) return;
    tab.df = df;
    tab.stale = false;
    if (tab.tv != null) {
      try {
        tab.tv.dataFrame = df;
        tab.tv.name = df.name || tab.paramName;
      } catch (e) {
        console.warn('Flow: output view refresh failed', e);
      }
    }
    else if (this.active === nodeId) {
      void this.ensureView(tab).then(() => {
        if (this.active === nodeId) this.cb.onActiveTabChanged(tab);
      });
    }
    this.refreshChip(tab);
  }

  markStale(nodeId: string): void {
    const tab = this.tabs.get(nodeId);
    if (!tab) return;
    tab.stale = true;
    this.refreshChip(tab);
  }

  clearValues(): void {
    for (const tab of this.tabs.values()) {
      if (tab.df != null) tab.stale = true;
      this.refreshChip(tab);
    }
  }

  activate(key: string): void {
    if (this.destroyed) return;
    if (key !== CANVAS_TAB && !this.tabs.has(key)) return;
    this.active = key;
    this.canvasChip.dataset.active = key === CANVAS_TAB ? 'true' : 'false';
    this.canvasPane.style.display = key === CANVAS_TAB ? '' : 'none';
    for (const [id, tab] of this.tabs) {
      tab.chip.dataset.active = id === key ? 'true' : 'false';
      tab.pane.style.display = id === key ? 'flex' : 'none';
    }
    if (key === CANVAS_TAB) {
      this.cb.onActiveTabChanged(null);
      return;
    }
    const tab = this.tabs.get(key)!;
    void this.ensureView(tab).then(() => {
      if (this.active === key) this.cb.onActiveTabChanged(tab);
    });
  }

  /** A live tv saves its real state; a never-opened tab keeps its loaded layout (must not lose it). */
  captureLayouts(): Record<string, string> {
    const res: Record<string, string> = {};
    for (const tab of this.tabs.values()) {
      if (tab.tv != null) {
        try {
          res[tab.paramName] = tab.tv.saveLayout().viewState;
          continue;
        } catch (e) {
          console.warn('Flow: layout capture failed for', tab.paramName, e);
        }
      }
      const kept = tab.pendingLayout ?? this.pendingByParam[tab.paramName];
      if (kept) res[tab.paramName] = kept;
    }
    return res;
  }

  setPendingLayouts(byParamName: Record<string, string>): void {
    this.pendingByParam = {...byParamName};
    for (const tab of this.tabs.values()) {
      const layout = byParamName[tab.paramName];
      if (layout) tab.pendingLayout = layout;
    }
  }

  destroy(): void {
    this.destroyed = true;
    for (const tab of this.tabs.values()) {
      try {
        tab.tv?.detach();
      } catch {/* already gone */}
      tab.pane.remove();
      tab.chip.remove();
    }
    this.tabs.clear();
    this.canvasChip.remove();
  }

  private buildTab(info: OutputTabInfo): OutputTab {
    const runBtn = ui.bigButton('Run', () => this.cb.runFlow());
    const emptyEl = ui.divV([
      ui.divText('This output has no value yet.', 'ff-output-view-empty-title'),
      ui.divText('Run the flow to see the table.', 'ff-output-view-empty-sub'),
      runBtn,
    ], 'ff-output-view-empty');
    setTid(emptyEl, 'output-view-empty');
    const pane = ui.div([emptyEl], 'ff-output-view-pane');
    setTid(pane, 'output-view-pane');
    pane.dataset.nodeId = info.nodeId;
    pane.style.display = 'none';

    const chip = this.buildChipEl(info.paramName, info.paramName);
    chip.dataset.nodeId = info.nodeId;
    chip.dataset.active = 'false';
    chip.onclick = (): void => this.activate(info.nodeId);

    const tab: OutputTab = {
      nodeId: info.nodeId,
      paramName: info.paramName,
      chip,
      chipLabel: chip.querySelector('.ff-view-tab-label') as HTMLElement,
      pane,
      emptyEl,
      tv: null,
      df: null,
      stale: false,
      pendingLayout: this.pendingByParam[info.paramName] ?? null,
      galleryConfigured: false,
    };
    chip.dataset.state = 'empty';
    ui.tooltip.bind(chip, () => this.chipTooltip(tab));
    return tab;
  }

  private buildChipEl(label: string, tidPart: string): HTMLElement {
    const dot = ui.div([], 'ff-view-tab-dot');
    const labelEl = ui.divText(label, 'ff-view-tab-label');
    const chip = ui.div([dot, labelEl], 'ff-view-tab');
    setTid(chip, 'view-tab', tidPart);
    chip.dataset.param = tidPart;
    return chip;
  }

  private chipTooltip(tab: OutputTab): string {
    if (tab.df == null)
      return `Flow output "${tab.paramName}" — not computed yet. Run the flow to see the table.`;
    const dims = `${tab.df.rowCount.toLocaleString()} × ${tab.df.columns.length}`;
    return tab.stale ?
      `Flow output "${tab.paramName}" — ${dims}. Out of date — run the flow to refresh.` :
      `Flow output "${tab.paramName}" — ${dims}`;
  }

  private refreshChip(tab: OutputTab): void {
    tab.chipLabel.textContent = tab.df?.name || tab.paramName;
    tab.chip.dataset.state = tab.df == null ? 'empty' : (tab.stale ? 'stale' : 'ready');
    setTid(tab.chip, 'view-tab', tab.paramName);
    tab.chip.dataset.param = tab.paramName;
  }

  private destroyTab(id: string, tab: OutputTab): void {
    this.tabs.delete(id);
    try {
      tab.tv?.detach();
    } catch {/* already gone */}
    tab.chip.remove();
    tab.pane.remove();
    if (this.active === id) this.activate(CANVAS_TAB);
  }

  /** Called only while the pane is visible and laid out — `_onAdded()` against a zero-size root corrupts docks. */
  private async ensureView(tab: OutputTab): Promise<void> {
    if (tab.df == null) return;
    if (tab.tv == null) {
      if (tab.creating) return;
      tab.creating = true;
      try {
        await this.createView(tab);
      } finally {
        tab.creating = false;
      }
      return;
    }
    (tab.tv as unknown as {_onAdded(): void})._onAdded(); // idempotent (MultiView pattern)
    this.applyPendingLayout(tab);
    this.nudgeResize(tab.tv);
  }

  private async createView(tab: OutputTab): Promise<void> {
    const df = tab.df!;
    try {
      await grok.data.detectSemanticTypes(df); // unregistered dfs skip auto-detection
    } catch (e) {
      console.warn('Flow: semantic type detection failed', e);
    }
    // The tab may have been destroyed while awaiting — mounting into its detached pane would leak the view.
    if (this.destroyed || tab.df !== df || this.tabs.get(tab.nodeId) !== tab) return;
    if (!df.name) df.name = tab.paramName;
    const tv = DG.TableView.create(df, false); // detached — no workspace pollution
    tv.name = df.name;
    tab.emptyEl.style.display = 'none';
    tv.root.style.width = '100%';
    tv.root.style.height = '100%';
    tab.pane.appendChild(tv.root);
    (tv as unknown as {_onAdded(): void})._onAdded();
    if (!tab.galleryConfigured) {
      try {
        await grok.functions.call('PowerPack:ConfigViewerGallery', {view: tv});
      } catch {/* PowerPack absent — the core gallery stays */}
      tab.galleryConfigured = true;
    }
    tab.tv = tv;
    // A newer value may have landed while the gallery call was in flight.
    if (tab.df !== df && tab.df != null) {
      try {
        tv.dataFrame = tab.df;
      } catch (e) {
        console.warn('Flow: output view refresh failed', e);
      }
    }
    this.applyPendingLayout(tab);
    this.refreshChip(tab);
  }

  private applyPendingLayout(tab: OutputTab): void {
    if (tab.tv == null || tab.pendingLayout == null) return;
    const state = tab.pendingLayout;
    tab.pendingLayout = null;
    try {
      tab.tv.loadLayout(DG.ViewLayout.fromViewState(state));
    } catch (e) {
      console.warn('Flow: saved layout could not be applied', e);
      grok.shell.warning(`Layout for output "${tab.paramName}" could not be applied`);
    }
  }

  /** Viewers added while hidden flush on resize; the js resize seam exists on DockView only. */
  private nudgeResize(tv: DG.TableView): void {
    const handle = (tv as unknown as {_handleResize?: () => void})._handleResize;
    if (typeof handle === 'function') {
      try {
        handle.call(tv);
        return;
      } catch {/* fall through */}
    }
    window.dispatchEvent(new Event('resize'));
  }
}
