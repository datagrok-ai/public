/** Main FuncFlow view — Datagrok ViewBase that hosts the Rete editor.
 *
 * Layout: [function-browser]  |  [Rete canvas container]
 *                                  + status bar at the bottom
 *                                  (property panel goes into the native context panel)
 *
 * The canvas container is a plain `<div>`; the Rete `AreaPlugin` builds the
 * scrollable/zoomable inner content inside it. Drag-and-drop from the
 * Datagrok browse tree is bound to this div via `ui.makeDroppable`. */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// @ts-ignore
import '../css/funcflow.css';

import {FlowEditor} from './rete/flow-editor';
import {FlowNode, isSetVarNode} from './rete/scheme';
import {FunctionBrowser, FF_DRAG_MIME} from './panel/function-browser';
import {PropertyPanel} from './panel/property-panel';
import {ColumnPicker} from './panel/column-picker';
import {FuncEditorLauncher} from './panel/func-editor-launcher';
import {
  registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs,
  createNode, FuncInfo,
} from './rete/node-factory';
import {validateGraph} from './compiler/validator';
import {emitScript} from './compiler/script-emitter';
import {emitCreationScript, emitCreationScriptsForTables} from './compiler/creation-script-emitter';
import {
  serializeFlow, deserializeFlow, downloadFlow, loadFlowFromFile,
} from './serialization/flow-serializer';
import {flowScriptText, parseFlowBody, FLOW_TAG, FLOW_LANGUAGE} from './serialization/flow-script-format';
import {SpacePicker} from './ui/space-picker';
import {buildFlowFromCreationScript} from './import/creation-script-importer';
import {FlowSettings, FuncFlowDocument} from './serialization/flow-schema';
import {ExecutionController} from './execution/execution-controller';
import {AutorunScheduler, AUTORUN_DEBOUNCE_MS, isAutorunByDefault} from './execution/autorun';
import {OutputPreviewPanel, OutputPanelState} from './execution/output-preview';
import {ValueSummary, NodeExecStatus} from './execution/execution-state';
import {OutputViewsManager, OutputTab, OutputTabInfo} from './views/output-views-manager';
import {buildPreview, setPreviewCellFocusHandler} from './execution/value-inspector';
import {SuggestionPane, FF_SUGGEST_MIME} from './panel/suggestion-pane';
import {
  collectSuggestContext, computeSuggestions, Suggestion, CellSignal,
} from './suggest/suggestion-engine';
import {_package} from './package';
import {setTid} from './utils/test-ids';
import {
  addPendingFile, getPendingFile, isPendingFileId, persistPendingFile, syncFlowFilePermissions,
} from './utils/uploaded-files';
import {GuideHost} from './guide/guide-model';
import {GuideRunner} from './guide/guide-runner';
import {createHelpButton, openGuideMenu} from './guide/guide-launcher';
import {TUTORIALS} from './guide/guide-content';
import {summarizeFlow} from './summary/summary-generator';
import {FlowAIContext} from './ai-tools';

/** Bundled starter flows (files in `files/`), surfaced on the Start panel so a
 *  scientist never faces a blank canvas. */
const FLOW_TEMPLATES: {label: string; file: string; desc: string}[] = [
  {label: 'Workflow demo', file: 'Workflow Demo.ffjson', desc: 'A sample multi-step data workflow.'},
  {label: 'Bio Molecules', file: 'Sequence demo.ffjson', desc: 'A Peptides conversion and calculation.'},
];

export class FuncFlowView extends DG.ViewBase {
  private flow!: FlowEditor;
  private functionBrowser!: FunctionBrowser;
  private propertyPanel!: PropertyPanel;
  private executionController!: ExecutionController;
  /** Toolbox Suggestions pane (bottom 25% of the browser). Public for tests. */
  suggestionPane!: SuggestionPane;
  /** The node the user last clicked — a suggestion-context signal that
   *  outlives a deselect ("the currently clicked node"). */
  private focusNodeId: string | null = null;
  /** The cell last clicked in the output preview (semType + value). */
  private previewCell: CellSignal | null = null;
  /** Bottom output panel — a pane of the view's vertical splitter (see
   *  `initUI`). Public so hosts (the creation-script dialog) and tests can
   *  reach it; created disabled when `options.outputPanel` is false. */
  outputPreview!: OutputPreviewPanel;
  /** Whether this view gets the bottom output panel (false for embedded
   *  hosts — the creation-script dialog). */
  private readonly outputPanelEnabled: boolean;
  /** Debounced rerun-on-change; toggled by the ribbon bolt icon, off by
   *  default. Created with the execution controller in `initEditor`. */
  private autorunScheduler: AutorunScheduler | null = null;
  private autorunIcon: HTMLElement | null = null;
  /** The ribbon Save button (non-creation mode only), greyed when there is
   *  nothing to save — an empty canvas or no change since the last save. */
  private saveButton: HTMLButtonElement | null = null;
  /** Serialized graph+settings at the last save / load, for detecting unsaved
   *  changes; null before the first baseline is recorded. */
  private savedSnapshot: string | null = null;
  /** Platform event subscriptions (file-drag overlay suppression, entity-shared
   *  permission sync) — disposed in `detach`. */
  private platformSubs: {unsubscribe(): void}[] = [];
  /** One-shot handler that pins this (preview) view on the first interaction so
   *  the toolbox appears; cleared once it fires (see `setupAutoPin`). */
  private autoPinHandler: (() => void) | null = null;
  private canvasContainer!: HTMLElement;
  private startPanel!: HTMLElement;
  private startBg!: HTMLElement;
  private startBgRaf = 0;
  private recentFlowsHost!: HTMLElement;
  private recentFlowsLoading = false;
  private helpButton!: HTMLElement;
  private readonly guideRunner = new GuideRunner();
  private statusBar!: HTMLElement;
  private nodeCountLabel!: HTMLElement;
  private linkCountLabel!: HTMLElement;
  private validationLabel!: HTMLElement;
  /** Internal output-view tabs (Canvas + one per table output), switched from
   *  the status bar. Public for tests. */
  outputViews!: OutputViewsManager;
  private viewTabStrip!: HTMLElement;
  /** Flow's own ribbon panels, kept by reference so they can be restored when
   *  switching back to the Canvas tab (`setRibbonPanels` MOVES elements). */
  private flowRibbonPanels: HTMLElement[][] = [];

  /** Desired initial minimap state, applied once the editor is created (the
   *  editor is built async). Hosts set this before the editor exists — e.g. the
   *  creation-script preview dialog opens with the minimap minimized. */
  private minimapCollapsed = false;

  /** When the view edits the creation scripts of specific tables (the
   *  `creationScriptEditor` entry point), these are those tables — keyed for
   *  the per-table split by their `.VariableName` tag. Presence adds a primary
   *  **Save** action to the ribbon that writes a creation script back to each. */
  private readonly tableInfos: DG.TableInfo[];

  private flowSettings: FlowSettings = {
    scriptName: 'MyFuncFlow',
    scriptDescription: 'Generated by FuncFlow',
    tags: ['funcflow'],
  };

  protected override afterPersist(): void {
    grok.shell.windows.showToolbox = false;
    grok.shell.windows.showBrowse = false;
    grok.shell.windows.showContextPanel = true;
    grok.shell.windows.showHelp = false;
    grok.shell.windows.showBrowse = true;
    grok.shell.windows.showToolbox = true;
  }

  /** The platform Script entity (language 'flow') this view edits, when opened
   *  from / saved to the server. Null for scratch flows — Save asks for a name. */
  private boundScript: DG.Script | null = null;

  /** The dashboard project this flow publishes into, persisted in the .ffjson:
   *  re-publishing updates that project instead of creating a new one per
   *  save. Cleared on New / Save As / the dialog's "publish as new" link. */
  private dashboardProjectId: string | null = null;

  /** Resolves once `initEditor()` has built the Rete editor (deferred a tick
   *  off the constructor). Load paths await this instead of retrying on
   *  timers, so opening an entity cannot race editor construction. */
  private editorReadyResolve!: () => void;
  private readonly editorReady = new Promise<void>((resolve) => {
    this.editorReadyResolve = resolve;
  });

  /** Watches the canvas for its first real layout when a fit-to-screen was
   *  requested before the view was attached (see `fitToScreen`). */
  private pendingFitObserver: ResizeObserver | null = null;

  /** @param tableInfos tables whose creation scripts this view edits — passing a
   *  non-empty array enables the **Save Creation Scripts** ribbon action.
   *  @param options `outputPanel: false` creates the view without the bottom
   *  output panel — for embedded hosts (the creation-script dialog), where a
   *  run-results pane makes no sense. `enableOutputPanel()` turns it on later
   *  (e.g. when that view is promoted into the real editor). */
  constructor(tableInfos: DG.TableInfo[] = [], options: {outputPanel?: boolean} = {}) {
    super();
    this.name = 'FuncFlow';
    this.aiDescription = 'Flow — Datagrok\'s visual pipeline editor: the user composes functions, ' +
      'queries, and scripts into an executable graph on a canvas. Act on it through the view ' +
      'functions (search list_view_functions with the word "flow"): listFlowNodes (call first to ' +
      'see the canvas), getFlowNodeDetails, findFlowNodeTypes, addFlowNode, connectFlowNodes, ' +
      'setFlowNodeInputs, selectFlowNode, runFlow; listFlowGuides / startFlowGuide launch ' +
      'interactive tutorials — offer one when the user asks how to do something here. When user asks about how to do something, ' +
      'Make sure to list guides (run the function listing the guides) to see if there are guides that cover the question, and if so and confidence is high, call that guide function';
    this.tableInfos = tableInfos;
    this.outputPanelEnabled = options.outputPanel !== false;

    registerBuiltinNodes();

    this.initUI();
    this.setupRibbon();
    this.setupStatusBar();

    // Browser lives in the platform toolbox window now — make sure both panels are visible.
    this.toolbox = this.functionBrowser.root;
    grok.shell.windows.showToolbox = true;
    grok.shell.windows.showContextPanel = true;

    setTimeout(() => {
      try {
        registerAllFunctions();
        this.functionBrowser.render();
      } catch (e) {
        console.warn('FuncFlow: error registering functions:', e);
      }
    }, 100);

    // Local-file drops: while this view is active, suppress the platform's
    // "drop to open" overlay so the canvas underneath receives the OS drag;
    // the drop handler in setupFileDropTarget turns each file into an
    // Uploaded File node. (Share-driven permission sync is global — the
    // package's flowShareSync autostart — so it works with no view open.)
    this.platformSubs.push(grok.events.onFileDragEnter.subscribe((ev) => {
      const target = ev.causedBy?.target;
      if (target !== null && target instanceof HTMLElement &&  (target === this.canvasContainer || this.canvasContainer.contains(target)))
        ev.preventDefault();
    }));
  }

  /** A view bound to a platform flow-script entity: loads its body and makes
   *  Save write back to the server. */
  static forScript(script: DG.Script): FuncFlowView {
    const view = new FuncFlowView();
    view.bindScript(script);
    return view;
  }

  private bindScript(script: DG.Script): void {
    this.boundScript = script;
    this.name = script.friendlyName || script.name || 'Flow';
    this.updatePath();
    try {
      const {doc} = parseFlowBody(script.script);
      // The entity is already on the server → its loaded state is the saved
      // baseline, so Save stays disabled until the user edits.
      void this.loadFromDoc(doc).then(() => this.markSaved());
    } catch (e) {
      grok.shell.error(`Cannot read flow "${this.name}": ${e instanceof Error ? e.message : e}`);
    }
  }

  /** Keeps the browser URL on the entity route (ScriptView's '/script/<id>'
   *  shape, which routes back into the visual editor). The Dart host may not
   *  be attached yet — core sets the same path at open time in that case. */
  private updatePath(): void {
    const id = this.boundScript?.id;
    if (!id) return;
    try {
      this.path = `/script/${id}`;
    } catch {/* view not yet attached to the platform shell */}
  }

  private initUI(): void {
    this.functionBrowser = new FunctionBrowser({
      onFunctionDoubleClick: (info: FuncInfo) => void this.addNodeByType(info.nodeTypeName),
      onBuiltinNodeDoubleClick: (typeName: string) => void this.addNodeByType(typeName),
      onFileDoubleClick: (file: DG.FileInfo) => void this.addOpenFileNode(file.fullPath),
      onLocalFilesPicked: (files: File[]) => void this.addUploadedFileNodes(files),
    });

    // Suggestions pane — the bottom 25% of the toolbox. Recomputes from the
    // full canvas context (selection, focus node, captured data, preview cell).
    this.suggestionPane = new SuggestionPane(
      async () => this.flow ?
        computeSuggestions(await collectSuggestContext(
          this.flow, this.executionController ?? null, this.focusNodeId, this.previewCell)) :
        [],
      (s) => void this.applySuggestion(s),
    );
    this.functionBrowser.root.appendChild(this.suggestionPane.root);
    // A cell clicked in the output preview is a context signal ("clicked a
    // Molecule value") — remember it and recompute.
    setPreviewCellFocusHandler((cell) => {
      this.previewCell = cell;
      this.suggestionPane.refresh();
    });

    // The canvas is a direct `.ui-box` child (the splitter pane), where core
    // css forces `overflow: auto !important` on EVERY child — `ui-div`-classed
    // ones via `div.ui-box > div.ui-div`, all others via a huge
    // `div.ui-box:not(…) > *:not(.ui-div):not(…)` rule that outranks anything
    // reasonable. Keep `ui.div` (the exempt, low-specificity branch) and win
    // it back with `div.ui-box > div.ui-div.funcflow-canvas-container
    // { overflow: hidden !important }` in funcflow.css — otherwise the
    // transformed Rete canvas grows giant scrollbars.
    this.canvasContainer = ui.div([], 'funcflow-canvas-container');
    setTid(this.canvasContainer, 'canvas');
    this.startPanel = this.buildStartPanel();
    this.canvasContainer.appendChild(this.startPanel);

    // Non-invasive floating help button (bottom-left) → tutorials & how-to menu.
    this.helpButton = createHelpButton((ev) => openGuideMenu(this.guideHost, this.guideRunner, ev));
    this.canvasContainer.appendChild(this.helpButton);

    this.nodeCountLabel = ui.divText('Nodes: 0');
    this.linkCountLabel = ui.divText('Links: 0');
    this.validationLabel = ui.divText('');
    setTid(this.nodeCountLabel, 'statusbar-nodes');
    setTid(this.linkCountLabel, 'statusbar-links');
    setTid(this.validationLabel, 'statusbar-validation');
    ui.tooltip.bind(this.nodeCountLabel, 'Total number of nodes in the graph');
    ui.tooltip.bind(this.linkCountLabel, 'Total number of connections between nodes');
    // The status bar leads with the view tabs (Canvas + one per table output,
    // Spotfire-style bottom pages); the counters sit right-aligned after them.
    this.viewTabStrip = ui.div([], 'ff-view-tabs');
    setTid(this.viewTabStrip, 'view-tabs');
    const statusLabels = ui.div(
      [this.nodeCountLabel, this.linkCountLabel, this.validationLabel], 'ff-status-labels');
    this.statusBar = ui.div([this.viewTabStrip, statusLabels], 'funcflow-status-bar');
    setTid(this.statusBar, 'statusbar');

    // Bottom output panel — a real pane of the view, not a dock-manager dock:
    // canvas and panel share a vertical splitter (`ui.splitV`), so the panel is
    // resizable via the divider, minimizes to its header strip, and can never
    // linger over other views. Created disabled for embedded hosts.
    this.outputPreview = new OutputPreviewPanel({enabled: this.outputPanelEnabled});
    const canvasBox = ui.box(this.canvasContainer);
    // The canvas absorbs all space the panel doesn't take — the panel's
    // explicit height (its pane is `flex: 0 0 auto`) is the single source of
    // truth, in every panel state and under `splitV`'s own resize handling.
    canvasBox.style.flex = '1 1 0';
    canvasBox.style.minHeight = '0';
    const split = ui.splitV([canvasBox, this.outputPreview.root], {style: {flex: '1 1 0', width: '100%', height: '100%'}}, true);
    // The divider only makes sense when there is something to resize.
    const divider = split.querySelector('.ui-split-v-divider') as HTMLElement | null;
    const syncDivider = (state: OutputPanelState): void => {
      if (divider) divider.style.display = state === 'expanded' ? '' : 'none';
    };
    // Opening the output preview minimizes the overview minimap so the two
    // don't crowd the same corner — fired only on the hidden → visible edge
    // (re-expanding the preview later, or the user reopening the minimap, is
    // left alone). This lives here rather than in OutputPreviewPanel because
    // the minimap belongs to the editor, not the preview.
    let lastPanelState: OutputPanelState = this.outputPreview.panelState;
    this.outputPreview.onStateChanged = (state) => {
      syncDivider(state);
      if (lastPanelState === 'hidden' && state !== 'hidden')
        this.flow?.setMinimapCollapsed(true);
      lastPanelState = state;
    };
    syncDivider(this.outputPreview.panelState);

    const mainLayout = setTid(ui.div([split], 'funcflow-root'), 'root');
    this.root.style.cssText = 'width:100%;height:100%;display:flex;flex-direction:column;';
    setTid(this.root, 'view');
    // The content host stacks the canvas layout and one pane per output-view
    // tab; panes are display-toggled (never unmounted — the hosted TableView
    // DOM must persist across tab switches).
    const contentHost = ui.div([mainLayout], 'ff-view-content');
    setTid(contentHost, 'view-content');
    this.root.appendChild(contentHost);
    this.root.appendChild(this.statusBar);
    mainLayout.style.flex = '1';
    mainLayout.style.overflow = 'hidden';

    this.outputViews = new OutputViewsManager(contentHost, mainLayout, this.viewTabStrip, {
      onActiveTabChanged: (tab) => this.onOutputTabChanged(tab),
      runFlow: () => this.runInstrumented(),
    });

    this.setupFileDropTarget();
    this.installPortContextMenu();
    this.setupAutoPin();

    setTimeout(() => this.initEditor(), 50);
  }

  // ---------- per-port "View Output" preview (KNIME pattern #2) ----------

  /** Right-click on an output socket → menu with "View output" if a runtime
   *  value is captured for that port. The handler runs in capture phase so
   *  it intercepts the contextmenu before AreaPlugin emits its `'contextmenu'`
   *  signal (which would route to the node menu). */
  private installPortContextMenu(): void {
    this.canvasContainer.addEventListener('contextmenu', (ev) => {
      const target = ev.target as HTMLElement | null;
      // Suppress root menu on input sockets — there's no output value to
      // show, and without this we'd otherwise leak the "Add annotation
      // here" item from the area-plugin's container-level emission.
      if (target?.closest('.ff-socket-row-input')) {
        ev.preventDefault();
        ev.stopPropagation();
        return;
      }
      // Match anywhere in an output row — dot, wrapper, or label — so the
      // user doesn't have to right-click the 12px socket dot precisely.
      const rowEl = target?.closest('.ff-socket-row-output') as HTMLElement | null;
      if (!rowEl) return;
      const nodeEl = rowEl.closest('.ff-node') as HTMLElement | null;
      if (!nodeEl) return;
      const nodeId = nodeEl.dataset.nodeId;
      if (!nodeId) return;
      const node = this.flow?.getNodeById(nodeId);
      if (!node) return;
      // Output key by row index — same convention as the suggestion menu.
      const rows = Array.from(nodeEl.querySelectorAll('.ff-socket-row-output')) as HTMLElement[];
      const idx = rows.indexOf(rowEl);
      const outputKey = Object.keys(node.outputs)[idx];
      if (!outputKey) return;

      ev.preventDefault();
      ev.stopPropagation();

      const state = this.executionController?.state.getNodeState(nodeId);
      const summary = state?.outputs?.[outputKey];

      const menu = DG.Menu.popup();
      if (summary && this.isPortPreviewable(summary))
        menu.item('View output', () => this.showPortPreview(rowEl, outputKey, summary));
      // Inspect-anywhere: run just the slice up to this node and preview it —
      // no full run, no Output node required.
      menu.item(summary ? 'Re-run up to here' : 'Run up to here & preview',
        () => this.previewNodeData(nodeId));
      menu.show({causedBy: ev});
    }, true);
  }

  /** Live subscription capturing a viewer's option changes into its node. */
  private viewerEditSub: {unsubscribe(): void} | undefined;

  /** Show a live viewer in the context panel and persist its option changes
   *  onto the node. `grok.shell.o = viewer` makes Datagrok render the viewer's
   *  full settings editor; we debounce `onPropertyValueChanged`, read
   *  `getOptions().look` (dropping the `#type` tag), and store it as the node's
   *  `viewerLook` so a re-run reproduces the exact look. */
  private editViewer(nodeId: string, viewer: unknown): void {
    const v = viewer as {
      root?: HTMLElement;
      getOptions?: () => {look?: unknown};
      onPropertyValueChanged?: unknown;
    };
    grok.shell.o = v as object;

    this.viewerEditSub?.unsubscribe();
    this.viewerEditSub = undefined;

    const capture = (): void => {
      try {
        const node = this.flow?.getNodeById(nodeId);
        if (!node) return;
        let look = v.getOptions?.().look as Record<string, unknown> | string | undefined;
        if (typeof look === 'string') look = JSON.parse(look) as Record<string, unknown>;
        if (!look || typeof look !== 'object') return;
        const clean: Record<string, unknown> = {};
        for (const [k, val] of Object.entries(look))
          if (k !== '#type') clean[k] = val;
        node.properties['viewerLook'] = clean;
        void this.flow.updateNode(node.id);
      } catch (e) {
        console.warn('FuncFlow: failed to capture viewer options', e);
      }
    };

    try {
      const obs = v.onPropertyValueChanged as {subscribe(cb: () => void): {unsubscribe(): void}} | undefined;
      if (obs)
        this.viewerEditSub = DG.debounce(obs as any, 300).subscribe(() => capture());
    } catch (e) {
      console.warn('FuncFlow: viewer onPropertyValueChanged unavailable', e);
    }
  }

  /** Run only the slice up to this node and open its data preview — the
   *  "inspect anywhere" entry point from an output-port right-click. */
  private previewNodeData(nodeId: string): void {
    this.executionController?.previewNodeData(nodeId, {
      name: this.flowSettings.scriptName,
      description: this.flowSettings.scriptDescription,
      tags: this.flowSettings.tags,
    });
  }

  /** Re-run just this node using upstream values captured from a prior run. */
  private rerunNode(nodeId: string): void {
    this.executionController?.rerunNode(nodeId, {
      name: this.flowSettings.scriptName,
      description: this.flowSettings.scriptDescription,
      tags: this.flowSettings.tags,
    });
  }

  private isPortPreviewable(summary: ValueSummary): boolean {
    if (summary.type === 'dataframe' && summary.clone) return true;
    if (summary.type === 'column' && Array.isArray(summary.sample) && summary.sample.length > 0) return true;
    if (summary.type === 'graphics' && typeof summary.value === 'string') return true;
    if ((summary.type === 'widget' || summary.type === 'viewer') && summary.value?.root instanceof Element) return true;
    return false;
  }

  private currentPortPopup: HTMLElement | null = null;

  /** Floating popover anchored next to a port's screen position. Reuses the
   *  same `buildPreview` renderer the docked panel uses. Auto-closes on
   *  Esc or click-outside. */
  private showPortPreview(anchorEl: HTMLElement, name: string, summary: ValueSummary): void {
    if (this.currentPortPopup) this.currentPortPopup.remove();
    const inner = buildPreview(name, summary);
    if (!inner) return;

    const popup = setTid(ui.div([inner], 'ff-port-preview'), 'port-preview');
    document.body.appendChild(popup);
    this.currentPortPopup = popup;

    const rect = anchorEl.getBoundingClientRect();
    const popupRect = popup.getBoundingClientRect();
    const vw = window.innerWidth; const vh = window.innerHeight;
    let left = rect.right + 12;
    if (left + popupRect.width > vw - 8) left = rect.left - popupRect.width - 12;
    if (left < 8) left = 8;
    let top = rect.top;
    if (top + popupRect.height > vh - 8) top = vh - popupRect.height - 8;
    if (top < 8) top = 8;
    popup.style.left = `${left}px`;
    popup.style.top = `${top}px`;

    const close = (): void => {
      popup.remove();
      this.currentPortPopup = null;
      document.removeEventListener('mousedown', onDoc, true);
      document.removeEventListener('keydown', onKey, true);
    };
    const onDoc = (ev: MouseEvent): void => {
      if (!popup.contains(ev.target as Node)) close();
    };
    const onKey = (ev: KeyboardEvent): void => {
      if (ev.key === 'Escape') close();
    };
    document.addEventListener('mousedown', onDoc, true);
    document.addEventListener('keydown', onKey, true);
  }

  private initEditor(): void {
    this.flow = new FlowEditor(this.canvasContainer, {
      onNodeSelected: (node: FlowNode) => {
        const execState = this.executionController?.state.getNodeState(node.id);
        this.propertyPanel.showNode(node, execState);
        grok.shell.o = this.propertyPanel.root;
        // Lazy: opens (or updates) the bottom-docked output panel only if
        // this node has captured runtime values from a prior run.
        this.executionController?.showOutputsForNode(node);
        // Clicking a node re-anchors the suggestion context; a stale preview
        // cell from another node's data no longer applies.
        this.focusNodeId = node.id;
        this.previewCell = null;
        this.suggestionPane?.refresh();
      },
      onNodeDeselected: () => {
        this.propertyPanel.clear();
        grok.shell.o = this.propertyPanel.root;
        this.suggestionPane?.refresh();
      },
      // Selection changes that bypass `nodepicked` — the marquee (its release
      // is swallowed before it can bubble to any container listener), Ctrl+A,
      // modifier-click removals, programmatic selects.
      onSelectionChanged: () => this.suggestionPane?.refresh(),
      onGraphChanged: () => {
        this.updateStatusBar();
        this.updateStartPanelVisibility();
        this.refreshNodeHints();
        this.flow?.refreshMinimap();
        this.updateSaveButtonState();
        this.suggestionPane?.refresh();
        this.outputViews?.syncTabs(this.tableOutputs());
      },
      // Classified edits: invalidate exactly the affected downstream cone, and
      // feed the same set to the autorun scheduler (a rerun of just that slice).
      onGraphEdited: (edit) => {
        const affected = this.executionController?.applyGraphEdit(edit) ?? new Set<string>();
        this.autorunScheduler?.onEdit(edit, affected);
        // `onGraphChanged` does NOT fire on parameter edits (only structural
        // add/remove/clear + annotations) — `notifyNodeParamsChanged` reports
        // via `onGraphEdited` only. So refresh Save state here too, else editing
        // a value on a saved flow would leave Save greyed.
        this.updateSaveButtonState();
        this.suggestionPane?.refresh();
        // paramName renames arrive as params-changed (no onGraphChanged) —
        // the tab set / labels must track them here too.
        this.outputViews?.syncTabs(this.tableOutputs());
      },
      onPreviewNode: (nodeId: string) => this.previewNodeData(nodeId),
      onRerunNode: (nodeId: string) => this.rerunNode(nodeId),
      canRerunNode: (nodeId: string) => this.executionController?.canRerunNode(nodeId) ?? false,
    });

    this.propertyPanel = new PropertyPanel(this.flow);
    this.executionController = new ExecutionController(this.flow, this.outputPreview);
    // Fresh captured values (node completed / invalidated) change what fits
    // next — recompute. Debounced in the pane, so per-node bursts are cheap.
    // Also fans out to the output-view tabs: a completed table output (or
    // SetVar terminal) refreshes its tab; an invalidated one goes stale.
    this.executionController.onNodeStateChanged = (nodeId) => {
      this.suggestionPane?.refresh();
      this.updateOutputViewValue(nodeId);
    };
    // Empty-canvas-click deselects happen inside rete with no callback — any
    // release on the canvas re-reads the selection. (Marquee releases never
    // get here — they are swallowed in `installRectSelect`; those report via
    // `onSelectionChanged` instead.)
    this.canvasContainer.addEventListener('pointerup', () => this.suggestionPane?.refresh());
    this.autorunScheduler = new AutorunScheduler((dirty, liveOnly) => {
      const settings = {
        name: this.flowSettings.scriptName,
        description: this.flowSettings.scriptDescription,
        tags: this.flowSettings.tags,
      };
      // Toggle off → the set holds only live-by-default nodes; run exactly
      // those (input-satisfied ones), never the rest of the canvas.
      return (liveOnly ?
        this.executionController?.runLiveNodes(dirty, settings) :
        this.executionController?.runAutorun(dirty, settings)) ?? 'skipped';
    },
    AUTORUN_DEBOUNCE_MS,
    // Live-by-default nodes (Open File, Add New Column, viewers, …) rerun on
    // change even with the ribbon toggle off.
    (id) => {
      const n = this.flow?.getNodeById(id);
      return n != null && isAutorunByDefault(n);
    });

    // Column inputs (column / column_list) get a picker dialog seeded by the
    // upstream table — running the flow up to that point on demand if needed.
    const columnPicker = new ColumnPicker(this.flow, this.executionController, () => ({
      name: this.flowSettings.scriptName,
      description: this.flowSettings.scriptDescription,
      tags: this.flowSettings.tags,
    }));
    this.propertyPanel.onPickColumns = (req) => void columnPicker.pick(req);

    // Functions with their own custom editor (an `editor:` meta or the explicit
    // allowlist) get an icon in the parameters pane header that opens that
    // editor seeded with the node's real upstream tables; the edited values are
    // written back into the panel, so re-render it on completion.
    const funcEditorLauncher = new FuncEditorLauncher(this.flow, this.executionController, () => ({
      name: this.flowSettings.scriptName,
      description: this.flowSettings.scriptDescription,
      tags: this.flowSettings.tags,
    }));
    this.propertyPanel.onEditFuncParams = (node) => {
      // No autorun while the editor dialog is open: the dialog intercepts the
      // global `d4-before-run-action` event, and a mid-dialog rerun of the same
      // function would be mistaken for the dialog's own run action — canceling
      // the rerun's call and resolving the round-trip with stale values. The
      // writeback below reports the edit, so the rerun happens right after.
      this.autorunScheduler?.hold();
      void funcEditorLauncher.open(node).then((applied) => {
        if (applied) {
          // The editor dialog wrote new values into `inputValues` — the same
          // invalidation as any parameter edit in the panel.
          this.flow.notifyNodeParamsChanged(node.id);
          this.propertyPanel.showNode(node, this.executionController?.state.getNodeState(node.id));
        }
      }).catch((e) => grok.shell.error(`Function editor failed: ${e instanceof Error ? e.message : e}`))
        .finally(() => this.autorunScheduler?.release());
    };
    this.executionController.onBreakpointHit = () => {
      grok.shell.info('Breakpoint hit — click Continue in the ribbon to resume');
    };
    this.executionController.onRunEnd = (success: boolean) => {
      // The save dialog's dashboard section refreshes when a run it started ends.
      this.saveDialogRunEnd?.(success);
      if (!success) return;
      grok.shell.info('Flow execution completed');
      // Land the user on the first output node — selecting it opens the
      // bottom-docked panel with the captured value (if any).
      this.autoSelectFirstOutputNode();
    };

    // "Edit settings" on a viewer preview → show the live viewer in the context
    // panel (Datagrok renders its full settings editor) and capture every change
    // back into the node's stored options, so a re-run reproduces the look.
    this.outputPreview.onEditViewer = (nodeRef, viewer) => this.editViewer(nodeRef.id, viewer);

    this.flow.setMinimapCollapsed(this.minimapCollapsed);
    this.updateStartPanelVisibility();
    // Baseline for the Save button: an empty scratch flow starts "unchanged"
    // (a load path records its own baseline once it finishes deserializing).
    this.markSaved();
    this.editorReadyResolve();
  }

  /** Re-render every node so the pre-run "Needs input" hint reflects the
   *  current wiring (connections don't otherwise re-render their endpoints).
   *  Cheap for the small graphs Flow targets; debounced via rAF. */
  private hintRaf = 0;
  private refreshNodeHints(): void {
    if (!this.flow || this.hintRaf) return;
    this.hintRaf = requestAnimationFrame(() => {
      this.hintRaf = 0;
      for (const n of this.flow.getNodes()) void this.flow.updateNode(n.id);
    });
  }

  /** Turn the bottom output panel on for a view created without it (the
   *  creation-script dialog) — e.g. when that view is promoted into the real
   *  editor via "Open In Editor". */
  enableOutputPanel(): void {
    this.outputPreview.setEnabled(true);
  }

  /** Set the overview minimap's collapsed state. Remembered and (re)applied when
   *  the editor finishes initializing, so it can be called before that. */
  setMinimapCollapsed(collapsed: boolean): void {
    this.minimapCollapsed = collapsed;
    this.flow?.setMinimapCollapsed(collapsed);
  }

  /** Find the first output node in the graph (preferring one that has
   *  captured runtime values) and programmatically select it. SetVar nodes
   *  count as outputs too (they compile to the same thing) — a flow whose
   *  terminals are SetVars (e.g. an imported creation script) lands on its
   *  first stored value just like one with Output nodes. */
  private autoSelectFirstOutputNode(): void {
    if (!this.flow) return;
    const outputs = this.flow.getNodes().filter((n) => n.dgNodeType === 'output' || isSetVarNode(n));
    if (outputs.length === 0) return;
    const withValue = outputs.find((n) => {
      const s = this.executionController?.state.getNodeState(n.id);
      return s && s.outputs && Object.keys(s.outputs).length > 0;
    });
    void this.flow.selectNode((withValue ?? outputs[0]).id);
  }

  // ---------- output-view tabs (Spotfire-style internal pages) ----------

  /** The flow's table-carrying terminals, one tab each: Table Outputs,
   *  dataframe-typed Value Outputs, and SetVar terminals whose stored value is
   *  a table (SetVar and Output are the same concept — Q2). */
  private tableOutputs(): OutputTabInfo[] {
    if (!this.flow) return [];
    const res: OutputTabInfo[] = [];
    for (const n of this.flow.getNodes()) {
      if (n.dgNodeType === 'output') {
        const isTable = n.dgTypeName === 'Outputs/Table Output' ||
          (n.dgTypeName === 'Outputs/Value Output' && n.properties['outputType'] === 'dataframe');
        if (isTable)
          res.push({nodeId: n.id, paramName: String(n.properties['paramName'] ?? '').trim() || n.label});
      }
      else if (isSetVarNode(n)) {
        const name = String(n.inputValues['variableName'] ?? '').trim();
        if (name === '') continue;
        const src = this.flow.getInputSource(n.id, 'value');
        const dgType = src ? (src.node.outputs[src.outputKey] as
          {socket?: {dgType?: string}} | undefined)?.socket?.dgType : undefined;
        if (dgType === 'dataframe') res.push({nodeId: n.id, paramName: name});
      }
    }
    return res;
  }

  /** `onNodeStateChanged` fan-out: a completed table output feeds its tab the
   *  run's cloned DataFrame (the same instance the output preview shows —
   *  shared selection for free); an invalidated one marks the tab stale. */
  private updateOutputViewValue(nodeId: string): void {
    if (!this.outputViews || !this.flow) return;
    const node = this.flow.getNodeById(nodeId);
    if (!node || (node.dgNodeType !== 'output' && !isSetVarNode(node))) return;
    const st = this.executionController?.state.getNodeState(nodeId);
    if (st?.status === NodeExecStatus.completed && st.outputs) {
      const summary = Object.values(st.outputs)
        .find((s) => s?.type === 'dataframe' && s.clone instanceof DG.DataFrame);
      if (summary) this.outputViews.setValue(nodeId, summary.clone as DG.DataFrame);
    }
    else if (st?.status === NodeExecStatus.stale)
      this.outputViews.markStale(nodeId);
  }

  /** Tab switched (or the active tab's TableView appeared): surface that
   *  view's ribbon panels and toolbox on this shell-attached view — the
   *  `DG.MultiView.currentView` recipe. Canvas (or an empty tab, which has no
   *  tools of its own) restores Flow's. The view's name and ribbon menu stay
   *  Flow's — the tab strip already names the table. */
  private onOutputTabChanged(tab: OutputTab | null): void {
    const tv = tab?.tv ?? null;
    try {
      if (tv == null) {
        if (this.flowRibbonPanels.length > 0) this.setRibbonPanels(this.flowRibbonPanels);
        this.toolbox = this.functionBrowser.root;
        return;
      }
      // Capture the tv's ribbon ONCE, as the inner (unwrapped) elements:
      // `setRibbonPanels` moves elements between wrappers, so re-reading
      // `tv.getRibbonPanels()` after a swap-and-restore returns empty husks.
      // Also drop items core inline-hides (the TableView's own Save button is
      // display:none'd for a table not in the workspace — copying it verbatim
      // leaves a vestigial empty panel), and lead with Flow's Save pill:
      // saving the flow is what persists this tab's layout.
      if (tab!.ribbonPanels == null) {
        const unwrap = (el: HTMLElement): HTMLElement =>
          el.classList.contains('d4-ribbon-item') && el.firstElementChild instanceof HTMLElement ?
            el.firstElementChild : el;
        tab!.ribbonPanels = tv.getRibbonPanels()
          .map((p) => p.map(unwrap).filter((el) => el.style.display !== 'none'))
          .filter((p) => p.length > 0);
      }
      const saveHost = this.flowRibbonPanels[0]?.[0];
      this.setRibbonPanels(saveHost != null ? [[saveHost], ...tab!.ribbonPanels] : tab!.ribbonPanels);
      this.toolbox = tv.toolbox;
    } catch (e) {
      console.warn('FuncFlow: ribbon/toolbox swap failed', e);
    }
  }

  /** Layouts of the output tabs for persistence in the `.ffjson`, keyed by
   *  paramName; undefined when there is nothing to save (keeps docs clean). */
  private captureOutputViews(): FuncFlowDocument['outputViews'] {
    const layouts = this.outputViews?.captureLayouts() ?? {};
    const names = Object.keys(layouts);
    if (names.length === 0) return undefined;
    const res: NonNullable<FuncFlowDocument['outputViews']> = {};
    for (const name of names) res[name] = {layout: layouts[name]};
    return res;
  }

  /** Accept drops of Datagrok files (→ OpenFile node), DG.Func (→ matching node),
   *  and HTML5 native drags from the function browser (carry typeName via the
   *  FF_DRAG_MIME data type). */
  private setupFileDropTarget(): void {
    ui.makeDroppable(this.canvasContainer, {
      acceptDrop: (drag: any) =>
        (drag instanceof DG.FileInfo && drag.isFile) || drag instanceof DG.Func,
      doDrop: (args) => {
        const drag = args.dragObject;
        // Place the node where it was dropped (like the native drag from the
        // function browser), not in the center. `dropEvent` carries the pointer.
        const ev = args.dropEvent;
        if (drag instanceof DG.Func) {void this.addFuncNode(drag, ev); return;}
        const fi = drag as DG.FileInfo;
        void this.addOpenFileNode(fi.fullPath, ev);
      },
    });

    // Native HTML5 drag/drop from the function browser and the Suggestions pane,
    // plus local (OS) files — the platform's drop overlay is suppressed while
    // this view is active (see the onFileDragEnter subscription).
    this.canvasContainer.addEventListener('dragover', (ev) => {
      if (!ev.dataTransfer) return;
      const types = Array.from(ev.dataTransfer.types);
      if (types.includes(FF_DRAG_MIME) || types.includes(FF_SUGGEST_MIME) || types.includes('Files')) {
        ev.preventDefault();
        ev.dataTransfer.dropEffect = 'copy';
      }
    });
    this.canvasContainer.addEventListener('drop', (ev) => {
      // Local files from the user's computer → one Uploaded File node each.
      const osFiles = ev.dataTransfer?.files;
      if (osFiles != null && osFiles.length > 0) {
        ev.preventDefault();
        ev.stopPropagation();
        void this.addUploadedFileNodes(Array.from(osFiles), ev);
        return;
      }
      // A dragged suggestion carries its full JSON payload — apply it (wiring,
      // prefill and all) at the drop point.
      const sug = ev.dataTransfer?.getData(FF_SUGGEST_MIME);
      if (sug) {
        ev.preventDefault();
        try {
          void this.applySuggestion(JSON.parse(sug) as Suggestion, {clientX: ev.clientX, clientY: ev.clientY});
        } catch {/* malformed payload — not ours, ignore */}
        return;
      }
      const typeName = ev.dataTransfer?.getData(FF_DRAG_MIME);
      if (!typeName) return;
      ev.preventDefault();
      void this.addNodeByTypeAt(typeName, ev.clientX, ev.clientY);
    });
  }

  /** Add a node and place it at the canvas position corresponding to the given
   *  screen-space pointer event. Used by drop handlers. */
  private async addNodeByTypeAt(typeName: string, clientX: number, clientY: number): Promise<FlowNode | null> {
    if (!this.flow) return null;
    const node = createNode(typeName);
    if (!node) {
      grok.shell.warning(`Unknown node type: ${typeName}`);
      return null;
    }
    const {x, y} = this.flow.screenToCanvas(clientX, clientY);
    await this.flow.addNodeAt(node, x, y);
    this.updateStatusBar();
    return node;
  }

  private async addNodeByType(typeName: string): Promise<FlowNode | null> {
    if (!this.flow) return null;
    const node = createNode(typeName);
    if (!node) {
      grok.shell.warning(`Unknown node type: ${typeName}`);
      return null;
    }
    await this.flow.addNodeAtCenter(node);
    this.updateStatusBar();
    return node;
  }

  /** Functions applicable to this view — collected by the AI assistant through
   *  `view.getFunctions()` (the Dart JsViewHost forwards the call here). */
  getFunctions(): DG.Func[] {
    return DG.Func.find({package: 'Flow', tags: ['flowViewFunction']});
  }

  /** Facade the registered Flow view functions (ai-tools.ts) use to act on this
   *  instance — they receive the generic view and reach it via `view.jsView`. */
  aiContext(): FlowAIContext | null {
    if (!this.flow)
      return null;
    return {
      flow: () => this.flow,
      execution: () => this.executionController ?? null,
      addNodeByType: (typeName: string) => this.addNodeByType(typeName),
      run: () => this.runInstrumented(),
      runGuide: (guide) => void this.guideRunner.run(guide, this.guideHost),
    };
  }

  /** Accept a toolbox suggestion: create the node — at the drop point when the
   *  suggestion was dragged onto the canvas (`at`), else next to its (first)
   *  source — prefill the suggested input values (column names, a clicked
   *  molecule), and wire every suggested connection — accepting "Join Tables"
   *  for two selected tables lands it connected to both. */
  private async applySuggestion(s: Suggestion, at?: {clientX: number; clientY: number}): Promise<void> {
    if (!this.flow) return;
    const node = createNode(s.typeName);
    if (!node) {
      grok.shell.warning(`Unknown node type: ${s.typeName}`);
      return;
    }
    const src = s.wire.length > 0 ? this.flow.getNodeById(s.wire[0].fromNodeId) : null;
    if (at) {
      const {x, y} = this.flow.screenToCanvas(at.clientX, at.clientY);
      await this.flow.addNodeAt(node, x, y);
    }
    else if (src)
      await this.flow.addNodeAt(node, src.pos.x + 340, src.pos.y + (s.wire.length > 1 ? 60 : 0));
    else
      await this.flow.addNodeAtCenter(node);

    for (const w of s.wire) {
      if (w.toInput && node.inputs[w.toInput])
        await this.flow.addConnectionByKeys(w.fromNodeId, w.fromOutputKey, node.id, w.toInput);
    }
    if (s.prefill && Object.keys(s.prefill).length > 0) {
      for (const [k, v] of Object.entries(s.prefill)) node.inputValues[k] = v;
      await this.flow.updateNode(node.id);
      // Report the programmatic writes like a panel edit — drives invalidation
      // and lets live-by-default nodes autorun with the prefilled value.
      this.flow.notifyNodeParamsChanged(node.id);
    }
    this.updateStatusBar();
    this.suggestionPane?.refresh();
  }

  /** Place a node of the given type at the drop pointer when one is provided
   *  (drag-and-drop), else in the center (double-click / programmatic add). */
  private addNodeByTypeAtDrop(typeName: string, dropEvent?: MouseEvent): Promise<FlowNode | null> {
    return dropEvent ?
      this.addNodeByTypeAt(typeName, dropEvent.clientX, dropEvent.clientY) :
      this.addNodeByType(typeName);
  }

  private async addOpenFileNode(filePath: string, dropEvent?: MouseEvent): Promise<void> {
    const typeName = this.findOpenFileNodeType();
    if (!typeName) {
      grok.shell.warning('OpenFile function not found in registered nodes');
      return;
    }
    const node = await this.addNodeByTypeAtDrop(typeName, dropEvent);
    if (node) {
      node.inputValues['fullPath'] = filePath;
      await this.flow.updateNode(node.id);
      // Report the programmatic param write like any panel edit — it drives
      // invalidation AND lets the live-by-default autorun load the file.
      this.flow.notifyNodeParamsChanged(node.id);
      grok.shell.info(`Added OpenFile node for: ${filePath}`);
    }
  }

  private async addFuncNode(func: DG.Func, dropEvent?: MouseEvent): Promise<void> {
    const info = getRegisteredFuncs().find((f) => f.func.name === func.name);
    if (!info) {
      grok.shell.warning(`Function "${func.name}" is not available as a node`);
      return;
    }
    if (await this.addNodeByTypeAtDrop(info.nodeTypeName, dropEvent))
      grok.shell.info(`Added node: ${func.name}`);
  }

  private findOpenFileNodeType(): string | null {
    for (const info of getRegisteredFuncs()) {
      if (info.func.name === 'OpenFile' || info.func.name === 'openFile')
        return info.nodeTypeName;
    }
    return null;
  }

  private findUploadedFileNodeType(): string | null {
    for (const info of getRegisteredFuncs()) {
      if (info.func.name === 'readUploadedFile')
        return info.nodeTypeName;
    }
    return null;
  }

  /** Local (OS) files dropped onto the canvas: register the bytes in the
   *  pending store (100 MB cap) and add an Uploaded File node per file. Bytes
   *  reach the server only when the flow is saved (`persistPendingUploads`);
   *  until then the node replays from memory. */
  private async addUploadedFileNodes(files: File[], dropEvent?: MouseEvent): Promise<void> {
    const typeName = this.findUploadedFileNodeType();
    if (!typeName) {
      grok.shell.warning('Uploaded File node is not available');
      return;
    }
    let offset = 0;
    for (const file of files) {
      let fileId: string;
      try {
        fileId = addPendingFile(file.name, new Uint8Array(await file.arrayBuffer()));
      } catch (e: any) {
        grok.shell.error(e?.message ?? String(e));
        continue;
      }
      const node = dropEvent ?
        await this.addNodeByTypeAt(typeName, dropEvent.clientX + offset, dropEvent.clientY + offset) :
        await this.addNodeByType(typeName);
      if (!node) continue;
      node.label = file.name;
      node.inputValues['fileId'] = fileId;
      node.inputValues['fileName'] = file.name;
      await this.flow.updateNode(node.id);
      // Report the programmatic param write like a panel edit — it drives
      // invalidation and lets the live-by-default autorun parse the file.
      this.flow.notifyNodeParamsChanged(node.id);
      offset += 40;
    }
    this.updateStatusBar();
  }

  /** Uploads every pending local file referenced by an Uploaded File node to
   *  the server's GUID-addressed file store and rewrites the node's `fileId`
   *  to the real entity id. Called before any serialization that outlives this
   *  session (entity save, creation-script save, .ffjson export). */
  private async persistPendingUploads(): Promise<void> {
    if (!this.flow) return;
    const upType = this.findUploadedFileNodeType();
    if (!upType) return;
    for (const node of this.flow.getNodes()) {
      if (node.dgTypeName !== upType) continue;
      const fileId = node.inputValues['fileId'];
      if (typeof fileId !== 'string' || !isPendingFileId(fileId)) continue;
      if (!getPendingFile(fileId)) {
        throw new Error(`Uploaded file "${node.inputValues['fileName'] ?? node.label}" is no longer ` +
          'available in this session — drop it onto the canvas again before saving');
      }
      const fi = await persistPendingFile(fileId);
      node.inputValues['fileId'] = fi.id;
      await this.flow.updateNode(node.id);
    }
  }

  /** Post-save permission sync: the freshly saved body references the real
   *  file ids, so delegate to the shared util (share-time syncs run globally
   *  through the flowShareSync autostart, view or no view). */
  private async syncUploadedFilePermissions(): Promise<void> {
    if (this.boundScript?.script)
      await syncFlowFilePermissions(this.boundScript);
  }

  // ---------- guide system (tutorials + how-to) ----------

  /** What the interactive guides need from this view. */
  private get guideHost(): GuideHost {
    return {
      getFlow: () => this.flow,
      showFunctionBrowser: () => {
        grok.shell.windows.showToolbox = true;
        try {
          this.functionBrowser.render();
        } catch {/* not ready yet */}
      },
      showToolboxTab: (name) => {
        try {
          this.functionBrowser.showTab(name);
        } catch {/* not ready yet */}
      },
      anchorEl: this.helpButton,
    };
  }

  // ---------- start panel (U1: never open empty) ----------

  /** A welcoming overlay shown over the empty canvas: pick a template, open a
   *  saved flow, import from a table, or start blank — instead of a blank page. */
  private buildStartPanel(): HTMLElement {
    const title = ui.divText('Start a flow', 'funcflow-start-title');
    const subtitle = ui.divText(
      'Build a data pipeline by chaining functions — no code required.',
      'funcflow-start-subtitle');

    const cards = FLOW_TEMPLATES.map((t) => {
      const card = ui.divV([
        ui.divText(t.label, 'funcflow-start-card-title'),
        ui.divText(t.desc, 'funcflow-start-card-desc'),
      ], 'funcflow-start-card');
      setTid(card, 'start-template', t.file.replace(/\.ffjson$/i, ''));
      card.onclick = (): void => void this.loadTemplate(t.file);
      ui.tooltip.bind(card, `Open the "${t.label}" template`);
      return card;
    });

    const blankCard = ui.divV([
      ui.divText('Blank canvas', 'funcflow-start-card-title'),
      ui.divText('Start from scratch.', 'funcflow-start-card-desc'),
    ], 'funcflow-start-card funcflow-start-card-blank');
    setTid(blankCard, 'start-blank');
    blankCard.onclick = (): void => this.hideStartPanel();
    cards.push(blankCard);

    // Primary call-to-action: launch the hands-on tour (the first tutorial),
    // which walks a newcomer through loading data and adding a column.
    const firstFlowBtn = ui.button('Create your first flow', () => {
      this.hideStartPanel();
      void this.guideRunner.run(TUTORIALS[0], this.guideHost);
    });
    firstFlowBtn.classList.add('funcflow-start-tour');
    ui.tooltip.bind(firstFlowBtn, `Hands-on: ${TUTORIALS[0].title}`);
    setTid(firstFlowBtn, 'start-first-flow');

    const openBtn = ui.button('Open a flow…', () => void this.openFromPlatform());
    setTid(openBtn, 'start-open');
    const actions = ui.divH([firstFlowBtn, openBtn], 'funcflow-start-actions');

    // Discovery hint, with an actionable link that launches the interface tour.
    const interfaceTour = TUTORIALS.find((t) => t.id === 'interface-tour');
    const tourLink = ui.link('take a tour of the interface', () => {
      this.hideStartPanel();
      if (interfaceTour) void this.guideRunner.run(interfaceTour, this.guideHost);
    }, 'Walk through every part of the UI — toolbox, ribbon, canvas, and context panel');
    setTid(tourLink, 'start-ui-tour');
    const hint = ui.div([], 'funcflow-start-hint');
    hint.appendChild(document.createTextNode('New here? Create your first flow above, or '));
    hint.appendChild(tourLink);
    hint.appendChild(document.createTextNode(
      '. You can also double-click a function in the list on the left, or drag a file onto the canvas.'));

    // The 10 most recently updated flows visible to the user, rendered as
    // one-line entity rows (same markup as the Browse panel). Populated
    // asynchronously; hidden entirely when the server has none.
    this.recentFlowsHost = setTid(ui.divV([], 'funcflow-start-recent'), 'start-recent');
    void this.refreshRecentFlows();

    const panel = ui.divV([
      title, subtitle,
      ui.divH(cards, 'funcflow-start-cards'),
      this.recentFlowsHost,
      actions, hint,
    ], 'funcflow-start-panel');
    setTid(panel, 'start-panel');
    return setTid(ui.div([this.buildStartBackground(), panel], 'funcflow-start-overlay'), 'start-overlay');
  }

  /** Fill the start panel's "Recent flows" section: the 10 most recently
   *  updated flow scripts the user can see (own + shared), newest first.
   *  Each row is the platform's one-line entity markup; clicking it opens
   *  the flow in this view. Server errors just leave the section empty. */
  private async refreshRecentFlows(): Promise<void> {
    if (this.recentFlowsLoading) return;
    this.recentFlowsLoading = true;
    try {
      const flows = await grok.dapi.scripts.list(
        {pageSize: 12, pageNumber: 1, filter: 'language="flow"', order: '!updatedOn'});
      this.recentFlowsHost.innerHTML = '';
      if (flows.length === 0) return;
      this.recentFlowsHost.appendChild(ui.divText('Recent flows', 'funcflow-start-recent-title'));
      const list = ui.divV([], 'funcflow-start-recent-list');
      for (const f of flows) {
        const row = ui.div([ui.render(f)], 'funcflow-start-recent-item');
        setTid(row, 'start-recent-item', f.friendlyName || f.name);
        // Capture-phase so the entity markup's own click (set current object)
        // doesn't swallow the open.
        row.addEventListener('click', (ev) => {
          ev.preventDefault();
          ev.stopPropagation();
          void this.openRecentFlow(f);
        }, true);
        list.appendChild(row);
      }
      this.recentFlowsHost.appendChild(list);
    } catch { /* server unreachable — no section */ }
    finally {
      this.recentFlowsLoading = false;
    }
  }

  private async openRecentFlow(f: DG.Script): Promise<void> {
    try {
      // Re-find by id so the body is guaranteed present, not a lean listing row.
      const full = await grok.dapi.scripts.find(f.id);
      this.bindScript((full as DG.Script | null) ?? f);
      this.hideStartPanel();
    } catch (e: any) {
      grok.shell.error(`Could not open "${f.friendlyName || f.name}": ${e?.message ?? e}`);
    }
  }

  /** Decorative animated backdrop host. The graph is drawn by
   *  `renderStartBackground` at the container's *actual* size (re-run on resize)
   *  so it fills any shape with no scaling distortion. */
  private buildStartBackground(): HTMLElement {
    this.startBg = ui.div([], 'funcflow-start-bg');
    setTid(this.startBg, 'start-bg');
    return this.startBg;
  }

  /** Lay a flowing-edges graph across the real container dimensions: nodes on a
   *  jittered grid sized to W×H, connected left-to-right with horizontal-flow
   *  beziers. Because the SVG viewBox equals the pixel size, dots stay round and
   *  curves keep their shape on portrait, square, or ultrawide alike. The jitter
   *  is deterministic, so resizing reflows smoothly instead of reshuffling. */
  private renderStartBackground(): void {
    if (!this.startBg) return;
    const W = this.canvasContainer.clientWidth || 1200;
    const H = this.canvasContainer.clientHeight || 800;
    const COLORS = ['#FF9100', '#42A5F5', '#66BB6A', '#AB47BC', '#26C6DA', '#EC407A'];

    // Stable pseudo-random in [0,1) from two integer seeds.
    const rnd = (i: number, j: number): number => {
      const x = Math.sin(i * 127.1 + j * 311.7) * 43758.5453;
      return x - Math.floor(x);
    };
    const cols = Math.min(8, Math.max(4, Math.round(W / 200)));
    const rows = Math.min(6, Math.max(3, Math.round(H / 150)));
    const cw = W / cols;
    const ch = H / rows;
    const nx = (c: number, r: number): number => (c + 0.5) * cw + (rnd(c, r) - 0.5) * cw * 0.55;
    const ny = (c: number, r: number): number => (r + 0.5) * ch + (rnd(c + 9, r + 4) - 0.5) * ch * 0.55;

    const paths: string[] = [];
    const dots: string[] = [];
    let k = 0;
    for (let c = 0; c < cols; c++) {
      for (let r = 0; r < rows; r++) {
        const ax = nx(c, r);
        const ay = ny(c, r);
        dots.push(`<circle cx="${ax.toFixed(1)}" cy="${ay.toFixed(1)}" r="5.5" ` +
          `fill="${COLORS[(c + r) % COLORS.length]}" stroke="#ffffff" stroke-width="2"/>`);
        if (c >= cols - 1) continue;
        // Right neighbor, plus an occasional diagonal branch — a flowing DAG.
        const targets = [r];
        if (rnd(c + 3, r + 7) > 0.5 && r + 1 < rows) targets.push(r + 1);
        if (rnd(c + 5, r + 2) > 0.62 && r - 1 >= 0) targets.push(r - 1);
        for (const tr of targets) {
          const bx = nx(c + 1, tr);
          const by = ny(c + 1, tr);
          const mx = (ax + bx) / 2;
          paths.push(`<path d="M${ax.toFixed(1)},${ay.toFixed(1)} ` +
            `C${mx.toFixed(1)},${ay.toFixed(1)} ${mx.toFixed(1)},${by.toFixed(1)} ${bx.toFixed(1)},${by.toFixed(1)}" ` +
            `stroke="${COLORS[k++ % COLORS.length]}"/>`);
        }
      }
    }
    this.startBg.innerHTML =
      `<svg viewBox="0 0 ${W} ${H}" width="100%" height="100%" xmlns="http://www.w3.org/2000/svg">` +
      paths.join('') + dots.join('') + `</svg>`;
  }

  /** Draw the backdrop sized to the canvas, on the next frame so the container
   *  is laid out. Done lazily each time the panel is shown (no resize tracking). */
  private drawStartBackgroundSoon(): void {
    if (this.startBgRaf) cancelAnimationFrame(this.startBgRaf);
    this.startBgRaf = requestAnimationFrame(() => {
      this.startBgRaf = 0;
      this.renderStartBackground();
    });
  }

  /** Show the start overlay (e.g. from the ribbon "Templates…" item). */
  private showStartPanel(): void {
    this.startPanel.style.display = 'flex';
    this.drawStartBackgroundSoon();
    void this.refreshRecentFlows();
  }

  private hideStartPanel(): void {
    this.startPanel.style.display = 'none';
  }

  /** Hide the overlay only while the canvas has content; show it on an empty one. */
  private updateStartPanelVisibility(): void {
    const wasHidden = this.startPanel.style.display === 'none';
    const empty = !this.flow || this.flow.getNodeCount() === 0;
    this.startPanel.style.display = empty ? 'flex' : 'none';
    if (empty) {
      this.drawStartBackgroundSoon();
      // Re-shown (New flow / everything deleted) → the recent list may be
      // stale (the user just saved a flow); reload it.
      if (wasHidden) void this.refreshRecentFlows();
    }
  }

  /** Load a bundled template flow from the package `files/` folder. */
  private async loadTemplate(file: string): Promise<void> {
    try {
      const json = await _package.files.readAsText(file);
      await this.loadFromJson(json);
      this.hideStartPanel();
      grok.shell.info(`Opened template: ${file.replace(/\.ffjson$/i, '')}`);
    } catch (e: any) {
      grok.shell.error(`Could not open template "${file}": ${e?.message ?? e}`);
    }
  }

  // ---------- ribbon ----------

  private setupRibbon(): void {
    // In creation-script mode Save writes creation scripts back to the tables,
    // and the platform-entity save/open options are hidden so the two save
    // targets cannot be confused.
    const creationMode = this.tableInfos.length > 0;

    // Menu is grouped around what a scientist wants to do — Flow / Run / Edit /
    // Arrange — with the script-and-code machinery tucked under "Advanced".
    let m = DG.Menu.create()
      .group('Flow')
      .item('New…', () => void this.newFlow());
    if (!creationMode)
      m = m.item('Open from platform…', () => void this.openFromPlatform());
    m = m.item('Save', () => void this.saveFlow());
    if (!creationMode)
      m = m.item('Save As…', () => void this.saveDialog({asNew: true}));
    this.ribbonMenu = m
      .separator()
      .item('Import .ffjson…', () => void this.openFlow())
      .item('Export .ffjson', () => void this.exportFfjson())
      .separator()
      .item('Templates…', () => this.showStartPanel())
      .item('Settings…', () => this.editSettings())
      .endGroup()
      .group('Run')
      .item('Run', () => this.runInstrumented())
      .item('Debug (stop at breakpoints)', () => this.debugInstrumented())
      .item('Continue', () => this.executionController?.continueBreakpoint())
      .item('Stop', () => this.executionController?.stopRun())
      .separator()
      .item('Clear run highlights', () => {
        this.executionController?.resetVisuals();
        this.outputViews?.clearValues();
      })
      .endGroup()
      .group('Edit')
      .item('Undo', () => void this.flow?.undo())
      .item('Redo', () => void this.flow?.redo())
      .endGroup()
      .group('Arrange')
      .item('Tidy up layout', () => this.cleanLayout())
      .item('Zoom to fit', () => void this.flow?.zoomToFit())
      .item('Zoom in', () => this.flow?.zoomIn())
      .item('Zoom out', () => this.flow?.zoomOut())
      .separator()
      .item('Show/hide function list', () => this.toggleToolbox())
      .endGroup()
      .group('Advanced')
      .item('Describe this flow…', () => this.describeFlow())
      .item('See the steps (generated script)…', () => this.generateAndPreview())
      .item('Copy script', () => this.copyScriptToClipboard())
      .item('Export as .js file', () => this.exportAsJs())
      .item('Check for problems', () => this.showValidation())
      .separator()
      .item('Run as plain script (no live view)', () => this.runScript())
      .separator()
      .item('Import from a table’s history…', () => this.importCreationScriptDialog())
      .item('Export as table-creation script…', () => this.compileToCreationScript())
      .endGroup();

    const ribbonIcon = (icon: string, action: () => void, tooltip: string, id: string): HTMLElement =>
      setTid(ui.iconFA(icon, action, tooltip), 'ribbon', id);

    // Autorun toggle: faded outline bolt when off (default), colored + filled
    // (font-weight 600) when on — see `.ff-autorun-toggle` in funcflow.css.
    // The tooltip is dynamic so it always names the current state.
    const autorunIcon = setTid(ui.iconFA('bolt', () => this.toggleAutorun()), 'ribbon', 'autorun');
    autorunIcon.classList.add('ff-autorun-toggle');
    ui.tooltip.bind(autorunIcon, () => this.autorunScheduler?.enabled ?
      'Autorun is on — the flow reruns the affected nodes after every change. Click to turn off.' :
      'Autorun is off — click to rerun the flow (only the affected nodes) automatically after every change. ' +
      'Live nodes (Open File, Add New Column, viewers) still rerun on change.');
    this.autorunIcon = autorunIcon;

    // Saving leads the ribbon; saveFlow routes to the right target (entity
    // update / Save As for never-saved flows / creation scripts).
    const saveButton = creationMode ?
      ui.bigButton('Save', () => void this.saveCreationScriptsDialog(), 'Review and save a creation script for each table') :
      ui.bigButton('Save', () => {if (this.saveAvailability().enabled) void this.saveFlow();}, '');
    saveButton.prepend(ui.iconFA('cloud-upload'));
    const saveEl = setTid(saveButton, 'ribbon', creationMode ? 'save-creation-scripts' : 'save');

    // The flow Save button is state-tracked: greyed and non-clickable
    // (`pointer-events:none` via `.ff-ribbon-btn-disabled`) when there's nothing
    // to save. A disabled button doesn't fire hover, so the dynamic "why"
    // tooltip lives on a wrapper that stays hoverable in both states.
    let saveHost: HTMLElement = saveEl;
    if (!creationMode) {
      this.saveButton = saveButton;
      saveHost = ui.div([saveEl], 'ff-save-btn-wrap');
      saveHost.style.marginRight = '8px';
      ui.tooltip.bind(saveHost, () => this.saveAvailability().tooltip);
    }
    const savePanel: HTMLElement[] = [saveHost];
    if (!creationMode)
      savePanel.push(ribbonIcon('folder-open', () => void this.openFromPlatform(), 'Open a flow', 'open'));
    const panels: HTMLElement[][] = [
      savePanel,
      [
        ribbonIcon('play', () => this.runInstrumented(), 'Run the flow', 'run'),
        ribbonIcon('bug', () => this.debugInstrumented(), 'Debug (stop at breakpoints)', 'debug'),
        ribbonIcon('forward', () => this.executionController?.continueBreakpoint(), 'Continue', 'continue'),
        ribbonIcon('stop', () => this.executionController?.stopRun(), 'Stop', 'stop'),
        autorunIcon,
      ],
      [
        ribbonIcon('eye', () => this.generateAndPreview(), 'See the steps (generated script)', 'view-script'),
      ],
      [
        ribbonIcon('undo', () => void this.flow?.undo(), 'Undo (Ctrl+Z)', 'undo'),
        ribbonIcon('redo', () => void this.flow?.redo(), 'Redo (Ctrl+Shift+Z)', 'redo'),
      ],
      [
        ribbonIcon('sitemap', () => this.cleanLayout(), 'Tidy up layout', 'layout'),
        ribbonIcon('search-plus', () => this.flow?.zoomIn(), 'Zoom in', 'zoom-in'),
        ribbonIcon('search-minus', () => this.flow?.zoomOut(), 'Zoom out', 'zoom-out'),
        ribbonIcon('compress-arrows-alt', () => void this.flow?.zoomToFit(), 'Zoom to fit', 'zoom-fit'),
        ribbonIcon('list-ul', () => this.toggleToolbox(), 'Show/hide function list', 'toggle-browser'),
      ],
      [
        ribbonIcon('graduation-cap', () => openGuideMenu(this.guideHost, this.guideRunner), 'Tutorials & help', 'help'),
      ],
    ];

    this.setRibbonPanels(panels);
    // Kept by reference: `setRibbonPanels` moves (reparents) elements, so these
    // same references restore Flow's ribbon when leaving an output-view tab.
    this.flowRibbonPanels = panels;
    this.updateSaveButtonState();
  }

  private setupStatusBar(): void {
    this.statusBarPanels = [this.statusBar as HTMLDivElement];
  }

  private toggleToolbox(): void {
    grok.shell.windows.showToolbox = !grok.shell.windows.showToolbox;
  }

  /** Flip autorun mode and reflect it on the ribbon icon (grey ↔ colored).
   *  Turning it ON immediately schedules everything that has no fresh result
   *  (a new flow runs entirely; a half-run one completes itself) — enabling
   *  autorun should not sit idle until the first edit. */
  private toggleAutorun(): void {
    if (!this.autorunScheduler) return;
    const on = this.autorunScheduler.toggle();
    this.autorunIcon?.classList.toggle('ff-autorun-on', on);
    if (on) {
      const pending = this.executionController?.pendingNodes() ?? new Set<string>();
      if (pending.size > 0) this.autorunScheduler.kick(pending);
    }
  }

  // ---------- Save button state ----------

  /** A stable serialization of the graph + settings, for detecting unsaved
   *  changes. Side-effect-free (unlike `entityBodyText`, which stamps a tag).
   *  `created` / `modified` / `author` are dropped — `serializeFlow` stamps
   *  fresh timestamps on every call, so keeping them would make two snapshots of
   *  the *same* graph always differ, and Save could never grey out. */
  private currentSnapshot(): string {
    if (!this.flow) return '';
    try {
      const doc = serializeFlow(this.flow, this.flowSettings) as unknown as Record<string, unknown>;
      delete doc.created;
      delete doc.modified;
      delete doc.author;
      return JSON.stringify(doc);
    } catch {
      return '';
    }
  }

  /** Record the current graph as the saved baseline (after a save, or after
   *  loading a flow that already lives on the server / a fresh empty flow). */
  private markSaved(): void {
    this.savedSnapshot = this.currentSnapshot();
    this.updateSaveButtonState();
  }

  /** Whether Save is available, and its tooltip. The button is the gateway to
   *  both saving the script AND publishing a dashboard, so it is enabled for
   *  any non-empty canvas; the tooltip reflects the dirty state. */
  private saveAvailability(): {enabled: boolean; tooltip: string} {
    if ((this.flow?.getNodeCount() ?? 0) === 0)
      return {enabled: false, tooltip: 'Nothing to save yet — the canvas is empty'};
    if (this.savedSnapshot !== null && this.currentSnapshot() === this.savedSnapshot)
      return {enabled: true, tooltip: 'No changes since the last save — open to publish a dashboard or save as new'};
    return {enabled: true, tooltip: 'Save this flow to the platform'};
  }

  /** Reflect save availability on the ribbon button — greyed + non-clickable
   *  when disabled (the wrapper carries the reason tooltip). A no-op in
   *  creation-script mode (no state-tracked button). */
  private updateSaveButtonState(): void {
    const btn = this.saveButton;
    if (!btn) return;
    btn.classList.toggle('ff-ribbon-btn-disabled', !this.saveAvailability().enabled);
  }

  // ---------- auto-pin the preview view ----------

  /** A flow app / flow script opens as an unpinned **preview** view (a Dart
   *  ViewBase concept): the toolbox stays hidden until the view is pinned. Pin
   *  it the moment the user does anything in it — a click on the canvas, an
   *  edit — so the toolbox appears. One-shot: the handler removes itself once it
   *  fires. Skipped for embedded hosts (the creation-script dialog), which are
   *  never the shell's current view. */
  private setupAutoPin(): void {
    if (!this.outputPanelEnabled) return;
    const opts: AddEventListenerOptions = {capture: true};
    const handler = (): void => {
      // Only when THIS view is the current one — never pin someone else's view.
      const cur = grok.shell.v;
      if (!this.dart || !cur || cur.dart !== this.dart) return;
      try {
        (grok.shell.v as DG.View).pin?.();
      } catch {/* not pinnable in this host — ignore */}
      this.teardownAutoPin();
    };
    this.autoPinHandler = handler;
    this.root.addEventListener('pointerdown', handler, opts);
    this.root.addEventListener('keydown', handler, opts);
  }

  private teardownAutoPin(): void {
    const h = this.autoPinHandler;
    if (!h) return;
    this.autoPinHandler = null;
    const opts: AddEventListenerOptions = {capture: true};
    this.root.removeEventListener('pointerdown', h, opts);
    this.root.removeEventListener('keydown', h, opts);
  }

  /** View closed — a pending debounced autorun must not fire into it. */
  detach(): void {
    this.autorunScheduler?.reset();
    this.teardownAutoPin();
    for (const sub of this.platformSubs) sub.unsubscribe();
    this.platformSubs = [];
    this.suggestionPane?.destroy();
    this.functionBrowser?.destroy();
    this.outputViews?.destroy();
    this.flow?.destroy();
    super.detach();
  }

  /** Re-arrange the existing graph with the importer's layered/banded layout. */
  private cleanLayout(): void {
    if (!this.flow || this.flow.getNodeCount() === 0) {
      grok.shell.info('Nothing to lay out');
      return;
    }
    void this.flow.autoLayout();
  }

  private updateStatusBar(): void {
    if (!this.flow) return;
    this.nodeCountLabel.textContent = `Nodes: ${this.flow.getNodeCount()}`;
    this.linkCountLabel.textContent = `Links: ${this.flow.getConnectionCount()}`;
  }

  // ---------- file actions ----------

  private async newFlow(): Promise<void> {
    await this.flow.clear();
    this.boundScript = null;
    this.dashboardProjectId = null;
    this.propertyPanel.clear();
    this.updateStatusBar();
    this.updateStartPanelVisibility();
    this.markSaved(); // empty canvas → nothing to save
    grok.shell.info('New flow created');
  }

  private async openFlow(): Promise<void> {
    const input = document.createElement('input');
    input.type = 'file';
    input.accept = '.ffjson,.json';
    input.addEventListener('change', async () => {
      const file = input.files?.[0];
      if (!file) return;
      try {
        const doc = await loadFlowFromFile(file);
        await this.loadFromDoc(doc);
        this.boundScript = null; // an imported file is a new, unsaved flow
        grok.shell.info(`Loaded flow: ${doc.name}`);
      } catch (e: any) {
        grok.shell.error(`Failed to load flow: ${e.message}`);
      }
    });
    input.click();
  }

  /** Save: in creation-script mode writes creation scripts back to the tables;
   *  otherwise always opens the combined save dialog — script name/space plus
   *  the dashboard section (run the flow → publish its result tables). */
  private async saveFlow(): Promise<void> {
    if (!this.flow) return;
    if (this.tableInfos.length > 0) {
      await this.saveCreationScriptsDialog();
      return;
    }
    await this.saveDialog();
  }

  /** Whether a script id resolves to a real, accessible server entity.
   *  `grok.dapi.scripts.find` throws (or resolves nullish) for a missing or
   *  inaccessible id — either way the flow isn't saved yet, so Save must ask
   *  for a name rather than silently update a non-existent entity. */
  private async scriptExistsOnServer(id: string): Promise<boolean> {
    try {
      return (await grok.dapi.scripts.find(id)) != null;
    } catch {
      return false;
    }
  }

  /** The `.flow` entity body for the current graph — the single writer, so the
   *  annotation header and the ffjson payload can never disagree. */
  private entityBodyText(): string {
    if (!this.flowSettings.tags.includes(FLOW_TAG))
      this.flowSettings.tags.push(FLOW_TAG);
    return flowScriptText(this.flow, this.flowSettings, {
      outputViews: this.captureOutputViews(),
      dashboard: this.dashboardProjectId != null ? {projectId: this.dashboardProjectId} : undefined,
    });
  }

  /** @param silent no "saved" balloon — for follow-up saves that only persist
   *  metadata (e.g. the dashboard binding right after a publish). */
  private async saveToServer(silent = false): Promise<void> {
    try {
      // Local files live in memory until this moment — persist them first so
      // the saved body references real, server-addressable file ids.
      await this.persistPendingUploads();
      const script = DG.Script.create(this.entityBodyText());
      if (this.boundScript?.id)
        script.id = this.boundScript.id;
      this.boundScript = await grok.dapi.scripts.save(script);
      this.name = this.flowSettings.scriptName;
      this.updatePath();
      this.markSaved(); // this graph is now the saved baseline
      void this.syncUploadedFilePermissions();
      if (!silent) grok.shell.info(`Flow "${this.flowSettings.scriptName}" saved`);
    } catch (e: any) {
      grok.shell.error(`Failed to save flow: ${e?.message ?? e}`);
    }
  }

  /** Fires while the save dialog is open and a run initiated from it ends —
   *  refreshes the dialog's dashboard section. */
  private saveDialogRunEnd: ((success: boolean) => void) | null = null;

  /** The combined save dialog: script name / description / space on top, the
   *  dashboard section below. A bound entity is updated in place (Save As
   *  forces a new one via `asNew`). The dashboard section is run-aware: no
   *  computed outputs → a Run button; outputs → the table list + a "Create
   *  dashboard" toggle that opens the core Save-project dialog after saving. */
  private async saveDialog(opts: {asNew?: boolean} = {}): Promise<void> {
    const nameInput = ui.input.string('Name', {value: this.flowSettings.scriptName,
      tooltipText: 'The flow is saved as a platform script entity under this name'});
    const descInput = ui.input.textArea('Description', {value: this.flowSettings.scriptDescription,
      tooltipText: 'Shown in galleries, previews and the context panel'});

    const noSpaceLabel = 'Saved as a plain script (not bound to a space)';
    let targetSpace: DG.Project | null = null;
    const spaceLabel = ui.divText(noSpaceLabel, 'funcflow-save-space-label');
    const clearLink = ui.link('clear', () => updateSpace(null),
      'Save as a plain script, not bound to a space');
    const warningDiv = ui.divText('', 'funcflow-save-name-warning');
    warningDiv.style.color = '#b26a00'; // matches the warnings strip in saveCreationScriptsDialog

    // Best-effort duplicate-name check: purely advisory — the server keeps
    // names unique within a space by suffixing, and plain scripts may repeat.
    let warnTimer: ReturnType<typeof setTimeout> | null = null;
    const refreshNameWarning = async () => {
      const name = nameInput.value.trim();
      warningDiv.textContent = '';
      if (name === '') return;
      const esc = name.replace(/"/g, '\\"');
      const scope = targetSpace ? ` and namespace = "${targetSpace.nqName}:"` : '';
      try {
        const clashes = (await grok.dapi.scripts
          .filter(`language = "${FLOW_LANGUAGE}" and friendlyName = "${esc}"${scope}`).list())
          .filter((s) => s.id !== this.boundScript?.id);
        if (clashes.length > 0) {
          warningDiv.textContent = targetSpace ?
            `A flow named "${name}" already exists in "${targetSpace.friendlyName}" — it will be saved under a unique name` :
            `A flow named "${name}" already exists`;
        }
      } catch {/* advisory only */}
    };
    const updateSpace = (space: DG.Project | null) => {
      targetSpace = space;
      spaceLabel.textContent = space ? `Space: ${space.friendlyName}` : noSpaceLabel;
      clearLink.style.display = space ? '' : 'none';
      void refreshNameWarning();
    };
    updateSpace(null);

    const pickerHost = ui.div([]);
    let picker: SpacePicker | null = null;
    const bindBtn = ui.button('Add to space…', async () => {
      if (picker == null) {
        picker = await SpacePicker.create();
        picker.onChanged = (space) => updateSpace(space);
        pickerHost.appendChild(picker.root);
      } else
        pickerHost.style.display = pickerHost.style.display === 'none' ? '' : 'none';
    });
    ui.tooltip.bind(bindBtn,
      'Choose a space (or subspace) to organize and share this flow; by default it is saved as a plain script in your namespace');
    ui.tooltip.bind(spaceLabel, 'Where this flow will live');

    // ---- dashboard section (run-aware) ----
    const computedTabs = (): OutputTab[] =>
      this.outputViews?.getTabs().filter((t) => t.df != null) ?? [];
    const dashHost = ui.divV([], 'ff-save-dash');
    const publishInput = ui.input.bool('Create dashboard', {value: true,
      tooltipText: 'After saving the flow, open the standard Save-project dialog seeded with the ' +
        'computed output tables and their layouts (data sync, sharing, upload)'});
    const refreshDash = (running = false): void => {
      ui.empty(dashHost);
      dashHost.appendChild(ui.divText('Dashboard', 'ff-save-dash-title'));
      const tabs = computedTabs();
      if (tabs.length > 0) {
        dashHost.appendChild(ui.divV(tabs.map((t) => ui.divText(
          `• ${t.df!.name || t.paramName} (${t.df!.rowCount.toLocaleString()} × ${t.df!.columns.length})`))));
        dashHost.appendChild(publishInput.root);
        if (this.dashboardProjectId != null) {
          const unbind = ui.link('publish as new', () => {
            this.dashboardProjectId = null;
            refreshDash();
          }, 'Forget the bound project — the next publish creates a new dashboard');
          const bound = ui.div([], 'ff-save-dash-hint');
          bound.appendChild(document.createTextNode('Updates the previously published dashboard — or '));
          bound.appendChild(unbind);
          dashHost.appendChild(bound);
        }
      }
      else if (running)
        dashHost.appendChild(ui.divText('Running the flow…', 'ff-save-dash-hint'));
      else {
        dashHost.appendChild(ui.divText(
          'Run the flow to publish its result tables as a dashboard.', 'ff-save-dash-hint'));
        dashHost.appendChild(ui.button('Run the flow', () => {
          refreshDash(true);
          this.runInstrumented();
        }));
      }
    };
    refreshDash();
    this.saveDialogRunEnd = () => refreshDash();

    const dlg = ui.dialog({title: 'Save Flow'})
      .add(ui.divV([nameInput.root, descInput.root, warningDiv,
        ui.divH([bindBtn, spaceLabel, clearLink], {style: {alignItems: 'center', gap: '8px'}}),
        pickerHost, dashHost]))
      .onOK(async () => {
        const name = nameInput.value.trim();
        if (name === '') { // reachable via Enter even while the Save button is disabled
          grok.shell.warning('Give the flow a name first');
          return;
        }
        this.flowSettings.scriptName = name;
        this.flowSettings.scriptDescription = descInput.value;
        // Update the bound entity in place; Save As (or a stale/deleted binding —
        // e.g. a template id `find` can't resolve) creates a new one. A new
        // script also gets a new dashboard.
        if (opts.asNew || !(this.boundScript?.id && await this.scriptExistsOnServer(this.boundScript.id))) {
          this.boundScript = null;
          this.dashboardProjectId = null;
        }
        await this.saveToServer();
        const saved = this.boundScript as DG.Script | null;
        if (targetSpace && saved?.id) {
          try {
            await grok.dapi.spaces.id(targetSpace.id).addEntity(saved.id, false);
            grok.shell.info(`Added to space "${targetSpace.friendlyName}"`);
          } catch (e: any) {
            grok.shell.error(`Could not add to space: ${e?.message ?? e}`);
          }
        }
        if (saved != null && publishInput.value === true && computedTabs().length > 0)
          await this.openDashboardDialog();
      });
    dlg.onClose.subscribe(() => this.saveDialogRunEnd = null);
    dlg.show({width: 500});
    // Validate before close: empty names never reach the OK handler.
    const okBtn = dlg.getButton('OK') as HTMLButtonElement | null;
    if (okBtn) okBtn.textContent = 'Save';
    const syncOk = () => {
      if (okBtn == null) return;
      const empty = nameInput.value.trim() === '';
      okBtn.disabled = empty;
      okBtn.classList.toggle('disabled', empty);
    };
    nameInput.onChanged.subscribe(() => {
      syncOk();
      if (warnTimer != null) clearTimeout(warnTimer);
      warnTimer = setTimeout(() => void refreshNameWarning(), 400);
    });
    syncOk();
    void refreshNameWarning();
  }

  /** Pick a flow entity from the server and open it in this view. */
  private async openFromPlatform(): Promise<void> {
    let flows: DG.Script[] = [];
    try {
      flows = await grok.dapi.scripts.filter('language = "flow"').list();
    } catch (e: any) {
      grok.shell.error(`Could not list flows: ${e?.message ?? e}`);
      return;
    }
    if (flows.length === 0) {
      grok.shell.info('No flows on this server yet — save one first');
      return;
    }
    const byLabel = new Map(flows.map((f) => [f.friendlyName || f.name, f] as const));
    const items = [...byLabel.keys()].sort((a, b) => a.localeCompare(b));
    const input = ui.input.choice('Flow', {value: items[0], items});
    ui.dialog({title: 'Open Flow'})
      .add(input.root)
      .onOK(async () => {
        const picked = byLabel.get(input.value ?? '');
        if (!picked) return;
        // Re-find by id so the body is guaranteed present, not a lean listing row.
        const full = await grok.dapi.scripts.find(picked.id);
        this.bindScript(full ?? picked);
      })
      .show();
  }

  // ---------- publish as dashboard (via the core Save-project dialog) ----------

  /** Input-node types that cannot be replayed from a literal default — a flow
   *  taking one of these cannot be data-synced (the creation script would have
   *  no value to pass), so its outputs publish as static snapshots. */
  private static readonly NON_SYNCABLE_INPUTS = new Set([
    'Inputs/Table Input', 'Inputs/Column Input', 'Inputs/Column List Input',
    'Inputs/File Input', 'Inputs/Blob Input', 'Inputs/Map Input', 'Inputs/Dynamic Input',
  ]);

  /** How many outputs the flow script declares (output nodes + SetVar
   *  terminals). Decides whether the creation script needs the output
   *  accessor — `.param` is only valid with more than one output. */
  private flowOutputCount(): number {
    let count = 0;
    for (const n of this.flow.getNodes()) {
      if (n.dgNodeType === 'output') count++;
      else if (isSetVarNode(n) && String(n.inputValues['variableName'] ?? '').trim() !== '' &&
          this.flow.getInputSource(n.id, 'value') != null) count++;
    }
    return count;
  }

  /** Stamp each computed output table with the producing call (`.script` +
   *  `.VariableName` df tags) so the core Save-project dialog offers data sync
   *  — identical producing calls dedup on project open, so the flow runs ONCE
   *  and every table binds its own output. Skipped (static snapshots) when the
   *  flow isn't saved or takes inputs without a literal default. */
  private stampCreationScripts(tabs: OutputTab[]): void {
    if (this.boundScript == null) return;
    const nonLiteral = this.flow.getNodes().some((n) =>
      n.dgNodeType === 'input' && FuncFlowView.NON_SYNCABLE_INPUTS.has(n.dgTypeName ?? ''));
    if (nonLiteral) return;
    let callStr: string;
    try {
      // Serialized by the platform itself (`prepare().toString()` — the same
      // source of truth the creation-script emitter uses). Optional inputs at
      // their defaults are omitted by the serializer.
      const params: Record<string, unknown> = {};
      for (const n of this.flow.getNodes()) {
        if (n.dgNodeType !== 'input') continue;
        const pname = String(n.properties['paramName'] ?? '').trim();
        const def = n.properties['defaultValue'];
        if (pname !== '' && def !== undefined && def !== '') params[pname] = def;
      }
      callStr = this.boundScript.prepare(params).toString();
    } catch (e) {
      console.warn('FuncFlow: could not serialize the producing call', e);
      return;
    }
    const accessor = this.flowOutputCount() > 1;
    let ts = Date.now();
    for (const tab of tabs) {
      const df = tab.df!;
      if (!df.name) df.name = tab.paramName;
      df.setTag(DG.Tags.VariableName, tab.paramName);
      df.setTag(DG.Tags.CreationScript,
        `${tab.paramName} = ${callStr}${accessor ? '.' + tab.paramName : ''} //{"timestamp": ${ts++}}`);
    }
  }

  /** The platform's standard Save-project dialog, seeded with the computed
   *  output tables and their tab views — data-sync toggles, dependency
   *  handling, layout linking, upload, and sharing all come from core
   *  (`DG.Project.showSaveDialog`). Tabs never opened this session ship their
   *  stored layout as a view state string (a layout saved with the flow
   *  applies even without visiting the tab). A previously published project is
   *  passed back so re-publishing UPDATES it; the binding persists in the
   *  flow (saved silently right after, since the script was just saved). */
  private async openDashboardDialog(): Promise<void> {
    const tabs = this.outputViews?.getTabs().filter((t) => t.df != null) ?? [];
    if (tabs.length === 0) return;
    this.stampCreationScripts(tabs);
    const showSaveDialog = (DG.Project as unknown as {
      showSaveDialog?: (o: object) => Promise<unknown>;
    }).showSaveDialog;
    if (showSaveDialog == null) {
      grok.shell.warning('This platform version cannot publish dashboards from Flow yet');
      return;
    }
    const layoutsByParam = this.outputViews.captureLayouts();
    try {
      const saved = await showSaveDialog.call(DG.Project, {
        tables: tabs.map((t) => t.df!),
        views: tabs.map((t) => t.tv),
        layouts: tabs.map((t) => t.tv != null ? null : layoutsByParam[t.paramName] ?? null),
        name: this.flowSettings.scriptName,
        description: this.flowSettings.scriptDescription,
        project: this.dashboardProjectId ?? undefined,
      }) as DG.Project | null;
      if (saved?.id && saved.id !== this.dashboardProjectId) {
        this.dashboardProjectId = saved.id;
        if (this.boundScript != null)
          await this.saveToServer(true); // persist the binding in the entity body
      }
    } catch (e) {
      grok.shell.error(`Publish failed: ${e instanceof Error ? e.message : e}`);
    }
  }

  /** Download the graph as a local `.ffjson` file (the pre-entity behavior). */
  private async exportFfjson(): Promise<void> {
    // Best-effort: an export with pending (in-memory) file ids would not be
    // portable — persist them first; on failure export anyway with a warning.
    try {
      await this.persistPendingUploads();
    } catch (e: any) {
      grok.shell.warning(`Uploaded files were not persisted — the export will not be portable. ${e?.message ?? e}`);
    }
    const doc = serializeFlow(this.flow, this.flowSettings);
    downloadFlow(doc);
    grok.shell.info('Flow exported as .ffjson');
  }

  private editSettings(): void {
    const nameInput = ui.input.string('Script Name', {value: this.flowSettings.scriptName});
    const descInput = ui.input.string('Description', {value: this.flowSettings.scriptDescription});
    const tagsInput = ui.input.string('Tags', {value: this.flowSettings.tags.join(', ')});

    ui.dialog({title: 'Flow Settings'})
      .add(nameInput).add(descInput).add(tagsInput)
      .onOK(() => {
        this.flowSettings.scriptName = nameInput.value;
        this.flowSettings.scriptDescription = descInput.value;
        this.flowSettings.tags = tagsInput.value.split(',').map((s: string) => s.trim()).filter(Boolean);
      })
      .show();
  }

  // ---------- run / debug ----------

  private runInstrumented(): void {
    this.executionController?.runInstrumented({
      name: this.flowSettings.scriptName,
      description: this.flowSettings.scriptDescription,
      tags: this.flowSettings.tags,
    });
  }

  private debugInstrumented(): void {
    this.executionController?.debugInstrumented({
      name: this.flowSettings.scriptName,
      description: this.flowSettings.scriptDescription,
      tags: this.flowSettings.tags,
    });
  }

  private generateScript(): string | null {
    const errors = validateGraph(this.flow);
    if (errors.some((e) => e.severity === 'error')) {
      const msgs = errors.filter((e) => e.severity === 'error').map((e) => e.message).join('\n');
      grok.shell.error('Validation errors:\n' + msgs);
      return null;
    }
    try {
      return emitScript(this.flow, {
        name: this.flowSettings.scriptName,
        description: this.flowSettings.scriptDescription,
        tags: this.flowSettings.tags,
      });
    } catch (e: any) {
      grok.shell.error(`Script generation failed: ${e.message}`);
      return null;
    }
  }

  private runScript(): void {
    const script = this.generateScript();
    if (!script) return;
    // Classic (non-instrumented) run — no per-node values are captured, so
    // the click-to-inspect docked panel doesn't apply here. Outputs are
    // surfaced by Datagrok's standard script-run dialog.
    const func = DG.Script.create(script);
    const fc = func.prepare();
    if (func.inputs.length === 0)
      void fc.call(undefined, undefined, {processed: true});
    else {
      fc.getEditor(false).then((e: HTMLElement) => {
        ui.dialog({title: func.friendlyName ?? func.name}).add(e).show().onOK(async () => {
          await fc.call(undefined, undefined, {processed: true});
        });
      });
    }
  }

  private generateAndPreview(): void {
    const script = this.generateScript();
    if (!script) return;
    const pre = document.createElement('pre');
    pre.className = 'funcflow-script-preview';
    pre.textContent = script;
    const d = ui.dialog({title: 'Generated Script'})
      .add(pre)
      .addButton('Copy to Clipboard', () => {
        navigator.clipboard.writeText(script);
        grok.shell.info('Script copied to clipboard');
      })
      .addButton('Export .js', () => this.downloadScriptAsJs(script))
      .addButton('Open in Script View', () => {
        const sv = DG.ScriptView.create(DG.Script.create(script));
        grok.shell.addView(sv);
        d.close();
      })
      .addButton('Run', () => {
        DG.Script.create(script).prepare().edit();
        d.close();
      })
      .show({width: 700, height: 500});
  }

  private copyScriptToClipboard(): void {
    const script = this.generateScript();
    if (!script) return;
    navigator.clipboard.writeText(script);
    grok.shell.info('Script copied to clipboard');
  }

  /** Compile the graph into a Datagrok **creation script** (the grok-language
   *  function-call cascade, the inverse of "Import Creation Script") and show it
   *  in a dialog with any warnings about nodes that have no creation-script form. */
  private compileToCreationScript(): void {
    if (!this.flow) return;
    let result;
    try {
      result = emitCreationScript(this.flow);
    } catch (e: any) {
      grok.shell.error(`Creation script generation failed: ${e.message}`);
      return;
    }
    const {script, warnings} = result;

    const blocks: HTMLElement[] = [];
    if (warnings.length > 0) {
      const list = ui.divV(warnings.map((m) => ui.divText(`• ${m}`)));
      list.style.color = '#b26a00';
      list.style.marginBottom = '8px';
      list.style.maxHeight = '120px';
      list.style.overflow = 'auto';
      blocks.push(ui.divText(`${warnings.length} warning(s) — these nodes have no creation-script form:`,
        {style: {fontWeight: 'bold', color: '#b26a00'}}));
      blocks.push(list);
    }
    const pre = document.createElement('pre');
    pre.className = 'funcflow-script-preview';
    pre.textContent = script || '// (nothing to emit)';
    blocks.push(pre);

    ui.dialog({title: 'Creation Script'})
      .add(ui.divV(blocks))
      .addButton('Copy to Clipboard', () => {
        navigator.clipboard.writeText(script);
        grok.shell.info('Creation script copied to clipboard');
      })
      .show({width: 720, height: 520});
  }

  /** Compile the graph into a **separate** creation script per edited table
   *  (split by the variable each line builds), let the user review each in a
   *  horizontal tab, and on Save write it back to the table via
   *  `TableInfo.saveCreationScript`. Only available when the view was opened with
   *  `tableInfos` (the `creationScriptEditor` entry point). */
  private async saveCreationScriptsDialog(): Promise<void> {
    if (!this.flow || this.tableInfos.length === 0) return;

    // Creation scripts outlive the session — pending local files must be on
    // the server before the emitted script can reference them.
    try {
      await this.persistPendingUploads();
    } catch (e: any) {
      grok.shell.error(e?.message ?? String(e));
      return;
    }

    // The variable name a table is referenced by in the script — its
    // `.VariableName` tag, matching the SetVar/anchor names the emitter splits on
    // (falls back to the table's display name).
    const varNames = this.tableInfos.map((ti) =>
      String(ti.tags[DG.Tags.VariableName] ?? '').trim() || ti.name);

    let result;
    try {
      result = emitCreationScriptsForTables(this.flow, varNames);
    } catch (e: any) {
      grok.shell.error(`Creation script generation failed: ${e.message}`);
      return;
    }
    const {tables, warnings} = result;

    // One horizontal tab per table, each showing its standalone creation script.
    const tabs = ui.tabControl();
    this.tableInfos.forEach((ti, i) => {
      tabs.addPane(ti.name, () => {
        const pre = document.createElement('pre');
        pre.className = 'funcflow-script-preview';
        pre.style.whiteSpace = 'pre-wrap';
        pre.style.wordBreak = 'break-word';
        pre.style.height = '360px';
        pre.style.margin = '0';
        pre.style.userSelect = 'text';
        // An empty script is legitimate (e.g. a table updated locally, with no
        // recorded creation steps) — show it plainly, not as a warning.
        pre.textContent = tables[i].script || '// No creation script for this table.';
        return pre;
      });
    });
    tabs.root.style.width = '100%';
    tabs.root.style.height = 'unset';

    const blocks: HTMLElement[] = [];
    if (warnings.length > 0) {
      const list = ui.divV(warnings.map((m) => ui.divText(`• ${m}`)));
      list.style.color = '#b26a00';
      list.style.maxHeight = '90px';
      list.style.overflow = 'auto';
      list.style.marginBottom = '8px';
      blocks.push(ui.divText(`${warnings.length} warning(s):`,
        {style: {fontWeight: 'bold', color: '#b26a00'}}));
      blocks.push(list);
    }
    blocks.push(tabs.root);

    const dlg = ui.dialog({title: 'Save Creation Scripts'})
      .add(ui.divV(blocks))
      .onOK(async () => {
        try {
          await Promise.all(this.tableInfos.map((ti, i) => ti.saveCreationScript(tables[i].script)));
          grok.shell.info(`Saved creation script for ${this.tableInfos.length} table(s)`);
        } catch (e: any) {
          grok.shell.error(`Failed to save creation scripts: ${e.message}`);
        }
      });
    dlg.getButton('OK').textContent = 'Save';
    dlg.show({width: 760, height: 560});
  }

  private exportAsJs(): void {
    const script = this.generateScript();
    if (script) this.downloadScriptAsJs(script);
  }

  private downloadScriptAsJs(script: string): void {
    const blob = new Blob([script], {type: 'text/javascript'});
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `${this.flowSettings.scriptName || 'script'}.js`;
    a.click();
    URL.revokeObjectURL(url);
  }

  /** Load a flow from a JSON string (file viewer entry point). */
  async loadFromJson(json: string): Promise<void> {
    await this.loadFromDoc(JSON.parse(json) as FuncFlowDocument);
  }

  /** Load a parsed flow document. Awaits editor construction, so it is safe
   *  to call right after the constructor (no timer race). */
  async loadFromDoc(doc: FuncFlowDocument): Promise<void> {
    await this.editorReady;
    await deserializeFlow(doc, this.flow);
    if (doc.metadata?.settings) this.flowSettings = doc.metadata.settings;
    // Output-view tabs: rebuild the tab set from the fresh graph, then stash
    // the saved layouts (keyed by paramName — node ids were just remapped);
    // each applies once its tab is activated with a value.
    this.outputViews?.syncTabs(this.tableOutputs());
    if (doc.outputViews) {
      this.outputViews?.setPendingLayouts(Object.fromEntries(
        Object.entries(doc.outputViews).map(([name, v]) => [name, v.layout])));
    }
    this.dashboardProjectId = doc.dashboard?.projectId ?? null;
    this.name = doc.name || 'FuncFlow';
    this.updateStatusBar();
    this.fitToScreen();
    // Live-by-default nodes of a freshly loaded flow (Open File, …) run at
    // once; anything with unsatisfied inputs is dropped at fire time.
    this.autorunScheduler?.kickLive(this.flow.getNodes().map((n) => n.id));
  }

  /** Fit the whole graph into the viewport, so an opened flow is always shown
   *  fitted. Load paths usually run before the view is attached to the shell
   *  (file viewers, `forScript`, the creation-script dialog), when the canvas
   *  is still 0×0 and zooming would target a degenerate viewport — in that
   *  case the fit is deferred to the canvas's first real layout. Never blocks
   *  the load itself: a view that is never shown just never fits. */
  public fitToScreen(): void {
    this.pendingFitObserver?.disconnect();
    this.pendingFitObserver = null;
    if (!this.flow || this.flow.getNodeCount() === 0) return;
    const el = this.canvasContainer;
    if (el.clientWidth > 0 && el.clientHeight > 0) {
      void this.flow.zoomToFit();
      return;
    }
    const observer = new ResizeObserver(() => {
      if (el.clientWidth === 0 || el.clientHeight === 0) return;
      observer.disconnect();
      if (this.pendingFitObserver === observer) this.pendingFitObserver = null;
      // The graph may have been cleared while the fit was pending.
      if (this.flow && this.flow.getNodeCount() > 0) void this.flow.zoomToFit();
    });
    this.pendingFitObserver = observer;
    observer.observe(el);
  }

  /** Rebuild the canvas from a table-creation script — the cascade of
   *  function calls Datagrok records for reproducibly-created tables
   *  (the script behind a project's data sync). Replaces the current graph. */
  async loadFromCreationScript(script: string): Promise<void> {
    await this.editorReady;
    try {
      await this.flow.clear();
      const result = await buildFlowFromCreationScript(this.flow, script);
      this.updateStatusBar();
      this.fitToScreen();
      for (const warning of result.warnings) grok.shell.warning(warning);
      grok.shell.info(`Flow imported: ${result.nodesAdded} nodes, ${result.connectionsAdded} connections`);
    } catch (e: any) {
      grok.shell.error(`Creation script import failed: ${e?.message ?? e}`);
    }
  }

  /** Dialog: paste a creation script (or prefill it from an open table that
   *  carries one) and rebuild the canvas from it. */
  private importCreationScriptDialog(): void {
    const scriptInput = ui.input.textArea('Script', {value: ''});
    scriptInput.input.style.minHeight = '180px';
    scriptInput.input.style.minWidth = '420px';
    (scriptInput.input as HTMLTextAreaElement).placeholder =
      'Mol1K = OpenFile("System:AppData/Chem/mol1K.csv")\nChem:addChemPropertiesColumns(Mol1K, "molecule", true, …)';

    const items: HTMLElement[] = [];
    const tables = grok.shell.tables.filter((t) => (t.getTag(DG.Tags.CreationScript) ?? '') !== '');
    if (tables.length > 0) {
      const tableInput = ui.input.choice<string>('From table', {
        items: tables.map((t) => t.name),
        onValueChanged: (name: string | null) => {
          const table = tables.find((t) => t.name === name);
          if (table) scriptInput.value = table.getTag(DG.Tags.CreationScript) ?? '';
        },
      });
      ui.tooltip.bind(tableInput.root, 'Prefill the script from an open table that has a creation script');
      items.push(tableInput.root);
    }
    items.push(scriptInput.root);

    ui.dialog({title: 'Import Creation Script'})
      .add(ui.divV(items))
      .onOK(() => {
        const script = scriptInput.value.trim();
        if (script === '') {
          grok.shell.warning('Creation script is empty');
          return;
        }
        void this.loadFromCreationScript(script);
      })
      .show({width: 560, height: 420});
  }

  /** Plain-language description of the whole flow: each disjoint pipeline shown
   *  as its own card with a top-to-bottom, numbered list of node captions (full,
   *  never truncated). */
  private describeFlow(): void {
    if (!this.flow) return;
    const summary = summarizeFlow(this.flow.getNodes(), this.flow.getConnections());

    if (summary.nodeCount === 0) {
      ui.dialog({title: 'Flow summary'})
        .add(ui.divText('This flow is empty — add some nodes first.'))
        .show({width: 460, height: 220});
      return;
    }

    const pipeWord = summary.pipelineCount === 1 ? 'pipeline' : 'independent pipelines';
    const header = ui.divText(
      `${summary.nodeCount} node${summary.nodeCount === 1 ? '' : 's'} in ` +
      `${summary.pipelineCount} ${pipeWord}`, 'ff-flow-summary-header');

    const cards = summary.pipelines.map((p, i) => {
      const steps = p.steps.map((s) => {
        const fromText = s.inputs.length === 0 ? '' :
          '← ' + s.inputs.map((inp) => inp.key ? `${inp.key} from step ${inp.from}` : `step ${inp.from}`).join(', ');
        const main = ui.divV([ui.divText(s.caption, 'ff-flow-step-text')], 'ff-flow-step-main');
        if (fromText) main.appendChild(ui.divText(fromText, 'ff-flow-step-from'));
        return ui.div([ui.divText(`${s.index}`, 'ff-flow-step-n'), main], 'ff-flow-step');
      });
      const title = summary.pipelineCount > 1 ?
        [ui.divText(`Pipeline ${i + 1}`, 'ff-flow-pipe-title')] : [];
      return ui.divV([...title, ui.divV(steps, 'ff-flow-steps')], 'ff-flow-pipe');
    });

    ui.dialog({title: 'Flow summary'})
      .add(ui.divV([header, ...cards], 'ff-flow-summary'))
      .show({width: 600, height: 460});
  }

  private showValidation(): void {
    const results = validateGraph(this.flow);
    if (results.length === 0) {
      grok.shell.info('Graph is valid!');
      return;
    }
    const items = results.map((r) => {
      const icon = r.severity === 'error' ? '!!' : '!';
      return ui.divText(`[${icon}] ${r.message}`);
    });
    ui.dialog({title: 'Validation Results'}).add(ui.divV(items)).show({width: 500, height: 400});
  }
}
