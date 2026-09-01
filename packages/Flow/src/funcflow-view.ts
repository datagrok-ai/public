/** Main FuncFlow view — the Datagrok ViewBase hosting the Rete canvas, ribbon, status bar, and panels. */

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
  createNode, ensureFuncNodeType, FuncInfo,
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
import {buildPreview, setPreviewCellFocusHandler, releasePreviewCellFocusHandler, hasRenderablePreview}
  from './execution/value-inspector';
import {SuggestionPane, FF_SUGGEST_MIME} from './panel/suggestion-pane';
import {
  collectSuggestContext, computeSuggestions, Suggestion, CellSignal,
} from './suggest/suggestion-engine';
import {_package} from './package';
import {setTid} from './utils/test-ids';
import {getFuncDisplayName} from './utils/dart-proxy-utils';
import {
  addPendingFile, getPendingFile, isPendingFileId, persistPendingFile, removePendingFile,
  syncFlowFilePermissions,
} from './utils/uploaded-files';
import {GuideHost} from './guide/guide-model';
import {GuideRunner} from './guide/guide-runner';
import {createHelpButton, openGuideMenu} from './guide/guide-launcher';
import {TUTORIALS} from './guide/guide-content';
import {summarizeFlow} from './summary/summary-generator';
import {FlowAIContext} from './ai-tools';

const FLOW_TEMPLATES: {label: string; file: string; desc: string}[] = [
  {label: 'Workflow demo', file: 'Workflow Demo.flow', desc: 'A sample multi-step data workflow.'},
  {label: 'Interactive viewers', file: 'Interactive Viewers.flow',
    desc: 'Live in-node viewers driven by a molecule sketcher.'},
];

export class FuncFlowView extends DG.ViewBase {
  private flow!: FlowEditor;
  private functionBrowser!: FunctionBrowser;
  private propertyPanel!: PropertyPanel;
  private executionController!: ExecutionController;
  /** Public for tests. */
  suggestionPane!: SuggestionPane;
  /** The node last clicked — outlives a deselect. */
  private focusNodeId: string | null = null;
  private previewCell: CellSignal | null = null;
  /** Public for hosts and tests; created disabled when `options.outputPanel` is false. */
  outputPreview!: OutputPreviewPanel;
  private readonly outputPanelEnabled: boolean;
  private autorunScheduler: AutorunScheduler | null = null;
  private autorunIcon: HTMLElement | null = null;
  private saveButton: HTMLButtonElement | null = null;
  private savedSnapshot: string | null = null;
  private platformSubs: {unsubscribe(): void}[] = [];
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
  /** Public for tests. */
  outputViews!: OutputViewsManager;
  private viewTabStrip!: HTMLElement;
  /** Kept by reference — `setRibbonPanels` MOVES elements. */
  private flowRibbonPanels: HTMLElement[][] = [];

  private minimapCollapsed = false;

  private readonly tableInfos: DG.TableInfo[];

  private flowSettings: FlowSettings = {
    scriptName: 'My flow',
    scriptDescription: '',
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

  private boundScript: DG.Script | null = null;

  private dashboardProjectId: string | null = null;

  /** Resolves once `initEditor()` built the editor; load paths await it instead of retrying on timers. */
  private editorReadyResolve!: () => void;
  private readonly editorReady = new Promise<void>((resolve) => {
    this.editorReadyResolve = resolve;
  });

  private pendingFitObserver: ResizeObserver | null = null;

  private detached = false;
  private initEditorTimer: ReturnType<typeof setTimeout> | null = null;
  /** Kept so detach() releases only this view's registration. */
  private previewCellHandler:
    ((cell: {semType: string | null; column: string; value: unknown}) => void) | null = null;
  /** Serializes canvas loads — two interleaved clear()+add sequences would merge both graphs. */
  private loadChain: Promise<void> = Promise.resolve();

  constructor(tableInfos: DG.TableInfo[] = [], options: {outputPanel?: boolean} = {}) {
    super();
    this.name = 'Flow';
    this.aiDescription = 'Flow — Datagrok\'s visual pipeline editor: the user composes functions, ' +
      'queries, and scripts into an executable graph on a canvas. Act on it through the view ' +
      'functions (search list_view_functions with the word "flow"): listFlowNodes (call first to ' +
      'see the canvas), getFlowNodeDetails, findFlowNodeTypes, addFlowNode, connectFlowNodes, ' +
      'setFlowNodeInputs, selectFlowNode, runFlow; listFlowGuides / startFlowGuide launch ' +
      'interactive tutorials — offer one when the user asks how to do something here. When user asks about how to do something, ' +
      'Make sure to list guides (run the function listing the guides) to see if there are guides that cover the question, and if so and confidence is high, call that guide function';
    this.tableInfos = tableInfos;
    this.outputPanelEnabled = options.outputPanel !== false;
    // Test-harness hook — lives on this view's root so it dies with the DOM (never a global).
    (this.root as unknown as {__ffView?: FuncFlowView}).__ffView = this;

    registerBuiltinNodes();

    this.initUI();
    this.setupRibbon();
    this.setupStatusBar();

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

    // Suppress the platform's "drop to open" overlay so the canvas receives OS file drags.
    this.platformSubs.push(grok.events.onFileDragEnter.subscribe((ev) => {
      const target = ev.causedBy?.target;
      if (target !== null && target instanceof HTMLElement &&  (target === this.canvasContainer || this.canvasContainer.contains(target)))
        ev.preventDefault();
    }));
  }

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
      void this.loadFromDoc(doc).then(() => this.markSaved())
        .catch((e) => grok.shell.error(`Cannot load flow "${this.name}": ${e instanceof Error ? e.message : e}`));
    } catch (e) {
      grok.shell.error(`Cannot read flow "${this.name}": ${e instanceof Error ? e.message : e}`);
    }
  }

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

    this.suggestionPane = new SuggestionPane(
      async () => this.flow ?
        computeSuggestions(await collectSuggestContext(
          this.flow, this.executionController ?? null, this.focusNodeId, this.previewCell)) :
        [],
      (s) => void this.applySuggestion(s),
    );
    this.functionBrowser.root.appendChild(this.suggestionPane.root);
    this.previewCellHandler = (cell): void => {
      this.previewCell = cell;
      this.suggestionPane.refresh();
    };
    setPreviewCellFocusHandler(this.previewCellHandler);

    // Core css forces `overflow: auto !important` on `.ui-box` children; funcflow.css wins it back
    // for this class with a higher-specificity `overflow: hidden !important`.
    this.canvasContainer = ui.div([], 'funcflow-canvas-container');
    setTid(this.canvasContainer, 'canvas');
    this.startPanel = this.buildStartPanel();
    this.canvasContainer.appendChild(this.startPanel);

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
    this.viewTabStrip = ui.div([], 'ff-view-tabs');
    setTid(this.viewTabStrip, 'view-tabs');
    const statusLabels = ui.div(
      [this.nodeCountLabel, this.linkCountLabel, this.validationLabel], 'ff-status-labels');
    this.statusBar = ui.div([this.viewTabStrip, statusLabels], 'funcflow-status-bar');
    setTid(this.statusBar, 'statusbar');

    this.outputPreview = new OutputPreviewPanel({enabled: this.outputPanelEnabled});
    const canvasBox = ui.box(this.canvasContainer);
    // The panel's explicit height (its pane is `flex: 0 0 auto`) is the single source of truth.
    canvasBox.style.flex = '1 1 0';
    canvasBox.style.minHeight = '0';
    const split = ui.splitV([canvasBox, this.outputPreview.root], {style: {flex: '1 1 0', width: '100%', height: '100%'}}, true);
    const divider = split.querySelector('.ui-split-v-divider') as HTMLElement | null;
    const syncDivider = (state: OutputPanelState): void => {
      if (divider) divider.style.display = state === 'expanded' ? '' : 'none';
    };
    // Minimize the minimap only on the hidden → visible edge — later re-expands are the user's choice.
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
    // Panes are display-toggled, never unmounted — hosted TableView DOM must persist.
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

    this.initEditorTimer = setTimeout(() => {
      this.initEditorTimer = null;
      if (!this.detached) this.initEditor();
    }, 50);
  }

  /** Capture phase — intercepts the contextmenu before AreaPlugin routes it to the node menu. */
  private installPortContextMenu(): void {
    this.canvasContainer.addEventListener('contextmenu', (ev) => {
      const target = ev.target as HTMLElement | null;
      // No output value on input sockets — also suppresses the area-plugin's "Add annotation here" item.
      if (target?.closest('.ff-socket-row-input')) {
        ev.preventDefault();
        ev.stopPropagation();
        return;
      }
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
      menu.item(summary ? 'Re-run up to here' : 'Run up to here & preview',
        () => this.previewNodeData(nodeId));
      menu.show({causedBy: ev});
    }, true);
  }

  private viewerEditSub: {unsubscribe(): void} | undefined;

  /** Shows the live viewer in the context panel and persists its option changes into the node's `viewerLook`. */
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

  private previewNodeData(nodeId: string): void {
    this.executionController?.previewNodeData(nodeId, {
      name: this.flowSettings.scriptName,
      description: this.flowSettings.scriptDescription,
      tags: this.flowSettings.tags,
    });
  }

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
        this.executionController?.showOutputsForNode(node);
        this.focusNodeId = node.id;
        this.previewCell = null;
        this.suggestionPane?.refresh();
      },
      onNodeDeselected: () => {
        this.propertyPanel.clear();
        grok.shell.o = this.propertyPanel.root;
        this.suggestionPane?.refresh();
      },
      // A re-click must restore a stale context (a tab switch replaced grok.shell.o, or an unshown preview).
      isNodeContextCurrent: (node: FlowNode) => {
        if (this.propertyPanel.shownNodeId !== node.id) return false;
        try {
          if (grok.shell.o !== this.propertyPanel.root) return false;
        } catch {/* current-object read failed — treat as stale and re-show */
          return false;
        }
        const preview = this.outputPreview;
        if (preview.isEnabled && preview.pinnedNodeId == null && preview.currentNodeId !== node.id) {
          const state = this.executionController?.state.getNodeState(node.id);
          if (state != null && hasRenderablePreview(state)) return false;
        }
        return true;
      },
      // Selection changes that bypass `nodepicked` — marquee, Ctrl+A, modifier-click removals.
      onSelectionChanged: () => this.suggestionPane?.refresh(),
      onGraphChanged: () => {
        this.updateStatusBar();
        this.updateStartPanelVisibility();
        this.refreshNodeHints();
        this.flow?.refreshMinimap();
        this.updateSaveButtonState();
        this.suggestionPane?.refresh();
        this.outputViews?.syncTabs(this.tableOutputs());
        this.propertyPanel.refreshShownNode();
        this.updateAutorunIndicator();
      },
      onGraphEdited: (edit) => {
        const affected = this.executionController?.applyGraphEdit(edit) ?? new Set<string>();
        this.autorunScheduler?.onEdit(edit, affected);
        // `onGraphChanged` does not fire on parameter edits — refresh Save state here too.
        this.updateSaveButtonState();
        this.suggestionPane?.refresh();
        // paramName renames arrive as params-changed — the tab set/labels must track them here.
        this.outputViews?.syncTabs(this.tableOutputs());
        if (edit.kind === 'params-changed') {
          this.propertyPanel.refreshShownNode();
          // params-changed doesn't fire onGraphChanged — recompute the "Needs input" hint here.
          this.refreshNodeHints();
        }
        this.updateAutorunIndicator();
      },
      onPreviewNode: (nodeId: string) => this.previewNodeData(nodeId),
      onRerunNode: (nodeId: string) => this.rerunNode(nodeId),
      // In-node preview: the node component mounts the captured live viewer/widget root.
      getInlinePreviewContent: (nodeId: string) =>
        this.executionController?.inlinePreviewRoot(nodeId) ?? null,
      isInlinePreviewPending: (nodeId: string) =>
        // eslint-disable-next-line @typescript-eslint/no-unsafe-return
        (this.executionController?.inlinePreviewPending(nodeId) ?? false),
      onInlinePreviewToggled: (nodeId: string) => {
        this.executionController?.syncInlinePreviewOwnership(nodeId);
        // The bottom panel must swap between the live root and the hosted-note.
        if (this.outputPreview.currentNodeId === nodeId) this.outputPreview.refresh();
      },
      canRerunNode: (nodeId: string) => this.executionController?.canRerunNode(nodeId) ?? false,
      // Only suggestions wired FROM the drag-source node apply to a socket drag.
      getSocketSuggestions: async (nodeId: string) => {
        if (!this.flow) return [];
        const ctx = await collectSuggestContext(this.flow, this.executionController ?? null, nodeId, null);
        return computeSuggestions(ctx)
          .filter((s) => s.wire.some((w) => w.fromNodeId === nodeId))
          .map((s) => ({typeName: s.typeName, reason: s.reason, prefill: s.prefill}));
      },
    });

    this.propertyPanel = new PropertyPanel(this.flow);
    this.executionController = new ExecutionController(this.flow, this.outputPreview);
    this.executionController.onNodeStateChanged = (nodeId) => {
      this.suggestionPane?.refresh();
      this.updateOutputViewValue(nodeId);
      this.refreshOpenFileTitle(nodeId);
      this.updateAutorunIndicator();
      this.propertyPanel.updateExecState(nodeId, this.executionController?.state.getNodeState(nodeId));
    };
    // Empty-canvas deselects happen inside rete with no callback — re-read the selection on any release.
    this.canvasContainer.addEventListener('pointerup', () => this.suggestionPane?.refresh());
    this.autorunScheduler = new AutorunScheduler((dirty, liveOnly) => {
      const settings = {
        name: this.flowSettings.scriptName,
        description: this.flowSettings.scriptDescription,
        tags: this.flowSettings.tags,
      };
      return (liveOnly ?
        this.executionController?.runLiveNodes(dirty, settings) :
        this.executionController?.runAutorun(dirty, settings)) ?? 'skipped';
    },
    AUTORUN_DEBOUNCE_MS,
    (id) => {
      const n = this.flow?.getNodeById(id);
      return n != null && isAutorunByDefault(n);
    });

    const columnPicker = new ColumnPicker(this.flow, this.executionController, () => ({
      name: this.flowSettings.scriptName,
      description: this.flowSettings.scriptDescription,
      tags: this.flowSettings.tags,
    }));
    this.propertyPanel.onPickColumns = (req) => void columnPicker.pick(req);
    // Captured-only read — a panel render must never kick off work.
    this.propertyPanel.getUpstreamTable = (sourceNodeId) =>
      this.executionController?.cloneForNode(sourceNodeId) ?? null;
    this.propertyPanel.runUpstreamNode = (sourceNodeId) =>
      this.executionController.produceTableForNode(sourceNodeId, {
        name: this.flowSettings.scriptName,
        description: this.flowSettings.scriptDescription,
        tags: this.flowSettings.tags,
      });

    const funcEditorLauncher = new FuncEditorLauncher(this.flow, this.executionController, () => ({
      name: this.flowSettings.scriptName,
      description: this.flowSettings.scriptDescription,
      tags: this.flowSettings.tags,
    }));
    this.propertyPanel.onEditFuncParams = (node) => {
      // Hold autorun — a mid-dialog rerun of the same function would be mistaken for the dialog's own run action.
      this.autorunScheduler?.hold();
      void funcEditorLauncher.open(node).then((applied) => {
        if (applied) {
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
      this.saveDialogRunEnd?.(success);
      if (!success) return;
      this.autoSelectFirstOutputNode();
    };

    this.outputPreview.onEditViewer = (nodeRef, viewer) => this.editViewer(nodeRef.id, viewer);

    this.flow.setMinimapCollapsed(this.minimapCollapsed);
    this.updateStartPanelVisibility();
    this.markSaved();
    this.editorReadyResolve();
  }

  private hintRaf = 0;
  private refreshNodeHints(): void {
    if (!this.flow || this.hintRaf) return;
    this.hintRaf = requestAnimationFrame(() => {
      this.hintRaf = 0;
      for (const n of this.flow.getNodes()) void this.flow.updateNode(n.id);
    });
  }

  enableOutputPanel(): void {
    this.outputPreview.setEnabled(true);
  }

  setMinimapCollapsed(collapsed: boolean): void {
    this.minimapCollapsed = collapsed;
    this.flow?.setMinimapCollapsed(collapsed);
  }

  /** SetVar terminals count as outputs — they compile to the same thing. */
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

  /** Table-carrying terminals, one tab each: Table Outputs, dataframe Value Outputs, SetVar-with-table. */
  private tableOutputs(): OutputTabInfo[] {
    if (!this.flow) return [];
    const res: OutputTabInfo[] = [];
    for (const n of this.flow.getOutputNodes()) {
      const isTable = n.dgTypeName === 'Outputs/Table Output' ||
        (n.dgTypeName === 'Outputs/Value Output' && n.properties['outputType'] === 'dataframe');
      if (isTable)
        res.push({nodeId: n.id, paramName: String(n.properties['paramName'] ?? '').trim() || n.label});
    }
    for (const n of this.flow.getNodes()) {
      if (isSetVarNode(n)) {
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

  /** Restamps an Open File node's title from its path — only after a completed run. */
  private refreshOpenFileTitle(nodeId: string): void {
    const node = this.flow?.getNodeById(nodeId);
    if (!node || (node.dgFuncName ?? '').toLowerCase() !== 'openfile') return;
    if (this.executionController?.state.getNodeState(nodeId)?.status !== NodeExecStatus.completed) return;
    const path = String(node.inputValues['fullPath'] ?? '');
    if (path === '') return;
    // Restamp only machine-stamped titles — a user-renamed node keeps its name.
    const base = (node.dgFunc ? getFuncDisplayName(node.dgFunc) : '') || 'Open File';
    if (node.label !== base && !node.label.startsWith(`${base}: `)) return;
    const fileName = path.split('/').pop() || path;
    const label = `${base}: ${fileName}`;
    if (node.label === label) return;
    node.label = label;
    void this.flow.updateNode(node.id);
  }

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

  /** Surfaces the active tab's ribbon/toolbox on this view (the `DG.MultiView.currentView` recipe). */
  private onOutputTabChanged(tab: OutputTab | null): void {
    const tv = tab?.tv ?? null;
    try {
      if (tv == null) {
        if (this.flowRibbonPanels.length > 0) this.setRibbonPanels(this.flowRibbonPanels);
        this.toolbox = this.functionBrowser.root;
        return;
      }
      // Capture ONCE, unwrapped — `setRibbonPanels` moves elements, so re-reading `getRibbonPanels()`
      // after a swap returns empty husks. Core-hidden items are dropped; Flow's Save pill leads.
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

  private captureOutputViews(): FuncFlowDocument['outputViews'] {
    const layouts = this.outputViews?.captureLayouts() ?? {};
    const names = Object.keys(layouts);
    if (names.length === 0) return undefined;
    const res: NonNullable<FuncFlowDocument['outputViews']> = {};
    for (const name of names) res[name] = {layout: layouts[name]};
    return res;
  }

  private setupFileDropTarget(): void {
    ui.makeDroppable(this.canvasContainer, {
      acceptDrop: (drag: any) =>
        (drag instanceof DG.FileInfo && drag.isFile) || drag instanceof DG.Func,
      doDrop: (args) => {
        const drag = args.dragObject;
        const ev = args.dropEvent;
        if (drag instanceof DG.Func) {void this.addFuncNode(drag, ev); return;}
        const fi = drag as DG.FileInfo;
        void this.addOpenFileNode(fi.fullPath, ev);
      },
    });

    this.canvasContainer.addEventListener('dragover', (ev) => {
      if (!ev.dataTransfer) return;
      const types = Array.from(ev.dataTransfer.types);
      if (types.includes(FF_DRAG_MIME) || types.includes(FF_SUGGEST_MIME) || types.includes('Files')) {
        ev.preventDefault();
        ev.dataTransfer.dropEffect = 'copy';
      }
    });
    this.canvasContainer.addEventListener('drop', (ev) => {
      const osFiles = ev.dataTransfer?.files;
      if (osFiles != null && osFiles.length > 0) {
        ev.preventDefault();
        ev.stopPropagation();
        void this.addUploadedFileNodes(Array.from(osFiles), ev);
        return;
      }
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

  /** Collected by the AI assistant through `view.getFunctions()` (the Dart JsViewHost forwards here). */
  getFunctions(): DG.Func[] {
    return DG.Func.find({package: 'Flow', tags: ['flowViewFunction']});
  }

  /** Facade the registered AI view functions use to act on this instance. */
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

  /** Creates the suggested node (at the drop point, else next to its source), prefills and wires it. */
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
      this.flow.notifyNodeParamsChanged(node.id);
    }
    this.updateStatusBar();
    this.suggestionPane?.refresh();
  }

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
      const fileName = filePath.split('/').pop() || filePath;
      node.label = `${node.label}: ${fileName}`;
      await this.flow.updateNode(node.id);
      this.flow.notifyNodeParamsChanged(node.id);
    }
  }

  private async addFuncNode(func: DG.Func, dropEvent?: MouseEvent): Promise<void> {
    // Explicit intent — a func outside the catalog allowlist registers on the fly.
    const info = getRegisteredFuncs().find((f) => f.func.name === func.name);
    let typeName: string;
    try {
      typeName = info?.nodeTypeName ?? ensureFuncNodeType(func);
    } catch {
      grok.shell.warning(`Function "${func.name}" is not available as a node`);
      return;
    }
    await this.addNodeByTypeAtDrop(typeName, dropEvent);
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

  /** OS files dropped on the canvas: bytes go to the pending store, reaching the server only on save. */
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
      this.flow.notifyNodeParamsChanged(node.id);
      offset += 40;
    }
    this.updateStatusBar();
  }

  /** Uploads pending files and rewrites node fileIds to real ids before any serialization that outlives the session. */
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

  private async syncUploadedFilePermissions(): Promise<void> {
    if (this.boundScript?.script)
      await syncFlowFilePermissions(this.boundScript);
  }

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
      hideStartPanel: () => this.hideStartPanel(),
      anchorEl: this.helpButton,
    };
  }

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
      setTid(card, 'start-template', t.file.replace(/\.flow$/i, ''));
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
        // Capture phase — the entity markup's own click would swallow the open.
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

  private buildStartBackground(): HTMLElement {
    this.startBg = ui.div([], 'funcflow-start-bg');
    setTid(this.startBg, 'start-bg');
    return this.startBg;
  }

  /** Deterministic jittered graph sized to the real container — resizing reflows instead of reshuffling. */
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

  private drawStartBackgroundSoon(): void {
    if (this.startBgRaf) cancelAnimationFrame(this.startBgRaf);
    this.startBgRaf = requestAnimationFrame(() => {
      this.startBgRaf = 0;
      this.renderStartBackground();
    });
  }

  private showStartPanel(): void {
    this.startPanel.style.display = 'flex';
    this.drawStartBackgroundSoon();
    void this.refreshRecentFlows();
  }

  private hideStartPanel(): void {
    this.startPanel.style.display = 'none';
  }

  /** Annotations count as content — the overlay must not cover (and swallow clicks on) a lone note. */
  private updateStartPanelVisibility(): void {
    const wasHidden = this.startPanel.style.display === 'none';
    const empty = (!this.flow ||
      (this.flow.getNodeCount() === 0 && this.flow.getAnnotations().length === 0)) &&
      !this.guideRunner.isRunning;
    this.startPanel.style.display = empty ? 'flex' : 'none';
    if (empty) {
      this.drawStartBackgroundSoon();
      if (wasHidden) void this.refreshRecentFlows();
    }
  }

  private async loadTemplate(file: string): Promise<void> {
    try {
      const json = await _package.files.readAsText(file);
      await this.loadFromJson(json);
      this.hideStartPanel();
      grok.shell.info(`Opened template: ${file.replace(/\.flow$/i, '')}`);
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
      .item('Import .flow…', () => void this.openFlow())
      .item('Export .flow', () => void this.exportFlowFile())
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
    ui.tooltip.bind(autorunIcon, () => {
      if (!this.autorunScheduler?.enabled) {
        return 'Autorun is off — click to rerun the flow (only the affected nodes) automatically after every change. ' +
          'Live nodes (Open File, Add New Column, viewers) still rerun on change.';
      }
      const blockers = this.executionController?.autorunBlockers() ?? [];
      if (blockers.length === 0)
        return 'Autorun is on — the flow reruns the affected nodes after every change. Click to turn off.';
      return ui.divV([
        ui.divText('Autorun is on, but waiting for:'),
        ...blockers.map((b) => ui.divText('• ' + b)),
        ui.divText('Set the values on the nodes (or in the context panel), or press Run to be asked for them.',
          'ff-autorun-tooltip-hint'),
      ]);
    });
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
    // The mode saves with the flow (metadata.settings.autorun) so it reopens live.
    if (on) this.flowSettings.autorun = true;
    else delete this.flowSettings.autorun;
    this.updateSaveButtonState();
    if (on) {
      const pending = this.executionController?.pendingNodes() ?? new Set<string>();
      if (pending.size > 0) this.autorunScheduler.kick(pending);
    }
    this.updateAutorunIndicator();
  }

  /** Drive the autorun toggle to a target state (applying a loaded flow's saved
   *  mode) — turning it on kicks pending nodes exactly like a ribbon click. */
  private applyAutorunSetting(on: boolean): void {
    if (!this.autorunScheduler || this.autorunScheduler.enabled === on) return;
    this.toggleAutorun();
  }

  /** Amber "waiting" badge on the bolt while autorun is on but can't run what's
   *  pending (an input node without a value, validation errors) — the dynamic
   *  tooltip names the exact blockers. Cleared the moment nothing blocks. */
  private updateAutorunIndicator(): void {
    if (!this.autorunIcon) return;
    const blocked = (this.autorunScheduler?.enabled ?? false) &&
      (this.executionController?.autorunBlockers() ?? []).length > 0;
    this.autorunIcon.classList.toggle('ff-autorun-blocked', blocked);
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

  /** Record the saved baseline (after a save, or after loading a flow that
   *  already lives on the server / a fresh empty flow). Pass the snapshot
   *  taken when the payload was built so mid-save edits stay "unsaved". */
  private markSaved(snapshot?: string): void {
    this.savedSnapshot = snapshot ?? this.currentSnapshot();
    this.updateSaveButtonState();
  }

  /** Whether Save is available, and its tooltip. The button is the gateway to
   *  both saving the script AND publishing a dashboard, so it stays enabled
   *  for any non-empty canvas — a freshly opened (unchanged) flow must still
   *  open the dialog to publish; the tooltip reflects the dirty state. */
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
      setTimeout(() => {
        const cur = grok.shell.v;
        if (!this.dart || !cur || cur.dart !== this.dart) return;
        try {
          grok.shell.v.pin?.();
        } catch (e) {
          console.error(e);
        }
      }, 100);
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

  /** View closed — release everything that outlives the view's DOM: pending
   *  timers/rAFs, the deferred editor init, the guide, the run (its event
   *  subscription and captured live values), panel editors, and the shell's
   *  current-object pointer if it is ours. */
  detach(): void {
    this.detached = true;
    if (this.initEditorTimer != null) {
      clearTimeout(this.initEditorTimer);
      this.initEditorTimer = null;
    }
    if (this.hintRaf) {
      cancelAnimationFrame(this.hintRaf);
      this.hintRaf = 0;
    }
    if (this.startBgRaf) {
      cancelAnimationFrame(this.startBgRaf);
      this.startBgRaf = 0;
    }
    this.autorunScheduler?.reset();
    this.guideRunner.stop();
    this.teardownAutoPin();
    for (const sub of this.platformSubs) sub.unsubscribe();
    this.platformSubs = [];
    this.viewerEditSub?.unsubscribe();
    this.viewerEditSub = undefined;
    this.pendingFitObserver?.disconnect();
    this.pendingFitObserver = null;
    this.currentPortPopup?.remove();
    this.currentPortPopup = null;
    // Un-mount the last-shown node's custom editors (their detach() releases
    // widget subscriptions), and stop the module-level preview-cell hook from
    // feeding this dead view.
    this.propertyPanel?.clear();
    releasePreviewCellFocusHandler(this.previewCellHandler);
    // The shell must not keep rendering a destroyed view's panel/viewer.
    try {
      if (grok.shell.o === this.propertyPanel?.root) grok.shell.o = null;
    } catch {/* shell not available in this host */}
    // Bytes of files dropped but never saved would otherwise sit in the
    // page-global pending registry forever.
    this.dropPendingUploads();
    this.suggestionPane?.destroy();
    this.functionBrowser?.destroy();
    this.outputViews?.destroy();
    // Dispose BEFORE the editor: clearing this flow's live-registry entries
    // iterates the editor's nodes.
    this.executionController?.dispose();
    this.flow?.destroy();
    super.detach();
  }

  /** Remove this flow's not-yet-persisted uploaded files from the pending
   *  registry — the view owned them, and they are gone with it. */
  private dropPendingUploads(): void {
    if (!this.flow) return;
    for (const node of this.flow.getNodes()) {
      const fileId = String(node.inputValues['fileId'] ?? '');
      if (fileId.startsWith('pending:')) removePendingFile(fileId);
    }
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
    DG.Utils.openFile({
      accept: '.flow',
      open: async (file) => {
        try {
          const doc = await loadFlowFromFile(file);
          await this.loadFromDoc(doc);
          this.boundScript = null; // an imported file is a new, unsaved flow
          grok.shell.info(`Loaded flow: ${doc.name}`);
        } catch (e: any) {
          grok.shell.error(`Failed to load flow: ${e.message}`);
        }
      },
    });
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
   *  annotation header and the flow JSON payload can never disagree. */
  private entityBodyText(): string {
    if (!this.flowSettings.tags.includes(FLOW_TAG))
      this.flowSettings.tags.push(FLOW_TAG);
    return flowScriptText(this.flow, this.flowSettings, {
      outputViews: this.captureOutputViews(),
      dashboard: this.dashboardProjectId != null ? {projectId: this.dashboardProjectId} : undefined,
    });
  }

  /** @param silent no "saved" balloon — for follow-up saves that only persist
   *  metadata (e.g. the dashboard binding right after a publish).
   *  @returns whether the save reached the server — callers must not proceed
   *  to share/publish on false. */
  private async saveToServer(silent = false): Promise<boolean> {
    try {
      // Local files live in memory until this moment — persist them first so
      // the saved body references real, server-addressable file ids.
      await this.persistPendingUploads();
      const script = DG.Script.create(this.entityBodyText());
      // Snapshot what is actually in the payload BEFORE the server round-trip:
      // an edit made while the save is in flight must keep Save enabled, not be
      // silently blessed as "saved" by a post-await snapshot.
      const savedSnapshot = this.currentSnapshot();
      if (this.boundScript?.id)
        script.id = this.boundScript.id;
      const saved = await grok.dapi.scripts.save(script);
      // The UPDATE path (save with an existing id) returns an entity WITHOUT
      // its namespace — a call serialized from it reads `FlowName()`, which
      // breaks for anyone else once shared. Re-fetch so the bound entity
      // always carries the server-qualified name (`user:FlowName`) — what
      // `stampCreationScripts` embeds into dashboard creation scripts.
      this.boundScript = (await grok.dapi.scripts.find(saved.id).catch(() => null)) ?? saved;
      this.name = this.flowSettings.scriptName;
      this.updatePath();
      this.markSaved(savedSnapshot); // the payload we sent is the saved baseline
      void this.syncUploadedFilePermissions();
      if (!silent) grok.shell.info(`Flow "${this.flowSettings.scriptName}" saved`);
      return true;
    } catch (e: any) {
      grok.shell.error(`Failed to save flow: ${e?.message ?? e}`);
      return false;
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
    descInput.root.style.marginBottom = '6px';
    // Names aren't prose — the browser's squiggle under "MyFuncFlow" reads as
    // an error right next to the real (amber) warning line.
    nameInput.input.setAttribute('spellcheck', 'false');
    descInput.input.setAttribute('spellcheck', 'false');

    let targetSpace: DG.Project | null = null;
    // Reserves its line (min-height in CSS) so appearing text never shifts the form.
    const warningDiv = ui.divText('', 'funcflow-save-name-warning');

    // Best-effort duplicate-name check: purely advisory — the server keeps
    // names unique within a space by suffixing, and plain scripts may repeat.
    let warnTimer: ReturnType<typeof setTimeout> | null = null;
    const findClashes = async (name: string, space: DG.Project | null): Promise<DG.Script[]> => {
      const esc = name.replace(/"/g, '\\"');
      const scope = space ? ` and namespace = "${space.nqName}:"` : '';
      return (await grok.dapi.scripts
        .filter(`language = "${FLOW_LANGUAGE}" and friendlyName = "${esc}"${scope}`).list())
        .filter((s) => s.id !== this.boundScript?.id);
    };
    const refreshNameWarning = async () => {
      const name = nameInput.value.trim();
      warningDiv.textContent = '';
      if (name === '') return;
      try {
        // Short enough for one line — the reserved line must not wrap and shift the form.
        if ((await findClashes(name, targetSpace)).length > 0) {
          warningDiv.textContent = targetSpace ?
            `Already used in "${targetSpace.friendlyName}" — it will get a unique name` :
            'This name is already in use — the flow will be saved separately';
        }
      } catch {/* advisory only */}
    };

    // ---- space binding: a regular form row (caption / value / actions) ----
    const noSpaceText = 'None — saved to your scripts';
    const spaceValue = ui.divText(noSpaceText, 'ff-save-space-value ff-save-space-none');
    ui.tooltip.bind(spaceValue, 'Where this flow will live');
    const pickerHost = ui.div([], 'ff-save-space-picker-host');
    pickerHost.style.display = 'none';
    let picker: SpacePicker | null = null;
    const togglePicker = async () => {
      if (picker == null) {
        picker = await SpacePicker.create();
        picker.onChanged = (space) => updateSpace(space);
        pickerHost.appendChild(picker.root);
        pickerHost.style.display = '';
      } else
        pickerHost.style.display = pickerHost.style.display === 'none' ? '' : 'none';
    };
    const chooseLink = ui.link('Choose…', () => void togglePicker(),
      'Choose a space (or subspace) to organize and share this flow; by default it is saved as a plain script in your namespace');
    // Clearing goes through the picker when it exists, so its highlighted row
    // and footer label ("New space…" vs "New subspace…") stay in sync.
    const clearIcon = ui.iconFA('times',
      () => picker != null ? picker.clearSelection() : updateSpace(null),
      'Remove from space — save as a plain script');
    clearIcon.classList.add('ff-save-space-clear');
    const spaceCaption = document.createElement('label');
    spaceCaption.className = 'ui-label ui-input-label';
    spaceCaption.textContent = 'Space';
    // ui-input-root makes the caption inherit the dialog form's label column
    // (right-aligned, 140px, grey) so the row lines up with Name/Description.
    const spaceRow = ui.divH([spaceCaption, spaceValue, chooseLink, clearIcon],
      {classes: 'ui-input-root ff-save-space-row', style: {alignItems: 'center', gap: '8px'}});

    const updateSpace = (space: DG.Project | null) => {
      targetSpace = space;
      spaceValue.textContent = space ? space.friendlyName : noSpaceText;
      spaceValue.classList.toggle('ff-save-space-none', space == null);
      chooseLink.textContent = space ? 'Change…' : 'Choose…';
      clearIcon.style.display = space ? '' : 'none';
      void refreshNameWarning();
    };
    updateSpace(null);

    // ---- dashboard section (run-aware) ----
    const computedTabs = (): OutputTab[] =>
      this.outputViews?.getTabs().filter((t) => t.df != null) ?? [];
    const dashHost = ui.divV([], 'ff-save-dash');
    const publishInput = ui.input.bool('Create dashboard', {value: true,
      tooltipText: 'After saving the flow, open the standard Save-project dialog seeded with the ' +
        'computed output tables and their layouts (data sync, sharing, upload)'});
    // 'new' is applied only at save time, so the choice stays reversible.
    let dashMode: 'update' | 'new' = 'update';
    // The publish-mode combo only makes sense while publishing is on.
    let modeRoot: HTMLElement | null = null;
    publishInput.onChanged.subscribe(() => {
      if (modeRoot != null) modeRoot.style.display = publishInput.value ? '' : 'none';
      syncOk();
    });
    let runningNow = false;
    let okBtn: HTMLButtonElement | null = null;
    const syncOk = () => {
      if (okBtn == null) return;
      const empty = nameInput.value.trim() === '';
      const disabled = empty || runningNow;
      okBtn.disabled = disabled;
      okBtn.classList.toggle('disabled', disabled);
      // Forewarn the second step: with "Create dashboard" on, OK leads into
      // the platform's Save-project dialog rather than straight to done.
      const publishes = publishInput.value === true && computedTabs().length > 0;
      okBtn.textContent = publishes ? 'Save & publish…' : 'Save';
      okBtn.title = empty ? 'Give the flow a name' : runningNow ? 'Waiting for the run to finish' :
        publishes ? 'Saves the flow, then opens the dashboard-project setup' : '';
    };
    const refreshDash = (state: {running?: boolean, failed?: boolean, ranOnce?: boolean} = {}): void => {
      runningNow = state.running === true;
      syncOk();
      ui.empty(dashHost);
      dashHost.appendChild(ui.divText('Dashboard', 'ff-save-dash-title'));
      const tabs = computedTabs();
      if (tabs.length > 0) {
        dashHost.appendChild(ui.divText('Result tables:', 'ff-save-dash-hint'));
        dashHost.appendChild(ui.divV(tabs.map((t) => ui.divH([
          ui.iconFA('table'),
          ui.divText(t.df!.name || t.paramName),
          ui.divText(`${t.df!.rowCount.toLocaleString()} × ${t.df!.columns.length}`, 'ff-save-dash-dims'),
        ], {classes: 'ff-save-dash-row', style: {alignItems: 'center', gap: '6px'}}))));
        dashHost.appendChild(publishInput.root);
        if (this.dashboardProjectId != null) {
          const modeInput = ui.input.choice('Publish', {
            value: dashMode === 'new' ? 'As a new dashboard' : 'Update the existing dashboard',
            items: ['Update the existing dashboard', 'As a new dashboard'],
            tooltipText: 'This flow already published a dashboard — update it in place, or leave it be and publish a new one',
            onValueChanged: (v) => dashMode = v === 'As a new dashboard' ? 'new' : 'update',
          });
          modeRoot = modeInput.root;
          modeRoot.style.display = publishInput.value ? '' : 'none';
          dashHost.appendChild(modeRoot);
        } else
          modeRoot = null;
      }
      else if (state.running) {
        dashHost.appendChild(ui.divH([ui.loader(), ui.divText('Running the flow…')],
          {classes: 'ff-save-dash-running', style: {alignItems: 'center', gap: '8px'}}));
      } else {
        // Three distinct idle states: never ran, ran and failed, ran fine but
        // produced no tables — the user who just clicked Run must not see the
        // section silently bounce back to the initial hint.
        const ranButNothing = state.failed === false && state.ranOnce === true;
        dashHost.appendChild(state.failed ?
          ui.divText('The run failed — check the highlighted step on the canvas.', 'ff-save-dash-failed') :
          ranButNothing ?
            ui.divText('The run produced no result tables — check the flow\'s output nodes.', 'ff-save-dash-failed') :
            ui.divText('Run the flow to publish its result tables as a dashboard.', 'ff-save-dash-hint'));
        const runBtn = ui.button(state.failed || ranButNothing ? 'Run again' : 'Run the flow', () => {
          refreshDash({running: true});
          this.runInstrumented();
        });
        runBtn.style.alignSelf = 'flex-start';
        runBtn.style.marginLeft = '0';
        dashHost.appendChild(runBtn);
      }
    };
    refreshDash();
    this.saveDialogRunEnd = (success) => refreshDash({failed: !success, ranOnce: true});

    const dlg = ui.dialog({title: 'Save Flow'})
      .add(ui.divV([nameInput.root, warningDiv, descInput.root,
        spaceRow, pickerHost, dashHost]))
      .onOK(async () => {
        const name = nameInput.value.trim();
        if (name === '') { // reachable via Enter even while the Save button is disabled
          grok.shell.warning('Give the flow a name first');
          return;
        }
        if (runningNow) {
          grok.shell.warning('Wait for the run to finish');
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
        // A failed save must not go on to share/publish a stale entity.
        if (!await this.saveToServer()) return;
        const saved = this.boundScript as DG.Script | null;
        if (targetSpace && saved?.id) {
          try {
            await grok.dapi.spaces.id(targetSpace.id).addEntity(saved.id, false);
            grok.shell.info(`Added to space "${targetSpace.friendlyName}"`);
          } catch (e: any) {
            grok.shell.error(`Could not add to space: ${e?.message ?? e}`);
          }
        }
        if (saved != null && publishInput.value === true && computedTabs().length > 0) {
          if (dashMode === 'new') this.dashboardProjectId = null;
          await this.openDashboardDialog();
        }
      });
    dlg.onClose.subscribe(() => {
      this.saveDialogRunEnd = null;
      if (warnTimer != null) clearTimeout(warnTimer);
    });
    dlg.show({width: 500});
    // Validate before close: empty names never reach the OK handler.
    okBtn = dlg.getButton('OK') as HTMLButtonElement | null;
    if (okBtn) okBtn.textContent = 'Save';
    nameInput.onChanged.subscribe(() => {
      syncOk();
      if (warnTimer != null) clearTimeout(warnTimer);
      warnTimer = setTimeout(() => void refreshNameWarning(), 400);
    });
    syncOk();
    void refreshNameWarning();
    // A fresh (unsaved) flow keeps the template name — pre-uniquify it so the
    // dialog doesn't open with an "already exists" warning the user didn't cause.
    if (this.boundScript == null && !opts.asNew) {
      const orig = nameInput.value;
      void (async () => {
        try {
          if ((await findClashes(orig.trim(), null)).length === 0 || nameInput.value !== orig) return;
          for (let i = 2; i <= 9; i++) {
            const cand = `${orig.trim()} ${i}`;
            if ((await findClashes(cand, null)).length === 0) {
              if (nameInput.value === orig) nameInput.value = cand;
              return;
            }
          }
        } catch {/* advisory only */}
      })();
    }
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

  /** Belt-and-braces before stamping creation scripts: a bound script whose
   *  `nqName` carries no namespace would serialize an unqualified call
   *  (`FlowName()` instead of `user:FlowName()`) that breaks when shared —
   *  re-fetch the entity from the server, which always returns it qualified. */
  private async ensureBoundScriptQualified(): Promise<void> {
    const s = this.boundScript;
    if (s?.id == null) return;
    let nq = '';
    try {
      nq = s.nqName ?? '';
    } catch {/* Dart proxy read can throw */}
    if (!nq.includes(':'))
      this.boundScript = (await grok.dapi.scripts.find(s.id).catch(() => null)) ?? s;
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
    await this.ensureBoundScriptQualified();
    this.stampCreationScripts(tabs);
    const layoutsByParam = this.outputViews.captureLayouts();
    try {
      const saved = await DG.Project.showSaveDialog({
        tables: tabs.map((t) => t.df!),
        views: tabs.map((t) => t.tv),
        layouts: tabs.map((t) => t.tv != null ? null : layoutsByParam[t.paramName] ?? null),
        name: this.flowSettings.scriptName,
        description: this.flowSettings.scriptDescription,
        project: this.dashboardProjectId ?? undefined,
      });
      if (saved?.id && saved.id !== this.dashboardProjectId) {
        this.dashboardProjectId = saved.id;
        if (this.boundScript != null)
          await this.saveToServer(true); // persist the binding in the entity body
      }
      // The end of a two-dialog journey must not be silence — the earlier
      // "saved" balloon fired before this dialog even opened.
      if (saved?.id)
        grok.shell.info(`Dashboard "${saved.friendlyName ?? this.flowSettings.scriptName}" published`);
    } catch (e) {
      grok.shell.error(`Publish failed: ${e instanceof Error ? e.message : e}`);
    }
  }

  /** Download the graph as a local `.flow` file (the pre-entity behavior). */
  private async exportFlowFile(): Promise<void> {
    // Best-effort: an export with pending (in-memory) file ids would not be
    // portable — persist them first; on failure export anyway with a warning.
    try {
      await this.persistPendingUploads();
    } catch (e: any) {
      grok.shell.warning(`Uploaded files were not persisted — the export will not be portable. ${e?.message ?? e}`);
    }
    const doc = serializeFlow(this.flow, this.flowSettings);
    downloadFlow(doc);
    grok.shell.info('Flow exported as .flow');
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
    const report = (e: unknown): void =>
      grok.shell.error(`Run failed: ${e instanceof Error ? e.message : e}`);
    if (func.inputs.length === 0)
      void fc.call(undefined, undefined, {processed: true}).catch(report);
    else {
      fc.getEditor(false).then((e: HTMLElement) => {
        ui.dialog({title: func.friendlyName ?? func.name}).add(e).show().onOK(() => {
          void fc.call(undefined, undefined, {processed: true}).catch(report);
        });
      }).catch(report);
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
      list.style.color = 'var(--orange-3, #805125)';
      list.style.marginBottom = '8px';
      list.style.maxHeight = '120px';
      list.style.overflow = 'auto';
      blocks.push(ui.divText(`${warnings.length} warning(s) — these nodes have no creation-script form:`,
        {style: {fontWeight: 'bold', color: 'var(--orange-3, #805125)'}}));
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
      list.style.color = 'var(--orange-3, #805125)';
      list.style.maxHeight = '90px';
      list.style.overflow = 'auto';
      list.style.marginBottom = '8px';
      blocks.push(ui.divText(`${warnings.length} warning(s):`,
        {style: {fontWeight: 'bold', color: 'var(--orange-3, #805125)'}}));
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
    DG.Utils.download(`${this.flowSettings.scriptName || 'script'}.js`, script, 'text/javascript');
  }

  /** Load a flow from a JSON string (file viewer entry point). */
  async loadFromJson(json: string): Promise<void> {
    // Accepts both `.flow` shapes: the annotated script body (header + JSON)
    // and the bare JSON document — parseFlowBody handles either.
    await this.loadFromDoc(parseFlowBody(json).doc);
  }

  /** Load a parsed flow document. Awaits editor construction, so it is safe
   *  to call right after the constructor (no timer race). Loads are
   *  serialized: deserialize is clear() + per-node adds, so two interleaved
   *  loads (a double-clicked recent-flow row) would merge both graphs. */
  async loadFromDoc(doc: FuncFlowDocument): Promise<void> {
    const run = this.loadChain.then(() => this.doLoadFromDoc(doc));
    this.loadChain = run.catch(() => {});
    return run;
  }

  private async doLoadFromDoc(doc: FuncFlowDocument): Promise<void> {
    await this.editorReady;
    await deserializeFlow(doc, this.flow);
    if (doc.metadata?.settings) this.flowSettings = doc.metadata.settings;
    // A flow saved with autorun on reopens live (and vice versa).
    this.applyAutorunSetting(this.flowSettings.autorun === true);
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
    // Same serialization as loadFromDoc — this path also clears + re-adds.
    const run = this.loadChain.then(async () => {
      await this.editorReady;
      try {
        await this.flow.clear();
        const result = await buildFlowFromCreationScript(this.flow, script);
        this.updateStatusBar();
        this.fitToScreen();
        for (const warning of result.warnings) grok.shell.warning(warning);
        // The imported graph arrives mostly un-run (files may autorun, the
        // rest sits idle) — say what to do next, not just what happened.
        grok.shell.info(`Imported ${result.nodesAdded} steps — press Run to compute them`);
      } catch (e: any) {
        grok.shell.error(`Creation script import failed: ${e?.message ?? e}`);
      }
    });
    this.loadChain = run.catch(() => {});
    return run;
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
        // Importing clears the canvas — and is not undoable. A non-empty flow
        // deserves a confirmation before it is silently replaced.
        const nodeCount = this.flow?.getNodes().length ?? 0;
        if (nodeCount > 0) {
          ui.dialog({title: 'Replace current flow?'})
            .add(ui.divText(`Importing replaces the current flow (${nodeCount} ` +
              `node${nodeCount === 1 ? '' : 's'}). This cannot be undone.`))
            .onOK(() => void this.loadFromCreationScript(script))
            .show();
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
