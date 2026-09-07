/** Orchestrates instrumented script runs: validation, event subscriptions,
 *  state tracking, visualization, output preview. */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {FlowEditor, GraphEdit} from '../rete/flow-editor';
import {ExecutionState, NodeExecStatus, ExecEvent} from './execution-state';
import {ExecutionVisualizer, errorSummary, normalizeErrorMessage} from './execution-visualizer';
import {OutputPreviewPanel} from './output-preview';
import {graphicsElement} from './value-inspector';
import {emitScript, ScriptSettings, EmitOptions} from '../compiler/script-emitter';
import {validateGraph} from '../compiler/validator';
import {sliceUpTo, sliceDownFrom} from '../compiler/graph-compiler';
import {ValueSummary} from './execution-state';
import {FlowNode, isExecKey, nodeMissingRequirements, inlinePreviewEnabled,
  INLINE_HOSTED_DATA_KEY} from '../rete/scheme';
import {resolveInputValue, inputBlockReason} from '../utils/input-values';

/** Grow a dirty node set upstream until every crossing connection comes from a
 *  captured (`hasLive`) output — the smallest slice runnable with `liveExternalInputs`. */
export function expandToLiveBoundary(
  flow: FlowEditor, dirty: Iterable<string>, hasLive: (nodeId: string, outputKey: string) => boolean,
): Set<string> {
  const slice = new Set(dirty);
  const connections = flow.getConnections();
  const stack = [...slice];
  while (stack.length > 0) {
    const id = stack.pop()!;
    for (const c of connections) {
      if (c.target !== id || slice.has(c.source)) continue;
      if (isExecKey(String(c.sourceOutput))) continue;
      if (hasLive(c.source, String(c.sourceOutput))) continue;
      slice.add(c.source);
      stack.push(c.source);
    }
  }
  return slice;
}

/** Count label for a wire — "1,204 × 8" for a table, "1,204" for a column. */
export function connectionCountLabel(summary?: ValueSummary | null): string | null {
  if (!summary) return null;
  const n = (v: number): string => v.toLocaleString('en-US');
  if (summary.type === 'dataframe' && typeof summary.rows === 'number') {
    return typeof summary.cols === 'number' ?
      `${n(summary.rows)} × ${n(summary.cols)}` : `${n(summary.rows)} rows`;
  }
  if (summary.type === 'column' && typeof summary.length === 'number') return n(summary.length);
  return null;
}

function summarizeOutputs(outputs?: Record<string, ValueSummary>): string | undefined {
  if (!outputs) return undefined;
  const n = (v: number): string => v.toLocaleString('en-US');
  for (const s of Object.values(outputs)) {
    if (s?.type === 'dataframe' && typeof s.rows === 'number' && typeof s.cols === 'number')
      return `${n(s.rows)} × ${n(s.cols)}`;
  }
  for (const s of Object.values(outputs)) {
    if (s?.type === 'column' && typeof s.length === 'number') return `${n(s.length)} values`;
  }
  return undefined;
}

export class ExecutionController {
  state: ExecutionState;
  private visualizer: ExecutionVisualizer;
  private flow: FlowEditor;
  private subscription: {unsubscribe(): void} | null = null;
  /** The in-flight run's call — kept so Stop / a superseding run can cancel it. */
  private currentCall: DG.FuncCall | null = null;
  /** The output panel's node when an invalidation closed it — the next autorun
   *  re-opens it (without stealing selection). */
  private autorunPreviewNodeId: string | null = null;
  /** On completion, open this node's data panel instead of the end-of-run UI. */
  private pendingPreviewNodeId: string | null = null;
  /** Headless-run completion callback; skips the preview / end-of-run UI. */
  private pendingOnComplete: (() => void) | null = null;
  private runNodeIds: Set<string> | null = null;
  private runFailedNodeId: string | null = null;
  private runIsBackground = false;
  outputPreview: OutputPreviewPanel;

  onBreakpointHit?: (nodeId: string) => void;
  onRunEnd?: (success: boolean) => void;
  onNodeStateChanged?: (nodeId: string) => void;

  constructor(flow: FlowEditor, outputPreview?: OutputPreviewPanel) {
    this.flow = flow;
    this.state = new ExecutionState();
    this.visualizer = new ExecutionVisualizer(flow);
    this.outputPreview = outputPreview ?? new OutputPreviewPanel();
    // On unpin, jump to the node selected while pinned — re-clicking an
    // already-selected node fires no selection event.
    this.outputPreview.onUnpinned = (): void => {
      const ids = this.flow.getSelectedNodeIds();
      if (ids.length !== 1) return;
      const node = this.flow.getNodeById(ids[0]);
      if (node && node.id !== this.outputPreview.currentNodeId) this.showOutputsForNode(node);
    };
  }

  runInstrumented(settings: ScriptSettings): void {
    this.executeFull(settings, false);
  }

  debugInstrumented(settings: ScriptSettings): void {
    this.executeFull(settings, true);
  }

  /** Full run minus any unready node (a required input/property unset) and its
   *  downstream cone; warns which nodes were skipped. */
  private executeFull(settings: ScriptSettings, debug: boolean): void {
    const errors = validateGraph(this.flow);
    if (errors.some((e) => e.severity === 'error')) {
      const msgs = errors.filter((e) => e.severity === 'error').map((e) => e.message).join('\n');
      grok.shell.error('Validation errors:\n' + msgs);
      // A run that never starts must still report an end — the save dialog
      // waits on onRunEnd.
      this.onRunEnd?.(false);
      return;
    }
    const {roots, cone} = this.invalidNodes();
    const runSet = this.flow.getNodes().map((n) => n.id).filter((id) => !cone.has(id));
    if (runSet.length === 0) {
      grok.shell.error('Nothing to run — every node is missing a required input (for a plot, connect a table).');
      this.onRunEnd?.(false);
      return;
    }
    if (roots.length > 0)
      grok.shell.warning(`Skipped ${cone.size} node(s) that aren't ready: ${roots.map((n) => n.label).join(', ')}`);
    this.executeInstrumented(settings, debug, {
      onlyNodeIds: cone.size > 0 ? new Set(runSet) : undefined,
      skipValidation: true,
    });
  }

  /** `roots` = unready nodes; `cone` = roots + all transitive successors. */
  private invalidNodes(): {roots: FlowNode[]; cone: Set<string>} {
    const roots: FlowNode[] = [];
    const cone = new Set<string>();
    for (const n of this.flow.getNodes()) {
      if (nodeMissingRequirements(n, (k) => this.flow.isInputConnected(n.id, k)).length === 0) continue;
      roots.push(n);
      for (const id of sliceDownFrom(this.flow, n.id)) cone.add(id);
    }
    return {roots, cone};
  }

  /** The nodes a full run would execute. */
  runnableNodes(): Set<string> {
    const {cone} = this.invalidNodes();
    const set = new Set<string>();
    for (const n of this.flow.getNodes()) if (!cone.has(n.id)) set.add(n.id);
    return set;
  }

  /** Run only the slice needed for `targetNodeId`, then open its data preview. */
  previewNodeData(targetNodeId: string, settings: ScriptSettings): void {
    const slice = sliceUpTo(this.flow, targetNodeId);
    this.executeInstrumented(settings, false, {onlyNodeIds: slice, focusNodeId: targetNodeId, skipValidation: true});
  }

  /** First captured dataframe clone among a node's outputs, or null. Only a
   *  fresh (completed, non-stale) result is reused — an edited upstream must
   *  recompute, not serve an outdated table. */
  cloneForNode(nodeId: string): DG.DataFrame | null {
    const st = this.state.getNodeState(nodeId);
    if (!st?.outputs || st.status !== NodeExecStatus.completed) return null;
    for (const s of Object.values(st.outputs))
      if (s.type === 'dataframe' && s.clone) return s.clone as DG.DataFrame;
    return null;
  }

  /** Headless run of the slice up to `sourceNodeId`, resolving with its first
   *  dataframe clone. `preserveState` keeps every other node's captured result
   *  intact so one column pick can't wipe another table's result. */
  produceTableForNode(sourceNodeId: string, settings: ScriptSettings): Promise<DG.DataFrame | null> {
    return new Promise((resolve) => {
      const slice = sliceUpTo(this.flow, sourceNodeId);
      this.executeInstrumented(settings, false, {
        onlyNodeIds: slice, skipValidation: true, preserveState: true,
        onComplete: () => resolve(this.cloneForNode(sourceNodeId)),
      });
    });
  }

  /** Elements built from captured graphics outputs, one per summary — React
   *  re-renders must re-attach the same element, never rebuild it. */
  private readonly graphicsPreviewEls = new WeakMap<ValueSummary, HTMLElement>();

  /** The in-node preview content for a node: the first captured live
   *  viewer/widget root, else an element built from a graphics output. Kept
   *  through `stale` — like the bottom panel, the node shows the last result. */
  inlinePreviewRoot(nodeId: string): HTMLElement | null {
    const live = this.liveObjectRoot(nodeId);
    if (live) return live;
    const outputs = this.state.getNodeState(nodeId)?.outputs;
    if (!outputs) return null;
    // Graphics arrive as data (SVG markup / base64 PNG), not a live object — the
    // node gets its own element and the bottom panel keeps its own copy too.
    for (const s of Object.values(outputs)) {
      if (s != null && s.type === 'graphics' && typeof s.value === 'string') {
        let el = this.graphicsPreviewEls.get(s);
        if (!el) {
          el = graphicsElement(s.value as string);
          el && (el.style.removeProperty('min-height'));
          this.graphicsPreviewEls.set(s, el);
        }
        return el;
      }
    }
    return null;
  }

  /** First captured live viewer/widget root — the only content the node and the
   *  bottom panel would otherwise fight over (graphics is copied, not shared). */
  private liveObjectRoot(nodeId: string): HTMLElement | null {
    const outputs = this.state.getNodeState(nodeId)?.outputs;
    if (!outputs) return null;
    for (const s of Object.values(outputs)) {
      if (s != null && (s.type === 'viewer' || s.type === 'widget') &&
          s.value?.root instanceof HTMLElement)
        return s.value.root as HTMLElement;
    }
    return null;
  }

  /** True while a run in progress still has this node ahead of it (part of the
   *  run set, neither completed nor errored yet) — drives the in-node preview's
   *  loader so an upstream computation never reads as a blank box. */
  inlinePreviewPending(nodeId: string): boolean {
    if (!this.state.isRunning) return false;
    if (this.runNodeIds && !this.runNodeIds.has(nodeId)) return false;
    const s = this.state.getNodeState(nodeId)?.status;
    return s !== NodeExecStatus.completed && s !== NodeExecStatus.errored;
  }

  /** Stamp/clear the hosted marker on the node's live root. A stamped root makes
   *  the bottom panel render a note instead of stealing the element the node
   *  preview is showing; stamping happens BEFORE any panel build so the outcome
   *  never depends on render order. */
  syncInlinePreviewOwnership(nodeId: string): void {
    const root = this.liveObjectRoot(nodeId);
    if (!root) return;
    const node = this.flow.getNodeById(nodeId);
    if (node && inlinePreviewEnabled(node) && !node.collapsed)
      root.dataset[INLINE_HOSTED_DATA_KEY] = 'true';
    else
      delete root.dataset[INLINE_HOSTED_DATA_KEY];
  }

  /** Whether the `__ff_stash` live-value registry holds this node's output. */
  hasLiveValue(nodeId: string, outputKey: string): boolean {
    const reg = (globalThis as {__ffFlowLive?: Record<string, Record<string, unknown>>}).__ffFlowLive;
    return !!(reg && reg[nodeId] && outputKey in reg[nodeId]);
  }

  /** The captured live value for a node's output from a prior run, or `undefined`. */
  liveValue(nodeId: string, outputKey: string): unknown {
    const reg = (globalThis as {__ffFlowLive?: Record<string, Record<string, unknown>>}).__ffFlowLive;
    return reg?.[nodeId]?.[outputKey];
  }

  private clearLiveRegistry(): void {
    // The registry is page-global and shared by all live Flow views — delete
    // only this flow's node ids, never wipe it wholesale.
    const reg = (globalThis as {__ffFlowLive?: Record<string, unknown>}).__ffFlowLive;
    if (!reg) return;
    for (const node of this.flow.getNodes()) delete reg[node.id];
  }

  /** Whether "Rerun this node only" should be offered: a compute node runnable
   *  without touching upstream. */
  canRerunNode(nodeId: string): boolean {
    const node = this.flow.getNodeById(nodeId);
    if (!node) return false;
    if (node.dgNodeType !== 'func' && node.dgNodeType !== 'utility') return false;
    const anyConnected = Object.keys(node.inputs)
      .some((k) => !isExecKey(k) && !!this.flow.getInputSource(nodeId, k));
    return anyConnected && this.readyForLiveRun(nodeId);
  }

  /** No missing requirements and every connected input fed by a captured live
   *  value; a node with no connected inputs (Open File) qualifies. */
  private readyForLiveRun(nodeId: string): boolean {
    const node = this.flow.getNodeById(nodeId);
    if (!node) return false;
    if (nodeMissingRequirements(node, (k) => this.flow.isInputConnected(nodeId, k)).length > 0) return false;
    for (const key of Object.keys(node.inputs)) {
      if (isExecKey(key)) continue;
      const src = this.flow.getInputSource(nodeId, key);
      if (src && !this.hasLiveValue(src.node.id, src.outputKey)) return false;
    }
    return true;
  }

  /** Re-run just this node — connected inputs resolve to `_ffLive(...)` registry
   *  reads, so nothing upstream re-executes. */
  rerunNode(nodeId: string, settings: ScriptSettings): void {
    this.executeInstrumented(settings, false, {
      onlyNodeIds: new Set([nodeId]), focusNodeId: nodeId,
      skipValidation: true, preserveState: true, liveExternalInputs: true,
    });
  }

  /** Every node with no fresh result plus its downstream — what switching
   *  autorun ON should immediately schedule. */
  pendingNodes(): Set<string> {
    const pending = new Set<string>();
    for (const n of this.flow.getNodes()) {
      if (pending.has(n.id)) continue;
      if (this.state.getNodeState(n.id)?.status === NodeExecStatus.completed) continue;
      for (const id of sliceDownFrom(this.flow, n.id)) pending.add(id);
    }
    return pending;
  }

  /** Live-node entry (autorun toggle OFF): runs ONLY the given ready nodes,
   *  leaving everything else — including their downstream — alone. */
  runLiveNodes(liveIds: Set<string>, settings: ScriptSettings): 'started' | 'busy' | 'skipped' {
    if (this.state.isRunning) return 'busy';
    // A mid-edit graph is often momentarily invalid; just wait for more edits.
    if (validateGraph(this.flow).some((e) => e.severity === 'error')) return 'skipped';

    const runSet = new Set<string>();
    for (const id of liveIds) {
      const node = this.flow.getNodeById(id);
      // An input node would emit `//input:` headers → a dialog; never live-run.
      if (!node || node.dgNodeType === 'input') continue;
      if (this.readyForLiveRun(id)) runSet.add(id);
    }
    if (runSet.size === 0) return 'skipped';

    const restorePreviewId = this.autorunPreviewNodeId;
    this.autorunPreviewNodeId = null;
    this.executeInstrumented(settings, false, {
      onlyNodeIds: runSet,
      liveExternalInputs: true,
      preserveState: true,
      skipValidation: true,
      onComplete: () => {
        if (!restorePreviewId) return;
        const node = this.flow.getNodeById(restorePreviewId);
        const state = this.state.getNodeState(restorePreviewId);
        if (node && state?.status === NodeExecStatus.completed)
          this.outputPreview.showForNode(node, state);
      },
    });
    return 'started';
  }

  runAutorun(dirty: Set<string>, settings: ScriptSettings): 'started' | 'busy' | 'skipped' {
    if (this.state.isRunning) return 'busy';
    // A mid-edit graph is often momentarily invalid; just wait for more edits.
    if (validateGraph(this.flow).some((e) => e.severity === 'error')) return 'skipped';

    let slice: Set<string> | null = null;
    if (dirty.size > 0 && this.state.nodeStates.size > 0) {
      const expanded = expandToLiveBoundary(this.flow, dirty,
        (nodeId, outputKey) => this.hasLiveValue(nodeId, outputKey));
      if (expanded.size < this.flow.getNodeCount()) slice = expanded;
    }

    const {cone: excluded} = this.invalidNodes();
    const base = slice ?? new Set(this.flow.getNodes().map((n) => n.id));
    const runSet = new Set([...base].filter((id) => !excluded.has(id)));
    if (runSet.size === 0) return 'skipped';
    const useSlice = slice !== null || excluded.size > 0 ? runSet : null;

    // An input node without a resolvable configured value would pop the
    // parameter dialog on every keystroke — skip; the bolt's blocked badge says why.
    const runsNode = (id: string): boolean => useSlice === null || useSlice.has(id);
    if (this.flow.getNodes().some((n) =>
      n.dgNodeType === 'input' && runsNode(n.id) && inputBlockReason(n) != null)) return 'skipped';

    const restorePreviewId = this.autorunPreviewNodeId;
    this.autorunPreviewNodeId = null;
    this.executeInstrumented(settings, false, {
      onlyNodeIds: useSlice ?? undefined,
      // Only an incremental slice reads live boundary values / preserves state;
      // a pruned full run is still a fresh from-scratch run.
      liveExternalInputs: slice !== null,
      preserveState: slice !== null,
      skipValidation: true,
      onComplete: () => {
        // Bring back the preview the invalidation closed — content only, no
        // selection change.
        if (!restorePreviewId) return;
        const node = this.flow.getNodeById(restorePreviewId);
        const state = this.state.getNodeState(restorePreviewId);
        if (node && state?.status === NodeExecStatus.completed)
          this.outputPreview.showForNode(node, state);
      },
    });
    return 'started';
  }

  /** Resolve the emitted script's parameters from the input nodes that declared
   *  them; `missing` params still need the parameter dialog. */
  private configuredInputValues(func: DG.Func): {values: Record<string, unknown>; missing: string[]} {
    const byParam = new Map<string, FlowNode>();
    for (const n of this.flow.getNodes()) {
      if (n.dgNodeType === 'input')
        byParam.set(String(n.properties['paramName'] ?? ''), n);
    }
    const values: Record<string, unknown> = {};
    const missing: string[] = [];
    for (const p of func.inputs) {
      const node = byParam.get(p.name);
      const r = node ? resolveInputValue(node) : {ok: false as const};
      if (r.ok && 'value' in r) values[p.name] = r.value;
      else missing.push(p.name);
    }
    return {values, missing};
  }

  /** Why autorun cannot currently run what's pending — drives the ribbon bolt's
   *  blocked badge and tooltip. */
  autorunBlockers(): string[] {
    const pending = this.pendingNodes();
    if (pending.size === 0) return [];
    const reasons = validateGraph(this.flow)
      .filter((e) => e.severity === 'error').map((e) => e.message);
    for (const n of this.flow.getNodes()) {
      if (n.dgNodeType !== 'input' || !pending.has(n.id)) continue;
      const r = inputBlockReason(n);
      if (r) reasons.push(r);
    }
    return reasons;
  }

  private executeInstrumented(
    settings: ScriptSettings, debug: boolean,
    opts?: {onlyNodeIds?: Set<string>; focusNodeId?: string; skipValidation?: boolean;
      onComplete?: () => void; preserveState?: boolean; liveExternalInputs?: boolean},
  ): void {
    if (!opts?.skipValidation) {
      const errors = validateGraph(this.flow);
      if (errors.some((e) => e.severity === 'error')) {
        const msgs = errors.filter((e) => e.severity === 'error').map((e) => e.message).join('\n');
        grok.shell.error('Validation errors:\n' + msgs);
        return;
      }
    }

    this.stopRun();

    const runId = crypto.randomUUID();
    if (opts?.preserveState) {
      this.state.runId = runId;
      this.state.isRunning = true;
    } else {
      // A pinned preview keeps its stale content in place instead of blinking:
      // node-start overlays the spinner, node-complete renders fresh.
      if (this.outputPreview.pinnedNodeId == null || this.outputPreview.currentNodeId == null)
        this.outputPreview.clear();
      this.state.startRun(runId);
      this.visualizer.beginRun();
      this.flow.clearConnectionLabels();
      this.clearLiveRegistry();
    }
    this.pendingPreviewNodeId = opts?.focusNodeId ?? null;
    this.pendingOnComplete = opts?.onComplete ?? null;
    this.runNodeIds = new Set(opts?.onlyNodeIds ?? this.flow.getNodes().map((n) => n.id));
    this.runFailedNodeId = null;
    this.runIsBackground = opts?.preserveState === true;

    const channel = `funcflow.exec.${runId}`;
    this.subscription = grok.events.onCustomEvent(channel).subscribe((event: ExecEvent) => {
      this.handleEvent(event);
    });

    const options: EmitOptions = {
      instrumented: true,
      runId,
      enableBreakpoints: debug,
      haltOnError: true,
      onlyNodeIds: opts?.onlyNodeIds,
      liveExternalInputs: opts?.liveExternalInputs,
    };

    try {
      const script = emitScript(this.flow, settings, options);
      const func = DG.Script.create(script);
      // Configured input-node values feed the prepared call directly; the
      // parameter dialog only appears for params still missing one.
      const {values, missing} = this.configuredInputValues(func);
      const fc = func.prepare(values);
      this.currentCall = fc;
      // The rejection is the same failure already surfaced via the node-error events.
      const startCall = (): void => {void fc.call(undefined, undefined, {processed: true}).catch(() => {})};
      if (func.inputs.length === 0 || missing.length === 0)
        startCall();
      else {
        fc.getEditor(false).then((e: HTMLElement) => {
          // Canceling the dialog must release the run state like a failed start
          // does — `isRunning` would otherwise stay true forever.
          let started = false;
          const dlg = ui.dialog({title: settings.name}).add(e);
          dlg.onOK(() => {
            started = true;
            startCall();
          });
          dlg.onClose.subscribe(() => {
            if (!started) this.abortPendingRun();
          });
          dlg.show();
        }).catch((err: unknown) => {
          grok.shell.error(`Could not build the run dialog: ${err instanceof Error ? err.message : err}`);
          this.abortPendingRun();
        });
      }
    } catch (e: any) {
      grok.shell.error(`Script generation failed: ${e.message}`);
      this.abortPendingRun();
    }
  }

  /** A run that was set up but never got going: release the run state and
   *  settle every waiter. */
  private abortPendingRun(): void {
    this.stopRun();
    this.pendingPreviewNodeId = null;
    if (this.pendingOnComplete) {
      const cb = this.pendingOnComplete;
      this.pendingOnComplete = null;
      cb();
    }
    else
      this.onRunEnd?.(false);
  }

  showOutputsForNode(node: {id: string; label: string}): void {
    const state = this.state.getNodeState(node.id);
    this.outputPreview.showForNode(node, state);
  }

  private handleEvent(event: ExecEvent): void {
    switch (event.type) {
    case 'run-start':
      break;
    case 'node-start':
      this.state.setNodeStatus(event.nodeId, NodeExecStatus.running, {startTime: event.timestamp});
      this.visualizer.highlightNode(event.nodeId, NodeExecStatus.running);
      if (this.outputPreview.pinnedNodeId === event.nodeId) this.outputPreview.markUpdating();
      this.onNodeStateChanged?.(event.nodeId);
      break;
    case 'node-complete':
      this.state.setNodeStatus(event.nodeId, NodeExecStatus.completed, {
        endTime: event.timestamp, outputs: event.outputs,
      });
      // Claim a fresh viewer/widget root for the in-node preview before any
      // panel below could mount it.
      this.syncInlinePreviewOwnership(event.nodeId);
      this.visualizer.highlightNode(event.nodeId, NodeExecStatus.completed, summarizeOutputs(event.outputs));
      this.labelOutgoingConnections(event.nodeId, event.outputs);
      this.onNodeStateChanged?.(event.nodeId);
      // A pinned preview tracks its node live: re-render the moment a fresh result lands.
      if (this.outputPreview.pinnedNodeId === event.nodeId) {
        const pinned = this.flow.getNodeById(event.nodeId);
        if (pinned) this.showOutputsForNode(pinned);
      } else if (this.outputPreview.pinnedNodeId == null) {
        // A fresh result for the sole selected node opens/updates its preview —
        // an autorun must not require an unselect/select round-trip to see it.
        const sel = this.flow.getSelectedNodeIds();
        if (sel.length === 1 && sel[0] === event.nodeId) {
          const node = this.flow.getNodeById(event.nodeId);
          if (node) this.showOutputsForNode(node);
        }
      }
      break;
    case 'node-error': {
      const msg = normalizeErrorMessage(event.error);
      if (this.runFailedNodeId == null) this.runFailedNodeId = event.nodeId;
      this.state.setNodeStatus(event.nodeId, NodeExecStatus.errored, {
        endTime: event.timestamp, error: msg, stack: event.stack,
      });
      this.visualizer.highlightNode(event.nodeId, NodeExecStatus.errored, errorSummary(msg));
      if (this.outputPreview.pinnedNodeId === event.nodeId) this.outputPreview.clearUpdating();
      this.onNodeStateChanged?.(event.nodeId);
      break;
    }
    case 'breakpoint-hit':
      this.state.setNodeStatus(event.nodeId, NodeExecStatus.running, {startTime: event.timestamp});
      this.visualizer.highlightNode(event.nodeId, NodeExecStatus.running);
      this.onBreakpointHit?.(event.nodeId);
      break;
    case 'run-complete':
      this.state.endRun();
      this.currentCall = null;
      // A halted run must explain itself on the nodes it never reached — a
      // silent grey branch gives the user no path to the actual cause.
      if (this.runFailedNodeId != null) {
        const failed = this.flow.getNodeById(this.runFailedNodeId);
        const failedLabel = failed?.label ?? 'a previous step';
        const skipped = [...(this.runNodeIds ?? [])].filter((id) => {
          const s = this.state.getNodeState(id)?.status;
          return id !== this.runFailedNodeId &&
            s !== NodeExecStatus.completed && s !== NodeExecStatus.errored;
        });
        this.visualizer.markSkipped(skipped, failedLabel);
        if (!this.runIsBackground)
          grok.shell.error(`Run stopped: "${failedLabel}" failed`);
        this.runFailedNodeId = null;
      }
      // Drop the subscription now, not at the next run — a view that never runs
      // again would keep it (and everything it captures) alive.
      if (this.subscription) {
        this.subscription.unsubscribe();
        this.subscription = null;
      }
      // A pinned node that never completed this run must not keep its spinner.
      this.outputPreview.clearUpdating();
      if (this.pendingOnComplete) {
        const cb = this.pendingOnComplete;
        this.pendingOnComplete = null;
        this.pendingPreviewNodeId = null;
        cb();
      } else if (this.pendingPreviewNodeId) {
        const id = this.pendingPreviewNodeId;
        this.pendingPreviewNodeId = null;
        const node = this.flow.getNodeById(id);
        if (node) {
          void this.flow.selectNode(id);
          this.showOutputsForNode(node);
        }
      } else {
        this.onRunEnd?.(event.success === true);
      }
      break;
    }
  }

  continueBreakpoint(): void {
    if (!this.state.runId) return;
    grok.events.fireCustomEvent(`funcflow.exec.${this.state.runId}.continue`, {type: 'continue'});
  }

  stopRun(): void {
    if (this.subscription) {
      this.subscription.unsubscribe();
      this.subscription = null;
    }
    // Mark the abandoned run so its still-executing script stops writing to the
    // live registry (`__ff_stash` checks this set) — late stashes would
    // otherwise overwrite the next run's values.
    if (this.state.isRunning && this.state.runId) {
      const g = globalThis as {__ffAbortedRuns?: Set<string>};
      const aborted = g.__ffAbortedRuns ?? (g.__ffAbortedRuns = new Set());
      aborted.add(this.state.runId);
      if (aborted.size > 100)
        aborted.delete(aborted.values().next().value!);
    }
    // Release a script paused at a breakpoint so its promise settles and the
    // closure (captured tables, the continue subscription) is freed.
    this.continueBreakpoint();
    if (this.currentCall) {
      // Best effort — `cancel()` may throw or return undefined; never let that
      // break the state release below.
      try {
        void Promise.resolve(this.currentCall.cancel()).catch(() => {});
      } catch {/* cancellation unsupported here — the run just finishes headless */}
      this.currentCall = null;
    }
    this.state.endRun();
  }

  /** Tag a completed node's outgoing data wires with the count flowing through them. */
  private labelOutgoingConnections(nodeId: string, outputs?: Record<string, ValueSummary>): void {
    if (!outputs) return;
    for (const c of this.flow.getConnections()) {
      if (c.source !== nodeId) continue;
      const key = String(c.sourceOutput);
      if (isExecKey(key)) continue;
      const label = connectionCountLabel(outputs[key]);
      if (label) this.flow.setConnectionLabel(c.id, label);
    }
  }

  /** Invalidate exactly what a classified graph edit can affect; returns the
   *  node ids whose results must be recomputed. */
  applyGraphEdit(edit: GraphEdit): Set<string> {
    switch (edit.kind) {
    case 'node-added':
      return new Set();
    case 'node-removed':
      this.forgetNode(edit.nodeId);
      return new Set();
    case 'connection-added':
    case 'connection-removed':
      return this.invalidateDownstream(edit.targetId);
    case 'params-changed':
      return this.invalidateDownstream(edit.nodeId);
    case 'cleared':
      this.resetVisuals();
      return new Set();
    }
  }

  /** Mark the node and its transitive successors "Out of date": state, visuals,
   *  outgoing wire labels, captured live values. */
  invalidateDownstream(rootId: string): Set<string> {
    const affected = sliceDownFrom(this.flow, rootId);
    this.state.markStale(affected);
    this.visualizer.markStale(affected);
    for (const c of this.flow.getConnections())
      if (affected.has(c.source)) this.flow.setConnectionLabel(c.id, null);
    const reg = (globalThis as {__ffFlowLive?: Record<string, unknown>}).__ffFlowLive;
    if (reg)
      for (const id of affected) delete reg[id];
    // Close the preview, remembering the node so autorun can bring it back
    // fresh. A PINNED preview keeps its stale content instead — hiding and
    // re-docking on every upstream edit reads as jumping.
    const previewId = this.outputPreview.currentNodeId;
    if (previewId !== null && affected.has(previewId) &&
        this.outputPreview.pinnedNodeId !== previewId) {
      this.outputPreview.clear();
      this.autorunPreviewNodeId = previewId;
    }
    for (const id of affected) this.onNodeStateChanged?.(id);
    return affected;
  }

  private forgetNode(nodeId: string): void {
    this.state.forgetNode(nodeId);
    this.visualizer.forgetNode(nodeId);
    const reg = (globalThis as {__ffFlowLive?: Record<string, unknown>}).__ffFlowLive;
    if (reg) delete reg[nodeId];
    if (this.outputPreview.currentNodeId === nodeId) this.outputPreview.clear();
    if (this.outputPreview.pinnedNodeId === nodeId) this.outputPreview.unpin();
    if (this.autorunPreviewNodeId === nodeId) this.autorunPreviewNodeId = null;
  }

  resetVisuals(): void {
    this.visualizer.resetAllNodes();
    this.flow.clearConnectionLabels();
    this.state.reset();
    this.outputPreview.clear();
    this.outputPreview.unpin();
    this.clearLiveRegistry();
  }

  dispose(): void {
    this.stopRun();
    this.visualizer.resetAllNodes();
    // The captured DataFrame clones on the page-global registry must not outlive the view.
    this.clearLiveRegistry();
  }
}
