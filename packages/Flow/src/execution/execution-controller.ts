/** Orchestrates instrumented script runs: validation, event subscriptions,
 *  state tracking, visualization, output preview. */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {FlowEditor, GraphEdit} from '../rete/flow-editor';
import {ExecutionState, NodeExecStatus, ExecEvent} from './execution-state';
import {ExecutionVisualizer} from './execution-visualizer';
import {OutputPreviewPanel} from './output-preview';
import {emitScript, ScriptSettings, EmitOptions} from '../compiler/script-emitter';
import {validateGraph} from '../compiler/validator';
import {sliceUpTo, sliceDownFrom} from '../compiler/graph-compiler';
import {ValueSummary} from './execution-state';
import {FlowNode, isExecKey, nodeMissingRequirements} from '../rete/scheme';
import {resolveInputValue, inputBlockReason} from '../utils/input-values';

/** Grow a dirty node set upstream until every connection crossing into it
 *  comes from a node whose output value is already captured (`hasLive`) — the
 *  smallest self-contained slice that can re-run with `liveExternalInputs`.
 *  Order (exec) edges never carry values, so they never force expansion. */
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

/** Short count label for a wire from a source output's summary — "1,204 × 8"
 *  (rows × columns) for a table, "1,204" for a column. Null when there's
 *  nothing countable. Exported for tests. */
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

/** A short, human data summary of a completed node's outputs for the status
 *  line under the node — the first table ("1,204 × 8") or column ("1,204 values")
 *  found, else nothing. Numbers are grouped with thousands separators. */
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
  /** The output panel's node when an invalidation closed it — the next autorun
   *  re-opens it (without stealing selection) so the preview stays live. */
  private autorunPreviewNodeId: string | null = null;
  /** When set, the next run is a single-node preview; on completion we open
   *  this node's data panel instead of the normal end-of-run behavior. */
  private pendingPreviewNodeId: string | null = null;
  /** When set, the next run is a headless slice (column picker); on completion
   *  we invoke this and skip the preview / end-of-run UI entirely. */
  private pendingOnComplete: (() => void) | null = null;
  outputPreview: OutputPreviewPanel;

  /** Called when execution stops at a breakpoint (so the view can prompt Continue). */
  onBreakpointHit?: (nodeId: string) => void;
  /** Called when the run ends, with overall success boolean. */
  onRunEnd?: (success: boolean) => void;
  /** Called on every per-node state change — used to refresh the property panel. */
  onNodeStateChanged?: (nodeId: string) => void;

  /** @param outputPreview the view-owned bottom output panel (a pane of the
   *  view's splitter). Defaults to a detached one for headless usage. */
  constructor(flow: FlowEditor, outputPreview?: OutputPreviewPanel) {
    this.flow = flow;
    this.state = new ExecutionState();
    this.visualizer = new ExecutionVisualizer(flow);
    this.outputPreview = outputPreview ?? new OutputPreviewPanel();
    // Unpinning via the header icon jumps to the node the user selected while
    // pinned — re-clicking an already-selected node fires no selection event,
    // so without this the preview would stay on the stale pinned content.
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

  /** Full run of the graph, minus any **unready** node (a required input or
   *  property unset — e.g. a plot with no table) and everything downstream of it
   *  (which would fail without its output). Blocks on real validation errors;
   *  otherwise runs the ready subgraph and warns which nodes it skipped. */
  private executeFull(settings: ScriptSettings, debug: boolean): void {
    const errors = validateGraph(this.flow);
    if (errors.some((e) => e.severity === 'error')) {
      const msgs = errors.filter((e) => e.severity === 'error').map((e) => e.message).join('\n');
      grok.shell.error('Validation errors:\n' + msgs);
      return;
    }
    const {roots, cone} = this.invalidNodes();
    const runSet = this.flow.getNodes().map((n) => n.id).filter((id) => !cone.has(id));
    if (runSet.length === 0) {
      grok.shell.error('Nothing to run — every node is missing a required input (for a plot, connect a table).');
      return;
    }
    if (roots.length > 0)
      grok.shell.warning(`Skipped ${cone.size} node(s) that aren't ready: ${roots.map((n) => n.label).join(', ')}`);
    this.executeInstrumented(settings, debug, {
      onlyNodeIds: cone.size > 0 ? new Set(runSet) : undefined,
      skipValidation: true,
    });
  }

  /** Nodes that can't run because a required input or property is unset, plus
   *  everything downstream of them. `roots` are the unready nodes themselves
   *  (for messaging); `cone` is roots + all their transitive successors — the
   *  full set the run must exclude so no emitted step references a dropped one. */
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

  /** The nodes a full run would execute — every node except the unready ones and
   *  their downstream cone. Pure; used by the run gate and tests. */
  runnableNodes(): Set<string> {
    const {cone} = this.invalidNodes();
    const set = new Set<string>();
    for (const n of this.flow.getNodes()) if (!cone.has(n.id)) set.add(n.id);
    return set;
  }

  /** Run only the slice needed to produce `targetNodeId`'s output (the node and
   *  all its upstream ancestors), then open that node's data preview. The
   *  "inspect anywhere" keystone — see any port's data without running the
   *  whole graph or wiring an Output node. */
  previewNodeData(targetNodeId: string, settings: ScriptSettings): void {
    const slice = sliceUpTo(this.flow, targetNodeId);
    this.executeInstrumented(settings, false, {onlyNodeIds: slice, focusNodeId: targetNodeId, skipValidation: true});
  }

  /** First captured dataframe clone among a node's outputs, or null. Lets the
   *  column picker read an upstream table that has already been run without
   *  re-running it. Only a *fresh* (completed, non-stale) result is reused — a
   *  graph edit marks nodes stale, so an edited upstream is recomputed rather
   *  than picked from an outdated table. (Heuristic "first dataframe": real
   *  catalog nodes that feed a table input emit a single table.) */
  cloneForNode(nodeId: string): DG.DataFrame | null {
    const st = this.state.getNodeState(nodeId);
    if (!st?.outputs || st.status !== NodeExecStatus.completed) return null;
    for (const s of Object.values(st.outputs))
      if (s.type === 'dataframe' && s.clone) return s.clone as DG.DataFrame;
    return null;
  }

  /** Run the slice up to `sourceNodeId` (and its ancestors) headlessly, then
   *  resolve with the first dataframe clone captured for that node — used by the
   *  column picker when an input's table is connected but not yet computed.
   *  Runs in `preserveState` mode so it only (re)computes this slice and leaves
   *  every other node's already-captured result intact: picking a key for one
   *  table doesn't wipe another table's result, and a later pick against the
   *  same table reuses the cached table instead of re-running it. */
  produceTableForNode(sourceNodeId: string, settings: ScriptSettings): Promise<DG.DataFrame | null> {
    return new Promise((resolve) => {
      const slice = sliceUpTo(this.flow, sourceNodeId);
      this.executeInstrumented(settings, false, {
        onlyNodeIds: slice, skipValidation: true, preserveState: true,
        onComplete: () => resolve(this.cloneForNode(sourceNodeId)),
      });
    });
  }

  /** Whether a captured value exists for a node's output (from a prior run) — the
   *  live-value registry the instrumented run populates via `__ff_stash`. */
  hasLiveValue(nodeId: string, outputKey: string): boolean {
    const reg = (globalThis as {__ffFlowLive?: Record<string, Record<string, unknown>>}).__ffFlowLive;
    return !!(reg && reg[nodeId] && outputKey in reg[nodeId]);
  }

  /** The captured live value for a node's output from a prior run, or
   *  `undefined` — used to seed a function-editor FuncCall with the values
   *  actually flowing through the node's connected inputs. */
  liveValue(nodeId: string, outputKey: string): unknown {
    const reg = (globalThis as {__ffFlowLive?: Record<string, Record<string, unknown>>}).__ffFlowLive;
    return reg?.[nodeId]?.[outputKey];
  }

  private clearLiveRegistry(): void {
    // The registry is page-global (the emitted script writes to it, outside any
    // view), and several Flow views can be live at once — never wipe it
    // wholesale. Node ids are unique per editor, so deleting this flow's ids
    // leaves other views' captured values (and their single-node rerun) intact.
    const reg = (globalThis as {__ffFlowLive?: Record<string, unknown>}).__ffFlowLive;
    if (!reg) return;
    for (const node of this.flow.getNodes()) delete reg[node.id];
  }

  /** Whether "Rerun this node only" should be offered: it's a compute node whose
   *  required inputs are all satisfied AND every connected input already has a
   *  captured upstream value (so it can run without re-running upstream). */
  canRerunNode(nodeId: string): boolean {
    const node = this.flow.getNodeById(nodeId);
    if (!node) return false;
    // Inputs/outputs have nothing to recompute on their own.
    if (node.dgNodeType !== 'func' && node.dgNodeType !== 'utility') return false;
    const anyConnected = Object.keys(node.inputs)
      .some((k) => !isExecKey(k) && !!this.flow.getInputSource(nodeId, k));
    return anyConnected && this.readyForLiveRun(nodeId);
  }

  /** Inputs satisfied for a run that must not touch upstream: no missing
   *  requirements, and every connected input is fed by a captured live value.
   *  Unlike {@link canRerunNode}, a node with no connected inputs qualifies —
   *  an Open File runs from its parameters alone. */
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

  /** Re-run just this node using upstream values captured from a prior run — its
   *  connected inputs resolve to `_ffLive(...)` registry reads, so nothing
   *  upstream re-executes. Preserves every other node's state; opens this node's
   *  preview on completion. */
  rerunNode(nodeId: string, settings: ScriptSettings): void {
    this.executeInstrumented(settings, false, {
      onlyNodeIds: new Set([nodeId]), focusNodeId: nodeId,
      skipValidation: true, preserveState: true, liveExternalInputs: true,
    });
  }

  /** Every node with no fresh result (never ran, stale, errored, mid-run) plus
   *  everything downstream of it — what switching autorun ON should
   *  immediately schedule, so enabling it on a new/half-run flow runs the
   *  missing part instead of idling until the first edit. */
  pendingNodes(): Set<string> {
    const pending = new Set<string>();
    for (const n of this.flow.getNodes()) {
      if (pending.has(n.id)) continue;
      if (this.state.getNodeState(n.id)?.status === NodeExecStatus.completed) continue;
      for (const id of sliceDownFrom(this.flow, n.id)) pending.add(id);
    }
    return pending;
  }

  /** Debounced autorun entry. Re-runs only the invalidated slice when its
   *  boundary can be fed from captured live values; falls back to a full run
   *  otherwise (nothing ran yet, upstream never completed, …). Everything is
   *  silent — no validation toasts, no run dialog, no selection stealing —
   *  because it fires after every edit. The outcome tells the scheduler
   *  whether to retry ('busy') or wait for the next edit ('skipped'). */
  /** Debounced live-node entry — the global autorun toggle is OFF and an edit
   *  touched live-by-default nodes (Open File, viewers, …). Runs ONLY the given
   *  nodes, and only those whose inputs are satisfied (see
   *  {@link readyForLiveRun}); everything else on the canvas — including the
   *  live nodes' own downstream — is left alone. Unready live nodes are
   *  silently dropped; if none remain, nothing runs. */
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

    // Never autorun an unready node (a plot with no table, a Select Column with
    // no column, …) or anything downstream of it — drop them from whatever we'd
    // run. If that leaves nothing, wait for the next edit.
    const {cone: excluded} = this.invalidNodes();
    const base = slice ?? new Set(this.flow.getNodes().map((n) => n.id));
    const runSet = new Set([...base].filter((id) => !excluded.has(id)));
    if (runSet.size === 0) return 'skipped';
    // An explicit run set is needed when we sliced OR when we pruned unready
    // nodes from a would-be full run; otherwise a plain full run (no filter).
    const useSlice = slice !== null || excluded.size > 0 ? runSet : null;

    // A run touching an input node emits `//input:` headers → a dialog on
    // every keystroke — UNLESS the node carries a resolvable configured value,
    // which `executeInstrumented` feeds to the prepared call silently. Skip
    // only for unconfigured/unresolvable ones; the bolt's blocked badge says
    // why. (A slice whose boundary covers the input nodes still autoruns.)
    const runsNode = (id: string): boolean => useSlice === null || useSlice.has(id);
    if (this.flow.getNodes().some((n) =>
      n.dgNodeType === 'input' && runsNode(n.id) && inputBlockReason(n) != null)) return 'skipped';

    const restorePreviewId = this.autorunPreviewNodeId;
    this.autorunPreviewNodeId = null;
    this.executeInstrumented(settings, false, {
      onlyNodeIds: useSlice ?? undefined,
      // A pruned full run (slice === null) is still a fresh from-scratch run of
      // the ready subgraph — only an incremental slice reads live boundary
      // values / preserves prior state.
      liveExternalInputs: slice !== null,
      preserveState: slice !== null,
      skipValidation: true,
      onComplete: () => {
        // Bring back the preview the invalidation closed — content only, no
        // selection change (the user may be mid-edit in the property panel).
        if (!restorePreviewId) return;
        const node = this.flow.getNodeById(restorePreviewId);
        const state = this.state.getNodeState(restorePreviewId);
        if (node && state?.status === NodeExecStatus.completed)
          this.outputPreview.showForNode(node, state);
      },
    });
    return 'started';
  }

  /** Resolve each of the emitted script's parameters from the input node that
   *  declared it (matched by `paramName`). `values` feeds `func.prepare`;
   *  `missing` lists params with no resolvable configured value — non-empty
   *  means an explicit run still needs the parameter dialog. */
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

  /** Why autorun cannot currently run what's pending — validation errors plus
   *  input nodes without a resolvable configured value. Empty when nothing is
   *  pending (autorun has nothing to do, so nothing blocks it). Drives the
   *  ribbon bolt's blocked badge and tooltip. */
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
      // Headless slice (column picker) or autorun slice: keep prior node
      // states, visuals, connection labels and the output panel — only the
      // slice's nodes get recomputed and re-highlighted; everything else
      // stays as it was.
      this.state.runId = runId;
      this.state.isRunning = true;
    } else {
      // A new full/slice run invalidates anything we were showing — and the
      // captured live values (they're recomputed by this run's `__ff_stash`).
      // A pinned preview keeps its stale content in place instead of blinking:
      // node-start overlays the spinner, node-complete renders fresh.
      if (this.outputPreview.pinnedNodeId == null || this.outputPreview.currentNodeId == null)
        this.outputPreview.clear();
      this.state.startRun(runId);
      this.visualizer.resetAllNodes();
      this.flow.clearConnectionLabels();
      this.clearLiveRegistry();
    }
    this.pendingPreviewNodeId = opts?.focusNodeId ?? null;
    this.pendingOnComplete = opts?.onComplete ?? null;

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
      // Values configured on the input nodes feed the prepared call directly —
      // the parameter dialog only appears for params still missing one, and
      // then comes prefilled with everything that was configured.
      const {values, missing} = this.configuredInputValues(func);
      const fc = func.prepare(values);
      // No auto-dock at completion — the user opens the panel by clicking a
      // completed node (see `showOutputsForNode`).
      if (func.inputs.length === 0 || missing.length === 0)
        void fc.call(undefined, undefined, {processed: true});
      else {
        fc.getEditor(false).then((e: HTMLElement) => {
          ui.dialog({title: settings.name}).add(e).show().onOK(async () => {
            await fc.call(undefined, undefined, {processed: true});
          });
        });
      }
    } catch (e: any) {
      grok.shell.error(`Script generation failed: ${e.message}`);
      this.stopRun();
      // A headless slice run (column picker) is awaiting completion — release it
      // so the caller's promise resolves (with null, since nothing was captured).
      if (this.pendingOnComplete) {
        const cb = this.pendingOnComplete;
        this.pendingOnComplete = null;
        cb();
      }
    }
  }

  /** Lazily open or update the docked output panel with this node's runtime
   *  values. No-op if the node has nothing captured. Called from the view's
   *  selection callback. */
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
      // The pinned node is recomputing — spinner over the kept stale content.
      if (this.outputPreview.pinnedNodeId === event.nodeId) this.outputPreview.markUpdating();
      this.onNodeStateChanged?.(event.nodeId);
      break;
    case 'node-complete':
      this.state.setNodeStatus(event.nodeId, NodeExecStatus.completed, {
        endTime: event.timestamp, outputs: event.outputs,
      });
      this.visualizer.highlightNode(event.nodeId, NodeExecStatus.completed, summarizeOutputs(event.outputs));
      this.labelOutgoingConnections(event.nodeId, event.outputs);
      this.onNodeStateChanged?.(event.nodeId);
      // A pinned preview tracks its node live: the moment a fresh result lands,
      // re-render it (selection can't bring it back — other nodes are gated).
      if (this.outputPreview.pinnedNodeId === event.nodeId) {
        const pinned = this.flow.getNodeById(event.nodeId);
        if (pinned) this.showOutputsForNode(pinned);
      }
      break;
    case 'node-error':
      this.state.setNodeStatus(event.nodeId, NodeExecStatus.errored, {
        endTime: event.timestamp, error: event.error, stack: event.stack,
      });
      this.visualizer.highlightNode(event.nodeId, NodeExecStatus.errored);
      // A failed pinned recompute: drop the spinner, keep the last good content.
      if (this.outputPreview.pinnedNodeId === event.nodeId) this.outputPreview.clearUpdating();
      this.onNodeStateChanged?.(event.nodeId);
      break;
    case 'breakpoint-hit':
      this.state.setNodeStatus(event.nodeId, NodeExecStatus.running, {startTime: event.timestamp});
      this.visualizer.highlightNode(event.nodeId, NodeExecStatus.running);
      this.onBreakpointHit?.(event.nodeId);
      break;
    case 'run-complete':
      this.state.endRun();
      // Safety: if the pinned node never completed this run (skipped/halted),
      // don't leave the spinner spinning over the kept content.
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
          void this.flow.selectNode(id);   // focuses the node + opens its panel via the view
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
    this.state.endRun();
  }

  /** After a node completes, tag its outgoing data wires with the row/value
   *  count flowing through them (Make/n8n-style on-edge counts). */
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

  /** React to a classified graph edit: invalidate exactly what the change can
   *  affect, and nothing else. Returns the set of node ids whose results must
   *  be recomputed — the autorun scheduler accumulates these.
   *
   *  - a fresh node has no wiring, so adding one invalidates nothing;
   *  - a connection change invalidates its *target* and everything downstream
   *    (the source's value is untouched — it still computed what it computed);
   *  - a parameter edit invalidates the edited node and everything downstream;
   *  - removing a node just drops its state (its connections' removal events
   *    have already invalidated the affected downstream nodes). */
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

  /** Mark the node and its transitive successors "Out of date": state, node
   *  visuals, outgoing wire labels, captured live values. Upstream nodes keep
   *  their completed results (and stay eligible for live-value reuse). */
  invalidateDownstream(rootId: string): Set<string> {
    const affected = sliceDownFrom(this.flow, rootId);
    this.state.markStale(affected);
    this.visualizer.markStale(affected);
    // Wire labels announce the value that flowed through — gone for any wire
    // *leaving* an invalidated node; wires feeding the cone from valid
    // upstream nodes keep theirs (that data is still what will flow in).
    for (const c of this.flow.getConnections())
      if (affected.has(c.source)) this.flow.setConnectionLabel(c.id, null);
    // Captured live values of invalidated nodes are lies now — drop them so
    // single-node rerun / column-picker reuse recompute instead.
    const reg = (globalThis as {__ffFlowLive?: Record<string, unknown>}).__ffFlowLive;
    if (reg)
      for (const id of affected) delete reg[id];
    // Stale values aren't worth previewing; close (and remember the node so an
    // autorun can bring the preview back once fresh values exist). A PINNED
    // preview keeps its stale content instead — hiding and re-docking on every
    // upstream edit reads as jumping; node-start overlays the recalculating
    // spinner and node-complete swaps in the fresh render in place.
    const previewId = this.outputPreview.currentNodeId;
    if (previewId !== null && affected.has(previewId) &&
        this.outputPreview.pinnedNodeId !== previewId) {
      this.outputPreview.clear();
      this.autorunPreviewNodeId = previewId;
    }
    for (const id of affected) this.onNodeStateChanged?.(id);
    return affected;
  }

  /** A node was deleted — drop every trace of it (state, visuals tracking,
   *  captured values, preview). */
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
  }
}
