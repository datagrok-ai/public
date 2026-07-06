/** Orchestrates instrumented script runs: validation, event subscriptions,
 *  state tracking, visualization, output preview. */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {FlowEditor} from '../rete/flow-editor';
import {ExecutionState, NodeExecStatus, ExecEvent} from './execution-state';
import {ExecutionVisualizer} from './execution-visualizer';
import {OutputPreviewPanel} from './output-preview';
import {emitScript, ScriptSettings, EmitOptions} from '../compiler/script-emitter';
import {validateGraph} from '../compiler/validator';
import {sliceUpTo} from '../compiler/graph-compiler';
import {ValueSummary} from './execution-state';
import {isExecKey, missingRequiredInputs} from '../rete/scheme';

/** Short count label for a wire from a source output's summary — "1,204 rows"
 *  for a table, "1,204" for a column. Null when there's nothing countable. */
function connectionCountLabel(summary?: ValueSummary): string | null {
  if (!summary) return null;
  const n = (v: number): string => v.toLocaleString('en-US');
  if (summary.type === 'dataframe' && typeof summary.rows === 'number') return `${n(summary.rows)} rows`;
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
  private graphVersion = 0;
  /** When set, the next run is a single-node preview; on completion we open
   *  this node's data panel instead of the normal end-of-run behavior. */
  private pendingPreviewNodeId: string | null = null;
  /** When set, the next run is a headless slice (column picker); on completion
   *  we invoke this and skip the preview / end-of-run UI entirely. */
  private pendingOnComplete: (() => void) | null = null;
  outputPreview: OutputPreviewPanel = new OutputPreviewPanel();

  /** Called when execution stops at a breakpoint (so the view can prompt Continue). */
  onBreakpointHit?: (nodeId: string) => void;
  /** Called when the run ends, with overall success boolean. */
  onRunEnd?: (success: boolean) => void;
  /** Called on every per-node state change — used to refresh the property panel. */
  onNodeStateChanged?: (nodeId: string) => void;

  constructor(flow: FlowEditor) {
    this.flow = flow;
    this.state = new ExecutionState();
    this.visualizer = new ExecutionVisualizer(flow);
  }

  runInstrumented(settings: ScriptSettings): void {
    this.executeInstrumented(settings, false);
  }

  debugInstrumented(settings: ScriptSettings): void {
    this.executeInstrumented(settings, true);
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
    (globalThis as {__ffFlowLive?: unknown}).__ffFlowLive = {};
  }

  /** Whether "Rerun this node only" should be offered: it's a compute node whose
   *  required inputs are all satisfied AND every connected input already has a
   *  captured upstream value (so it can run without re-running upstream). */
  canRerunNode(nodeId: string): boolean {
    const node = this.flow.getNodeById(nodeId);
    if (!node) return false;
    // Inputs/outputs have nothing to recompute on their own.
    if (node.dgNodeType !== 'func' && node.dgNodeType !== 'utility') return false;
    if (missingRequiredInputs(node, (k) => this.flow.isInputConnected(nodeId, k)).length > 0) return false;
    let anyConnected = false;
    for (const key of Object.keys(node.inputs)) {
      if (isExecKey(key)) continue;
      const src = this.flow.getInputSource(nodeId, key);
      if (!src) continue;
      anyConnected = true;
      if (!this.hasLiveValue(src.node.id, src.outputKey)) return false;
    }
    return anyConnected;
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
      // Headless slice (column picker): keep prior node states, visuals,
      // connection labels and the output panel — only the slice's nodes get
      // recomputed and re-highlighted; everything else stays as it was.
      this.state.runId = runId;
      this.state.isRunning = true;
      this.state.graphVersionAtRun = this.graphVersion;
    } else {
      // A new full/slice run invalidates anything we were showing — and the
      // captured live values (they're recomputed by this run's `__ff_stash`).
      this.outputPreview.close();
      this.state.startRun(runId, this.graphVersion);
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
      const fc = func.prepare();
      // No auto-dock at completion — the user opens the panel by clicking a
      // completed node (see `showOutputsForNode`).
      if (func.inputs.length === 0)
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
      this.onNodeStateChanged?.(event.nodeId);
      break;
    case 'node-complete':
      this.state.setNodeStatus(event.nodeId, NodeExecStatus.completed, {
        endTime: event.timestamp, outputs: event.outputs,
      });
      this.visualizer.highlightNode(event.nodeId, NodeExecStatus.completed, summarizeOutputs(event.outputs));
      this.labelOutgoingConnections(event.nodeId, event.outputs);
      this.onNodeStateChanged?.(event.nodeId);
      break;
    case 'node-error':
      this.state.setNodeStatus(event.nodeId, NodeExecStatus.errored, {
        endTime: event.timestamp, error: event.error, stack: event.stack,
      });
      this.visualizer.highlightNode(event.nodeId, NodeExecStatus.errored);
      this.onNodeStateChanged?.(event.nodeId);
      break;
    case 'breakpoint-hit':
      this.state.setNodeStatus(event.nodeId, NodeExecStatus.running, {startTime: event.timestamp});
      this.visualizer.highlightNode(event.nodeId, NodeExecStatus.running);
      this.onBreakpointHit?.(event.nodeId);
      break;
    case 'run-complete':
      this.state.endRun();
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

  /** Bumped every time the graph changes; if results are stale, mark them so. */
  onGraphChanged(): void {
    this.graphVersion++;
    if (this.state.nodeStates.size > 0 && this.state.isStale(this.graphVersion)) {
      this.state.markAllStale();
      this.visualizer.markAllStale();
      this.flow.clearConnectionLabels();
      // Stale values aren't worth previewing; close to avoid stale impressions.
      this.outputPreview.close();
      // A structural change invalidates captured values — disable single-node
      // re-run (which reads them) until a fresh run repopulates the registry.
      this.clearLiveRegistry();
    }
  }

  resetVisuals(): void {
    this.visualizer.resetAllNodes();
    this.flow.clearConnectionLabels();
    this.state.reset();
    this.outputPreview.close();
    this.clearLiveRegistry();
  }

  dispose(): void {
    this.stopRun();
    this.visualizer.resetAllNodes();
  }
}
