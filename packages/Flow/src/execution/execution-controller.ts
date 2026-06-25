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
import {isExecKey} from '../rete/scheme';

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

  private executeInstrumented(
    settings: ScriptSettings, debug: boolean,
    opts?: {onlyNodeIds?: Set<string>; focusNodeId?: string; skipValidation?: boolean},
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
    // A new run invalidates anything we were showing.
    this.outputPreview.close();

    const runId = crypto.randomUUID();
    this.state.startRun(runId, this.graphVersion);
    this.visualizer.resetAllNodes();
    this.flow.clearConnectionLabels();
    this.pendingPreviewNodeId = opts?.focusNodeId ?? null;

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
      if (this.pendingPreviewNodeId) {
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
    }
  }

  resetVisuals(): void {
    this.visualizer.resetAllNodes();
    this.flow.clearConnectionLabels();
    this.state.reset();
    this.outputPreview.close();
  }

  dispose(): void {
    this.stopRun();
    this.visualizer.resetAllNodes();
  }
}
