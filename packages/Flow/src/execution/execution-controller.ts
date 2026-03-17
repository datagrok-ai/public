/** Orchestrates instrumented script runs: event subscriptions, state, visualization */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {LGraph} from 'litegraph.js';
import {ExecutionState, NodeExecStatus, ExecEvent} from './execution-state';
import {ExecutionVisualizer} from './execution-visualizer';
import {CanvasController} from '../canvas/canvas-controller';
import {GraphManager} from '../canvas/graph-manager';
import {emitScript, ScriptSettings, EmitOptions} from '../compiler/script-emitter';
import {validateGraph} from '../compiler/validator';

export class ExecutionController {
  state: ExecutionState;
  private visualizer: ExecutionVisualizer;
  private graphManager: GraphManager;
  private subscription: any = null;
  private graphVersion: number = 0;

  /** Callback fired when a breakpoint is hit (so the view can show Continue button) */
  onBreakpointHit?: (nodeId: number) => void;
  /** Callback fired when run ends */
  onRunEnd?: (success: boolean) => void;
  /** Callback fired on any node state change (for property panel refresh) */
  onNodeStateChanged?: (nodeId: number) => void;

  constructor(canvasController: CanvasController, graphManager: GraphManager) {
    this.state = new ExecutionState();
    this.visualizer = new ExecutionVisualizer(canvasController, graphManager);
    this.graphManager = graphManager;
  }

  /** Run the graph with instrumented script (no breakpoints) */
  runInstrumented(graph: LGraph, settings: ScriptSettings): void {
    this.executeInstrumented(graph, settings, false);
  }

  /** Run the graph in debug mode (breakpoints enabled) */
  debugInstrumented(graph: LGraph, settings: ScriptSettings): void {
    this.executeInstrumented(graph, settings, true);
  }

  private executeInstrumented(graph: LGraph, settings: ScriptSettings, debug: boolean): void {
    // Validate first
    const errors = validateGraph(graph);
    if (errors.some((e) => e.severity === 'error')) {
      const msgs = errors.filter((e) => e.severity === 'error').map((e) => e.message).join('\n');
      grok.shell.error('Validation errors:\n' + msgs);
      return;
    }

    // Stop any previous run
    this.stopRun();

    // Generate run ID and reset state
    const runId = crypto.randomUUID();
    this.state.startRun(runId, this.graphVersion);
    this.visualizer.resetAllNodes();

    // Subscribe to events BEFORE launching the script
    const channel = `funcflow.exec.${runId}`;
    this.subscription = grok.events.onCustomEvent(channel).subscribe((event: ExecEvent) => {
      this.handleEvent(event);
    });

    // Generate instrumented script
    const options: EmitOptions = {
      instrumented: true,
      runId,
      enableBreakpoints: debug,
      haltOnError: true,
    };

    try {
      const script = emitScript(graph, settings, options);
      DG.Script.create(script).prepare().edit();
    } catch (e: any) {
      grok.shell.error(`Script generation failed: ${e.message}`);
      this.stopRun();
    }
  }

  private handleEvent(event: ExecEvent): void {
    switch (event.type) {
    case 'run-start':
      // Already handled in executeInstrumented
      break;

    case 'node-start':
      this.state.setNodeStatus(event.nodeId, NodeExecStatus.running, {
        startTime: event.timestamp,
      });
      this.visualizer.highlightNode(event.nodeId, NodeExecStatus.running);
      this.onNodeStateChanged?.(event.nodeId);
      break;

    case 'node-complete':
      this.state.setNodeStatus(event.nodeId, NodeExecStatus.completed, {
        endTime: event.timestamp,
        outputs: event.outputs,
      });
      this.visualizer.highlightNode(event.nodeId, NodeExecStatus.completed);
      this.onNodeStateChanged?.(event.nodeId);
      break;

    case 'node-error':
      this.state.setNodeStatus(event.nodeId, NodeExecStatus.errored, {
        endTime: event.timestamp,
        error: event.error,
        stack: event.stack,
      });
      this.visualizer.highlightNode(event.nodeId, NodeExecStatus.errored);
      this.onNodeStateChanged?.(event.nodeId);
      break;

    case 'breakpoint-hit':
      this.state.setNodeStatus(event.nodeId, NodeExecStatus.running, {
        startTime: event.timestamp,
      });
      this.visualizer.highlightNode(event.nodeId, NodeExecStatus.running);
      this.onBreakpointHit?.(event.nodeId);
      break;

    case 'run-complete':
      this.state.endRun();
      this.onRunEnd?.(event.success === true);
      break;
    }
  }

  /** Fire continue event to resume from breakpoint */
  continueBreakpoint(): void {
    if (!this.state.runId) return;
    grok.events.fireCustomEvent(`funcflow.exec.${this.state.runId}.continue`, {type: 'continue'});
  }

  /** Stop listening for events and mark run as ended */
  stopRun(): void {
    if (this.subscription) {
      this.subscription.unsubscribe();
      this.subscription = null;
    }
    this.state.endRun();
  }

  /** Called when the graph structure changes — invalidates previous results */
  onGraphChanged(): void {
    this.graphVersion++;
    if (this.state.nodeStates.size > 0 && this.state.isStale(this.graphVersion)) {
      this.state.markAllStale();
      this.visualizer.markAllStale();
    }
  }

  /** Reset all visual state */
  resetVisuals(): void {
    this.visualizer.resetAllNodes();
    this.state.reset();
  }

  dispose(): void {
    this.stopRun();
    this.visualizer.resetAllNodes();
  }
}
