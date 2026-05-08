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

export class ExecutionController {
  state: ExecutionState;
  private visualizer: ExecutionVisualizer;
  private flow: FlowEditor;
  private subscription: {unsubscribe(): void} | null = null;
  private graphVersion = 0;
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

  private executeInstrumented(settings: ScriptSettings, debug: boolean): void {
    const errors = validateGraph(this.flow);
    if (errors.some((e) => e.severity === 'error')) {
      const msgs = errors.filter((e) => e.severity === 'error').map((e) => e.message).join('\n');
      grok.shell.error('Validation errors:\n' + msgs);
      return;
    }

    this.stopRun();

    const runId = crypto.randomUUID();
    this.state.startRun(runId, this.graphVersion);
    this.visualizer.resetAllNodes();

    const channel = `funcflow.exec.${runId}`;
    this.subscription = grok.events.onCustomEvent(channel).subscribe((event: ExecEvent) => {
      this.handleEvent(event);
    });

    const options: EmitOptions = {
      instrumented: true,
      runId,
      enableBreakpoints: debug,
      haltOnError: true,
    };

    try {
      const script = emitScript(this.flow, settings, options);
      const typeHints = this.getOutputTypeHints();
      const func = DG.Script.create(script);
      const fc = func.prepare();
      const onComplete = (outputs: Record<string, any>) => {
        this.outputPreview.showOutputs(outputs, typeHints);
      };

      if (func.inputs.length === 0)
        fc.call(undefined, undefined, {processed: true}).then(() => onComplete(fc.outputs));
      else {
        fc.getEditor(false).then((e: HTMLElement) => {
          ui.dialog({title: settings.name}).add(e).show().onOK(async () => {
            await fc.call(undefined, undefined, {processed: true});
            onComplete(fc.outputs);
          });
        });
      }
    } catch (e: any) {
      grok.shell.error(`Script generation failed: ${e.message}`);
      this.stopRun();
    }
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
      this.visualizer.highlightNode(event.nodeId, NodeExecStatus.completed);
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
      this.onRunEnd?.(event.success === true);
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

  /** Bumped every time the graph changes; if results are stale, mark them so. */
  onGraphChanged(): void {
    this.graphVersion++;
    if (this.state.nodeStates.size > 0 && this.state.isStale(this.graphVersion)) {
      this.state.markAllStale();
      this.visualizer.markAllStale();
    }
  }

  resetVisuals(): void {
    this.visualizer.resetAllNodes();
    this.state.reset();
  }

  /** Map output paramName → declared DG type (for output-preview classification). */
  getOutputTypeHints(): Record<string, string> {
    const hints: Record<string, string> = {};
    for (const node of this.flow.getNodes()) {
      if (node.dgNodeType !== 'output') continue;
      const paramName = node.properties['paramName'] as string | undefined;
      const outputType = (node.properties['outputType'] as string | undefined) ?? node.dgOutputType;
      if (paramName && outputType) hints[paramName] = outputType;
    }
    return hints;
  }

  dispose(): void {
    this.stopRun();
    this.visualizer.resetAllNodes();
  }
}
