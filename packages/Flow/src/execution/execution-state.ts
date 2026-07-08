/** Execution state tracking for instrumented Flow runs.
 *
 * Node IDs are strings (Rete uses UUID-style strings, not LiteGraph's integer
 * IDs). The instrumented script emits events keyed by these strings. */

export enum NodeExecStatus {
  idle = 'idle',
  running = 'running',
  completed = 'completed',
  errored = 'errored',
  stale = 'stale',
}

export interface ValueSummary {
  type: 'dataframe' | 'column' | 'primitive' | 'object' | 'null' | 'graphics' | 'widget' | 'viewer';
  [key: string]: any;
}

export interface NodeExecState {
  status: NodeExecStatus;
  startTime?: number;
  endTime?: number;
  outputs?: Record<string, ValueSummary>;
  error?: string;
  stack?: string;
}

export interface ExecEvent {
  type: 'run-start' | 'node-start' | 'node-complete' | 'node-error' |
        'breakpoint-hit' | 'run-complete';
  nodeId: string;
  timestamp: number;
  outputs?: Record<string, ValueSummary>;
  error?: string;
  stack?: string;
  success?: boolean;
}

export class ExecutionState {
  runId: string = '';
  nodeStates: Map<string, NodeExecState> = new Map();
  isRunning: boolean = false;

  reset(): void {
    this.runId = '';
    this.nodeStates.clear();
    this.isRunning = false;
  }

  startRun(runId: string): void {
    this.reset();
    this.runId = runId;
    this.isRunning = true;
  }

  endRun(): void {
    this.isRunning = false;
  }

  setNodeStatus(nodeId: string, status: NodeExecStatus, data?: Partial<NodeExecState>): void {
    const existing = this.nodeStates.get(nodeId) ?? {status: NodeExecStatus.idle};
    this.nodeStates.set(nodeId, {...existing, status, ...data});
  }

  getNodeState(nodeId: string): NodeExecState | undefined {
    return this.nodeStates.get(nodeId);
  }

  /** Mark the given completed/errored nodes stale (a graph edit invalidated
   *  them); nodes outside the set — and idle/running ones — are untouched. */
  markStale(ids: Iterable<string>): void {
    for (const id of ids) {
      const state = this.nodeStates.get(id);
      if (state && (state.status === NodeExecStatus.completed || state.status === NodeExecStatus.errored))
        this.nodeStates.set(id, {...state, status: NodeExecStatus.stale});
    }
  }

  /** Drop a removed node's state entirely (the node no longer exists). */
  forgetNode(id: string): void {
    this.nodeStates.delete(id);
  }
}
