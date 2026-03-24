/** Execution state tracking for instrumented Flow runs */

export enum NodeExecStatus {
  idle = 'idle',
  running = 'running',
  completed = 'completed',
  errored = 'errored',
  stale = 'stale',
}

export interface ValueSummary {
  type: 'dataframe' | 'column' | 'primitive' | 'object' | 'null' | 'graphics';
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

/** Execution event payloads fired by instrumented scripts */
export interface ExecEvent {
  type: 'run-start' | 'node-start' | 'node-complete' | 'node-error' |
        'breakpoint-hit' | 'run-complete';
  nodeId: number;
  timestamp: number;
  outputs?: Record<string, ValueSummary>;
  error?: string;
  stack?: string;
  success?: boolean;
}

/** Tracks execution state for all nodes in a single run */
export class ExecutionState {
  runId: string = '';
  nodeStates: Map<number, NodeExecState> = new Map();
  isRunning: boolean = false;
  graphVersionAtRun: number = 0;

  reset(): void {
    this.runId = '';
    this.nodeStates.clear();
    this.isRunning = false;
  }

  startRun(runId: string, graphVersion: number): void {
    this.reset();
    this.runId = runId;
    this.isRunning = true;
    this.graphVersionAtRun = graphVersion;
  }

  endRun(): void {
    this.isRunning = false;
  }

  setNodeStatus(nodeId: number, status: NodeExecStatus, data?: Partial<NodeExecState>): void {
    const existing = this.nodeStates.get(nodeId) || {status: NodeExecStatus.idle};
    this.nodeStates.set(nodeId, {...existing, status, ...data});
  }

  getNodeState(nodeId: number): NodeExecState | undefined {
    return this.nodeStates.get(nodeId);
  }

  /** Mark all completed/errored nodes as stale (graph changed since run) */
  markAllStale(): void {
    for (const [nodeId, state] of this.nodeStates) {
      if (state.status === NodeExecStatus.completed || state.status === NodeExecStatus.errored)
        this.nodeStates.set(nodeId, {...state, status: NodeExecStatus.stale});
    }
  }

  isStale(currentGraphVersion: number): boolean {
    return this.graphVersionAtRun !== currentGraphVersion;
  }
}
