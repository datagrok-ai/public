/** Maps execution state to node visuals (`dgStatus` → `data-status` → CSS)
 *  and connection styling via `FlowEditor.setConnectionStatus`. */

import {FlowEditor} from '../rete/flow-editor';
import {FlowNode} from '../rete/scheme';
import {NodeExecStatus} from './execution-state';

interface FlowNodeWithStatus extends FlowNode {
  dgStatus?: NodeExecStatus;
}

/** Plain-language label shown under the node title. */
export function statusLabel(status: NodeExecStatus, detail?: string): string {
  switch (status) {
  case NodeExecStatus.running:   return 'Running…';
  case NodeExecStatus.completed: return detail ? `Done · ${detail}` : 'Done';
  case NodeExecStatus.errored:   return detail ? `Error — ${detail}` : 'Error';
  case NodeExecStatus.stale:     return 'Out of date';
  default:                       return '';
  }
}

/** The message's first sentence, capped to fit a node's one-line status. */
export function errorSummary(message: string | undefined): string {
  const first = String(message ?? '').split('\n')[0].split(/(?<=\.)\s/)[0].trim();
  return first.length > 120 ? first.slice(0, 117) + '…' : first;
}

/** Strip the platform's exception wrapper and map server-side file paths back
 *  to the user's `System:AppData/…` form. */
export function normalizeErrorMessage(message: string | undefined): string {
  let m = String(message ?? '').trim();
  const wrapper = /^Operation caused an exception\s*\(([\s\S]*)\)$/.exec(m);
  if (wrapper) m = wrapper[1].trim();
  // A failed file read races into the d42 deserializer, whose parse error carries no path.
  if (/outside the bounds of the DataView/i.test(m))
    m = 'The server response was not a table — check the file path and that the file is a readable format';
  return m.replace(/['"]?(?:\/[^\s'"]+)+\/packages\/data\/([^\s'"]+)['"]?/g, 'System:AppData/$1');
}

export class ExecutionVisualizer {
  private flow: FlowEditor;
  private trackedNodes = new Set<string>();

  constructor(flow: FlowEditor) {
    this.flow = flow;
  }

  highlightNode(nodeId: string, status: NodeExecStatus, detail?: string): void {
    const node = this.flow.getNodeById(nodeId) as FlowNodeWithStatus | undefined;
    if (!node) return;
    node.dgStatus = status;
    node.statusText = statusLabel(status, detail);
    this.trackedNodes.add(nodeId);
    void this.flow.updateNode(nodeId);
    this.propagateToConnections(nodeId, status);
  }

  private propagateToConnections(nodeId: string, status: NodeExecStatus): void {
    const incoming = this.flow.getConnections().filter((c) => c.target === nodeId);
    let connStatus: 'idle' | 'active' | 'completed' | 'errored' | 'stale';
    switch (status) {
    case NodeExecStatus.running:   connStatus = 'active'; break;
    case NodeExecStatus.completed: connStatus = 'completed'; break;
    case NodeExecStatus.errored:   connStatus = 'errored'; break;
    case NodeExecStatus.stale:     connStatus = 'stale'; break;
    default:                       connStatus = 'idle';
    }
    for (const c of incoming) this.flow.setConnectionStatus(c.id, connStatus);
  }

  /** Completed nodes flip to stale but KEEP their status text — resetting to
   *  idle here made every card flash blank between runs. */
  beginRun(): void {
    for (const id of this.trackedNodes) {
      const node = this.flow.getNodeById(id) as FlowNodeWithStatus | undefined;
      if (!node || node.dgStatus !== NodeExecStatus.completed) continue;
      node.dgStatus = NodeExecStatus.stale;
      void this.flow.updateNode(id);
    }
    this.flow.resetConnectionStatuses();
  }

  resetAllNodes(): void {
    for (const id of this.trackedNodes) {
      const node = this.flow.getNodeById(id) as FlowNodeWithStatus | undefined;
      if (node) {
        node.dgStatus = NodeExecStatus.idle;
        node.statusText = '';
        void this.flow.updateNode(id);
      }
    }
    this.trackedNodes.clear();
    this.flow.resetConnectionStatuses();
  }

  markStale(ids: Iterable<string>): void {
    for (const id of ids) {
      if (!this.trackedNodes.has(id)) continue;
      const node = this.flow.getNodeById(id) as FlowNodeWithStatus | undefined;
      if (!node || node.dgStatus === NodeExecStatus.idle || node.dgStatus === undefined) continue;
      node.dgStatus = NodeExecStatus.stale;
      node.statusText = statusLabel(NodeExecStatus.stale);
      void this.flow.updateNode(id);
      this.propagateToConnections(id, NodeExecStatus.stale);
    }
  }

  /** Nodes a halted run never reached: stale visuals + a line naming the failure. */
  markSkipped(ids: Iterable<string>, failedLabel: string): void {
    for (const id of ids) {
      const node = this.flow.getNodeById(id) as FlowNodeWithStatus | undefined;
      if (!node) continue;
      node.dgStatus = NodeExecStatus.stale;
      node.statusText = `Skipped — "${failedLabel}" failed`;
      this.trackedNodes.add(id);
      void this.flow.updateNode(id);
    }
  }

  forgetNode(id: string): void {
    this.trackedNodes.delete(id);
  }
}
