/** Maps execution state to node DOM via the `dgStatus` field on FlowNode and
 *  to connection styling via `FlowEditor.setConnectionStatus`.
 *
 *  Node side: the React node component reads `dgStatus` and writes it as a
 *  `data-status` attribute on the rendered `<div class="ff-node">`. CSS handles
 *  the visuals (status circle, pulse animation, body tint).
 *
 *  Connection side: incoming edges of a running node go to `active` (animated
 *  marching dashes); when the node completes the same edges flip to
 *  `completed` (steady, source-colored); on error they go `errored` (red).
 *  This lets a viewer see the data-flow front advance through the graph. */

import {FlowEditor} from '../rete/flow-editor';
import {FlowNode} from '../rete/scheme';
import {NodeExecStatus} from './execution-state';

interface FlowNodeWithStatus extends FlowNode {
  dgStatus?: NodeExecStatus;
}

/** Plain-language label shown under the node title. `detail` is an optional
 *  short data summary for a completed node (e.g. "1,204 × 8"). */
export function statusLabel(status: NodeExecStatus, detail?: string): string {
  switch (status) {
  case NodeExecStatus.running:   return 'Running…';
  case NodeExecStatus.completed: return detail ? `Done · ${detail}` : 'Done';
  case NodeExecStatus.errored:   return 'Error';
  case NodeExecStatus.stale:     return 'Out of date';
  default:                       return '';
  }
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

  /** Mirror a node's status onto its incoming connections (the edges that
   *  delivered data into this step). */
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

  /** Flip only the given nodes to "Out of date" (and their incoming edges to
   *  stale). Nodes outside the set keep their completed/errored visuals — a
   *  graph edit invalidates its downstream cone, not the whole canvas. */
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

  /** Stop tracking a removed node (its element is gone; nothing to repaint). */
  forgetNode(id: string): void {
    this.trackedNodes.delete(id);
  }
}
