/** Utility helpers around `FlowEditor` data. With Rete the editor exposes
 *  `getNodes()` / `getConnections()` directly, so this file is mostly a thin
 *  shim kept for parity with the old codebase shape. */

import {FlowEditor} from '../rete/flow-editor';
import {FlowNode} from '../rete/scheme';

export function getGraphNodes(flow: FlowEditor): FlowNode[] {
  return flow.getNodes();
}
