import {LGraphNode, IWidget} from 'litegraph.js';
import * as DG from 'datagrok-api/dg';

/** Extended LGraphNode with FuncFlow-specific properties */
export interface FuncFlowNode extends LGraphNode {
  dgNodeType?: string;
  dgOutputType?: string;
  dgFunc?: DG.Func;
  dgFuncName?: string;
  dgRole?: string | null;
  inputWidgets?: Record<string, IWidget>;
  /** Number of pass-through output slots at the start (func nodes only) */
  _passthroughCount?: number;
}

/** Helper to access FuncFlow properties on a generic LGraphNode */
export function asFuncFlowNode(node: LGraphNode): FuncFlowNode {
  return node as FuncFlowNode;
}
