import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/** Interface for hierarchical data structure returned by parseNewick (from phylocanvas.gl Newick.parse_newick) */
export interface NodeType {
  name: string;
  branch_length?: number;
  children?: NodeType[];
}

export interface NodeCuttedType extends NodeType {
  cuttedLeafNameList: string[];
  cuttedChildren?: NodeType[];
}

export type TreeLeafDict = { [nodeName: string]: NodeType };
export type DataNodeDict = { [nodeName: string]: number };
export type NodeNameCallback = (nodeName: string) => void;

export enum TAGS {
  NEWICK = '.newick',
  NEWICK_JSON = '.newick_json',
  OBJECTS = '.newick_objects',
  NEWICK_LEAF_COL_NAME = '.newickLeafColumn',
}

export function isLeaf(node: NodeType) {
  return !node.children || node.children.length == 0;
}
