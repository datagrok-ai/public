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

export type ClusterMatrix = {
  mergeRow1:Int32Array;
  mergeRow2:Int32Array;
  heightsResult:Float32Array;
}

export type TreeLeafDict = { [nodeName: string]: NodeType };
export type DataNodeDict = { [nodeName: string]: number };
export type NodeNameCallback = (nodeName: string) => void;
