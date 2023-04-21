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

export enum TreeColorNames {
  Main = 'Main',
  Light = 'Light',
  Current = 'Current',
  MouseOver = 'MouseOver',
  Selection = 'Selection',
}

export enum DistanceMetric{
  Euclidean = 'euclidean',
  Manhattan = 'manhattan',
}

export enum LinkageMethod{
  Single = 'single',
  Complete = 'complete',
  Average = 'average',
  Weighted = 'weighted',
  Centroid = 'centroid',
  Median = 'median',
  Ward = 'ward',
}


export const TreeDefaultPalette: { [name: string]: number } = {
  [TreeColorNames.Main]: DG.Color.categoricalPalette[12],
  [TreeColorNames.Light]: DG.Color.categoricalPalette[13],
  [TreeColorNames.Current]: DG.Color.currentRow,
  [TreeColorNames.MouseOver]: DG.Color.mouseOverRows,
  [TreeColorNames.Selection]: DG.Color.selectedRows,
};
