import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// https://gitlab.com/cgps/phylocanvas/phylocanvas.gl

import {NodeType} from './index';

export const TreeTypes = {
  Circular: 'cr',
  Diagonal: 'dg',
  Hierarchical: 'hr',
  Radial: 'rd',
  Rectangular: 'rc',
};

export const Shapes = {
  Chevron: 'chevron',
  ChevronInverted: 'chevron-inverted',
  ChevronLeft: 'chevron-left',
  ChevronRight: 'chevron-right',
  Circle: 'circle',
  Cross: 'cross',
  Diamond: 'diamond',
  Dot: 'dot',
  DoubleChevron: 'double-chevron',
  DoubleChevronInverted: 'double-chevron-inverted',
  DoubleChevronLeft: 'double-chevron-left',
  DoubleChevronRight: 'double-chevron-right',
  Heptagon: 'heptagon',
  HeptagonInverted: 'heptagon-inverted',
  Heptagram: 'heptagram',
  HeptagramInverted: 'heptagram-inverted',
  Hexagon: 'hexagon',
  Hexagram: 'hexagram',
  Octagon: 'octagon',
  Octagram: 'octagram',
  Pentagon: 'pentagon',
  PentagonInverted: 'pentagon-inverted',
  Pentagram: 'pentagram',
  PentagramInverted: 'pentagram-inverted',
  Plus: 'plus',
  Square: 'square',
  Star: 'pentagram',
  StarInverted: 'pentagram-inverted',
  Tetragram: 'tetragram',
  Triangle: 'triangle',
  TriangleInverted: 'triangle-inverted',
  TriangleLeft: 'triangle-left',
  TriangleRight: 'triangle-right',
  Wye: 'wye',
  WyeInverted: 'wye-inverted',
};

/** Represents a single tree node */
export type PhylocanvasTreeNode = NodeType & {
  branchLength: number,
  children: PhylocanvasTreeNode[],
  parent?: PhylocanvasTreeNode,
  id: string,
  isCollapsed: boolean,
  isHidden: boolean,
  isLeaf: boolean,
  name?: string,
  postIndex: number,
  preIndex: number,
  totalLeaves: number,
  totalNodes: number,
  totalSubtreeLength: number,
  visibleLeaves: number,
}

/** Based on phylocanvas.gl */
export function parseNewick(s: string): NodeType {
  const ancestors = [];
  let tree: NodeType = {name: ''};
  const tokens = s.split(/\s*(;|\(|\)|,|:)\s*/);
  for (let i = 0; i < tokens.length; i++) {
    const token = tokens[i];
    switch (token) {
    case '(': { // new children
      const subtree: NodeType = {name: ''};
      tree.children = [subtree];
      ancestors.push(tree);
      tree = subtree;
    }
      break;
    case ',': { // another branch
      const subtree: NodeType = {name: ''};
      ancestors[ancestors.length - 1].children!.push(subtree);
      tree = subtree;
    }
      break;
    case ')': // optional name next
      tree = ancestors.pop()!;
      break;
    case ':': // optional length next
      break;
    default:
      const x = tokens[i - 1];
      if (x == ')' || x == '(' || x == ',')
        tree.name = token;
      else if (x == ':')
        tree.branch_length = parseFloat(token);
    }
  }
  return tree;
}
