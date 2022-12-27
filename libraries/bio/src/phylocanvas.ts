// https://gitlab.com/cgps/phylocanvas/phylocanvas.gl

import {NodeType} from './types/index';

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

/*

// https://en.wikipedia.org/wiki/Tree_traversal#Post-order
function getPostorderTraversal<TNode extends NodeType>(rootNode: TNode): TNode[] {
  const nodes: TNode[] = [];
  const queue: TNode[] = [rootNode];

  while (queue.length) {
    const node: TNode = queue.pop()!;
    if (node.children)
      Array.prototype.push.apply(queue, node.children);
    nodes.push(node);
  }

  return nodes.reverse();
}

// https://en.wikipedia.org/wiki/Tree_traversal#Pre-order
function getPreorderTraversal<TNode extends NodeType>(rootNode: TNode): TNode[] {
  const nodes: TNode[] = [];
  const queue: TNode[] = [rootNode];

  while (queue.length) {
    const node: TNode = queue.shift()!;
    nodes.push(node);
    if (node.children)
      Array.prototype.unshift.apply(queue, node.children);
  }

  return nodes;
}

export function treeTraversal(rootNode: PhylocanvasTreeNode, {trimQuotes = true} = {}) {
  performance.mark('getPostorderTraversal');
  const postorderTraversal = getPostorderTraversal(rootNode);
  performance.measure('    getPostorderTraversal', 'getPostorderTraversal');
  performance.mark('getPreorderTraversal');
  const preorderTraversal = getPreorderTraversal(rootNode);
  performance.measure('    getPreorderTraversal', 'getPreorderTraversal');

  // Detect cladograms
  const isCladogram = postorderTraversal.every((x) => (x.branchLength || x.branch_length || 0) === 0);
  if (isCladogram) {
    rootNode.branchLength = 0;
    for (let nodeIndex = 0; nodeIndex < preorderTraversal.length; nodeIndex++) {
      const node = preorderTraversal[nodeIndex];
      if (node.children) {
        for (const child of node.children)
          child.branchLength = node.branchLength + 1;
      }
    }
  }

  performance.mark('bottom-up traversal');
  // bottom-up traversal starting from leaves to root
  for (let nodeIndex = 0; nodeIndex < postorderTraversal.length; nodeIndex++) {
    const node: PhylocanvasTreeNode = postorderTraversal[nodeIndex];
    node.postIndex = nodeIndex;
    node.isLeaf = !Array.isArray(node.children);
    node.branchLength = Math.abs(node.branchLength || node.branch_length || 0);
    delete node.branch_length;
    if (node.isLeaf && typeof node.name === 'string') {
      if (trimQuotes)
        node.id = node.name.trim().replace(/^['"]|['"]$/g, '');
      else
        node.id = node.name;

      // @ts-ignore
      delete node.name;
    }
    node.totalNodes = 1;
    node.totalLeaves = 1;
    node.totalSubtreeLength = 0;
    if (!node.isLeaf) {
      node.totalNodes = 1;
      node.totalLeaves = 0;
      let totalSubtreeLength = 0;
      for (const child of node.children) {
        node.totalNodes += child.totalNodes;
        node.totalLeaves += child.totalLeaves;
        if (child.totalSubtreeLength + child.branchLength > totalSubtreeLength)
          totalSubtreeLength = child.totalSubtreeLength + child.branchLength;
        child.parent = node;
      }
      node.totalSubtreeLength = totalSubtreeLength;
    }
  }
  performance.measure('    bottom-up traversal', 'bottom-up traversal');

  performance.mark('top-down traversal');
  const nodeById: { [id: string]: PhylocanvasTreeNode } = {};
  // top-down traversal starting from root to leaves
  for (let nodeIndex = 0; nodeIndex < preorderTraversal.length; nodeIndex++) {
    const node = preorderTraversal[nodeIndex];
    node.preIndex = nodeIndex;
    if (!node.id)
      node.id = nodeIndex.toString();
    nodeById[node.id] = node;
    node.visibleLeaves = node.totalLeaves;
    node.isCollapsed = false;
    node.isHidden = false;
  }
  performance.measure('    top-down traversal', 'top-down traversal');

  return {
    nodeById,
    rootNode,
    postorderTraversal,
    preorderTraversal,
  };
}

/**/

export function parseNewick(s: string): any {
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
