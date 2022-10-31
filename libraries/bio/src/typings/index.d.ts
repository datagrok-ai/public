declare module '@phylocanvas/phylocanvas.gl' {
  import {Deck} from '@deck.gl/core/typed';
  import {NodeType} from '../types';

  export class PhylocanvasGL {
    get deck(): Deck;

    get view(): HTMLDivElement;

    constructor(element: HTMLElement, props: { [propName: string]: any });

    render(): void;

    resume(): void;

    setProps(props: { [propName: string]: any }): void;

    selectNode: (nodeOrId: any, append: boolean = false) => void;

    destroy(): void;

    getBranchScale(...arguments): number;
  }

  export const TreeTypes;
  export const Shapes: { [key: string]: string };

  module Newick {
    function parse_newick(newick: string): NodeType;
  }

  module Utils {
    function treeTraversal(root: PhylocanvasTreeNode);
  }

  // TODO: add test for these properties existing.
  /**
   * Represents a single tree node.
   */
  export interface PhylocanvasTreeNode {
    branchLength: 0;
    children: PhylocanvasTreeNode[];
    id: string;
    isCollapsed: boolean;
    isHidden: boolean;
    isLeaf: boolean;
    name: string;
    postIndex: number;
    preIndex: number;
    totalLeaves: number;
    totalNodes: number;
    totalSubtreeLength: number;
    visibleLeaves: number;
  }
}
