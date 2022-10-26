// import {NodeType} from '@phylocanvas/phylocanvas.gl';

import {Observable} from 'rxjs';

export interface INode<TNode extends INode<TNode>> {
  name: string;
  branch_length?: number;
  children: TNode[];

  clone(): TNode;
}

export function isLeaf<TNode extends INode<TNode>>(node: TNode) {
  return !node.children || node.children.length == 0;
}


export class Node implements INode<Node> {
  name: string;
  children: Node[];
  branch_length?: number;

  constructor(name: string, branch_length?: number, children: Node[] = [],) {
    this.name = name;
    this.branch_length = branch_length;
    this.children = children;
  }

  /** Shallow copy, copies children list but children itself remains reference to the same */
  clone() {
    return new Node(this.name, this.branch_length, [...this.children],);
  }
}

export type Monomer = {
  at: { [R: string]: string },
  id: string,
  m: string,
  n: string,
  na: string,
  rs: number;
}

//expected types: HELM_AA, HELM_BASE, HELM_CHEM, HELM_LINKER, HELM_SUGAR
export interface IMonomerLib {
  get(monomerType: string, monomerName: string): Monomer | null;

  // TODO:
  get onChanged(): Observable<any>;
}