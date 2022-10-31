// import {NodeType} from '@phylocanvas/phylocanvas.gl';

import {Observable} from 'rxjs';


/* Interface for hierarchical data structure returned by Newick.parse_newick */
export interface NodeType {
  name: string;
  branch_length?: number;
  children?: NodeType[];
}

export function isLeaf(node: NodeType) {
  return !node.children || node.children.length == 0;
}

export interface NodeCuttedType extends NodeType {
  cuttedLeafNameList: string[];
  cuttedChildren?: NodeType[];
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