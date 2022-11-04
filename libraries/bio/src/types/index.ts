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
  symbol: string,
  name: string,
  naturalAnalog: string,
  molfile: string,
  polymerType: string,
  monomerType: string,
  rgroups: {capGroupSmiles: string, alternateId: string, capGroupName: string, label: string }[],
  data: {[property: string]: any}
};

export interface IMonomerLib {
  getMonomer(monomerType: string, monomerName: string): Monomer | null;
  getMonomersByType(type: string): {[symbol: string]: string} | null;
  getTypes(): string[];
  update(monomers: { [type: string]: { [name: string]: Monomer } }): void;
  get onChanged(): Observable<any>;
}
