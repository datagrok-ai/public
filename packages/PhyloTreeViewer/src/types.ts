import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';

import {PhylotreeNode} from 'phylotree';
import * as Phylotree from 'phylotree';


export class ExtNode extends bio.Node {
  cuttedLeafNameList?: string[];

  constructor(name: string, branch_length?: number, cuttedLeafNameList?: string[]) {
    super(name, branch_length, []);
    this.cuttedLeafNameList = cuttedLeafNameList;
  }
}


// export class TreeNode<TChild extends PhylotreeNode, PhylocanvasNode>
//   implements PhylotreeNode, PhylocanvasNode {
//
//   _name: string;
//   _children: TChild[];
//   _branchLength: number;
//
//   get name(): string { return this._name; }
//
//   get children(): TreeNode[] { return this._children; }
//
//
//   constructor(name: string, branchLength: number, children: TreeNode[]) {
//
//     this._name = name;
//     this._branchLength = branchLength;
//     this._children = children;
//   }
// }