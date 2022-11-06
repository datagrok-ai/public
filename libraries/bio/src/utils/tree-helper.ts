import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {NodeType} from '../types';


export interface ITreeHelper {
  toNewick(node: NodeType | null): string;

  getLeafList(node: NodeType): NodeType[];

  getNodeList(node: NodeType): NodeType[];

  treeFilterByLeaves(node: NodeType, leaves: { [name: string]: any }): NodeType | null;

  treeCutAsLeaves(node: NodeType, cutHeight: number, currentHeight?: number): NodeType[];

  treeCutAsTree(node: NodeType, cutHeight: number, currentHeight?: number): NodeType | null;

  setGridOrder(node: NodeType, grid: DG.Grid, leafColName: string): [NodeType, string[]];

  cutTreeToGrid(node:NodeType, cutHeight: number, dataDf: DG.DataFrame, leafColName: string, clusterColName: string): void;
}

export async function getTreeHelper(): Promise<ITreeHelper> {
  const funcList = DG.Func.find({package: 'PhyloTreeViewer', name: 'getTreeHelper'});
  if (funcList.length === 0)
    throw new Error('Package "PhyloTreeViewer"" must be installed for TreeHelper.');

  const res: ITreeHelper = (await funcList[0].prepare().call()).getOutputParamValue() as ITreeHelper;
  return res;
}