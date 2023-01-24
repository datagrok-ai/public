import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {DataLoader} from './data-loader';
import {TreeDataFrame, MlbDataFrame} from '../types/dataframe';
import {cleanMlbNewick} from '../mlb-tree';
import {getVId} from './tree-stats';
import {ProgressIndicator, TaskBarProgressIndicator} from 'datagrok-api/dg';
import {getTreeHelper, ITreeHelper} from '@datagrok-libraries/bio/src/trees/tree-helper';
import {NodeType} from '@datagrok-libraries/bio/src/trees';
import {parseNewick} from '@datagrok-libraries/bio/src/trees/phylocanvas';

/** Opens TableView with dataframe containing VR ids presented in tree but not in MLB */
export async function checkVrForTreeUI(dl: DataLoader, pi: ProgressIndicator): Promise<void> {
  const th: ITreeHelper = await getTreeHelper();

  const res: { ag: string[], treeClone: string[], treeVId: string[], treeTree: string[], vId: string[] } =
    {ag: [], treeClone: [], treeVId: [], treeTree: [], vId: []};

  const allMlbDf: DG.DataFrame = await dl.loadMlbDf();
  const allVIds: Set<string> = new Set<string>(allMlbDf.getCol('v_id').toList());

  const agRowCount: number = dl.antigens.rowCount;
  const agRowStart: number = 0;
  const agRowStop: number = 100;
  for (let agRowI: number = agRowStart; agRowI < agRowStop; agRowI++) {
    const agName = dl.antigens.antigenCol.get(agRowI);

    const agMlbDf: MlbDataFrame = await dl.getMlbByAntigen(agName);
    try {
      // antigen VRs at MLB
      const agMlbVIds: { [vId: string]: number } = {};
      const agVIds: { [vId: string]: string } = {};
      const agMlbRowCount: number = agMlbDf.rowCount;
      for (let agMlbRowI = 0; agMlbRowI < agMlbRowCount; agMlbRowI++) {
        const vId: string = agMlbDf.vIdCol.get(agMlbRowI);
        agMlbVIds[vId] = agMlbRowI;
      }

      // antigen VRs at trees (igphyml data)
      const agTreeDf: TreeDataFrame = await dl.getTreeByAntigen(agName);
      const agTreeRowCount: number = agTreeDf.rowCount;
      for (let agTreeRowI = 0; agTreeRowI < agTreeRowCount; agTreeRowI++) {
        const treeClone = agTreeDf.cloneCol.get(agTreeRowI);
        const treeTree = agTreeDf.treeCol.get(agTreeRowI);
        const treeNewick = cleanMlbNewick(agTreeDf.treeCol.get(agTreeRowI));
        const treeRoot: NodeType = parseNewick(treeNewick);
        const treeLeafList: NodeType[] = th.getLeafList(treeRoot);
        const treeVIdList: string[] = treeLeafList.map((l) => getVId(l.name));
        for (const treeVId of treeVIdList) {
          if (!(treeVId in agMlbVIds)) {
            res.ag.push(agName);
            res.treeClone.push(treeClone);
            res.treeVId.push(treeVId);
            res.treeTree.push(treeTree);
            res.vId.push((allVIds.has(treeVId)) ? treeVId : '');
          }
        }
      }
    } finally {
      const progress: number = Math.floor(100 * agRowI / agRowCount);
      pi.update(progress, `check VR for tree, ${progress}%`);
      console.debug(`MLB: checkVrForTreeUI() antigen: ${agName}, progress: ${progress}%, res.length: ${res.ag.length}`);
    }
  }

  const resDf: DG.DataFrame = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('antigen', res.ag),
    DG.Column.fromStrings('v id', res.treeVId),
    DG.Column.fromStrings('CLONE', res.treeClone),
    DG.Column.fromStrings('TREE', res.treeTree),
    DG.Column.fromStrings('vId', res.vId)]);
  grok.shell.addTableView(resDf);
}
