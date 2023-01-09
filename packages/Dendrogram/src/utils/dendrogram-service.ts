import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {GridNeighbor} from '@datagrok-libraries/gridext/src/ui/GridNeighbor';
import {IDendrogramService, NodeType} from '@datagrok-libraries/bio';
import {TreeCutOptions} from '@datagrok-libraries/bio';
import {injectTreeForGridUI2} from '../viewers/inject-tree-for-grid2';


export class DendrogramService implements IDendrogramService {
  injectTreeForGrid(
    grid: DG.Grid, treeRoot: NodeType, leafColName?: string, neighborWidth?: number, cut?: TreeCutOptions
  ): GridNeighbor {
    const neighborWidthVal = neighborWidth ?? 100;
    return injectTreeForGridUI2(grid, treeRoot, leafColName, neighborWidthVal, cut);
  }
}
