import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {GridNeighbor} from '@datagrok-libraries/gridext/src/ui/GridNeighbor';
import {NodeType} from '.';

export type TreeCutOptions = {
  min: number, max: number, dataDf: DG.DataFrame,
  clusterDf: DG.DataFrame, clusterColName: string,
}

export interface IDendrogramService {
  /** Inject Dendrogram tree to {@see grid}. Requires Dendrogram package. */
  injectTreeForGrid(
    grid: DG.Grid, treeRoot: NodeType, leafColName?: string, neighborWidth?: number, cut?: TreeCutOptions
  ): GridNeighbor;
}

export async function getDendrogramService(): Promise<IDendrogramService> {
  return await grok.functions.call('Dendrogram:getDendrogramService', {});
}
