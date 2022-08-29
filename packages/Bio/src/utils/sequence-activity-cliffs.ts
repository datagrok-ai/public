import {ITooltipAndPanelParams} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {getSimilarityFromDistance} from '@datagrok-libraries/utils/src/similarity-metrics';
import {AvailableMetrics} from '@datagrok-libraries/ml/src/typed-metrics';

export async function sequenceGetSimilarities(col: DG.Column, seq: string): Promise<DG.Column | null> {
  const stringArray = col.toList();
  const distances = new Array(stringArray.length).fill(0.0);
  for (let i = 0; i < stringArray.length; ++i)
    distances[i] = stringArray[i] ? getSimilarityFromDistance(AvailableMetrics['String']['Levenshtein'](stringArray[i], seq)) : 0;
  return DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'distances', distances);
}

export function drawSequences(params: ITooltipAndPanelParams) {
  params.line.mols.forEach((mol: number, index: number) => {
    ui.empty(params.hosts[index]);
    params.hosts[index].append(ui.divText(params.seqCol.get(mol)));
  });
}
