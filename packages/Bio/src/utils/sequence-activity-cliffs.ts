import {IDrawTooltipParams} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {getSimilarityFromDistance} from '@datagrok-libraries/utils/src/similarity-metrics';
import {AvailableMetrics} from '@datagrok-libraries/ml/src/typed-metrics';

export async function sequenceGetSimilarities(col: DG.Column, seq: string): Promise<DG.Column | null> {
  const stringArray = col.toList();
  const distances = new Array(stringArray.length).fill(0.0);
  for (let i = 0; i < stringArray.length; ++i)
    distances[i] = getSimilarityFromDistance(AvailableMetrics['String']['Levenshtein'](stringArray[i], seq));
  return DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'distances', distances);
}

export function drawTooltip(params: IDrawTooltipParams) {
  params.tooltips[params.line.id] = ui.divH([]);
  const columnNames = ui.divV([
    ui.divText('sequense'),
    ui.divText(params.activity.name),
  ]);
  columnNames.style.fontWeight = 'bold';
  columnNames.style.display = 'flex';
  columnNames.style.justifyContent = 'space-between';
  params.tooltips[params.line.id].append(columnNames);
  params.line.mols.forEach((mol: number) => {
    const seq = ui.divText(params.df.get(params.seqCol.name, mol));
    const activity = ui.divText(params.df.get(params.activity.name, mol).toFixed(2));
    activity.style.display = 'flex';
    activity.style.justifyContent = 'left';
    activity.style.paddingLeft = '30px';
    params.tooltips[params.line.id].append(ui.divV([
      seq,
      activity,
    ], {style: {paddingLeft: '5px'}}));
  });
}
