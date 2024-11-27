import * as DG from 'datagrok-api/dg';
import {MMP_NAMES} from './mmp-constants';
import {PaletteCodes} from './palette';
import * as grok from 'datagrok-api/grok';

export function getMmpTrellisPlot(fGrid: DG.Grid, activityMeanNames: Array<string>,
  palette: PaletteCodes): DG.Viewer {
  const schemes = new Array<any>(activityMeanNames.length);
  for (let i = 0; i < activityMeanNames.length; i++)
    schemes[i] = [palette.numerical[i]];
  const aggregations = activityMeanNames.map((_) => DG.STATS.MED);

  const tp = DG.Viewer.fromType(DG.VIEWER.TRELLIS_PLOT, fGrid.table, {
    xColumnNames: [MMP_NAMES.FROM],
    yColumnNames: [MMP_NAMES.TO],
    viewerType: 'Summary',
    innerViewerLook: {
      columnNames: activityMeanNames,
      aggregations: aggregations,
      visualizationType: 'bars',
      colorColumnName: MMP_NAMES.COLOR,
      colorAggrType: DG.STATS.MED,
      colorSchemes: schemes,
    },
  });

  return tp;
}
