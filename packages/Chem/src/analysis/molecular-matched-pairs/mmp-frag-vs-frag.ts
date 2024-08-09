import * as DG from 'datagrok-api/dg';
import {MMP_NAMES} from './mmp-constants';
import {PaletteCodes} from './mmp-mol-rendering';

export function getMmpTrellisPlot(allPairsGrid: DG.Grid, activityMeanNames: Array<string>,
  palette: PaletteCodes): DG.Viewer {
  const schemes = new Array<any>(activityMeanNames.length);
  for (let i = 0; i < activityMeanNames.length; i++)
    schemes[i] = [palette.numerical[i]];

  const tp = DG.Viewer.fromType(DG.VIEWER.TRELLIS_PLOT, allPairsGrid.table, {
    xColumnNames: [allPairsGrid.table.columns.byIndex(0).name],
    yColumnNames: [allPairsGrid.table.columns.byIndex(1).name],
    viewerType: 'Summary',
    innerViewerLook: {
      columnNames: activityMeanNames,
      aggrColumnName: '',
      aggregation: DG.STATS.MED,
      aggregations: [DG.STATS.MED, DG.STATS.MED, DG.STATS.MED],
      visualizationType: 'bars',
      colorColumnName: MMP_NAMES.COLOR,
      colorAggrType: DG.STATS.MED,
      colorSchemes: schemes,
    },
  });

  return tp;
}
