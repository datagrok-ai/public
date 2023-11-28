import * as DG from 'datagrok-api/dg';
import {MMP_COLOR} from './mmp-constants';

export function getMmpTrellisPlot(allPairsGrid: DG.Grid, activityMeanNames: Array<string>): DG.Viewer {
  const tp = DG.Viewer.fromType(DG.VIEWER.TRELLIS_PLOT, allPairsGrid.table, {
    xColumnNames: [allPairsGrid.table.columns.byIndex(0).name],
    yColumnNames: [allPairsGrid.table.columns.byIndex(1).name],
    viewerType: 'Summary',
    innerViewerLook: {
      columnNames: activityMeanNames,
      aggregations: [DG.STATS.MED, DG.STATS.MED],
      visualizationType: 'bars',
      colorColumnName: MMP_COLOR,
      colorAggrType: DG.STATS.MED,
      colorSchemes: [[DG.Color.categoricalPalette]],
      //invertColorScheme: true,
    },
  });

  return tp;
}
