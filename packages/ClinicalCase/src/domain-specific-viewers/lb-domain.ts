import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {LAB_TEST, LAB_TEST_CAT, PLANNED_TRT_ARM, ACT_TRT_ARM} from '../constants/columns-constants';
import {calculateLBBaselineColumns} from '../data-preparation/data-preparation';

export function createLbHeatmap(lbDomain: DG.DataFrame): DG.Viewer {
  // Check if required columns exist
  const labTestCol = lbDomain.col(LAB_TEST);
  const labCatCol = lbDomain.col(LAB_TEST_CAT);
  const armCol = lbDomain.col(PLANNED_TRT_ARM) || lbDomain.col(ACT_TRT_ARM);
  const pctChgCol = lbDomain.col('LB_PCT_CHG');

  if (!labTestCol || !armCol)
    return;


  // Ensure LB_PCT_CHG column exists
  if (!pctChgCol)
    calculateLBBaselineColumns(lbDomain);

  const finalPctChgCol = lbDomain.col('LB_PCT_CHG');
  if (!finalPctChgCol)
    return;


  // Create pivot dataframe
  // Group by: lab test + lab category (if exists)
  const groupByCols = labCatCol ? [LAB_TEST_CAT, LAB_TEST] : [LAB_TEST];

  // Create pivoted dataframe: rows = lab test + category, columns = treatment arms, values = mean % change
  const pivotedDf = lbDomain
    .groupBy(groupByCols)
    .pivot(armCol.name)
    .avg('LB_PCT_CHG')
    .aggregate();

  // Clean up column names (remove " avg(LB_PCT_CHG)" suffix)
  const cols = pivotedDf.columns.names();
  for (const colName of cols) {
    if (colName.endsWith(' avg(LB_PCT_CHG)'))
        pivotedDf.col(colName)!.name = colName.replace(' avg(LB_PCT_CHG)', '');
  }

  // Create heatmap viewer
  const heatmapViewer = DG.Viewer.fromType(DG.VIEWER.HEAT_MAP, pivotedDf);

  grok.data.linkTables(lbDomain, pivotedDf,
    groupByCols, groupByCols,
    [DG.SYNC_TYPE.FILTER_TO_FILTER]);

  return heatmapViewer;
}
