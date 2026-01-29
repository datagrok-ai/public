import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {studies} from '../utils/app-utils';
import {CDISC_STANDARD} from '../utils/types';
import {createVisitDayStrCol} from '../data-preparation/data-preparation';
import {DOMAIN, PLANNED_TRT_ARM, SEX, SUBJECT_ID, VISIT_DAY_STR} from '../constants/columns-constants';

export function createMeasurementProfileTableView(studyId: string): any {
  const isSend = studies[studyId].config.standard === CDISC_STANDARD.SEND;
  let resDf: DG.DataFrame | null = null;

  //first creating a long dataframe with results from all applicable domains
  studies[studyId].domains.all().forEach((it) => {
    let requiredColumns = [SUBJECT_ID, DOMAIN, `${it.name.toUpperCase()}TEST`,
      `${it.name.toUpperCase()}STRESN`, `${it.name.toUpperCase()}DY`];
    if (requiredColumns.every((colName) => it.columns.names().includes(colName))) {
      if (isSend) {
        createVisitDayStrCol(it);
        requiredColumns = requiredColumns.concat([VISIT_DAY_STR]);
      }
      // if there is a category column, create temporary test column name, containing both category and test name
      const catColName = `${it.name.toUpperCase()}CAT`;
      const containsCatCol = it.col(catColName);
      const fullTestNameCol = 'full_test_name';
      if (containsCatCol) {
        it.columns.addNewString(fullTestNameCol)
          .init((i) => `${it.get(`${it.name.toUpperCase()}TEST`, i)} (${it.get(catColName, i)})`);
        requiredColumns.push(fullTestNameCol);
        requiredColumns = requiredColumns.filter((i) => i !== `${it.name.toUpperCase()}TEST`);
      }

      const df = it.clone(null, requiredColumns);
      df.getCol(containsCatCol ? fullTestNameCol : `${it.name.toUpperCase()}TEST`).name = 'test';
      df.getCol(`${it.name.toUpperCase()}STRESN`).name = 'result';
      df.getCol(`${it.name.toUpperCase()}DY`).name = 'visit_day';

      // remove temporary full test name column
      if (containsCatCol)
        it.columns.remove(fullTestNameCol);

      if (!resDf)
        resDf = df;
      else
        resDf.append(df, true);
    }
  });

  //add dm domain fields for further filtering
  const columnsFromDm = studies[studyId].domains.dm.columns.names().filter((it) => it !== SUBJECT_ID);
  grok.data.joinTables(resDf, studies[studyId].domains.dm, [SUBJECT_ID], [SUBJECT_ID], null,
    columnsFromDm, DG.JOIN_TYPE.LEFT, true);

  resDf.name = 'Measurements';
  resDf.rows.filter((row) => row.test == resDf.get('test', 0));

  let trellisPlot: DG.Viewer | null = null;
  // onTableViewAdded function to add viewers
  const onTableViewAdded = async (tableView: DG.TableView) => {
    trellisPlot = await DG.Viewer.fromType(DG.VIEWER.TRELLIS_PLOT, tableView.dataFrame, {
      xColumnNames: [SEX],
      yColumnNames: ['test'],
      viewerType: DG.VIEWER.LINE_CHART,
      innerViewerLook: {
        showXAxis: true,
        aggrType: DG.STATS.AVG,
        splitColumnNames: [PLANNED_TRT_ARM],
      },
    });

    tableView.dockManager.dock(trellisPlot, DG.DOCK_TYPE.FILL);

    const fg = tableView.getFiltersGroup({createDefaultFilters: false});
    fg.updateOrAdd({
      type: DG.FILTER_TYPE.CATEGORICAL,
      column: 'test',
      columnName: 'test',
    });
  };

  return {df: resDf, onTableViewAddedFunc: onTableViewAdded};
}
