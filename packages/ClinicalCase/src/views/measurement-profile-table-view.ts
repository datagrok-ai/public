import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {studies} from '../utils/app-utils';
import {CDISC_STANDARD} from '../utils/types';
import {createVisitDayStrCol} from '../data-preparation/data-preparation';
import {DOMAIN, PLANNED_TRT_ARM, SEX, SUBJECT_ID, VISIT_DAY_STR} from '../constants/columns-constants';
import {createSubjectProfileView} from './subject-profile-view';

export function createMeasurementProfileTableView(studyId: string): any {
  const isSend = studies[studyId].config.standard === CDISC_STANDARD.SEND;
  let resDf: DG.DataFrame | null = null;

  //first creating a long dataframe with results from all applicable domains
  studies[studyId].domains.all().forEach((it) => {
    let requiredColumns = [SUBJECT_ID, DOMAIN, `${it.name.toUpperCase()}TEST`,
      `${it.name.toUpperCase()}STRESN`, `${it.name.toUpperCase()}STRESC`,
      `${it.name.toUpperCase()}STRESU`, `${it.name.toUpperCase()}DY`];
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
          .init((i) => `${it.get(`${it.name.toUpperCase()}TEST`, i)} ${
            it.get(catColName, i) ? `(${it.get(catColName, i)})` : ``}`);
        requiredColumns.push(fullTestNameCol);
        requiredColumns = requiredColumns.filter((i) => i !== `${it.name.toUpperCase()}TEST`);
      } else //create category column for compatibility in case it is missing
        it.columns.addNewString(catColName);
      requiredColumns.push(`${it.name.toUpperCase()}CAT`);

      const df = it.clone(null, requiredColumns);
      df.getCol(containsCatCol ? fullTestNameCol : `${it.name.toUpperCase()}TEST`).name = 'test';
      df.getCol(`${it.name.toUpperCase()}CAT`).name = 'category';
      df.getCol(`${it.name.toUpperCase()}STRESN`).name = 'result';
      df.getCol(`${it.name.toUpperCase()}DY`).name = 'visit_day';
      df.getCol(`${it.name.toUpperCase()}STRESC`).name = 'string_result';
      df.getCol(`${it.name.toUpperCase()}STRESU`).name = 'units';

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

  //set correct type to visit day column
  const intVisitDay = resDf.col('visit_day')!.convertTo( DG.TYPE.INT);
  resDf.columns.remove('visit_day');
  resDf.columns.add(intVisitDay);

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

    let subjectView: DG.DockNode | null = null;
    tableView.dataFrame.onCurrentRowChanged.subscribe((_) => {
      if (tableView.dataFrame.currentRowIdx === -1)
        return;
      const look = trellisPlot.getOptions().look;
      if (look.viewerType === DG.VIEWER.LINE_CHART && look.innerViewerLook.splitColumnNames.includes(SUBJECT_ID)) {
        const v = createSubjectProfileView(studyId, tableView.dataFrame,
          tableView.dataFrame.get(SUBJECT_ID, tableView.dataFrame.currentRowIdx));

        if (subjectView)
          (tableView as DG.TableView).dockManager.close(subjectView);
        subjectView = (tableView as DG.TableView).dockManager.dock(v, DG.DOCK_TYPE.DOWN);
      }
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
