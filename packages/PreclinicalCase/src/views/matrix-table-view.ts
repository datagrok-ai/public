import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {studies} from '../utils/app-utils';
import {createVisitDayStrCol} from '../data-preparation/data-preparation';
import {DOMAIN, PLANNED_TRT_ARM, SUBJECT_ID, VISIT_DAY_STR} from '../constants/columns-constants';
import {addDomainFilters, restoreBrowsePanelOnRemoval} from '../utils/utils';
import {StudyTableViewParams} from '../utils/views-creation-utils';

export function createMatrixTableView(studyId: string): StudyTableViewParams {
  let resDf: DG.DataFrame | null = null;
  let pivotedDF = DG.DataFrame.create();

  studies[studyId].domains.all().forEach((it) => {
    let requiredColumns = [SUBJECT_ID, DOMAIN, `${it.name.toUpperCase()}TEST`,
      `${it.name.toUpperCase()}STRESN`, `${it.name.toUpperCase()}DY`];
    if (requiredColumns.every((colName) => it.columns.names().includes(colName))) {
      createVisitDayStrCol(it);
      requiredColumns = requiredColumns.concat([VISIT_DAY_STR]);

      const df = it.clone(null, requiredColumns);
      df.getCol(`${it.name.toUpperCase()}TEST`).name = 'test';
      df.getCol(`${it.name.toUpperCase()}STRESN`).name = 'res_num';
      df.getCol(`${it.name.toUpperCase()}DY`).name = 'visit_day';
      if (!resDf)
        resDf = df;
      else
        resDf.append(df, true);
    }
  });

  if (resDf) {
    const matrixDataframe = (resDf as DG.DataFrame)
      .groupBy([SUBJECT_ID, 'visit_day', VISIT_DAY_STR])
      .pivot('test')
      .avg('res_num')
      .aggregate();

    const cols = matrixDataframe.columns.names();
    for (const colName of cols) {
      if (colName.endsWith(' avg(res_num)'))
        matrixDataframe.col(colName)!.name = colName.replace(' avg(res_num)', '');
    }

    const columnsFromDm = studies[studyId].domains.dm!.columns.names().filter((it) => it !== SUBJECT_ID);
    grok.data.joinTables(matrixDataframe, studies[studyId].domains.dm!, [SUBJECT_ID], [SUBJECT_ID], null,
      columnsFromDm, DG.JOIN_TYPE.LEFT, true);
    pivotedDF = matrixDataframe;
    pivotedDF.name = 'Matrix dataframe';
  }

  const onTableViewAdded = async (tableView: DG.TableView) => {
    const matrixViewer = await DG.Viewer.fromType(DG.VIEWER.MATRIX_PLOT, pivotedDF, {cellPlotType: 'Scatter plot'});
    matrixViewer.setOptions({
      xColumnNames: matrixViewer.getOptions().look.xColumnNames
        .filter((it: string) => ![SUBJECT_ID, 'visit_day', VISIT_DAY_STR].includes(it)),
      yColumnNames: matrixViewer.getOptions().look.yColumnNames
        .filter((it: string) => ![SUBJECT_ID, 'visit_day', VISIT_DAY_STR].includes(it)),
    });
    if (pivotedDF.columns.names().includes(PLANNED_TRT_ARM)) {
      const l = matrixViewer.getOptions().look.innerViewerLook;
      l.colorColumnName = PLANNED_TRT_ARM;
      matrixViewer.setOptions({
        innerViewerLook: l,
      });
    }
    tableView.dockManager.dock(matrixViewer, DG.DOCK_TYPE.FILL);
    addDomainFilters(tableView, studyId);
    restoreBrowsePanelOnRemoval();
  };

  return {df: pivotedDF, onTableViewAddedFunc: onTableViewAdded};
}
