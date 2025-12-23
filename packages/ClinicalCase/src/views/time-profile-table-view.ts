import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {studies} from '../utils/app-utils';
import {CDISC_STANDARD} from '../utils/types';
import {createVisitDayStrCol} from '../data-preparation/data-preparation';
import {DOMAIN, SUBJECT_ID, VISIT_DAY_STR} from '../constants/columns-constants';

export function createTimeProfileTableView(studyId: string): any {
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

  //add dm domain fields fro further filtering
  const columnsFromDm = studies[studyId].domains.dm.columns.names().filter((it) => it !== SUBJECT_ID);
  grok.data.joinTables(resDf, studies[studyId].domains.dm, [SUBJECT_ID], [SUBJECT_ID], null,
    columnsFromDm, DG.JOIN_TYPE.LEFT, true);

  // onTableViewAdded function to add viewers
  const onTableViewAdded = async (tableView: DG.TableView) => {
    const lineChart = await DG.Viewer.fromType(DG.VIEWER.LINE_CHART, resDf);
    lineChart.setOptions({yColumnNames: ['res_num'], yAggrTypes: ['avg'], whiskersType: 'Med | Q1, Q3'});

    tableView.addViewer(lineChart);
    const fg = tableView.getFiltersGroup();
    fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL,
      column: 'test',
      selected: [resDf.col('test').categories[0]]});
  };


  return {df: resDf, onTableViewAddedFunc: onTableViewAdded};
}
