import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {studies} from '../utils/app-utils';
import {CDISC_STANDARD} from '../utils/types';
import {createVisitDayStrCol} from '../data-preparation/data-preparation';
import {DOMAIN, PLANNED_TRT_ARM, SUBJECT_ID, VISIT_DAY_STR} from '../constants/columns-constants';
import {addDomainFilters} from '../utils/utils';

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

      const df = it.clone(null, requiredColumns);
      df.getCol(`${it.name.toUpperCase()}TEST`).name = 'test';
      df.getCol(`${it.name.toUpperCase()}STRESN`).name = 'result';
      df.getCol(`${it.name.toUpperCase()}DY`).name = 'visit_day';
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
  const tests = resDf.col('test')!.categories;
  let selectedTest = tests[0];
  resDf.rows.match({test: selectedTest}).filter();
  const dfFortableView = resDf.clone(resDf.filter);
  const testChoiceInput = ui.input.choice('Test', {
    items: tests,
    value: selectedTest,
    onValueChanged: () => {
      selectedTest = testChoiceInput.value;
      resDf.rows.match({test: selectedTest}).filter();
      grok.shell.tv.dataFrame = resDf.clone(resDf.filter);
      //need to set boxPlot options since they are reset after table is changed
      boxPlot.setOptions({value: 'result', category1: PLANNED_TRT_ARM});
    },
  });
  let boxPlot: DG.Viewer | null = null;
  // onTableViewAdded function to add viewers
  const onTableViewAdded = async (tableView: DG.TableView) => {
    const ribbons = tableView.getRibbonPanels();
    ribbons.push([testChoiceInput.root]);
    tableView.setRibbonPanels(ribbons);

    const lineChart = await DG.Viewer.fromType(DG.VIEWER.LINE_CHART, tableView.dataFrame);
    lineChart.setOptions({yColumnNames: ['result'], yAggrTypes: ['avg'], whiskersType: 'Med | Q1, Q3'});
    //tableView.addViewer(lineChart);

    boxPlot = await DG.Viewer.fromType(DG.VIEWER.BOX_PLOT, tableView.dataFrame);
    boxPlot.setOptions({value: 'result', category1: PLANNED_TRT_ARM});
    //tableView.addViewer(boxPlot);

    tableView.dockManager.dock(boxPlot, DG.DOCK_TYPE.FILL);
    tableView.dockManager.dock(lineChart, DG.DOCK_TYPE.FILL);

    addDomainFilters(tableView, studyId);
  };

  return {df: dfFortableView, onTableViewAddedFunc: onTableViewAdded};
}
