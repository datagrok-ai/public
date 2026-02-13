import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {studies} from '../utils/app-utils';
import {createAllMeasurementsDf, createVisitDayStrCol} from '../data-preparation/data-preparation';
import {DOMAIN, PLANNED_TRT_ARM, SEX, SUBJECT_ID, VISIT_DAY_STR} from '../constants/columns-constants';
import {createSubjectProfileView} from './subject-profile-view';
import {StudyTableViewParams} from '../utils/views-creation-utils';
import {restoreBrowsePanelOnRemoval} from '../utils/utils';

export function createMeasurementProfileTableView(studyId: string): StudyTableViewParams {
  let resDf = createAllMeasurementsDf(studyId);

  if (!resDf)
    return {df: DG.DataFrame.create()};

  const columnsFromDm = studies[studyId].domains.dm!.columns.names().filter((it) => it !== SUBJECT_ID);
  grok.data.joinTables(resDf, studies[studyId].domains.dm!, [SUBJECT_ID], [SUBJECT_ID], null,
    columnsFromDm, DG.JOIN_TYPE.LEFT, true);

  const intVisitDay = (resDf as DG.DataFrame).col('visit_day')!.convertTo( DG.TYPE.INT);
  (resDf as DG.DataFrame).columns.remove('visit_day');
  (resDf as DG.DataFrame).columns.add(intVisitDay);

  (resDf as DG.DataFrame).name = 'Measurements';
  (resDf as DG.DataFrame).rows.filter((row) => row.test == (resDf as DG.DataFrame).get('test', 0));

  let trellisPlot: DG.Viewer | null = null;
  const onTableViewAdded = async (tableView: DG.TableView) => {
    const splitBy = ui.input.choice('Split by', {
      items: ['Treatment arm', 'Subject'],
      value: 'Treatment arm',
      onValueChanged: () => {
        trellisPlot?.setOptions({
          innerViewerLook: {
            splitColumnNames: splitBy.value === 'Treatment arm' ? [PLANNED_TRT_ARM] : [SUBJECT_ID],
          },
        });
      },
    });

    tableView.setRibbonPanels([[splitBy.root]]);

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
      const look = trellisPlot?.getOptions().look;
      if (look?.viewerType === DG.VIEWER.LINE_CHART && look.innerViewerLook.splitColumnNames.includes(SUBJECT_ID)) {
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
    restoreBrowsePanelOnRemoval();
  };

  return {df: resDf, onTableViewAddedFunc: onTableViewAdded};
}
