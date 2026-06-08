import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {studies} from '../utils/app-utils';
import {SUBJECT_ID, DOMAIN} from '../constants/columns-constants';

export function createSubjectProfileView(studyId: string, measurementsDf: DG.DataFrame, subjectId?: string):
HTMLDivElement {
  const div = ui.div('', {style: {overflowY: 'scroll'}});

  const study = studies[studyId];

  const subjectIds = study.domains.dm.col(SUBJECT_ID)!.categories;
  let selectedSubjectId = subjectId ?? (subjectIds.length > 0 ? subjectIds[0] : null);

  const subjectIdInput = ui.input.choice('Subject ID', {
    value: selectedSubjectId!,
    items: subjectIds,
    onValueChanged: (value) => {
      selectedSubjectId = value;
      updateMeasurements();
    },
  });
  subjectIdInput.root.style.marginLeft = '12px';

  div.appendChild(subjectIdInput.root);

  const measurementsContainer = ui.divV([], {style: {width: '100%', height: '100%'}});

  const acc = ui.accordion();

  acc.addPane('Measurements', () => measurementsContainer, true);

  div.appendChild(acc.root);

  function updateMeasurements() {
    ui.empty(measurementsContainer);

    const pivotedGrid = collectMeasurements(selectedSubjectId, measurementsDf);

    measurementsContainer.appendChild(pivotedGrid.root);
  }

  // Initial update
  updateMeasurements();

  return div;
}

function filterRequiredColumns(df: DG.DataFrame, cols: string[]) {
  for (const colName of df.columns.names()) {
    if (!cols.includes(colName))
      df.columns.remove(colName);
  }
}

function collectMeasurements(selectedSubjectId?: string, measurementsDf?: DG.DataFrame): DG.Grid {
  if (!selectedSubjectId || !measurementsDf)
    return DG.DataFrame.create().plot.grid();

  // Filter by selected subject and leave only required columns
  const filteredDf: DG.DataFrame =
      measurementsDf.rows.match({[SUBJECT_ID]: selectedSubjectId}).toDataFrame();

  if (filteredDf.rowCount === 0)
    return DG.DataFrame.create().plot.grid();

  const colNames = [DOMAIN, 'category', 'test', 'visit_day', 'string_result', 'units'];
  filterRequiredColumns(filteredDf, colNames);

  const allVisitDays = filteredDf.col('visit_day').categories.map((it) => parseInt(it));
  const correctColumnsOrder = [DOMAIN, 'category', 'test', 'units'].concat(allVisitDays.map((it) => it.toString()));

  const pivotedDf = filteredDf
    .groupBy([DOMAIN, 'category', 'test', 'units'])
    .pivot('visit_day')
    .first('string_result', '')
    .aggregate();

  const grid = pivotedDf.plot.grid();
  grid.columns.setOrder(correctColumnsOrder);

  //set grid styles
  grid.root.style.width = '100%';
  for (const col of allVisitDays)
    grid.col(col.toString()).width = 35;
  grid.col('test').width = 150;
  grid.setOptions({allowColumnMenu: true});

  return grid;
}
