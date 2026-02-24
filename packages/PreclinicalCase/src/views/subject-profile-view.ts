import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {studies} from '../utils/app-utils';
import {SUBJECT_ID, DOMAIN} from '../constants/columns-constants';

export function createSubjectProfileView(studyId: string, measurementsDf: DG.DataFrame, subjectId?: string):
HTMLDivElement {
  const div = ui.div('', {style: {overflowY: 'scroll'}});

  const study = studies[studyId];
  const currentStudyId = studyId;

  const subjectIds = study.domains.dm!.col(SUBJECT_ID)!.categories;
  let selectedSubjectId = subjectId ?? (subjectIds.length > 0 ? subjectIds[0] : null);

  const subjectIdInput = ui.input.choice('Subject ID', {
    value: selectedSubjectId!,
    items: subjectIds,
    onValueChanged: (value) => {
      selectedSubjectId = value;
      updateMeasurements();
      updateClinicalObservations();
    },
  });
  subjectIdInput.root.style.marginLeft = '12px';

  div.appendChild(subjectIdInput.root);

  const measurementsContainer = ui.divV([], {style: {width: '100%', height: '100%'}});
  const observationsContainer = ui.divV([], {style: {width: '100%', height: '100%'}});

  const acc = ui.accordion();

  acc.addPane('Measurements', () => measurementsContainer, true);
  acc.addPane('Clinical Observations', () => observationsContainer, true);

  div.appendChild(acc.root);

  function updateMeasurements() {
    ui.empty(measurementsContainer);

    const pivotedGrid = collectMeasurements(selectedSubjectId!, measurementsDf);

    measurementsContainer.appendChild(pivotedGrid.root);
  }

  function updateClinicalObservations() {
    ui.empty(observationsContainer);

    const pivotedGrid = collectClinicalObservations(currentStudyId, selectedSubjectId!);

    observationsContainer.appendChild(pivotedGrid.root);
  }

  updateMeasurements();
  updateClinicalObservations();

  return div;
}

function collectClinicalObservations(studyId: string, selectedSubjectId?: string): DG.Grid {
  if (!selectedSubjectId || !studies[studyId].domains.cl)
    return DG.DataFrame.create().plot.grid();

  const filteredDf: DG.DataFrame =
    studies[studyId].domains.cl.rows.match({[SUBJECT_ID]: selectedSubjectId}).toDataFrame();

  if (filteredDf.rowCount === 0)
    return DG.DataFrame.create().plot.grid();

  const colNames = ['CLCAT', 'CLSCAT', 'CLSTRESC', 'CLLOC', 'CLDY'];
  const existingCols: string[] = [];
  for (const colName of colNames) {
    if (filteredDf.col(colName))
      existingCols.push(colName);
  }

  if (!['CLSTRESC', 'CLDY'].every((it) => existingCols.includes(it)))
    return DG.DataFrame.create().plot.grid();

  filterRequiredColumns(filteredDf, existingCols);

  const intVisitDay = filteredDf.col('CLDY')!.convertTo( DG.TYPE.INT);
  filteredDf.columns.remove('CLDY');
  filteredDf.columns.add(intVisitDay);

  const allVisitDays = filteredDf.col('CLDY')!.categories.map((it) => parseInt(it));
  const colNamesDict: {[key: string]: string} = {
    'CLCAT': 'category',
    'CLSCAT': 'subcategory',
    'CLLOC': 'location',
    'CLSTRESC': 'result',
  };

  const pivotedDf = filteredDf
    .groupBy(existingCols.filter((it) => it !== 'CLDY'))
    .pivot('CLDY')
    .count('')
    .aggregate();

  Object.keys(colNamesDict).forEach((key) => {
    if (pivotedDf.col(key))
      pivotedDf.col(key)!.setTag('friendlyName', colNamesDict[key]);
  });

  const correctColumnsOrder = ['category', 'subcategory', 'location', 'result']
    .filter((it) => pivotedDf.columns.names().includes(it))
    .concat(allVisitDays.map((it) => it.toString()));

  const grid = pivotedDf.plot.grid();
  grid.columns.setOrder(correctColumnsOrder);

  grid.root.style.width = '100%';
  for (const col of allVisitDays)
    grid.col(col.toString())!.width = 35;
  grid.setOptions({allowColumnMenu: true});

  return grid;
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

  const filteredDf: DG.DataFrame =
      measurementsDf.rows.match({[SUBJECT_ID]: selectedSubjectId}).toDataFrame();

  if (filteredDf.rowCount === 0)
    return DG.DataFrame.create().plot.grid();

  const colNames = [DOMAIN, 'category', 'test', 'visit_day', 'string_result', 'units'];
  filterRequiredColumns(filteredDf, colNames);

  const allVisitDays = filteredDf.col('visit_day')!.categories.map((it) => parseInt(it));
  const correctColumnsOrder = [DOMAIN, 'category', 'test', 'units'].concat(allVisitDays.map((it) => it.toString()));

  const pivotedDf = filteredDf
    .groupBy([DOMAIN, 'category', 'test', 'units'])
    .pivot('visit_day')
    .first('string_result', '')
    .aggregate();

  const grid = pivotedDf.plot.grid();
  grid.columns.setOrder(correctColumnsOrder);

  grid.root.style.width = '100%';
  for (const col of allVisitDays)
    grid.col(col.toString())!.width = 35;
  grid.col('test')!.width = 150;
  grid.setOptions({allowColumnMenu: true});

  return grid;
}
