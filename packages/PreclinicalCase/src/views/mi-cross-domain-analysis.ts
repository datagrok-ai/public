import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {studies} from '../utils/app-utils';
import {createVisitDayStrCol} from '../data-preparation/data-preparation';
import {DOMAIN, MIANTREG, MIDIR, MIDY, MIDY_STR, MILAT, MISEV, MISPEC, MISTRESC, PLANNED_TRT_ARM, SUBJECT_ID,
  VISIT_DAY_STR} from '../constants/columns-constants';
import {awaitCheck} from '@datagrok-libraries/test/src/test';
import { StudyTableViewParams } from '../utils/views-creation-utils';

export function createMICrossDomainView(studyId: string): StudyTableViewParams {
  let resDf: DG.DataFrame | null = null;

  studies[studyId].domains.all().forEach((it) => {
    let requiredColumns = [SUBJECT_ID, DOMAIN, `${it.name.toUpperCase()}TEST`,
      `${it.name.toUpperCase()}STRESN`, `${it.name.toUpperCase()}STRESC`,
      `${it.name.toUpperCase()}STRESU`, `${it.name.toUpperCase()}DY`];
    if (requiredColumns.every((colName) => it.columns.names().includes(colName))) {
      createVisitDayStrCol(it);
      requiredColumns = requiredColumns.concat([VISIT_DAY_STR]);

      const catColName = `${it.name.toUpperCase()}CAT`;
      const containsCatCol = it.col(catColName);
      const fullTestNameCol = 'full_test_name';
      if (containsCatCol) {
        it.columns.addNewString(fullTestNameCol)
          .init((i) => `${it.get(`${it.name.toUpperCase()}TEST`, i)} ${
            it.get(catColName, i) ? `(${it.get(catColName, i)})` : ``}`);
        requiredColumns.push(fullTestNameCol);
        requiredColumns = requiredColumns.filter((i) => i !== `${it.name.toUpperCase()}TEST`);
      } else
        it.columns.addNewString(catColName);
      requiredColumns.push(`${it.name.toUpperCase()}CAT`);

      const df = it.clone(null, requiredColumns);
      df.getCol(containsCatCol ? fullTestNameCol : `${it.name.toUpperCase()}TEST`).name = 'test';
      df.getCol(`${it.name.toUpperCase()}CAT`).name = 'category';
      df.getCol(`${it.name.toUpperCase()}STRESN`).name = 'result';
      df.getCol(`${it.name.toUpperCase()}DY`).name = 'visit_day';
      df.getCol(`${it.name.toUpperCase()}STRESC`).name = 'string_result';
      df.getCol(`${it.name.toUpperCase()}STRESU`).name = 'units';

      if (containsCatCol)
        it.columns.remove(fullTestNameCol);

      if (!resDf)
        resDf = df;
      else
        resDf.append(df, true);
    }
  });

  const columnsFromDm = studies[studyId].domains.dm!.columns.names().filter((it) => it !== SUBJECT_ID);
  grok.data.joinTables(resDf!, studies[studyId].domains.dm!, [SUBJECT_ID], [SUBJECT_ID], null,
    columnsFromDm, DG.JOIN_TYPE.LEFT, true);

  const intVisitDay = resDf!.col('visit_day')!.convertTo( DG.TYPE.INT);
  resDf!.columns.remove('visit_day');
  resDf!.columns.add(intVisitDay);

  const miDf = studies[studyId].domains.mi;
  if (!miDf)
   return {df: DG.DataFrame.create()};
  if (!miDf.col(MIDY_STR)) {
    miDf.columns.addNewString(MIDY_STR).init((i) => {
      const val = miDf.get(MIDY, i);
      return val ? val.toString() : '';
    });
  };

  const createFindingsSummaryTable = () => {
    const subjectByArmCount = studies[studyId].domains.mi!.groupBy([PLANNED_TRT_ARM])
      .uniqueCount(SUBJECT_ID, 'total_subjects').aggregate();
    const listOfArms: string[] = subjectByArmCount.col(PLANNED_TRT_ARM)!.toList();
    const listOfTotalSubjNums: number[] = subjectByArmCount.col('total_subjects')!.toList();

    const subjByFinding = miDf.groupBy([PLANNED_TRT_ARM, MISPEC, MISTRESC])
      .uniqueCount(SUBJECT_ID, 'subjects_with_finding').aggregate();

    const pivotedCount = subjByFinding.groupBy([MISPEC, MISTRESC])
      .pivot(PLANNED_TRT_ARM).first('subjects_with_finding').aggregate();

    const result = DG.DataFrame.create(pivotedCount.rowCount);
    result.columns.addNewString(MISPEC).init((i) => pivotedCount.get(MISPEC, i));
    result.columns.addNewString(MISTRESC).init((i) => pivotedCount.get(MISTRESC, i));

    for (let a = 0; a < listOfArms.length; a++) {
      const arm = listOfArms[a];
      const total = listOfTotalSubjNums[a];
      const srcCountCol = pivotedCount.columns.names()
        .find((c) => c.startsWith(arm) && c.includes('first(subjects_with_finding)'));
      if (!srcCountCol)
        continue;

      const countColName = `N ${arm} (n=${total})`;
      const pctColName = `% ${arm} (n=${total})`;
      result.columns.addNewInt(countColName).init((i) => {
        const v = pivotedCount.get(srcCountCol, i);
        return v == null || v === '' ? 0 : Number(v);
      });
      const pctCol = result.columns.addNewFloat(pctColName).init((i) => {
        const v = pivotedCount.get(srcCountCol, i);
        const n = v == null || v === '' ? 0 : Number(v);
        return total > 0 ? (n / total) * 100 : 0;
      });
    }

    result.name = 'findings summary';
    grok.data.linkTables(miDf, result, [MISPEC, MISTRESC], [MISPEC, MISTRESC], [DG.SYNC_TYPE.FILTER_TO_FILTER], true);
    return result;
  };

  const miSplitCategories: {[key: string]: string} = {
    [MISPEC]: 'Sample type',
    [MISTRESC]: 'Standardized result',
    [MIANTREG]: 'Anatomical Region',
    [MILAT]: 'Laterality within subject',
    [MIDIR]: 'Directionality within subject',
    [MISEV]: 'Severity',
    [PLANNED_TRT_ARM]: 'Treatment arm',
  };

  const existingMiSplitCats: string[] = [];
  Object.keys(miSplitCategories).forEach((key) => {
    if (miDf.col(key)) {
      existingMiSplitCats.push(key);
      miDf.col(key)!.setTag('friendlyName', miSplitCategories[key]);
    }
  });
  const catsForSeverityBoxplot = [MISTRESC, MIANTREG, MILAT, MIDIR];
  grok.data.joinTables(resDf!, miDf, [SUBJECT_ID], [SUBJECT_ID], null,
    [MIDY_STR], DG.JOIN_TYPE.LEFT, true);

  let latestFilterFromMi = DG.BitSet.create(resDf!.rowCount, () => true);
  const tests = resDf!.col('test')!.categories;
  let selectedTest = tests[0];
  const testChoiceInput = ui.input.choice('Test', {
    items: tests,
    value: selectedTest,
    onValueChanged: () => {
      selectedTest = testChoiceInput.value!;
      updateBoxPlots();
    },
  });

  const createBoxPlot = async (df: DG.DataFrame, cat1: string, title: string) => {
    const boxPlot = await DG.Viewer.fromType(DG.VIEWER.BOX_PLOT, df, {
      category1: cat1,
      category2: PLANNED_TRT_ARM,
      showStatistics: false,
      showCategorySelector: false,
      showColorSelector: false,
      showValueSelector: false,
      title: title,
    });
    boxPlot.root.style.width = '100%';
    boxPlot.root.style.height = '100%';
    return boxPlot;
  };

  const NO_FINDING = 'No Finding';

  const updateBoxPlots = async () => {
    ui.empty(boxplotBySeverityDiv);
    ui.empty(boxplotByDayDiv);
    const filteredMiDf = miDf.clone(miDf.filter);
    resDf!.rows.match({test: selectedTest}).filter();
    const testOnlyFilter = resDf!.filter.clone();
    const filteredMeasurements = resDf!.clone(resDf!.filter.and(latestFilterFromMi));
    const miDayBoxplot = await createBoxPlot(filteredMeasurements, MIDY_STR,
      `${selectedTest} before Microscopic finding day`);
    boxplotByDayDiv.append(miDayBoxplot.root);
    const catsToSelect: string[] = [];
    for (const cat of catsForSeverityBoxplot) {
      if (filteredMiDf.col(cat) && filteredMiDf.col(cat)!.categories.length > 1)
        catsToSelect.push(cat);
    }
    if (catsToSelect.length) {
      const catsFriendlyNames = catsToSelect.map((cat) => miSplitCategories[cat]);
      boxplotBySeverityDiv.append(
        ui.divText(`Select values for ${catsFriendlyNames.join(', ')} to create severity boxplot`,
          {style: {textAlign: 'center'}}));
      return;
    }
    // Include ALL subjects for the selected test — not just those with MI findings.
    // Subjects without the selected finding get "No Finding", providing the reference
    // distribution needed to assess correlation between severity and lab values.
    const allSubjectMeasurements = resDf!.clone(testOnlyFilter);
    const hasMiSev = filteredMiDf.col(MISEV) !== null;
    const joinCol = hasMiSev ? MISEV : MISTRESC;
    grok.data.joinTables(allSubjectMeasurements, filteredMiDf, [SUBJECT_ID], [SUBJECT_ID], null,
      [joinCol], DG.JOIN_TYPE.LEFT, true);
    const joinedCol = allSubjectMeasurements.col(joinCol);
    if (!joinedCol)
      return;
    allSubjectMeasurements.columns.addNewString('severity').init((i) => {
      if (joinedCol.isNone(i) || joinedCol.get(i) === '')
        return NO_FINDING;
      return hasMiSev ? joinedCol.get(i) : 'Finding';
    });
    const severityBoxPlot = await createBoxPlot(allSubjectMeasurements, 'severity',
      `${selectedTest} vs Microscopic finding severity`);
    boxplotBySeverityDiv.append(severityBoxPlot.root);
  };

  const boxplotByDayDiv = ui.div('', {style: {width: '100%', height: '100%'}});
  const boxplotBySeverityDiv = ui.div('', {style: {width: '100%', height: '100%', alignContent: 'center'}});
  const boxplotsDiv = ui.divV([
    boxplotByDayDiv,
    boxplotBySeverityDiv,
  ], 'cross-domain-viewers-div');

  grok.data.linkTables(miDf, resDf!, [SUBJECT_ID], [SUBJECT_ID], [DG.SYNC_TYPE.FILTER_TO_FILTER], true);

  const summaryResult = createFindingsSummaryTable();
  const summaryFindingsGrid = summaryResult.plot.grid({title: 'Findings frequency table'});

  const onTableViewAdded = async (tableView: DG.TableView) => {
    tableView.setRibbonPanels([[testChoiceInput.root]]);

    await awaitCheck(() => tableView.grid !== null, '', 1000);
    tableView.dockManager.dock(boxplotsDiv, DG.DOCK_TYPE.FILL);
    tableView.dockManager.dock(summaryFindingsGrid, DG.DOCK_TYPE.DOWN, undefined, 'Findings frequency table', 0.2);

    const fg = tableView.getFiltersGroup({createDefaultFilters: false});
    for (const col of existingMiSplitCats)
      fg.add({type: DG.FILTER_TYPE.CATEGORICAL, column: col, sortColumn: col === MISPEC ? 'total' : 'filtered', sortDirection: 'false'});

    tableView.dataFrame.onFilterChanged.subscribe(() => {
      latestFilterFromMi = resDf!.filter.clone();
      updateBoxPlots();
    });

    updateBoxPlots();
    //workaround to switch back to browse from toolbox
    (grok.shell.browsePanel.root.closest('.dock-container')?.querySelector('[name = "view-handle: Browse"]') as HTMLElement)?.click();
  };

  return {df: miDf, onTableViewAddedFunc: onTableViewAdded};
}
