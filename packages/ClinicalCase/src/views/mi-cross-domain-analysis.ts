import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {studies} from '../utils/app-utils';
import {CDISC_STANDARD} from '../utils/types';
import {createVisitDayStrCol} from '../data-preparation/data-preparation';
import {DOMAIN, MIANTREG, MIDIR, MIDY, MIDY_STR, MILAT, MISEV, MISPEC, MISTRESC, PLANNED_TRT_ARM, SUBJECT_ID,
  VISIT_DAY_STR} from '../constants/columns-constants';
import {awaitCheck} from '@datagrok-libraries/test/src/test';

export function createMICrossDomainView(studyId: string): any {
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

  const miDf = studies[studyId].domains.mi;
  //in case there is no MI day categorical column - create it
  if (!miDf.col(MIDY_STR)) {
    miDf.columns.addNewString(MIDY_STR).init((i) => {
      const val = miDf.get(MIDY, i);
      return val ? val.toString() : '';
    });
  };

  const createFindingsSummaryTable = () => {
    //get total subjectsCountByTreatmentArm
    const subjectByArmCount = studies[studyId].domains.mi.groupBy([PLANNED_TRT_ARM])
      .count('total_subjects').aggregate();
    const listOfArms = subjectByArmCount.col(PLANNED_TRT_ARM).toList();
    const listOfTotalSubjNums = subjectByArmCount.col('total_subjects').toList();
    //join with microscopic findings domain
    if (!miDf.col('total_subjects')) {
      grok.data.joinTables(miDf, subjectByArmCount, [PLANNED_TRT_ARM], [PLANNED_TRT_ARM],
        undefined, ['total_subjects'], DG.JOIN_TYPE.LEFT, true);
      miDf.name = 'mi';
    }
    //count number of each finding occurence for each treatment group
    const subjByFinding = miDf.groupBy(['ARM', 'MISPEC', 'MISTRESC', 'total_subjects'])
      .count('subjects_with_finding').aggregate();
    //calculate percentage
    subjByFinding.columns.addNewString('percent')
      .init((i) => {
        const subjWithFinding = subjByFinding.col('subjects_with_finding').get(i);
        const percent = (subjWithFinding / subjByFinding.col('total_subjects').get(i)) * 100;
        return `${subjWithFinding} (${percent.toFixed(2)}%)`;
      });
    //final pivoted by treatment arm
    const pivoted = subjByFinding.groupBy(['MISPEC', 'MISTRESC']).pivot('ARM').first('percent').aggregate();
    pivoted.columns.names().forEach((col) => {
      if (col.endsWith(' first(percent)')) {
        let newColName = col.replace(' first(percent)', '');
        const idxOfArm = listOfArms.findIndex((it) => it === newColName);
        if (idxOfArm !== -1)
          newColName = `${newColName} (n = ${listOfTotalSubjNums[idxOfArm]})`;
        pivoted.col(col)!.name = newColName;
      }
    });
    pivoted.name = 'findings sumary';
    grok.data.linkTables(miDf, pivoted, [MISPEC, MISTRESC], [MISPEC, MISTRESC], [DG.SYNC_TYPE.FILTER_TO_FILTER], true);

    return pivoted;
  };

  const miSplitCategories = {
    [MISPEC]: 'Sample type',
    [MISTRESC]: 'Standardized result',
    [MIANTREG]: 'Anatomical Region',
    [MILAT]: 'Laterality within subject',
    [MIDIR]: 'Directionality within subject',
    [MISEV]: 'Severity',
    [PLANNED_TRT_ARM]: 'Treatment arm',
  };

  const existingMiSplitCats: string[] = [];
  //set friendly name for proper filters naming
  Object.keys(miSplitCategories).forEach((key) => {
    if (miDf.col(key)) {
      existingMiSplitCats.push(key);
      miDf.col(key).setTag('friendlyName', miSplitCategories[key]);
    }
  });
  const catsForSeverityBoxplot = [MISTRESC, MIANTREG, MILAT, MIDIR];
  //adding MIDY_STR to measurements for further boxplot creation
  grok.data.joinTables(resDf, miDf, [SUBJECT_ID], [SUBJECT_ID], null,
    [MIDY_STR], DG.JOIN_TYPE.LEFT, true);

  let latestFilterFromMi = DG.BitSet.create(resDf.rowCount, () => true);
  const tests = resDf.col('test')!.categories;
  let selectedTest = tests[0];
  const testChoiceInput = ui.input.choice('Test', {
    items: tests,
    value: selectedTest,
    onValueChanged: () => {
      selectedTest = testChoiceInput.value;
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

  const updateBoxPlots = async () => {
    ui.empty(boxplotBySeverityDiv);
    ui.empty(boxplotByDayDiv);
    const filteredMiDf = miDf.clone(miDf.filter);
    resDf.rows.match({test: selectedTest}).filter();
    //combine filters from dropdown and MI domain
    const filteredMeasurements = resDf.clone(resDf.filter.and(latestFilterFromMi));
    const miDayBoxplot = await createBoxPlot(filteredMeasurements, MIDY_STR,
      `${selectedTest} before Microscopic finding day`);
    boxplotByDayDiv.append(miDayBoxplot.root);
    //go through all categories equired for creation of severity boxplot to make sure unique value for each subject
    const catsToSelect: string[] = [];
    for (const cat of catsForSeverityBoxplot) {
      if (filteredMiDf.col(cat)!.categories.length > 1)
        catsToSelect.push(cat);
    }
    if (catsToSelect.length) {
      const catsFriendlyNames = catsToSelect.map((cat) => miSplitCategories[cat]);
      boxplotBySeverityDiv.append(
        ui.divText(`Select values for ${catsFriendlyNames.join(', ')} to create severity boxplot`,
          {style: {textAlign: 'center'}}));
      return;
    }
    grok.data.joinTables(filteredMeasurements, filteredMiDf, [SUBJECT_ID], [SUBJECT_ID], null,
      [MISEV], DG.JOIN_TYPE.LEFT, true);
    const severityBoxPlot = await createBoxPlot(filteredMeasurements, MISEV,
      `${selectedTest} vs Microscopic finding severity`);
    boxplotBySeverityDiv.append(severityBoxPlot.root);
  };

  const boxplotByDayDiv = ui.div('', {style: {width: '100%', height: '100%'}});
  const boxplotBySeverityDiv = ui.div('', {style: {width: '100%', height: '100%', alignContent: 'center'}});
  const boxplotsDiv = ui.divV([
    boxplotByDayDiv,
    boxplotBySeverityDiv,
  ], 'cross-domain-viewers-div');

  grok.data.linkTables(miDf, resDf, [SUBJECT_ID], [SUBJECT_ID], [DG.SYNC_TYPE.FILTER_TO_FILTER], true);

  const summaryFindingsGrid = createFindingsSummaryTable().plot.grid({title: 'Findings frequency table'});

  const onTableViewAdded = async (tableView: DG.TableView) => {
    const ribbons = tableView.getRibbonPanels();
    ribbons.push([testChoiceInput.root]);
    tableView.setRibbonPanels(ribbons);

    await awaitCheck(() => tableView.grid !== null, '', 1000);
    tableView.dockManager.dock(boxplotsDiv, DG.DOCK_TYPE.FILL);
    tableView.dockManager.dock(summaryFindingsGrid, DG.DOCK_TYPE.DOWN, undefined, 'Findings frequency table', 0.2);

    const fg = tableView.getFiltersGroup({createDefaultFilters: false});
    for (const col of existingMiSplitCats)
      fg.add({type: DG.FILTER_TYPE.CATEGORICAL, column: col});

    tableView.dataFrame.onFilterChanged.subscribe(() => {
      latestFilterFromMi = resDf.filter.clone();
      updateBoxPlots();
    });

    updateBoxPlots();
  };

  return {df: miDf, onTableViewAddedFunc: onTableViewAdded};
}
