import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {LAB_TEST, LAB_TEST_CAT, LAB_DAY, LAB_RES_N, PLANNED_TRT_ARM, ACT_TRT_ARM} from '../constants/columns-constants';
import {calculateLBBaselineColumns} from '../data-preparation/data-preparation';

// Common types and interfaces
interface CategoryTestData {
  categories: string[];
  allTests: string[];
  labTestCol: DG.Column;
  labCatCol: DG.Column | null;
}

interface SelectedValues {
  category: string | null;
  test: string | null;
}

// Common functions for both linechart and boxplots
function getCategoryTestData(lbDomain: DG.DataFrame, cache: Map<string, CategoryTestData>): CategoryTestData | null {
  // Create a unique key for this dataframe (using name and row count for uniqueness)
  const cacheKey = `${lbDomain.name}_${lbDomain.rowCount}`;

  // Check if we have cached data
  if (cache.has(cacheKey))
    return cache.get(cacheKey)!;

  const labTestCol = lbDomain.col(LAB_TEST);
  const labCatCol = lbDomain.col(LAB_TEST_CAT);

  if (!labTestCol)
    return null;

  const categoryTestData: CategoryTestData = {
    categories: labCatCol ? labCatCol.categories : [],
    allTests: labTestCol.categories,
    labTestCol,
    labCatCol,
  };

  // Cache the result
  cache.set(cacheKey, categoryTestData);

  return categoryTestData;
}

function getTestsForCategory(
  lbDomain: DG.DataFrame,
  category: string | null,
  allTests: string[],
  labCatCol: DG.Column | null): string[] {
  if (!labCatCol || !category)
    return allTests;

  const tests = [];
  for (const test of allTests) {
    const matchDf = lbDomain.rows.match({[LAB_TEST_CAT]: category, [LAB_TEST]: test}).toDataFrame();
    if (matchDf.rowCount)
      tests.push(test);
  }
  return tests;
}

function buildFilterString(category: string | null, test: string | null, labCatCol: DG.Column | null): string {
  let filterString = '';
  if (labCatCol && category)
    filterString = `\${${LAB_TEST_CAT}} == '${category}'`;
  if (test) {
    if (filterString)
      filterString += ` && \${${LAB_TEST}} == '${test}'`;
    else
      filterString = `\${${LAB_TEST}} == '${test}'`;
  }
  return filterString;
}

interface DropdownControls {
  container: HTMLElement;
  categoryDropdown: DG.InputBase | null;
  testDropdown: DG.InputBase;
  selectedCategory: {value: string | null};
  selectedTest: {value: string | null};
}

function createDropdowns(
  lbDomain: DG.DataFrame,
  categoryTestData: CategoryTestData,
  initialValues: SelectedValues,
  onValueChanged: () => void): DropdownControls {
  const {categories, allTests, labCatCol} = categoryTestData;
  const selectedCategory = {value: initialValues.category};
  const selectedTest = {value: initialValues.test};

  // Create dropdowns container
  const dropdownsContainer = ui.divH([], {style: {gap: '10px', padding: '10px', alignItems: 'center'}});

  // Function to update test dropdown items when category changes
  const updateTestDropdownItems = () => {
    const testsForCategory = labCatCol && selectedCategory.value ?
      getTestsForCategory(lbDomain, selectedCategory.value, allTests, labCatCol) : allTests;
    (testDropdown as any).items = testsForCategory;
    if (testsForCategory.length > 0 && (!selectedTest.value || !testsForCategory.includes(selectedTest.value))) {
      selectedTest.value = testsForCategory[0];
      testDropdown.value = selectedTest.value;
    } else if (testsForCategory.length === 0)
      selectedTest.value = null;
    onValueChanged();
  };

  // Create test dropdown
  const testDropdown = ui.input.choice('Lab Test', {
    value: selectedTest.value!,
    items: labCatCol && selectedCategory.value ?
      getTestsForCategory(lbDomain, selectedCategory.value, allTests, labCatCol) : allTests,
    onValueChanged: () => {
      selectedTest.value = testDropdown.value;
      onValueChanged();
    },
  });

  // Create category dropdown if LAB_TEST_CAT exists
  let categoryDropdown: DG.InputBase | null = null;
  if (labCatCol && categories.length > 0) {
    categoryDropdown = ui.input.choice('Lab Category', {
      value: selectedCategory.value!,
      items: categories,
      onValueChanged: () => {
        selectedCategory.value = categoryDropdown!.value;
        updateTestDropdownItems();
      },
    });
    dropdownsContainer.appendChild(categoryDropdown.root);
  }

  dropdownsContainer.appendChild(testDropdown.root);

  return {
    container: dropdownsContainer,
    categoryDropdown,
    testDropdown,
    selectedCategory,
    selectedTest,
  };
}

function initializeDefaultValues(
  lbDomain: DG.DataFrame,
  categoryTestData: CategoryTestData,
  findMaxTest?: () => SelectedValues): SelectedValues {
  const {categories, allTests, labCatCol} = categoryTestData;
  let selectedCategory: string | null = null;
  let selectedTest: string | null = null;

  if (findMaxTest) {
    // For linechart - find test with maximum absolute mean % change
    const maxTestResult = findMaxTest();
    if (maxTestResult.test) {
      selectedTest = maxTestResult.test;
      selectedCategory = maxTestResult.category;
    }
  }

  // Fallback to first available if no valid data found
  if (!selectedTest) {
    if (labCatCol && categories.length > 0) {
      selectedCategory = categories[0];
      const testsForCategory = getTestsForCategory(lbDomain, selectedCategory, allTests, labCatCol);
      selectedTest = testsForCategory.length > 0 ? testsForCategory[0] : null;
    } else
      selectedTest = allTests.length > 0 ? allTests[0] : null;
  }

  return {category: selectedCategory, test: selectedTest};
}

export function getLbContextPane(lbDomain: DG.DataFrame): HTMLElement {
  // Create accordion
  const acc = ui.accordion();

  acc.addPane('LB linechart by visit day', () => createLbLineChart(lbDomain), true);
  acc.addPane('LB boxplots by treatment arm', () => createLbBoxPlots(lbDomain), true);

  return acc.root;
}

// function createHeatmapPane(lbDomain: DG.DataFrame): HTMLElement {
//   const container = ui.div([], {style: {width: '100%', height: '100%'}});

//   // Check if required columns exist
//   const labTestCol = lbDomain.col(LAB_TEST);
//   const labCatCol = lbDomain.col(LAB_TEST_CAT);
//   const armCol = lbDomain.col(PLANNED_TRT_ARM) || lbDomain.col(ACT_TRT_ARM);
//   const pctChgCol = lbDomain.col('LB_PCT_CHG');

//   if (!labTestCol || !armCol) {
//     container.appendChild(ui.divText('Required columns not found in LB domain'));
//     return container;
//   }

//   // Ensure LB_PCT_CHG column exists
//   if (!pctChgCol)
//     calculateLBBaselineColumns(lbDomain);

//   const finalPctChgCol = lbDomain.col('LB_PCT_CHG');
//   if (!finalPctChgCol) {
//     container.appendChild(ui.divText('Could not calculate LB_PCT_CHG column'));
//     return container;
//   }

//   // Create pivot dataframe
//   // Group by: lab test + lab category (if exists)
//   const groupByCols = labCatCol ? [LAB_TEST_CAT, LAB_TEST] : [LAB_TEST];

//   // Create pivoted dataframe: rows = lab test + category, columns = treatment arms, values = mean % change
//   let pivotedDf: DG.DataFrame;

//   try {
//     pivotedDf = lbDomain
//       .groupBy(groupByCols)
//       .pivot(armCol.name)
//       .avg('LB_PCT_CHG')
//       .aggregate();

//     // Clean up column names (remove " avg(LB_PCT_CHG)" suffix)
//     const cols = pivotedDf.columns.names();
//     for (const colName of cols) {
//       if (colName.endsWith(' avg(LB_PCT_CHG)'))
//         pivotedDf.col(colName)!.name = colName.replace(' avg(LB_PCT_CHG)', '');
//     }

//     // Create heatmap viewer
//     const heatmapViewer = DG.Viewer.fromType(DG.VIEWER.HEAT_MAP, pivotedDf);
//     heatmapViewer.root.style.width = '100%';
//     container.appendChild(heatmapViewer.root);
//   } catch (error: any) {
//     console.error('Error creating pivot dataframe:', error);
//     container.appendChild(ui.divText(`Error creating pivot dataframe: ${error.message}`));
//   }
//   return container;
// }

function createLbLineChart(lbDomain: DG.DataFrame): HTMLElement {
  const categoryTestDataCache = new Map<string, CategoryTestData>();

  const container = ui.divV([], {style: {width: '100%', height: '100%'}});

  // Check if required columns exist
  const labDayCol = lbDomain.col(LAB_DAY);
  const armCol = lbDomain.col(PLANNED_TRT_ARM) || lbDomain.col(ACT_TRT_ARM);
  const pctChgCol = lbDomain.col('LB_PCT_CHG');

  if (!labDayCol || !armCol) {
    container.appendChild(ui.divText('Required columns not found in LB domain'));
    return container;
  }

  // Ensure LB_PCT_CHG column exists
  if (!pctChgCol)
    calculateLBBaselineColumns(lbDomain);

  const finalPctChgCol = lbDomain.col('LB_PCT_CHG');
  if (!finalPctChgCol) {
    container.appendChild(ui.divText('Could not calculate LB_PCT_CHG column'));
    return container;
  }

  // Get category and test data
  const categoryTestData = getCategoryTestData(lbDomain, categoryTestDataCache);
  if (!categoryTestData) {
    container.appendChild(ui.divText('Required columns not found in LB domain'));
    return container;
  }

  const {labCatCol} = categoryTestData;

  // Find test with maximum absolute mean % change from baseline
  const findTestWithMaxMeanPctChg = (): SelectedValues => {
    const groupByCols = labCatCol ? [LAB_TEST_CAT, LAB_TEST] : [LAB_TEST];

    // Create a temporary dataframe with absolute values of LB_PCT_CHG
    const tempDf = lbDomain.clone();
    const absPctChgCol = tempDf.columns.addNewFloat('ABS_LB_PCT_CHG');
    const pctChgCol = tempDf.col('LB_PCT_CHG');
    absPctChgCol.init((i) => {
      if (pctChgCol.isNone(i))
        return null;
      const value = pctChgCol.get(i);
      return value !== null && !isNaN(value) ? Math.abs(value) : null;
    });

    // Calculate mean of absolute values
    const meanDf = tempDf
      .groupBy(groupByCols)
      .avg('ABS_LB_PCT_CHG')
      .aggregate();

    const avgColName = 'avg(ABS_LB_PCT_CHG)';
    const avgCol = meanDf.col(avgColName);
    if (!avgCol || meanDf.rowCount === 0)
      return {category: null, test: null};

    let maxAbsMean = -Infinity;
    let maxCategory: string | null = null;
    let maxTest: string | null = null;

    for (let i = 0; i < meanDf.rowCount; i++) {
      if (avgCol.isNone(i))
        continue;

      const meanValue = avgCol.get(i);
      if (meanValue !== null && !isNaN(meanValue) && meanValue > maxAbsMean) {
        maxAbsMean = meanValue;
        maxTest = meanDf.get(LAB_TEST, i);
        if (labCatCol)
          maxCategory = meanDf.get(LAB_TEST_CAT, i);
      }
    }

    return {category: maxCategory, test: maxTest};
  };

  // Initialize selected values - find test with maximum mean % change
  const initialValues = initializeDefaultValues(lbDomain, categoryTestData, findTestWithMaxMeanPctChg);

  // Create line chart (will be initialized later)
  let lineChart: DG.Viewer | null = null;

  // Function to update line chart filter
  const updateLineChartFilter = () => {
    if (!lineChart)
      return;

    const filterString = buildFilterString(selectedCategory.value, selectedTest.value, labCatCol);
    if (filterString) {
      lineChart.setOptions({
        filter: filterString,
      });
    }
  };

  // Create dropdowns using common function
  const dropdownControls = createDropdowns(
    lbDomain,
    categoryTestData,
    initialValues,
    updateLineChartFilter);
  const {container: dropdownsContainer, selectedCategory, selectedTest} = dropdownControls;

  container.appendChild(dropdownsContainer);

  // Create line chart viewer
  try {
    lineChart = DG.Viewer.fromType(DG.VIEWER.LINE_CHART, lbDomain, {
      xColumnName: LAB_DAY,
      yColumnNames: ['LB_PCT_CHG'],
      splitColumnName: armCol.name,
      rowSource: 'All',
    });
    updateLineChartFilter();
    lineChart.root.style.width = '100%';
    container.appendChild(lineChart.root);
  } catch (error: any) {
    console.error('Error creating line chart:', error);
    container.appendChild(ui.divText(`Error creating line chart: ${error.message}`));
  }

  return container;
}

function createLbBoxPlots(lbDomain: DG.DataFrame): HTMLElement {
  const categoryTestDataCache = new Map<string, CategoryTestData>();

  const container = ui.divV([], {style: {width: '100%', height: '100%'}});

  // Check if required columns exist
  const armCol = lbDomain.col(PLANNED_TRT_ARM) || lbDomain.col(ACT_TRT_ARM);
  const labResNCol = lbDomain.col(LAB_RES_N);

  if (!armCol || !labResNCol) {
    container.appendChild(ui.divText('Required columns not found in LB domain'));
    return container;
  }

  // Get category and test data
  const categoryTestData = getCategoryTestData(lbDomain, categoryTestDataCache);
  if (!categoryTestData) {
    container.appendChild(ui.divText('Required columns not found in LB domain'));
    return container;
  }

  // Initialize selected values - use first available
  const initialValues = initializeDefaultValues(lbDomain, categoryTestData);

  // Container for boxplot viewer
  const boxplotContainer = ui.div([], {style: {width: '100%', height: '400px'}});
  let boxPlot: DG.Viewer | null = null;

  // Function to create filtered dataframe and boxplot
  const updateBoxPlot = () => {
    if (!selectedTest.value)
      return;

    const {labCatCol} = categoryTestData;

    // Build filter object
    const filterObj: any = {};
    if (labCatCol && selectedCategory.value)
      filterObj[LAB_TEST_CAT] = selectedCategory.value;
    if (selectedTest.value)
      filterObj[LAB_TEST] = selectedTest.value;

    // Create filtered dataframe
    const filteredDf = lbDomain.rows.match(filterObj).toDataFrame();

    // Remove old boxplot if exists
    if (boxPlot) {
      boxPlot.root.remove();
      boxPlot = null;
    }

    // Create new boxplot
    try {
      boxPlot = DG.Viewer.fromType(DG.VIEWER.BOX_PLOT, filteredDf, {
        category1: armCol.name,
        value: LAB_RES_N,
      });

      boxPlot.root.style.width = '100%';
      boxplotContainer.appendChild(boxPlot.root);
    } catch (error: any) {
      console.error('Error creating boxplot:', error);
      boxplotContainer.appendChild(ui.divText(`Error creating boxplot: ${error.message}`));
    }
  };

  // Create dropdowns using common function
  const dropdownControls = createDropdowns(
    lbDomain,
    categoryTestData,
    initialValues,
    updateBoxPlot);
  const {container: dropdownsContainer, selectedCategory, selectedTest} = dropdownControls;

  container.appendChild(dropdownsContainer);
  container.appendChild(boxplotContainer);

  // Create initial boxplot
  updateBoxPlot();

  return container;
}
