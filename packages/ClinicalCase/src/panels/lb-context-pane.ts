import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {LAB_TEST, LAB_TEST_CAT, LAB_DAY, PLANNED_TRT_ARM, ACT_TRT_ARM} from '../constants/columns-constants';
import {calculateLBBaselineColumns} from '../data-preparation/data-preparation';

export function getLbContextPane(lbDomain: DG.DataFrame): HTMLElement {
  // Create accordion
  const acc = ui.accordion();

  // Add Heatmap pane
  acc.addPane('LB linechart by visit day', () => createLbLineChart(lbDomain), true);

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
  const container = ui.divV([], {style: {width: '100%', height: '100%'}});

  // Check if required columns exist
  const labTestCol = lbDomain.col(LAB_TEST);
  const labCatCol = lbDomain.col(LAB_TEST_CAT);
  const labDayCol = lbDomain.col(LAB_DAY);
  const armCol = lbDomain.col(PLANNED_TRT_ARM) || lbDomain.col(ACT_TRT_ARM);
  const pctChgCol = lbDomain.col('LB_PCT_CHG');

  if (!labTestCol || !labDayCol || !armCol) {
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

  // Get unique categories and tests
  const categories = labCatCol ? labCatCol.categories : [];
  const allTests = labTestCol.categories;
  let selectedCategory: string | null = null;
  let selectedTest: string | null = null;

  // Filter tests by category if category column exists
  const getTestsForCategory = (category: string | null): string[] => {
    if (!labCatCol || !category)
      return allTests;

    const cats = [];
    for (const test of allTests) {
      const matchDf = lbDomain.rows.match({[LAB_TEST_CAT]: category, [LAB_TEST]: test}).toDataFrame();
      if (matchDf.rowCount)
        cats.push(test);
    }
    return cats;
  };

  // Find test with maximum absolute mean % change from baseline
  const findTestWithMaxMeanPctChg = (): {category: string | null, test: string | null} => {
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
  const maxTestResult = findTestWithMaxMeanPctChg();
  if (maxTestResult.test) {
    selectedTest = maxTestResult.test;
    selectedCategory = maxTestResult.category;
  } else {
    // Fallback to first available if no valid data found
    if (labCatCol && categories.length > 0) {
      selectedCategory = categories[0];
      const testsForCategory = getTestsForCategory(selectedCategory);
      selectedTest = testsForCategory.length > 0 ? testsForCategory[0] : null;
    } else
      selectedTest = allTests.length > 0 ? allTests[0] : null;
  }

  // Create dropdowns container
  const dropdownsContainer = ui.divH([], {style: {gap: '10px', padding: '10px', alignItems: 'center'}});

  // Create line chart (will be initialized later)
  let lineChart: DG.Viewer | null = null;

  // Function to build filter string
  const buildFilterString = (): string => {
    let filterString = '';
    if (labCatCol && selectedCategory)
      filterString = `\${${LAB_TEST_CAT}} == '${selectedCategory}'`;
    if (selectedTest) {
      if (filterString)
        filterString += ` && \${${LAB_TEST}} == '${selectedTest}'`;
      else
        filterString = `\${${LAB_TEST}} == '${selectedTest}'`;
    }
    return filterString;
  };

  // Function to update line chart filter
  const updateLineChartFilter = () => {
    if (!lineChart)
      return;

    const filterString = buildFilterString();
    if (filterString) {
      lineChart.setOptions({
        filter: buildFilterString() || null,
      });
    }
  };

  // Create test dropdown
  const testDropdown = ui.input.choice('Lab Test', {
    value: selectedTest!,
    items: labCatCol && selectedCategory ? getTestsForCategory(selectedCategory) : allTests,
    onValueChanged: () => {
      selectedTest = testDropdown.value;
      updateLineChartFilter();
    },
  });

  // Function to update test dropdown items when category changes
  const updateTestDropdownItems = () => {
    const testsForCategory = labCatCol && selectedCategory ? getTestsForCategory(selectedCategory) : allTests;
    (testDropdown as any).items = testsForCategory;
    if (testsForCategory.length > 0 && (!selectedTest || !testsForCategory.includes(selectedTest!))) {
      selectedTest = testsForCategory[0];
      testDropdown.value = selectedTest;
    } else if (testsForCategory.length === 0)
      selectedTest = null;
  };

  // Create category dropdown if LAB_TEST_CAT exists
  let categoryDropdown: DG.InputBase | null = null;
  if (labCatCol && categories.length > 0) {
    categoryDropdown = ui.input.choice('Lab Category', {
      value: selectedCategory!,
      items: categories,
      onValueChanged: () => {
        selectedCategory = categoryDropdown!.value;
        updateTestDropdownItems();
        updateLineChartFilter();
      },
    });
    dropdownsContainer.appendChild(categoryDropdown.root);
  }

  dropdownsContainer.appendChild(testDropdown.root);

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
