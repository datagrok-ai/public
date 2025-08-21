
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { ComplexCondition, Operators, QueryBuilder, QueryBuilderLayout, SUGGESTIONS_FUNCTION } from '@datagrok-libraries/utils/src/query-builder/query-builder';
import { convertComplexConditionToSignalsSearchQuery, SignalsSearchQuery } from './signals-search-query';
import { materialsCondition } from './compounds';
import { getRevvityLibraries } from './libraries';
import { getDefaultProperties, REVVITY_FIELD_TO_PROP_TYPE_MAPPING } from './properties';
import { getTerms } from './package';
import { createPath, initSearchQuery, resetInitSearchQuery } from './view-utils';

export const viewQueryBuilders: { [key: string]: QueryBuilder } = {};

export async function runSearchQuery(libId: string, compoundType: string,
  queryBuilderCondition: ComplexCondition): Promise<DG.DataFrame> {
  const condition: ComplexCondition = {
    logicalOperator: Operators.Logical.and,
    conditions: [
      materialsCondition //now we have only materials loaded, so adding condition to search through materials
    ]
  }
  condition.conditions.push(
    {
      field: "assetTypeEid",
      operator: Operators.EQ,
      value: libId
    }
  );
  condition.conditions.push(
    {
      field: "type",
      operator: Operators.EQ,
      value: compoundType
    }
  );
  condition.conditions.push(queryBuilderCondition);
  const signalsQuery: SignalsSearchQuery = convertComplexConditionToSignalsSearchQuery(condition);
  console.log(signalsQuery);
  const resultDf = await grok.functions.call('RevvitySignalsLink:searchEntitiesWithStructures', {
    query: JSON.stringify(signalsQuery),
    params: '{}'
  });
  return resultDf;
}


export async function initializeFilters(tv: DG.TableView, filtersDiv: HTMLDivElement, libName: string, compoundType: string) {
  const libs = await getRevvityLibraries();
  const selectedLib = libs.filter((l) => l.name === libName);
  if (selectedLib.length) {
    //create filters icon
    const svgElement = document.createElement('div');
    svgElement.innerHTML = `
<svg width="24" height="24" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
  <!-- Funnel Body -->
  <path d="M4 4H15L10 10V17L8 19V10L4 4Z" stroke="currentColor" stroke-width="1" stroke-linecap="round" stroke-linejoin="round"/>
  
  <!-- The arrow group is rotated -45 degrees around point (16, 13) -->
  <g transform="rotate(-45, 16, 13)">
    <!-- Long Horizontal Arrow Line (starts at the new x=16, y=13) -->
    <path d="M16 13H24" stroke="currentColor" stroke-width="1" stroke-linecap="round"/>
    <!-- Arrowhead (adjusted for the new starting point) -->
    <path d="M21 10L24 13L21 16" stroke="currentColor" stroke-width="1" stroke-linecap="round" stroke-linejoin="round"/>
  </g>
</svg>`;
    svgElement.className = 'filters-button-icon';
    svgElement.onclick = () => {
      tv.dockManager.dock(filtersDiv, 'left', null, 'Filters', 0.2);
    };
    ui.tooltip.bind(svgElement, 'Add Revvity Signals filters');
    const filtersButton = ui.div(svgElement);
    tv.setRibbonPanels([[filtersButton]]);
    //create filters panel
    tv.dockManager.dock(filtersDiv, 'left', null, 'Filters', 0.2);
    ui.setUpdateIndicator(filtersDiv, true, 'Loading filters...');
    await initializeQueryBuilder(selectedLib[0].id, compoundType);

    const runSearchButton = ui.bigButton('Search', async () => {
      runSearch(tv, selectedLib[0].id, compoundType, libName);
    });

    ui.setUpdateIndicator(filtersDiv, false);
    filtersDiv.append(viewQueryBuilders[`${selectedLib[0].id}|${compoundType}`].root);
    filtersDiv.append(ui.div(runSearchButton, { style: { paddingLeft: '4px' } }));

    ui.onSizeChanged(filtersDiv).subscribe(() => {
      updateQueryBuilderLayout(filtersDiv.clientWidth, selectedLib[0].id, compoundType);
    });

    if (initSearchQuery)
      runSearch(tv, selectedLib[0].id, compoundType, libName);
  }
}

export async function runSearch(tv: DG.TableView, libId: string, compoundType: string, libName: string) {
  viewQueryBuilders[`${libId}|${compoundType}`]!.saveConditionToHistory();
  ui.setUpdateIndicator(tv.grid.root, true, 'Searching...');
  const condition = viewQueryBuilders[`${libId}|${compoundType}`]!.condition;
  const resultDf = await runSearchQuery(libId, compoundType, condition);
  tv.dataFrame = resultDf;
  tv.path = createPath([libName, compoundType, 'search', JSON.stringify(condition)]);
  ui.setUpdateIndicator(tv.grid.root, false);
  if (initSearchQuery)
    resetInitSearchQuery();
}


export async function initializeQueryBuilder(libId: string, compoundType: string) {
  const filterFields = getDefaultProperties();
  const tagsStr = await grok.functions.call('RevvitySignalsLink:getTags', {
    assetTypeId: libId,
    type: compoundType
  });
  const tags: { [key: string]: string } = JSON.parse(tagsStr);
  Object.keys(tags).forEach((tagName) => {
    const propOptions: { [key: string]: any } = {
      name: tagName,
      type: REVVITY_FIELD_TO_PROP_TYPE_MAPPING[tags[tagName]],
    };
    const nameArr = tagName.split('.');
    if (nameArr.length > 1)
      propOptions.friendlyName = nameArr[1];
    const prop = DG.Property.fromOptions(propOptions);
    prop.options[SUGGESTIONS_FUNCTION] = async (text: string) => {
      const termsStr = await getTerms(tagName, compoundType, libId, true);
      const terms: string[] = JSON.parse(termsStr);
      return terms.filter((it) => it.toLowerCase().includes(text.toLowerCase()));
    }
    filterFields.push(prop);
  });

  if (!viewQueryBuilders[`${libId}|${compoundType}`])
    viewQueryBuilders[`${libId}|${compoundType}`] = new QueryBuilder(filterFields, initSearchQuery, QueryBuilderLayout.Narrow, `Revvity Signals|${libId}|${compoundType}`);
}


function updateQueryBuilderLayout(width: number, libId: string, compoundType: string) {
  if (!viewQueryBuilders[`${libId}|${compoundType}`]) return;

  // Switch to narrow layout if width is less than 300px, otherwise use standard
  const newLayout = width < 300 ? QueryBuilderLayout.Narrow : QueryBuilderLayout.Standard;

  if (viewQueryBuilders[`${libId}|${compoundType}`].getLayout() !== newLayout) {
    viewQueryBuilders[`${libId}|${compoundType}`].setLayout(newLayout);
  }
};