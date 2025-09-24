
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { ComplexCondition, Operators, QueryBuilder, QueryBuilderLayout, SUGGESTIONS_FUNCTION } from '@datagrok-libraries/utils/src/query-builder/query-builder';
import { convertComplexConditionToSignalsSearchQuery, SignalsSearchQuery } from './signals-search-query';
import { materialsCondition } from './compounds';
import { createLibsObjectForStatistics, getRevvityLibraries, RevvityLibrary } from './libraries';
import { getDefaultProperties, REVVITY_FIELD_TO_PROP_TYPE_MAPPING } from './properties';
import { getTerms } from './package';
import { createPath, createViewForExpandabelNode, createViewFromPreDefinedQuery, openRevvityNode, setUserColumnsStyle } from './view-utils';
import { getAppHeader, getCompoundTypeByViewName, getViewNameByCompoundType } from './utils';

export const SAVED_SEARCH_STORAGE = 'RevvitySignalsLinkSavedSearch'

export async function runSearchQuery(libId: string, compoundType: string,
  queryBuilderCondition: ComplexCondition): Promise<DG.DataFrame> {
  const condition: ComplexCondition = {
    logicalOperator: Operators.Logical.and,
    conditions: [
      materialsCondition //now we have only materials section in our revvity instance, so adding condition to search through materials
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


export async function initializeFilters(tv: DG.TableView, filtersDiv: HTMLDivElement, libName: string,
  compoundType: string, initialSearchQuery?: ComplexCondition, changePath?: boolean) {
  const libs = await getRevvityLibraries();
  const selectedLib = libs.filter((l) => l.name === libName);
  if (selectedLib.length) {
    //create filters icon
    const externalFilterIcon = document.createElement('div');
    externalFilterIcon.innerHTML = `
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
    externalFilterIcon.className = 'revvity-filters-button-icon';
    externalFilterIcon.onclick = () => {
      tv.dockManager.dock(filtersDiv, 'left', null, 'Filters', 0.2);
      externalFilterIcon.classList.remove('revvity-filters-button-icon-show');
    };
    ui.tooltip.bind(externalFilterIcon, 'Add Revvity Signals filters');
    
    //save search icon
    const saveSearchIcon = ui.icons.save(() => {
      saveSearchQuery(selectedLib[0].name, compoundType, qb)
    }, 'Save current search query');


    //load search icon
    const loadSearchIcon = ui.iconFA('folder-open', () => {
      loadSearchQuery(selectedLib[0].name, compoundType, qb);
    }, 'Load saved search query');
    
    const filtersButton = ui.div(externalFilterIcon);
    const saveButton = ui.div(saveSearchIcon);
    const loadButton = ui.div(loadSearchIcon);
    tv.setRibbonPanels([[filtersButton]]);
    //create filters panel
    tv.dockManager.dock(filtersDiv, 'left', null, 'Filters', 0.2);
    tv.dockManager.onClosed.subscribe((el: any) => {
      externalFilterIcon.classList.add('revvity-filters-button-icon-show');
    })
    ui.setUpdateIndicator(filtersDiv, true, 'Loading filters...');
    const qb = await initializeQueryBuilder(selectedLib[0].id, compoundType, initialSearchQuery);
    qb.validationError.subscribe((isInvalid: boolean) => {
      runSearchButton.disabled = isInvalid;
    });

    const runSearchButton = ui.button('Search', async () => {
      runSearch(qb, tv, selectedLib[0].id, compoundType, libName, true);
    });

    //set initial validation status
    runSearchButton.disabled = qb.invalid;

    ui.setUpdateIndicator(filtersDiv, false);
    filtersDiv.append(ui.divH([saveButton, loadButton], 'revvity-saved-searches-icons-div'));
    filtersDiv.append(qb.root);
    filtersDiv.append(ui.div(runSearchButton, { style: { paddingLeft: '4px' } }));

    ui.onSizeChanged(filtersDiv).subscribe(() => {
      updateQueryBuilderLayout(qb, filtersDiv.clientWidth, selectedLib[0].id, compoundType);
    });

    if (initialSearchQuery)
      runSearch(qb, tv, selectedLib[0].id, compoundType, libName, changePath);
  }
}

export async function runSearch(qb: QueryBuilder, tv: DG.TableView, libId: string, compoundType: string, libName: string, changePath?: boolean) {
  qb.saveConditionToHistory();
  ui.setUpdateIndicator(tv.grid.root, true, 'Searching...');
  const condition = qb.condition;
  const resultDf = await runSearchQuery(libId, compoundType, condition);
  tv.dataFrame = resultDf;
  setUserColumnsStyle(tv);
  if (changePath)
    tv.path = createPath([libName, getViewNameByCompoundType(compoundType), 'search', JSON.stringify(condition)]);
  ui.setUpdateIndicator(tv.grid.root, false);
}


export async function initializeQueryBuilder(libId: string, compoundType: string,
  initialSearchQuery?: ComplexCondition): Promise<QueryBuilder> {
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

  return new QueryBuilder(filterFields, initialSearchQuery, QueryBuilderLayout.Narrow, `Revvity Signals|${libId}|${compoundType}`, true);

}


function updateQueryBuilderLayout(qb: QueryBuilder, width: number, libId: string, compoundType: string) {
  if (!qb) return;

  // Switch to narrow layout if width is less than 300px, otherwise use standard
  const newLayout = width < 300 ? QueryBuilderLayout.Narrow : QueryBuilderLayout.Standard;

  if (qb.getLayout() !== newLayout) {
    qb.setLayout(newLayout);
  }
};

// Function to save search query
async function saveSearchQuery(libName: string, compoundType: string, queryBuilder: QueryBuilder) {

  const storageKey = `${libName}|${compoundType}`;
  const savedSearchesStr = grok.userSettings.getValue(SAVED_SEARCH_STORAGE, storageKey) || '{}';
  const savedSearches: { [key: string]: string } = JSON.parse(savedSearchesStr);
  
  // Generate default name
  let defaultName = 'new search';
  let counter = 1;
  while (savedSearches[defaultName]) {
    defaultName = `new search (${counter})`;
    counter++;
  }

  const nameInput = ui.input.string('Search Name', { value: defaultName });
  const dialog = ui.dialog('Save Search Query')
    .add(nameInput)
    .onOK(async () => {
      const searchName = nameInput.value;
      if (!searchName.trim()) {
        grok.shell.error('Search name cannot be empty');
        return;
      }

      // Save the current query condition
      savedSearches[searchName] = JSON.stringify(queryBuilder.condition);
      grok.userSettings.add(SAVED_SEARCH_STORAGE, storageKey, JSON.stringify(savedSearches));
      
      grok.shell.info(`Search query "${searchName}" saved successfully`);
    });

  dialog.show();
}

// Function to load search query
async function loadSearchQuery(libName: string, compoundType: string, queryBuilder: QueryBuilder) {
  const storageKey = `${libName}|${compoundType}`;
  const savedSearchesStr = grok.userSettings.getValue(SAVED_SEARCH_STORAGE, storageKey) || '{}';
  const savedSearches: { [key: string]: string } = JSON.parse(savedSearchesStr);
  
  const searchNames = Object.keys(savedSearches);
  if (searchNames.length === 0) {
    grok.shell.info('No saved searches found for this view');
    return;
  }

  const searchNamesDf = DG.DataFrame.fromColumns([DG.Column.fromList('string', 'name', searchNames)]);

  // Create typeahead input with saved search names
  const searchInput = ui.typeAhead('Search name', {
    source: {
      local: searchNames.map(name => ({ label: name, value: name }))
    },
    minLength: 0,
    limit: 10,
    highlight: true,
    debounceRemote: 100,
    preventSubmit: true,
  });

  searchInput.onChanged.subscribe(() => {
      searchNamesDf.rows.filter((row) => (row.name as string).includes(searchInput.value));
  });
  
  const dialog = ui.dialog('Load Search Query')
    .add(ui.divV([searchInput, searchNamesDf.plot.grid().root]))
    .onOK(async () => {
      if (searchNamesDf.currentRowIdx === -1) {
        grok.shell.error('Select saved search in the grid');
        return;
      }
      const selectedSearchName = searchNamesDf.get('name', searchNamesDf.currentRowIdx);
      // Load the saved query condition into the query builder
      try {
        const condition = JSON.parse(savedSearches[selectedSearchName]);
        queryBuilder.loadCondition(condition);
      } catch (e) {
        grok.shell.error('Failed to parse saved search query');
      }
    });

  dialog.show();
}


export async function createSavedSearchesSatistics(statsDiv: HTMLElement, libName?: string, entityType?: string) {

  const searchesDiv = ui.divV([getAppHeader()]);
  statsDiv.append(searchesDiv);
  const tableDiv = ui.div('', { style: { position: 'relative', paddingTop: '15px' } });
  searchesDiv.append(tableDiv);
  ui.setUpdateIndicator(tableDiv, true, 'Loading saved searches...');

  getRevvityLibraries().then(async (libs: RevvityLibrary[]) => {
    const libObjForTable: any[] = [];
    for (const lib of libs) {
      if (libName && lib.name !== libName)
        continue;
      for (const libType of lib.types) {
        if (entityType && libType.name !== entityType)
          continue;
        const storageKey = `${lib.name}|${libType.name}`;
        const savedSearchesStr = grok.userSettings.getValue(SAVED_SEARCH_STORAGE, storageKey) || '{}';
        const savedSearches: { [key: string]: string } = JSON.parse(savedSearchesStr);
        const searchNames = Object.keys(savedSearches);
        for (const search of searchNames)
          libObjForTable.push({ libName: lib.name, libType: libType.name, search: search });
      }
    }

    const statsElement = await createSavedSearchStats(libObjForTable, libName, entityType);
    tableDiv.append(statsElement);
    ui.setUpdateIndicator(tableDiv, false);
  });
}

async function createSavedSearchStats(savedSearchesForTable: any[], libName?: string, libType?: string): Promise<HTMLElement> {
  const outputArr: string[] = ['Library', 'Type', 'Search'];
  let libObjForTable: any[] = [];
  if (!savedSearchesForTable.length) {
    libObjForTable = await createLibsObjectForStatistics(libName);
    if (libType)
      libObjForTable = libObjForTable.filter((it) => it.libType === libType);
  }
  let header = 'Saved Searches';

  if (libName) {
    outputArr.splice(0, 1);
    header += ` > ${libName}`;
  }

  if (libType) {
    outputArr.splice(0, 1);
    const entityTypeName = getViewNameByCompoundType(libType);
    header += ` > ${entityTypeName.charAt(0).toUpperCase()}${entityTypeName.slice(1)}`;
  }

  const createStatsTable = () => {
    const table = ui.table(!savedSearchesForTable.length ? libObjForTable : savedSearchesForTable, (search) => {
      const arr = [
        search.libName ?? '',
        search.libType ? search.libType.charAt(0).toUpperCase() + search.libType.slice(1) : '',
        ui.link(savedSearchesForTable.length ? search.search : 'New search', () => {
          const node = grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps').getOrCreateGroup('Chem').getOrCreateGroup('Revvity Signals');
          node.expanded = true;
          if (!savedSearchesForTable.length)
            openRevvityNode(node, [search.libName, getViewNameByCompoundType(search.libType)], search.libType, search.libName, search.libType);
          else
            openRevvityNode(node, ['saved searches', search.libName, getViewNameByCompoundType(search.libType)], search.search, search.libName, search.libType, undefined, true);
        }),
      ];
      if (libName)
        arr.splice(0, 1);

      if (libType)
        arr.splice(0, 1);

      return arr;
    },
      outputArr);
    return table;
  }

  const statsDiv = ui.div([ui.h1(header, {style: {paddingLeft: '10px'}})]);

  if (!savedSearchesForTable.length) {
    let infoStr = `No searches have been saved `;
    if (libName) {
      infoStr += `for ${libName}`;
      if (libType) {
        infoStr += ` > ${libType.charAt(0).toUpperCase() + libType.slice(1)}`;
      }
    } 
    statsDiv.append(ui.info(infoStr));
  }

  const table = createStatsTable();
  if (outputArr.length === 1)
      table.classList.add('revvity-signals-statistics-hide-header');
  statsDiv.append(table);
  return statsDiv;
}