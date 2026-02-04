/* Do not change these import lines to match external modules in webpack configuration */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {
  QueryBuilder,
  QueryBuilderLayout,
  ComplexCondition,
  Operators,
  SUGGESTIONS_FUNCTION,
  ConditionRegistry,
  MoleculeConditionEditor,
} from '@datagrok-libraries/utils/src/query-builder/query-builder';
import {MolTrackDockerService} from '../services/moltrack-docker-service';
import {
  MolTrackSearchQuery,
  MolTrackFilter,
  MolTrackComplexCondition,
  MolTrackProperty,
  MolTrackEntityType,
  molTrackSearchMapping,
  MolTrackSearchAggregation,
  MolTrackSearch,
  MolTrackSearchHistoryItem,
} from '../utils/types';
import {
  EXCLUDE_SEARCH_FIELDS,
  EXCLUDE_SEARCH_OUTPUT_FIELDS,
  Scope,
  STRUCTURE_FIELDS,
  STRUCTURE_SEARCH_FIELD,
  STRING_AGGREGATIONS,
  NUMERIC_AGGREGATIONS,
  MOLTRACK_ENTITY_TYPE,
  MOLTRACK_ENDPOINT,
  MOLTRACK_ENTITY_LEVEL,
  MOLTRACK_IS_STATIC_FIELD,
  PROP_NUM_TYPES,
  MOL_COL_NAME,
  CORPORATE_COMPOUND_ID_COL_NAME,
  SEARCH_NODE,
  SAVED_SEARCHES_NODE,
  excludedScopes,
  MAX_USER_SETTINGS_STRING_LENGTH,
  CHUNK_KEY_SUFFIX,
} from '../utils/constants';
import {funcs} from '../package-api';
import dayjs, {Dayjs} from 'dayjs';
import {_package} from '../package';
import {Subject} from 'rxjs';
import {awaitCheck} from '@datagrok-libraries/test/src/test';
import {createPathFromArr, getAppHeader, getStatisticsWidget} from '../utils/view-utils';
import {applyMolTrackLayout, saveMolTrackLayout} from '../utils/layout';

export type MolTrackSearchFields = {
  direct?: DG.Property[],
  dynamic?: DG.Property[],
}

const FIRST_COL_NAMES = [
  `${Scope.COMPOUNDS}.${MOL_COL_NAME}`,
  `${Scope.COMPOUNDS}.details.${CORPORATE_COMPOUND_ID_COL_NAME}`];
export const SAVED_SEARCH_STORAGE = 'MolTrackSavedSearch';

export let molTrackSearchFieldsArr: DG.Property[] | null = null;
export let openedSearchView: DG.ViewBase | null = null;

export function getSavedSearches(entityLevel: Scope): Record<string, string> {
  const allSavedSearches = grok.userSettings.get(SAVED_SEARCH_STORAGE) || {};
  const savedSearches: Record<string, string> = {};
  const prefix = `${entityLevel}_`;

  for (const key of Object.keys(allSavedSearches)) {
    if (!key.startsWith(prefix)) continue;

    const parts = key.split('_');
    const searchName = parts[1];

    if (!savedSearches[searchName])
      savedSearches[searchName] = getFullFilterStr(entityLevel, searchName);
  }

  return savedSearches;
}

async function saveSearch(entityLevel: Scope, savedSearch: MolTrackSearch) {
  const savedSearches: { [key: string]: string } = getSavedSearches(entityLevel);

  // Generate default name
  let defaultName = 'new search';
  let counter = 1;
  while (savedSearches[defaultName]) {
    defaultName = `new search (${counter})`;
    counter++;
  }

  const nameInput = ui.input.string('Search Name', {value: defaultName});
  nameInput.addValidator(
    (search: string) => savedSearches[search] != null ? 'A saved search with this name already exists' : null);

  const dialog = ui.dialog('Save Search Query')
    .add(nameInput)
    .onOK(async () => {
      const searchName = nameInput.value;
      if (!searchName.trim()) {
        grok.shell.error('Search name cannot be empty');
        return;
      }

      //add corresponding saved search to browse view tree
      const savedSearchesNode = grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps')
        .getOrCreateGroup('Chem').getOrCreateGroup('MolTrack').getOrCreateGroup(SAVED_SEARCHES_NODE)
        .getOrCreateGroup(entityLevel.charAt(0).toUpperCase() + entityLevel.slice(1));
      const newSearchNode = savedSearchesNode.item(searchName);
      newSearchNode.onSelected.subscribe(async () => {
        createSearchView(searchName, entityLevel, JSON.parse(savedSearches[searchName]), true);
      });

      // Save the current query condition
      const savedSearchStr = JSON.stringify(savedSearch);
      savedSearches[searchName] = savedSearchStr;
      saveSearchChunked(searchName, entityLevel, savedSearchStr);
      grok.shell.info(`Search query "${searchName}" saved successfully`);
    });

  dialog.show();
}

function saveSearchChunked(searchName: string, entityLevel: string, savedSearch: string) {
  const baseKey = `${entityLevel}_${searchName}`;
  const chunkSize = MAX_USER_SETTINGS_STRING_LENGTH;
  const chunksCount = Math.ceil(savedSearch.length / chunkSize);

  for (let i = 0; i < chunksCount; i++) {
    const start = i * chunkSize;
    const end = start + chunkSize;
    const chunkData = savedSearch.substring(start, end);

    const key = i === 0 ? baseKey : `${baseKey}_${i}${CHUNK_KEY_SUFFIX}`;
    const nextKey = i < chunksCount - 1 ? `${baseKey}_${i + 1}${CHUNK_KEY_SUFFIX}` : '';
    const value = JSON.stringify({data: chunkData, next: nextKey});

    grok.userSettings.add(SAVED_SEARCH_STORAGE, key, value);
  }
}

function deleteSearchChunked(entityLevel: string, searchName: string) {
  const baseKey = `${entityLevel}_${searchName}`;
  let key = baseKey;

  while (key) {
    const stored = grok.userSettings.getValue(SAVED_SEARCH_STORAGE, key);
    if (!stored) break;

    const {next} = JSON.parse(stored);
    grok.userSettings.delete(SAVED_SEARCH_STORAGE, key);
    key = next;
  }
}

export function getFullFilterStr(entityLevel: string, searchName: string): string {
  const baseKey = `${entityLevel}_${searchName}`;
  let key = baseKey;
  let fullStr = '';

  while (key) {
    const stored = grok.userSettings.getValue(SAVED_SEARCH_STORAGE, key);
    if (!stored) break;

    const {data, next} = JSON.parse(stored);
    fullStr += data;
    key = next;
  }

  return fullStr;
}

function openSavedSearch(entityLevel: Scope, queryBuilder: QueryBuilder, outputsFieldsInput: DG.InputBase,
  aggrContainer: HTMLElement, menuFieldsForAggr: DG.Property[], aggregations: MolTrackSearchAggregation[],
  validationErrorSubj: Subject<string>): void {
  const savedSearches: { [key: string]: string } = getSavedSearches(entityLevel);
  const searchNames = Object.keys(savedSearches);
  if (searchNames.length === 0) {
    grok.shell.info('No saved searches found for this view');
    return;
  }

  const searchNamesDf = DG.DataFrame.fromColumns([DG.Column.fromList('string', 'name', searchNames)]);
  searchNamesDf.columns.addNewString('delete');

  // Create typeahead input with saved search names
  const searchInput = ui.typeAhead('Search name', {
    source: {
      local: searchNames.map((name) => ({label: name, value: name})),
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

  const savedSearchesGrid = searchNamesDf.plot.grid({allowColumnMenu: false, allowSorting: false});
  savedSearchesGrid.onCellDoubleClick.subscribe((gc: DG.GridCell) => {
    if (gc.isTableCell && gc.gridColumn.name === 'name') {
      const searchName = gc.grid.dataFrame.col('name')?.get(gc.gridRow);
      const search: MolTrackSearch = JSON.parse(savedSearches[searchName]);
      loadSearchQuery(search, queryBuilder, outputsFieldsInput, aggrContainer, menuFieldsForAggr,
        aggregations, validationErrorSubj);
      dialog.close();
    }
  });
  const deleteSearchCol = savedSearchesGrid.columns.byName('delete');

  if (deleteSearchCol) {
    deleteSearchCol.cellType = 'html';

    const savedSearchesNode = grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps')
      .getOrCreateGroup('Chem').getOrCreateGroup('MolTrack').getOrCreateGroup(SAVED_SEARCHES_NODE)
      .getOrCreateGroup(entityLevel.charAt(0).toUpperCase() + entityLevel.slice(1));
    savedSearchesGrid.onCellPrepare(function(gc) {
      if (gc.isTableCell && gc.gridColumn.name === 'delete') {
        const removeIcon = ui.div(ui.icons.delete(() => {
          const searchName = gc.grid.dataFrame.col('name')?.get(gc.gridRow);
          //remove row from table
          searchNamesDf.rows.removeAt(gc.gridRow, 1);
          //remove search from tree
          const items = savedSearchesNode.items.filter((it) => it.text === searchName);
          if (items.length)
            items[0].remove();
          //save changes to storage
          if (savedSearches[searchName])
            delete savedSearches[searchName];
          deleteSearchChunked(entityLevel, searchName);
        }), 'moltrack-saved-search-delete');
        gc.style.element = removeIcon;
      }
    });
  }

  const dialog = ui.dialog('Open saved search')
    .add(ui.divV([searchInput, savedSearchesGrid.root]))
    .onOK(async () => {
      if (searchNamesDf.currentRowIdx === -1) {
        grok.shell.error('Select saved search in the grid');
        return;
      }
      const selectedSearchName = searchNamesDf.get('name', searchNamesDf.currentRowIdx);
      // Load the saved query condition into the query builder
      try {
        const savedSearch: MolTrackSearch = JSON.parse(savedSearches[selectedSearchName]);
        loadSearchQuery(savedSearch, queryBuilder, outputsFieldsInput, aggrContainer, menuFieldsForAggr,
          aggregations, validationErrorSubj);
      } catch (e) {
        grok.shell.error('Failed to load saved search query');
      }
    });

  dialog.show();
}

export function loadSearchQuery(savedSearch: MolTrackSearch, queryBuilder: QueryBuilder,
  outputsFieldsInput: DG.InputBase, aggrContainer: HTMLElement, menuFieldsForAggr: DG.Property[],
  aggregations: MolTrackSearchAggregation[],
  validationErrorSubj: Subject<string>) {
  queryBuilder.loadCondition(savedSearch.condition);
  savedSearch.condition = queryBuilder.condition;
  outputsFieldsInput.value = savedSearch.outputCols
    .map((it) => DG.Column.fromType((it.type as string === 'uuid' ? 'string' : it.type) as any, it.name));
  ui.empty(aggrContainer);
  aggregations = [];
  savedSearch.aggregations?.forEach((aggr) => {
    const aggregationRow = createAggregationRow(aggrContainer, menuFieldsForAggr, aggregations,
      validationErrorSubj, aggr);
    aggrContainer.append(aggregationRow);
  });
}

export function saveSearchHistory(key: string, value: string) {
  let collection: MolTrackSearchHistoryItem[] = getSearchHistory(key);
  const newItem = {date: new Date(), value: value};
  if (!collection.length)
    localStorage.setItem(key, JSON.stringify([newItem]));
  else {
    collection = collection.filter((it) => it.value !== value);
    localStorage.setItem(key, JSON.stringify([newItem, ...collection.slice(0, 9)]));
  }
}

export function getSearchHistory(key: string): MolTrackSearchHistoryItem[] {
  const stringItems = localStorage.getItem(key);
  return JSON.parse(!stringItems ? '[]' : stringItems);
}

export async function loadSearchFields() {
  if (!molTrackSearchFieldsArr) {
    try {
      molTrackSearchFieldsArr = await createSearchFileds();
    } catch (e: any) {
      molTrackSearchFieldsArr = null;
      throw e;
    }
  }
}

export async function createSearchPanel(tv: DG.TableView, entityLevel: Scope, initialQuery?: MolTrackSearch,
  changePath?: boolean) {
  const filtersDiv = ui.div([]);
  let queryBuilder: QueryBuilder | undefined = undefined;

  const updateQueryBuilderLayout = (width: number) => {
    if (!queryBuilder) return;

    // Switch to narrow layout if width is less than 300px, otherwise use standard
    const newLayout = width < 300 ? QueryBuilderLayout.Narrow : QueryBuilderLayout.Standard;

    if (queryBuilder.getLayout() !== newLayout)
      queryBuilder.setLayout(newLayout);
  };

  ui.onSizeChanged(filtersDiv).subscribe(() => {
    updateQueryBuilderLayout(filtersDiv.clientWidth);
  });

  const initializeQueryBuilder = async (entityLevel: Scope) => {
    const validationErrorSubj = new Subject<string>();
    tv.dockManager.dock(filtersDiv, 'left', null, 'Filters', 0.2);
    ui.setUpdateIndicator(filtersDiv, true, 'Loading filters...');
    try {
      await loadSearchFields();
      if (!molTrackSearchFieldsArr) {
        grok.shell.error(`Filter fields haven't been loaded`);
        return;
      }
      const filterFields: DG.Property[] = molTrackSearchFieldsArr
        .filter((it) => it.options[MOLTRACK_ENTITY_LEVEL] === entityLevel && !EXCLUDE_SEARCH_FIELDS.includes(it.name));
      if (!filterFields.length) {
        grok.shell.warning(`No search fields found for ${entityLevel}`);
        return;
      }
      const entityType = filterFields[0].options[MOLTRACK_ENTITY_TYPE];
      const endpoint = filterFields[0].options[MOLTRACK_ENDPOINT];

      queryBuilder = new QueryBuilder(filterFields, undefined, QueryBuilderLayout.Narrow);
      queryBuilder.validationError.subscribe((error) => {
        runSearchButton.disabled = error;
        saveSearchIcon.classList.toggle('moltrack-disable-icon', error);
      });
      const df = DG.DataFrame.create();
      const outputFields = molTrackSearchFieldsArr
        .filter((it) => it.options[MOLTRACK_ENTITY_LEVEL] === entityLevel &&
          !EXCLUDE_SEARCH_OUTPUT_FIELDS.includes(it.name) && !STRUCTURE_FIELDS.includes(it.name));
      const defaultFields = outputFields.map((it) => it.name);
      outputFields.forEach((it) => df.columns
        .add(DG.Column.fromType((it.type as string === 'uuid' ? 'string' : it.type) as any, it.name)));
      const outputFieldsInput = ui.input.columns('Output', {
        table: df,
        value: df.columns.toList(),
        checked: defaultFields,
        nullable: false,
        onValueChanged(value, input) {
          const isEmptyList = value.length == 0;
          runSearchButton.classList.toggle('dim', isEmptyList);
          runSearchButton.disabled = isEmptyList;
        },
      });
      outputFieldsInput.root.style.paddingLeft = '12px';

      validationErrorSubj.subscribe((val: string) => {
        runSearchButton.disabled = val !== '';
      });

      const menuFieldsForAggr = molTrackSearchFieldsArr
        .filter((it) => it.type === 'string' || PROP_NUM_TYPES.includes(it.type) &&
      (!it.semType || it.semType && !Object.values(DG.SEMTYPE).includes(it.semType)));

      const aggregations: MolTrackSearchAggregation[] = [];
      // Create aggregations section
      const aggregationsDiv = ui.divV([]);
      const aggregationsHeader = ui.divH([
        ui.divText('Aggregations'),
        ui.icons.add(() => {
          const aggregationRow = createAggregationRow(aggregationsContainer, menuFieldsForAggr, aggregations,
            validationErrorSubj);
          aggregationsContainer.append(aggregationRow);
        }),
      ], 'moltrack-search-aggr-header-div');

      const aggregationsContainer = ui.divV([]);
      aggregationsDiv.append(aggregationsHeader);
      aggregationsDiv.append(aggregationsContainer);

      const runSearchButton = ui.button('Search', async () => {
        if (queryBuilder) {
          const search: MolTrackSearch = {
            condition: queryBuilder?.condition,
            outputCols: outputFieldsInput.value!.map((col) => { return {name: col.name, type: col.type}; }),
            aggregations: aggregations,
          };
          await runSearch(tv, search, entityLevel, entityType, endpoint, true);
        }
      });

      //save search icon
      const saveSearchIcon = ui.icons.save(() => {
        if (queryBuilder) {
          const savedSearch: MolTrackSearch = {
            condition: queryBuilder?.condition,
            outputCols: outputFieldsInput.value!.map((col) => { return {name: col.name, type: col.type}; }),
            aggregations: aggregations,
          };
          saveSearch(entityLevel, savedSearch);
        }
      }, 'Save current search query');

      queryBuilder?.structureChanged.subscribe(async () => {
        if (!queryBuilder || !outputFieldsInput.value) return;

        const fields = outputFieldsInput.value.map((col) => col.name);
        const molTrackQuery = convertQueryBuilderConditionToMolTrackQuery(
          queryBuilder.condition,
          entityLevel,
          fields,
          aggregations ?? [],
          entityType
        );

        const isValid = await MolTrackDockerService.validateSearchQuery(molTrackQuery);
        queryBuilder.validationError.next(!isValid);
      });

      //load search icon
      const openSearchIcon = ui.iconFA('folder-open', () => {
        if (queryBuilder) {
          openSavedSearch(entityLevel, queryBuilder, outputFieldsInput, aggregationsContainer,
            menuFieldsForAggr, aggregations, validationErrorSubj);
        };
      }, 'Load saved search query');

      //search history icon
      const historyIcon = ui.iconFA('history', () => {
        const items: MolTrackSearchHistoryItem[] = getSearchHistory(`${_package.name}|${entityLevel}`);
        const menu = DG.Menu.popup();
        menu.items(items.map((it) => ui.tools.click(
          ui.divH([
            ui.divText(it.date.toString(), {style: {color: '#7990A5'}}),
            ui.divText(molTrackSerachToString(JSON.parse(it.value) as MolTrackSearch, queryBuilder),
              'moltrack-serch-history-item'),
          ], {style: {gap: '5px'}}),
          () => {
            loadSearchQuery(JSON.parse(it.value), queryBuilder!, outputFieldsInput, aggregationsContainer,
              menuFieldsForAggr, aggregations, validationErrorSubj);
          })), () => {});
        menu.show();
      }, 'Query History');

      filtersDiv.append(ui.divH([historyIcon, saveSearchIcon, openSearchIcon], 'moltrack-saved-searches-icons-div'));
      filtersDiv.append(queryBuilder.root);
      filtersDiv.append(outputFieldsInput.root);
      filtersDiv.append(aggregationsDiv);
      filtersDiv.append(ui.div(runSearchButton, {style: {paddingLeft: '4px', paddingTop: '10px'}}));
      createFiltersIcon(tv, filtersDiv);
      if (initialQuery) {
        if (initialQuery.condition.conditions.length) {
          loadSearchQuery(initialQuery, queryBuilder, outputFieldsInput, aggregationsContainer, menuFieldsForAggr,
            aggregations, validationErrorSubj);
        }
        await runSearch(tv, initialQuery, entityLevel, entityType, endpoint, !!changePath);
      }
    } catch (e: any) {
      grok.shell.error(`Error loading filters: ${e?.message ?? e}`);
    } finally {
      ui.setUpdateIndicator(filtersDiv, false);
    }
  };

  initializeQueryBuilder(entityLevel);
}

export function molTrackSerachToString(search: MolTrackSearch, queryBuilder?: QueryBuilder): string {
  if (!queryBuilder)
    return '';
  const parts: string[] = [];

  // Add condition string
  if (search.condition) {
    const conditionStr = queryBuilder.conditionToString(search.condition);
    if (conditionStr && conditionStr !== 'No conditions')
      parts.push(`Conditions: ${conditionStr}`);
  }

  // Add output columns
  if (search.outputCols && search.outputCols.length > 0) {
    const outputFields = search.outputCols.map((col) => col.name).join(', ');
    parts.push(`Output: ${outputFields}`);
  }

  // Add aggregations
  if (search.aggregations && search.aggregations.length > 0) {
    const aggrStrings = search.aggregations.map((aggr) => `${aggr.operation}(${aggr.field})`);
    parts.push(`Aggregations: ${aggrStrings.join(', ')}`);
  }

  return parts.length > 0 ? parts.join(' | ') : 'Empty search';
}

export async function createEmptyMolTrackQuery(scope: Scope): Promise<MolTrackSearch> {
  await loadSearchFields();
  let outputFields: { name: string, type: string }[] = [];
  if (!molTrackSearchFieldsArr) {
    grok.shell.warning(`Properties haven't been loaded for ${scope}`);
    outputFields = [{name: CORPORATE_COMPOUND_ID_COL_NAME, type: DG.TYPE.STRING}];
  } else {
    outputFields = (molTrackSearchFieldsArr
      .filter((it) => it.options[MOLTRACK_ENTITY_LEVEL] === scope &&
        !EXCLUDE_SEARCH_OUTPUT_FIELDS.includes(it.name))).map((it) => { return {name: it.name, type: it.type}; });
  }

  const search: MolTrackSearch = {
    condition: {
      logicalOperator: Operators.Logical.and,
      conditions: [],
    },
    outputCols: outputFields,
  };
  return search;
}

export async function runSearch(tv: DG.TableView, search: MolTrackSearch, entityLevel: Scope,
  entityType: MolTrackEntityType, endpoint: string, changePath: boolean) {
  try {
    ui.setUpdateIndicator(tv.grid.root, true, 'Searching...');
    if (changePath)
      tv.path = createPathFromArr(tv, [SEARCH_NODE, entityLevel, JSON.stringify(search)]);
    if (!search.outputCols.length)
      throw new Error(`At least one input field should be selected`);
    const outputFields = search.outputCols!
      .map((it) => `${entityLevel}.${isDynamicField(it.name, entityType) ? 'details.' : ''}${it.name}`);
    const molTrackQuery = convertQueryBuilderConditionToMolTrackQuery(search!.condition,
      entityLevel, outputFields, search.aggregations ?? [], entityType);
    const result = await grok.functions.call('MolTrack:search', {
      query: JSON.stringify(molTrackQuery),
      entityEndpoint: endpoint,
    }); ;
    const resultDf = DG.DataFrame.fromObjects(result.data);
    if (resultDf) {
      await loadSearchFields();
      const filterFields: DG.Property[] = molTrackSearchFieldsArr ?
        molTrackSearchFieldsArr.filter((it) => it.options[MOLTRACK_ENTITY_LEVEL] === entityLevel) : [];
      resultDf?.columns.names().forEach(async (colName: string) => {
        const propIdx = filterFields
          .findIndex((it) => it.name.toLowerCase() ===
            colName.replace(`${entityLevel.toLowerCase()}.`, '').replace(`details.`, '').toLowerCase());
        const friendlyName = propIdx !== -1 ?
          filterFields[propIdx].friendlyName ?? filterFields[propIdx].name :
          colName.replace(`${entityLevel.toLowerCase()}.`, '').replace(`details.`, '');
        resultDf.col(colName)!.setTag('friendlyName', friendlyName);
      });
      //save history in case of successful search
      saveSearchHistory(`${_package.name}|${entityLevel}`, JSON.stringify(search));
    }
    updateView(tv, resultDf ?? DG.DataFrame.create(), entityLevel, tv.name);
  } catch (e: any) {
    grok.shell.error(e?.message ?? e);
    tv.dataFrame = DG.DataFrame.create();
  } finally {
    ui.setUpdateIndicator(tv.grid.root, false);
  }
}

export function createFiltersIcon(tv: DG.TableView, filtersDiv: HTMLDivElement) {
  //create filters icon
  const externalFilterIcon = document.createElement('div');
  externalFilterIcon.innerHTML = `
<svg width="24" height="24" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
  <!-- Funnel Body -->
  <path d="M4 4H15L10 10V17L8 19V10L4 4Z" stroke="currentColor" stroke-width="1" 
        stroke-linecap="round" stroke-linejoin="round"/>
  
  <!-- The arrow group is rotated -45 degrees around point (16, 13) -->
  <g transform="rotate(-45, 16, 13)">
    <!-- Long Horizontal Arrow Line (starts at the new x=16, y=13) -->
    <path d="M16 13H24" stroke="currentColor" stroke-width="1" stroke-linecap="round"/>
    <!-- Arrowhead (adjusted for the new starting point) -->
    <path d="M21 10L24 13L21 16" stroke="currentColor" stroke-width="1" stroke-linecap="round" stroke-linejoin="round"/>
  </g>
</svg>`;
  externalFilterIcon.className = 'moltrack-filters-button-icon';
  externalFilterIcon.onclick = () => {
    tv.dockManager.dock(filtersDiv, 'left', null, 'Filters', 0.2);
    externalFilterIcon.classList.remove('moltrack-filters-button-icon-show');
  };
  tv.dockManager.onClosed.subscribe((el: any) => {
    externalFilterIcon.classList.add('moltrack-filters-button-icon-show');
  });
  ui.tooltip.bind(externalFilterIcon, 'Add MolTrack filters');
  const filtersButton = ui.div(externalFilterIcon);
  tv.setRibbonPanels([[filtersButton]], false);
}

export async function createSearchFileds(): Promise<DG.Property[]> {
  const propArr: DG.Property[] = [];
  const promises = [
    funcs.fetchDirectSchema(),
    funcs.fetchSchema(),
  ];
  const results = await Promise.all(promises);
  const staticFields: MolTrackProperty[] = JSON.parse(results[0]);
  const dynamicFields: MolTrackProperty[] = JSON.parse(results[1]);

  const addProperties = (props: MolTrackProperty[], isStatic: boolean) => {
    for (const prop of props) {
      const isStructureField = STRUCTURE_FIELDS.includes(prop.name);
      const structureAlreadyExists = propArr.filter((it) => it.options[MOLTRACK_ENTITY_TYPE] === prop.entity_type &&
            it.name === STRUCTURE_SEARCH_FIELD).length > 0;
      if (isStructureField && structureAlreadyExists)
        continue;
      const propOptions = isStructureField ?
        {
          name: STRUCTURE_SEARCH_FIELD,
          type: 'string',
          semType: DG.SEMTYPE.MOLECULE,
          friendlyName: prop.friendly_name,
        } :
        {
          name: prop.name,
          type: prop.value_type === 'uuid' ? 'string' : prop.value_type,
          semType: prop.semantic_type?.name,
          friendlyName: prop.friendly_name,
        };
      const dgProp = DG.Property.fromOptions(propOptions);
      dgProp.options[MOLTRACK_ENTITY_TYPE] = prop.entity_type;
      dgProp.options[MOLTRACK_IS_STATIC_FIELD] = isStatic;
      const propMapping = molTrackSearchMapping[prop.entity_type];
      if (propMapping) {
        dgProp.options[MOLTRACK_ENDPOINT] = propMapping.searchEndpoint;
        dgProp.options[MOLTRACK_ENTITY_LEVEL] = propMapping.level;
      }

      propArr.push(dgProp);
      const registeredSemType = dgProp.semType && Object.values(DG.SEMTYPE).includes(dgProp.semType);
      if (prop.value_type === 'string' && propMapping.level && !registeredSemType) {
        const fieldName = `${propMapping.level}${isStatic ? '' : '.details'}.${prop.name}`;
        const getPropSuggestions = async (text: string): Promise<string[]> => {
          const query: MolTrackSearchQuery = {
            level: propMapping.level,
            filter: {
              field: fieldName,
              operator: Operators.CONTAINS,
              value: text,
            },
            output: [fieldName],
          };
          const suggestions: string[] = [];
          try {
            const res = await MolTrackDockerService.search(query, propMapping.searchEndpoint);
            if (res.data.length)
              return res.data.map((it) => it[fieldName.toLowerCase()]);
          } catch (e: any) {
            grok.shell.error(e);
          }
          return suggestions;
        };
        dgProp.options[SUGGESTIONS_FUNCTION] = getPropSuggestions;
      }
    }
  };

  addProperties(staticFields, true);
  addProperties(dynamicFields, false);
  //Put boolean fields to the end of the list not to show boolean input by default
  const sortProps = (props: DG.Property[]) => {
    return props.sort((a, b) => {
      if (a.type === DG.TYPE.BOOL && b.type !== DG.TYPE.BOOL) return 1;
      if (a.type !== DG.TYPE.BOOL && b.type === DG.TYPE.BOOL) return -1;
      return 0;
    });
  };
  sortProps(propArr);
  return propArr;
}


export function convertQueryBuilderConditionToMolTrackQuery(
  cond: ComplexCondition,
  level: Scope,
  output: string[],
  aggregations: MolTrackSearchAggregation[],
  type: MolTrackEntityType,
): MolTrackSearchQuery {
  // If query builder condition contains only one condition in array,
  // create MolTrackSearchQuery where filter is a MolTrackSimpleCondition
  const query: MolTrackSearchQuery = {
    level: level,
    output: output,
  };
  if (cond.conditions.length) {
    query.filter = cond.conditions.length === 1 && isSimpleCondition(cond.conditions[0]) ?
      convertSimpleCondition(cond.conditions[0], level, type) :
      convertComplexCondition(cond, level, type);
  }
  if (aggregations.length)
    query.aggregations = aggregations;

  return query;
}

function convertComplexCondition(cond: ComplexCondition, level: Scope,
  type: MolTrackEntityType): MolTrackComplexCondition {
  if (cond.conditions.length === 0) {
    return {
      operator: Operators.Logical.and,
      conditions: [],
    };
  }

  const convertedConditions: MolTrackFilter[] = cond.conditions.map((condition) => {
    if (isSimpleCondition(condition))
      return convertSimpleCondition(condition, level, type);
    else
      return convertComplexCondition(condition as ComplexCondition, level, type);
  });

  return {
    operator: cond.logicalOperator,
    conditions: convertedConditions,
  };
}

function isDynamicField(field: string, type: MolTrackEntityType): boolean {
  let isDynamicProp = false;
  if (type) {
    const dynamicPropIdx = molTrackSearchFieldsArr!
      .findIndex((it) => it.name === field && it.options[MOLTRACK_ENTITY_TYPE] === type &&
        !it.options[MOLTRACK_IS_STATIC_FIELD]);
    isDynamicProp = dynamicPropIdx !== -1;
  }
  return isDynamicProp;
}

function convertSimpleCondition(cond: any, level: Scope, type: MolTrackEntityType): MolTrackFilter {
  const isDynamicProp = isDynamicField(cond.field, type);
  // Handle chemical similarity search special case
  if (cond.operator === Operators.IS_SIMILAR && cond.value &&
    typeof cond.value === 'object' && 'molecule' in cond.value) {
    return {
      field: `${level}${isDynamicProp ? '.details' : ''}.${cond.field}`,
      operator: Operators.IS_SIMILAR,
      value: cond.value.molecule,
      threshold: cond.value.threshold || null,
    };
  }

  // Handle regular simple conditions - operators are similar, no mapping needed
  return {
    field: `${level}${isDynamicProp ? '.details' : ''}.${cond.field}`,
    operator: cond.operator,
    value: cond.value ? cond.value instanceof dayjs ? convertDateToString(cond.value) : cond.value : null,
    threshold: cond.threshold || null,
  };
}

function isSimpleCondition(condition: any): boolean {
  return 'field' in condition && 'operator' in condition && 'value' in condition;
}

function convertDateToString(date: Dayjs) {
  const yyyy = date.year();
  const mm = String(date.month() + 1).padStart(2, '0');
  const dd = String(date.date()).padStart(2, '0');
  return `${yyyy}-${mm}-${dd}`;
}

function createAggregationRow(container: HTMLElement, menuFieldsForAggr: DG.Property[],
  aggregations: MolTrackSearchAggregation[], validationErrorSubj: Subject<string>,
  defaultAggr?: MolTrackSearchAggregation): HTMLElement {
  const row = ui.divH([], 'moltrack-search-aggr-row');

  const newAggr: MolTrackSearchAggregation = defaultAggr ?? {field: '', operation: ''};
  aggregations.push(newAggr);

  const sendValidationError = (errorMessage: string) => {
    validationErrorSubj.next(errorMessage);
    fieldTooltip = errorMessage;
    newAggr.field = '';
    Array.from(fieldInput.root.children).forEach((it) => it.classList.add('moltrack-invalid-aggr-field'));
  };

  const getFieldType = (): DG.TYPE | null => {
    const error = `Field ${fieldInput.value} not found. Click to select from list of available fields`;
    const splittedName = fieldInput.value.split('.');
    if (splittedName.length < 2) {
      sendValidationError(error);
      return null;
    }
    const prop = menuFieldsForAggr
      .filter((it) => it.options[MOLTRACK_ENTITY_LEVEL] === splittedName[0] &&
        it.name === (splittedName.length > 2 && splittedName[1] === 'details' ? splittedName[2] : splittedName[1]));
    if (!prop.length) {
      sendValidationError(error);
      return null;
    }
    newAggr.field = `${splittedName[0]}.${isDynamicField(prop[0].name,
      prop[0].options[MOLTRACK_ENTITY_TYPE]) ? 'details.' : ''}${splittedName[1]}`;
    const fieldType = prop[0].type;
    return fieldType;
  };

  const changeFieldValue = () => {
    const fieldType = getFieldType();
    updateAggregationOptions(aggregationInputDiv, fieldType, newAggr);
    Array.from(fieldInput.root.children).forEach((it) => it.classList.remove('moltrack-invalid-aggr-field'));
    //check if aggregations list contains invalid rows
    const containInvalidRows = aggregations.filter((it) => it.field === '').length > 0 ? 'Invalid aggregations' : '';
    validationErrorSubj.next(containInvalidRows);
  };

  const nonErrorTooltip = 'Click the field to see options';
  let fieldTooltip = nonErrorTooltip;
  const fieldInput = ui.input.string('', {
    value: defaultAggr?.field ?? '',
    nullable: false,
    onValueChanged: () => {
      changeFieldValue();
    },
  });
  fieldInput.classList.add('moltrack-search-aggr-field');
  ui.tooltip.bind(fieldInput.input, () => fieldTooltip);

  const aggregationInputDiv = ui.div();

  const removeButton = ui.icons.delete(() => {
    const idx = aggregations.findIndex((it) => it.field === newAggr.field && it.operation === newAggr.operation);
    if (idx !== -1)
      aggregations.splice(idx, 1);
    container.removeChild(row);
    const containInvalidRows = aggregations.filter((it) => it.field === '').length > 0 ? 'Invalid aggregations' : '';
    validationErrorSubj.next(containInvalidRows);
  });

  if (menuFieldsForAggr.length > 0) {
    const menu = DG.Menu.popup();
    if (menuFieldsForAggr.length) {
      menuFieldsForAggr.forEach((field) => {
        menu.item(`${field.options[MOLTRACK_ENTITY_LEVEL]}.${field.name}`, () => {
          fieldInput.value = `${field.options[MOLTRACK_ENTITY_LEVEL]}.${field.name}`;
        });
      });
      fieldInput.root.onclick = () => menu.show({element: fieldInput.root, y: fieldInput.root.offsetHeight});
    }
  }

  row.append(fieldInput.root);
  row.append(aggregationInputDiv);
  row.append(removeButton);

  //when added a new row, field input is empty, so disable SEARCH button
  if (!fieldInput.validate())
    validationErrorSubj.next(`Invalid aggregation field: ${fieldInput.value}`);

  if (defaultAggr)
    changeFieldValue();

  return row;
}

function updateAggregationOptions(aggregationInputDiv: any, fieldType: any, aggr: MolTrackSearchAggregation,
  defaultVal?: string) {
  const options = PROP_NUM_TYPES.includes(fieldType) ? NUMERIC_AGGREGATIONS : STRING_AGGREGATIONS;
  ui.empty(aggregationInputDiv);
  aggr.operation = options[0];
  const input = ui.input.choice('', {
    value: defaultVal ?? options[0],
    items: options,
    nullable: false,
    onValueChanged: () => {
      aggr.operation = input.value!;
    },
  });
  aggregationInputDiv.append(input.root);
}


export function createSearchNode(appNode:DG.TreeViewGroup, scope: Scope, initialQuery?: MolTrackSearch) {
  const formattedScope = scope
    .toLowerCase()
    .replace(/_/g, ' ')
    .replace(/\b\w/g, (c) => c.toUpperCase());

  appNode.getOrCreateGroup(SEARCH_NODE).item(formattedScope).onSelected.subscribe(async () => {
    createSearchView(formattedScope, scope, initialQuery);
  });
}

export async function createSearchView(viewName: string, scope: Scope, initialQuery?: MolTrackSearch,
  isSavedSearch?: boolean): Promise<DG.TableView> {
  if (openedSearchView?.name === viewName) {
    grok.shell.v = openedSearchView;
    return grok.shell.tv;
  }
  const df = DG.DataFrame.create();
  openedSearchView?.close();
  const tv = DG.TableView.create(df, false);
  openedSearchView = tv;
  openedSearchView.name = viewName;
  grok.shell.addPreview(openedSearchView);
  const initPath = [
    isSavedSearch ? SAVED_SEARCHES_NODE : SEARCH_NODE,
    isSavedSearch ? scope : viewName,
  ];
  if (isSavedSearch)
    initPath.push(viewName);
  else if (initialQuery)
    initPath.push(JSON.stringify(initialQuery));
  openedSearchView.path = createPathFromArr(openedSearchView, initPath);
  //if there is no initial search query - retrieve all data to initially show in the view
  if (!initialQuery)
    initialQuery = await createEmptyMolTrackQuery(scope);
  createSearchPanel(tv, scope as Scope, initialQuery, !isSavedSearch && initialQuery.condition.conditions.length > 0);
  return openedSearchView as DG.TableView;
}

export async function updateView(tv: DG.TableView, df: DG.DataFrame, scope: Scope, viewName?: string) {
  await grok.data.detectSemanticTypes(df);
  reorderSearchResColummns(df);
  tv.dataFrame = df;
  if (viewName)
    tv.dataFrame.name = viewName;

  applyMolTrackLayout(tv.grid, scope);
  const columns = tv.dataFrame.columns.names();
  if (columns.length)
    tv.dataFrame.currentCell = tv.dataFrame.cell(-1, tv.dataFrame.col(columns[0])!.name);

  /** Ensure the dataframe is available before saving the MolTrack layout. */
  const waitForGrid = (tv: DG.TableView, callback: (grid: DG.Grid) => void, interval = 100) => {
    const handle = setInterval(() => {
      if (tv.grid.dataFrame) {
        clearInterval(handle);
        callback(tv.grid);
      }
    }, interval);
  };

  waitForGrid(tv, (grid) => {
    grid.onPropertyValueChanged.subscribe(() => {
      saveMolTrackLayout(grid, scope);
    });

    tv.dataFrame.onMetadataChanged.subscribe(() => {
      saveMolTrackLayout(grid, scope);
    });
  });
}

export async function reorderSearchResColummns(df: DG.DataFrame) {
  const colNames = df.columns.names();
  const newColOrder = [];
  for (const colName of FIRST_COL_NAMES) {
    const index = colNames.indexOf(colName);
    if (index > -1) {
      colNames.splice(index, 1);
      newColOrder.push(colName);
    }
  }
  df.columns.setOrder(newColOrder.concat(colNames));
}

export async function openMolTrackSearchNode(nodesToExpand: string[]) {
  const treeNode = grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps')
    .getOrCreateGroup('Chem').getOrCreateGroup('MolTrack');
  treeNode.expanded = true;
  for (const nodeName of nodesToExpand) {
    await awaitCheck(() => treeNode.items.
      find((node) => node.text.toLowerCase() === nodeName.toLowerCase()) !== undefined,
    `${nodeName} haven't been loaded in 10 seconds`, 10000);
    const node = treeNode.items
      .find((node) => node.text.toLowerCase() === nodeName.toLowerCase()) as DG.TreeViewGroup;
    node.expanded = true;
  }
}

export async function handleSearchURL(url: string): Promise<DG.ViewBase> {
  if (url.startsWith('/'))
    url = url.slice(1);
  const componentsArr = url.split('/');
  if (componentsArr.length && (componentsArr[0] === SEARCH_NODE || componentsArr[0] === SAVED_SEARCHES_NODE)) {
    let idx = 0;
    const nodesToExpand = [];
    nodesToExpand.push(componentsArr[0]);
    idx++;
    if (componentsArr.length === idx) {
      openMolTrackSearchNode(nodesToExpand);
      if (nodesToExpand[0].toLowerCase() === SEARCH_NODE.toLowerCase())
        return await createSearchExpandableNode([SEARCH_NODE], () => getStatisticsWidget(createSearchView), false);
      if (nodesToExpand[0].toLowerCase() === SAVED_SEARCHES_NODE.toLowerCase()) {
        return await createSearchExpandableNode(
          [SAVED_SEARCHES_NODE],
          () => createSavedSearchesSatistics(undefined),
          false
        );
      }
    }
    if (componentsArr.length > idx) {
      const scope = componentsArr[1];
      if (componentsArr[0].toLowerCase() === SAVED_SEARCHES_NODE.toLowerCase()) {
        const scope = componentsArr[1];
        nodesToExpand.push(scope);
        idx++;
        if (componentsArr.length === idx) {
          openMolTrackSearchNode(nodesToExpand);
          if ((Object.values(Scope) as string[]).includes(scope.toLowerCase())) {
            return await createSearchExpandableNode([SAVED_SEARCHES_NODE, scope],
              () => createSavedSearchesSatistics(scope as Scope), false);
          } else {
            grok.shell.error(`Entity ${scope} doesn't exist in Moltrack`);
            return await createSearchExpandableNode([SAVED_SEARCHES_NODE],
              () => createSavedSearchesSatistics(undefined), false);
          }
        }
      } else {
        let initialSearch = undefined;
        //check for initial query in case of search url
        if (componentsArr.length > idx + 1)
          initialSearch = JSON.parse(componentsArr[2]);

        openMolTrackSearchNode(nodesToExpand);
        return await createSearchView(scope, scope.toLowerCase() as Scope, initialSearch, false);
      }
      //in case of saved search
      if (componentsArr.length > idx) {
        const savedSearchName = componentsArr[idx];
        const savedSearches = getSavedSearches(scope as Scope);
        if (savedSearches[savedSearchName]) {
          openMolTrackSearchNode(nodesToExpand);
          return await createSearchView(savedSearchName, scope.toLowerCase() as Scope,
            JSON.parse(savedSearches[savedSearchName]), true);
        } else {
          grok.shell.error(`Search ${savedSearchName} not found for ${scope}`);
          createSearchView(scope, scope.toLowerCase() as Scope, undefined, false);
          return await createSearchView(scope, scope.toLowerCase() as Scope, undefined, false);
        }
      }
    }
  }
  throw Error(`incorrect search path`);
}

export async function createSearchExpandableNode(viewpath: string[],
  getElement: () => Promise<HTMLElement>,
  addPreview: boolean = true,
): Promise<DG.ViewBase> {
  openedSearchView?.close();
  const header = getAppHeader();
  const contentDiv = ui.div('', 'moltrack-search-stats-div');
  const view = DG.View.fromRoot(ui.divV([header, contentDiv]));
  openedSearchView = addPreview ? grok.shell.addPreview(view) : view;

  const hasName = viewpath.length > 0;
  const name = hasName ? viewpath[viewpath.length - 1] : '';
  if (hasName) {
    openedSearchView.name = name;
    openedSearchView.path = createPathFromArr(openedSearchView, viewpath);
  }
  ui.setUpdateIndicator(contentDiv, true, `Loading ${name}...`);

  getElement()
    .then((res) => contentDiv.append(res))
    .catch((e) => grok.shell.error(e))
    .finally(() => ui.setUpdateIndicator(contentDiv, false));
  return openedSearchView;
}

export async function createSavedSearchesSatistics(scope?: Scope): Promise<HTMLElement> {
  const objForTable: any[] = [];
  const scopes = scope ? [scope] : Object.values(Scope).filter((scope) => !excludedScopes.includes(scope));
  const collectSearchesForScope = (scope: Scope) => {
    const savedSearches: { [key: string]: string } = getSavedSearches(scope);
    const searchNames = Object.keys(savedSearches);
    for (const search of searchNames)
      objForTable.push({scope: scope, search: search, query: JSON.parse(savedSearches[search])});
  };

  scopes.forEach((it) => collectSearchesForScope(it));

  const statsElement = await createSavedSearchStats(objForTable, scope);
  return statsElement;
}

async function createSavedSearchStats(savedSearchesForTable: any[], scope?: string): Promise<HTMLElement> {
  const outputArr: string[] = ['Scope', 'Search'];


  let scopesObjForTable: any[] = [];
  if (!savedSearchesForTable.length) {
    scopesObjForTable = Object.values(Scope)
      .filter((scope) => !excludedScopes.includes(scope)).map((it) => { return {scope: it}; });
    if (scope)
      scopesObjForTable = scopesObjForTable.filter((it) => it.scope === scope);
  }

  let header = SAVED_SEARCHES_NODE;

  if (scope) {
    outputArr.splice(0, 1);
    header += ` > ${scope.charAt(0).toUpperCase() + scope.slice(1)}`;
  }

  const createStatsTable = () => {
    const table = ui.table(!savedSearchesForTable.length ? scopesObjForTable : savedSearchesForTable, (search) => {
      const arr = [
        search.scope ? search.scope.charAt(0).toUpperCase() + search.scope.slice(1) : '',
        ui.link(savedSearchesForTable.length ? search.search : 'New search', () => {
          const node = grok.shell.browsePanel.mainTree
            .getOrCreateGroup('Apps').getOrCreateGroup('Chem').getOrCreateGroup('MolTrack');
          node.expanded = true;
          if (!savedSearchesForTable.length)
            createSearchView(search.scope.charAt(0).toUpperCase() + search.scope.slice(1), search.scope);
          else
            createSearchView(search.search, search.scope, search.query, true);
        }),
      ];

      if (scope)
        arr.splice(0, 1);

      return arr;
    },
    outputArr);
    return table;
  };

  const statsDiv = ui.div([ui.h1(header, {style: {paddingLeft: '10px'}})]);

  if (!savedSearchesForTable.length) {
    let infoStr = `No searches have been saved `;
    if (scope)
      infoStr += `for ${scope}`;
    statsDiv.append(ui.info(infoStr));
  }

  const table = createStatsTable();
  if (outputArr.length === 1)
    table.classList.add('moltrack-statistics-hide-header');
  statsDiv.append(table);
  return statsDiv;
}


//register operators for molecule semType, applicable for MolTrack
ConditionRegistry.getInstance()
  .registerSemTypeOperators(DG.SEMTYPE.MOLECULE, [Operators.IS_SIMILAR, 'has substructure', 'is substructure of']);
ConditionRegistry.getInstance()
  .registerEditor(DG.TYPE.STRING, DG.SEMTYPE.MOLECULE, 'has substructure', MoleculeConditionEditor);
ConditionRegistry.getInstance()
  .registerEditor(DG.TYPE.STRING, DG.SEMTYPE.MOLECULE, 'is substructure of', MoleculeConditionEditor);
