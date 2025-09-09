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
import { MolTrackDockerService } from './moltrack-docker-service';
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
} from './types';
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
} from './constants';
import { funcs } from '../package-api';
import dayjs, { Dayjs } from 'dayjs';
import { _package } from '../package';
import { Subject } from 'rxjs';
import { createPathFromArr } from './utils';
import { awaitCheck } from '@datagrok-libraries/utils/src/test';

export type MolTrackSearchFields = {
  direct?: DG.Property[],
  dynamic?: DG.Property[],
}

const FIRST_COL_NAMES = [
  `${Scope.COMPOUNDS}.${MOL_COL_NAME}`,
  `${Scope.COMPOUNDS}.details.${CORPORATE_COMPOUND_ID_COL_NAME}`];
export const SAVED_SEARCH_STORAGE = 'MolTrackSavedSearch';

export let molTrackSearchFieldsArr: DG.Property[] | null = null;
let openedSearchView: DG.TableView | null = null;

export function getSavedSearches(entityLevel: Scope): {[key: string]: string} {
  const savedSearchesStr = grok.userSettings.getValue(SAVED_SEARCH_STORAGE, entityLevel) || '{}';
  const savedSearches: { [key: string]: string } = JSON.parse(savedSearchesStr);
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
      savedSearches[searchName] = JSON.stringify(savedSearch);
      grok.userSettings.add(SAVED_SEARCH_STORAGE, entityLevel, JSON.stringify(savedSearches));
      grok.shell.info(`Search query "${searchName}" saved successfully`);
    });

  dialog.show();
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

  // Create typeahead input with saved search names
  const searchInput = ui.typeAhead('Search name', {
    source: {
      local: searchNames.map((name) => ({ label: name, value: name })),
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

  const dialog = ui.dialog('Open saved search')
    .add(ui.divV([searchInput, searchNamesDf.plot.grid().root]))
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
  const newItem = { date: new Date(), value: value };
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
      if (!molTrackSearchFieldsArr) {
        try {
          molTrackSearchFieldsArr = await createSearchFileds();
        } catch (e: any) {
          molTrackSearchFieldsArr = null;
          throw e;
        }
      }
      const filterFields: DG.Property[] = molTrackSearchFieldsArr
        .filter((it) => it.options[MOLTRACK_ENTITY_LEVEL] === entityLevel);
      if (!filterFields.length) {
        grok.shell.warning(`No search fields found for ${entityLevel}`);
        return;
      }
      const entityType = filterFields[0].options[MOLTRACK_ENTITY_TYPE];
      const endpoint = filterFields[0].options[MOLTRACK_ENDPOINT];

      queryBuilder = new QueryBuilder(filterFields, undefined, QueryBuilderLayout.Narrow);
      queryBuilder.validationError.subscribe((error) => {
        runSearchButton.disabled = error;
      });
      const df = DG.DataFrame.create();
      const outputFields = filterFields
        .filter((it) => !EXCLUDE_SEARCH_OUTPUT_FIELDS.includes(it.name) && it.name !== STRUCTURE_SEARCH_FIELD);
      const defaultFields = STRUCTURE_FIELDS.concat(filterFields.map((it) => it.name));
      outputFields.forEach((it) => df.columns
        .add(DG.Column.fromType((it.type as string === 'uuid' ? 'string' : it.type) as any, it.name)));
      //adding structure fields to return
      STRUCTURE_FIELDS.forEach((it) => df.columns.addNewString(it));
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
            outputCols: outputFieldsInput.value!.map((col) => {return {name: col.name, type: col.type};}),
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
            outputCols: outputFieldsInput.value!.map((col) => {return {name: col.name, type: col.type};}),
            aggregations: aggregations,
          };
          saveSearch(entityLevel, savedSearch);
        }
      }, 'Save current search query');


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
      filtersDiv.append(ui.div(runSearchButton, { style: { paddingLeft: '4px', paddingTop: '10px' } }));
      createFiltersIcon(tv, filtersDiv);
      if (initialQuery) {
        loadSearchQuery(initialQuery, queryBuilder, outputFieldsInput, aggregationsContainer, menuFieldsForAggr,
          aggregations, validationErrorSubj);
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

export async function runSearch(tv: DG.TableView, search: MolTrackSearch, entityLevel: Scope,
  entityType: MolTrackEntityType, endpoint: string, changePath: boolean) {
  try {
    ui.setUpdateIndicator(tv.grid.root, true, 'Searching...');
    if (changePath)
      tv.path = createPathFromArr([SEARCH_NODE, entityLevel, JSON.stringify(search)]);
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
      resultDf?.columns.names().forEach((colName: string) => {
        const friendlyName =
          resultDf.col(colName)!.name.replace(`${entityLevel.toLowerCase()}.`, '').replace(`details.`, '');
        resultDf.col(colName)!.setTag('friendlyName', friendlyName);
      });
      reorderSearchResColummns(resultDf);
      //save history in case of successful search
      saveSearchHistory(`${_package.name}|${entityLevel}`, JSON.stringify(search));
    }
    tv.dataFrame = resultDf ?? DG.DataFrame.create();
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
  tv.setRibbonPanels([[filtersButton]]);
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
      if (EXCLUDE_SEARCH_FIELDS.includes(prop.name))
        continue;
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
        }:
        {
          name: prop.name,
          type: prop.value_type,
          semType: prop.semantic_type?.name,
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
    filter: cond.conditions.length === 1 && isSimpleCondition(cond.conditions[0]) ?
      convertSimpleCondition(cond.conditions[0], level, type):
      convertComplexCondition(cond, level, type),
  };
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
  if (cond.operator === 'is_similar' && cond.value && typeof cond.value === 'object' && 'molecule' in cond.value) {
    return {
      field: `${level}${isDynamicProp ? '.details': ''}.${cond.field}`,
      operator: Operators.IS_SIMILAR,
      value: cond.value.molecule,
      threshold: cond.value.threshold || null,
    };
  }

  // Handle regular simple conditions - operators are similar, no mapping needed
  return {
    field: `${level}${isDynamicProp ? '.details': ''}.${cond.field}`,
    operator: cond.operator,
    value: cond.value instanceof dayjs ? convertDateToString(cond.value) : cond.value,
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
      fieldInput.root.onclick = () => menu.show({ element: fieldInput.root, y: fieldInput.root.offsetHeight });
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


export function createSearchNode(appNode:DG.TreeViewGroup, scope: string, initialQuery?: MolTrackSearch) {
  const formattedScope = scope
    .toLowerCase()
    .replace(/_/g, ' ')
    .replace(/\b\w/g, (c) => c.toUpperCase());

  appNode.getOrCreateGroup(SEARCH_NODE).item(formattedScope).onSelected.subscribe(async () => {
    createSearchView(formattedScope, scope, initialQuery);
  });
}

export async function createSearchView(viewName: string, scope: string, initialQuery?: MolTrackSearch,
  isSavedSearch?: boolean): Promise<DG.TableView> {
  if (openedSearchView?.name === viewName && openedSearchView.dockNode.parent) {
    grok.shell.v = openedSearchView;
    return grok.shell.tv;
  }
  const df = DG.DataFrame.create();
  openedSearchView?.close();
  openedSearchView = grok.shell.addTablePreview(df);
  openedSearchView.name = viewName;
  const initPath = [
    isSavedSearch ? SAVED_SEARCHES_NODE : SEARCH_NODE,
    isSavedSearch ? scope : viewName,
  ];
  if (isSavedSearch)
    initPath.push(viewName);
  else if (initialQuery)
    initPath.push(JSON.stringify(initialQuery));
  openedSearchView.path = createPathFromArr(initPath);
  //if there is no initial search query - retrieve all data to initially show in the view
  if (!initialQuery) {
    ui.setUpdateIndicator(openedSearchView.root, true, `Loading ${viewName.toLowerCase()}...`);
    try {
      const data: DG.DataFrame = await grok.functions.call('MolTrack:retrieveEntity', { scope });
      openedSearchView.dataFrame = data;
      openedSearchView.dataFrame.name = viewName;
    } finally {
      ui.setUpdateIndicator(openedSearchView.root, false);
    }
  }
  createSearchPanel(openedSearchView!, scope as Scope, initialQuery, !isSavedSearch);
  return openedSearchView;
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

export async function handleSearchURL(url: string): Promise<DG.TableView> {
  if (url.startsWith('/'))
    url = url.slice(1);
  const componentsArr = url.split('/');
  if (componentsArr.length && (componentsArr[0] === SEARCH_NODE || componentsArr[0] === SAVED_SEARCHES_NODE)) {
    const tv = grok.shell.addTablePreview(DG.DataFrame.create());
    ui.setUpdateIndicator(tv.root, true, 'Loading Revvity Signals...');
    let idx = 0;
    const nodesToExpand = [];
    nodesToExpand.push(componentsArr[0]);
    idx++;
    if (componentsArr.length > idx) {
      const scope = componentsArr[1];
      if (componentsArr[0] === SAVED_SEARCHES_NODE) {
        nodesToExpand.push(componentsArr[1]);
        idx++;
      } else {
        let initialSearch = undefined;
        //check for initial query in case of search url
        if (componentsArr.length > idx + 1)
          initialSearch = JSON.parse(componentsArr[2]);

        openMolTrackSearchNode(nodesToExpand);
        return await createSearchView(scope, scope.toLowerCase(), initialSearch, false);
      }
      //in case of saved search
      if (componentsArr.length > idx) {
        const savedSearchName = componentsArr[idx];
        const savedSearches = getSavedSearches(scope as Scope);
        if (savedSearches[savedSearchName]) {
          openMolTrackSearchNode(nodesToExpand);
          return await createSearchView(savedSearchName, scope.toLowerCase(),
            JSON.parse(savedSearches[savedSearchName]), true);
        } else {
          grok.shell.error(`Search ${savedSearchName} not found for ${scope}`);
          return await createSearchView(scope, scope.toLowerCase(), undefined, false);
        }
      }
    }
  }
  throw Error(`incorrect search path`);
}

//register operators for molecule semType, applicable for MolTrack
ConditionRegistry.getInstance()
  .registerSemTypeOperators(DG.SEMTYPE.MOLECULE, [Operators.IS_SIMILAR, 'has substructure', 'is substructure of']);
ConditionRegistry.getInstance()
  .registerEditor(DG.TYPE.STRING, DG.SEMTYPE.MOLECULE, 'has substructure', MoleculeConditionEditor);
ConditionRegistry.getInstance()
  .registerEditor(DG.TYPE.STRING, DG.SEMTYPE.MOLECULE, 'is substructure of', MoleculeConditionEditor);
