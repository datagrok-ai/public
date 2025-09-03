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
} from './constants';
import { funcs } from '../package-api';
import dayjs, { Dayjs } from 'dayjs';
import { _package } from '../package';
import { Subject } from 'rxjs';

export type MolTrackSearchFields = {
  direct?: DG.Property[],
  dynamic?: DG.Property[],
}

export let molTrackSearchFieldsArr: DG.Property[] | null = null;
let openedSearchView: DG.TableView | null = null;

export async function createSearchPanel(tv: DG.TableView, entityLevel: Scope) {
  const filtersDiv = ui.div([]);
  let queryBuilder: QueryBuilder | null = null;

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
    const validationErrorSubj = new Subject<boolean>();
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

      queryBuilder = new QueryBuilder(filterFields, undefined, QueryBuilderLayout.Narrow,
        `${_package.name}|${entityLevel}`);
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

      validationErrorSubj.subscribe((val: boolean) => runSearchButton.disabled = val);

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

      const runSearchButton = ui.bigButton('Search', async () => {
        await runSearch(tv, outputFieldsInput.value!, aggregations, entityLevel, queryBuilder!, entityType, endpoint);
      });
      filtersDiv.append(queryBuilder.root);
      filtersDiv.append(outputFieldsInput.root);
      filtersDiv.append(aggregationsDiv);
      filtersDiv.append(ui.div(runSearchButton, { style: { paddingLeft: '4px', paddingTop: '10px' } }));
      createFiltersIcon(tv, filtersDiv);
    } catch (e: any) {
      grok.shell.error(`Error loading filters: ${e?.message ?? e}`);
    } finally {
      ui.setUpdateIndicator(filtersDiv, false);
    }
  };

  initializeQueryBuilder(entityLevel);
}

export async function runSearch(tv: DG.TableView, outputFieldsList: DG.Column[],
  aggregations: MolTrackSearchAggregation[], entityLevel: Scope,
  queryBuilder: QueryBuilder, entityType: MolTrackEntityType, endpoint: string) {
  try {
    ui.setUpdateIndicator(tv.grid.root, true, 'Searching...');
    if (!outputFieldsList.length)
      throw new Error(`At least one input field should be selected`);
    const outputFields = outputFieldsList!
      .map((it) => `${entityLevel}.${isDynamicField(it.name, entityType) ? 'details.' : ''}${it.name}`);
    const molTrackQuery = convertQueryBuilderConditionToMolTrackQuery(queryBuilder!.condition,
      entityLevel, outputFields, aggregations, entityType);
    queryBuilder?.saveConditionToHistory();
    const result = await MolTrackDockerService.search(molTrackQuery, endpoint);
    const resultDf = DG.DataFrame.fromObjects(result.data);
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
  aggregations: MolTrackSearchAggregation[], validationErrorSubj: Subject<boolean>): HTMLElement {
  const row = ui.divH([], 'moltrack-search-aggr-row');

  const newAggr: MolTrackSearchAggregation = {field: '', operation: ''};
  aggregations.push(newAggr);

  const fieldValidationError = (errorMessage: string) => {
    validationErrorSubj.next(true);
    fieldTooltip = errorMessage;
    newAggr.field = '';
    Array.from(fieldInput.root.children).forEach((it) => it.classList.add('moltrack-invalid-aggr-field'));
  };

  const nonErrorTooltip = 'Click the field to see options';
  let fieldTooltip = nonErrorTooltip;
  const fieldInput = ui.input.string('', {
    value: '',
    nullable: false,
    onValueChanged: () => {
      const error = `Field ${fieldInput.value} not found. Click to select from list of available fields`;
      const splittedName = fieldInput.value.split('.');
      if (splittedName.length < 2) {
        fieldValidationError(error);
        return;
      }
      const prop = menuFieldsForAggr
        .filter((it) => it.options[MOLTRACK_ENTITY_LEVEL] === splittedName[0] && it.name === splittedName[1]);
      if (!prop.length) {
        fieldValidationError(error);
        return;
      }
      newAggr.field = `${splittedName[0]}.${isDynamicField(prop[0].name,
        prop[0].options[MOLTRACK_ENTITY_TYPE]) ? 'details.' : ''}${splittedName[1]}`;
      const fieldType = prop[0].type;
      updateAggregationOptions(aggregationInputDiv, fieldType, newAggr);
      Array.from(fieldInput.root.children).forEach((it) => it.classList.remove('moltrack-invalid-aggr-field'));
      const containInvalidRows = aggregations.filter((it) => it.field === '').length > 0;
      validationErrorSubj.next(containInvalidRows);
      fieldTooltip = nonErrorTooltip;
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
    const containInvalidRows = aggregations.filter((it) => it.field === '').length > 0;
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
  validationErrorSubj.next(true);

  return row;
}

function updateAggregationOptions(aggregationInputDiv: any, fieldType: any, aggr: MolTrackSearchAggregation) {
  const options = PROP_NUM_TYPES.includes(fieldType) ? NUMERIC_AGGREGATIONS : STRING_AGGREGATIONS;
  ui.empty(aggregationInputDiv);
  aggr.operation = options[0];
  const input = ui.input.choice('', {
    value: options[0],
    items: options,
    nullable: false,
    onValueChanged: () => {
      aggr.operation = input.value!;
    },
  });
  aggregationInputDiv.append(input.root);
}


export function createSearchNode(appNode:DG.TreeViewGroup, scope: string) {
  const formattedScope = scope
    .toLowerCase()
    .replace(/_/g, ' ')
    .replace(/\b\w/g, (c) => c.toUpperCase());

  appNode.getOrCreateGroup('Search').item(formattedScope).onSelected.subscribe(async () => {
    if (openedSearchView?.name === formattedScope)
      return;
    const df = DG.DataFrame.create();
    openedSearchView?.close();
    openedSearchView = grok.shell.addTablePreview(df);
    openedSearchView.name = formattedScope;

    createSearchPanel(openedSearchView!, scope as Scope);
  });
}

//register operators for molecule semType, applicable for MolTrack
ConditionRegistry.getInstance()
  .registerSemTypeOperators(DG.SEMTYPE.MOLECULE, [Operators.IS_SIMILAR, 'has substructure', 'is substructure of']);
ConditionRegistry.getInstance()
  .registerEditor(DG.TYPE.STRING, DG.SEMTYPE.MOLECULE, 'has substructure', MoleculeConditionEditor);
ConditionRegistry.getInstance()
  .registerEditor(DG.TYPE.STRING, DG.SEMTYPE.MOLECULE, 'is substructure of', MoleculeConditionEditor);
