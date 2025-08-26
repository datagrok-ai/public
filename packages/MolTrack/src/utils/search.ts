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
  searchTypeMapping,
  MolTrackProperty,
  MolTrackEntityType,
  MolTrackSearchMapping,
} from './types';
import { EXCLUDE_SEARCH_FIELDS, EXCLUDE_SEARCH_OUTPUT_FIELDS, Scope } from './constants';
import { funcs } from '../package-api';
import dayjs, { Dayjs } from 'dayjs';
import { _package } from '../package';

export type MolTrackSearchFields = {
  direct?: DG.Property[],
  dynamic?: DG.Property[],
}

export let molTrackSearchFields: {[key: string]: MolTrackSearchFields} | null = null;

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
    ui.setUpdateIndicator(filtersDiv, true, 'Loading filters...');
    const entityTypeArr = searchTypeMapping.filter((it) => it.level === entityLevel);
    const entityType = entityTypeArr.length ? entityTypeArr[0].propEntityType : undefined;
    const endpointLevelArr = searchTypeMapping.filter((it) => it.level === entityLevel);
    const endpointLevel = endpointLevelArr.length ? endpointLevelArr[0] : undefined;
    const filterFields: DG.Property[] = await getSearchFiledsForEntity(entityType);
    queryBuilder = new QueryBuilder(filterFields, undefined, QueryBuilderLayout.Narrow,
      `${_package.name}|${entityLevel}`);
    const df = DG.DataFrame.create();
    const outputFields = filterFields.filter((it) => !EXCLUDE_SEARCH_OUTPUT_FIELDS.includes(it.name));
    outputFields.forEach((it) => df.columns
      .add(DG.Column.fromType((it.type as string === 'uuid' ? 'string' : it.type) as any, it.name)));
    const outputFieldsInput = ui.input.columns('Output', {table: df});
    outputFieldsInput.root.style.paddingLeft = '15px';
    const runSearchButton = ui.bigButton('Search', async () => {
      try {
        ui.setUpdateIndicator(tv.grid.root, true, 'Searching...');
        const outputFields = outputFieldsInput.value!
          .map((it) => `${entityLevel}.${isDynamicField(it.name, entityType) ? 'details.' : ''}${it.name}`);
        const molTrackQuery = convertQueryBuilderConditionToMolTrackQuery(queryBuilder!.condition,
          entityLevel, outputFields, entityType);
        if (!endpointLevel) {
          ui.setUpdateIndicator(tv.grid.root, false);
          throw new Error(`No search endpoint for ${entityLevel}`);
        }
        queryBuilder?.saveConditionToHistory();
        const result = await MolTrackDockerService.search(molTrackQuery, endpointLevel.searchEndpoint);
        const resultDf = DG.DataFrame.fromObjects(result.data);
        tv.dataFrame = resultDf ?? DG.DataFrame.create();
      } catch (e: any) {
        grok.shell.error(e?.message ?? e);
        tv.dataFrame = DG.DataFrame.create();
      } finally {
        ui.setUpdateIndicator(tv.grid.root, false);
      }
    });
    ui.setUpdateIndicator(filtersDiv, false);
    filtersDiv.append(queryBuilder.root);
    filtersDiv.append(outputFieldsInput.root);
    filtersDiv.append(ui.div(runSearchButton, { style: { paddingLeft: '4px', paddingTop: '10px' } }));
    tv.dockManager.dock(filtersDiv, 'left', null, 'Filters', 0.2);
    createFiltersIcon(tv, filtersDiv);
  };

  initializeQueryBuilder(entityLevel);
}

export function createFiltersIcon(tv: DG.TableView, filtersDiv: HTMLDivElement) {
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

export async function getSearchFiledsForEntity(entityType?: MolTrackEntityType): Promise<DG.Property[]> {
  if (!entityType)
    return [];
  if (!molTrackSearchFields) {
    molTrackSearchFields = {};
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
        const propOptions: { [key: string]: any } = {
          name: prop.name,
          type: prop.value_type,
        };
        const dgProp = DG.Property.fromOptions(propOptions);
        if (!molTrackSearchFields![prop.entity_type]) {
          molTrackSearchFields![prop.entity_type] =
            { direct: isStatic ? [dgProp] : [], dynamic: isStatic ? [] : [dgProp] };
        } else {
          isStatic ? molTrackSearchFields![prop.entity_type].direct!.push(dgProp) :
            molTrackSearchFields![prop.entity_type].dynamic!.push(dgProp);
        }
        if (dgProp.name === 'canonical_smiles')
          dgProp.semType = DG.SEMTYPE.MOLECULE;
        //for string values - implement suggestions
        const entityLevelArr = searchTypeMapping.filter((it) => it.propEntityType === prop.entity_type);
        const entityLevel = entityLevelArr.length ? entityLevelArr[0] : undefined;
        if (prop.value_type === 'string' && entityLevel && dgProp.semType !== DG.SEMTYPE.MOLECULE) {
          const fieldName = `${entityLevel.level}${isStatic ? '' : '.details'}.${prop.name}`;
          const getPropSuggestions = async (text: string): Promise<string[]> => {
            const query: MolTrackSearchQuery = {
              level: entityLevel.level,
              filter: {
                field: fieldName,
                operator: Operators.CONTAINS,
                value: text,
              },
              output: [fieldName],
            };
            const suggestions: string[] = [];
            try {
              const res = await MolTrackDockerService.search(query, entityLevel.searchEndpoint);
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
  }
  return (molTrackSearchFields[entityType].direct ?? []).concat(molTrackSearchFields[entityType].dynamic ?? []);
}

export function convertQueryBuilderConditionToMolTrackQuery(
  cond: ComplexCondition,
  level: Scope,
  output: string[],
  type?: MolTrackEntityType,
): MolTrackSearchQuery {
  // If query builder condition contains only one condition in array,
  // create MolTrackSearchQuery where filter is a MolTrackSimpleCondition
  if (cond.conditions.length === 1 && isSimpleCondition(cond.conditions[0])) {
    const filter = convertSimpleCondition(cond.conditions[0], level, type);

    return {
      level: level,
      output: output,
      filter: filter,
    };
  }

  const filter = convertComplexCondition(cond, level, type);

  return {
    level: level,
    output: output,
    filter: filter,
  };
}

function convertComplexCondition(cond: ComplexCondition, level: Scope,
  type?: MolTrackEntityType): MolTrackComplexCondition {
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

function isDynamicField(field: string, type?: MolTrackEntityType): boolean {
  let isDynamicProp = false;
  if (type) {
    const dynamicPropIdx = molTrackSearchFields![type].dynamic?.findIndex((it) => it.name === field);
    isDynamicProp = dynamicPropIdx !== -1;
  }
  return isDynamicProp;
}

function convertSimpleCondition(cond: any, level: Scope, type?: MolTrackEntityType): MolTrackFilter {
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


//register operators for molecule semType, applicable for MolTrack
ConditionRegistry.getInstance()
  .registerSemTypeOperators(DG.SEMTYPE.MOLECULE, [Operators.IS_SIMILAR, 'has substructure', 'is substructure of']);
ConditionRegistry.getInstance()
  .registerEditor(DG.TYPE.STRING, DG.SEMTYPE.MOLECULE, 'has substructure', MoleculeConditionEditor);
ConditionRegistry.getInstance()
  .registerEditor(DG.TYPE.STRING, DG.SEMTYPE.MOLECULE, 'is substructure of', MoleculeConditionEditor);
