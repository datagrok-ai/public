/* Do not change these import lines to match external modules in webpack configuration */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {
  QueryBuilder,
  QueryBuilderLayout,
  ComplexCondition,
  Operators,
} from '@datagrok-libraries/utils/src/query-builder/query-builder';
import { MolTrackDockerService } from './moltrack-docker-service';
import {
  MolTrackSearchQuery,
  MolTrackFilter,
  MolTrackComplexCondition,
  searchTypeMapping,
} from './types';
import { Scope } from './constants';

export let molTrackEntitiesProps: {[key: string]: DG.Property[]} | null = null;

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
    const filterFields: DG.Property[] = await getPropertiesForEntity(entityLevel);
    queryBuilder = new QueryBuilder(filterFields, undefined, QueryBuilderLayout.Narrow);
    const outputFieldsInput = ui.input.list('Output', {});
    outputFieldsInput.root.style.paddingLeft = '15px';
    const runSearchButton = ui.bigButton('Search', async () => {
      ui.setUpdateIndicator(tv.grid.root, true, 'Searching...');
      const molTrackQuery = convertQueryBuilderConditionToMolTrackQuery(queryBuilder!.condition,
        entityLevel, outputFieldsInput.value!.map((it) => `${entityLevel}.${it}`));
      const endpointLevel = searchTypeMapping.filter((it) => it.level === entityLevel);
      if (!endpointLevel.length) {
        ui.setUpdateIndicator(tv.grid.root, false);
        throw new Error(`No search endpoint for ${entityLevel}`);
      }
      const result = await MolTrackDockerService.search(molTrackQuery, endpointLevel[0].searchEndpoint);
      const resultDf = DG.DataFrame.fromObjects(result.data);
      tv.dataFrame = resultDf ?? DG.DataFrame.create();
      ui.setUpdateIndicator(tv.grid.root, false);
    });
    ui.setUpdateIndicator(filtersDiv, false);
    filtersDiv.append(queryBuilder.root);
    filtersDiv.append(outputFieldsInput.root);
    filtersDiv.append(ui.div(runSearchButton, { style: { paddingLeft: '4px' } }));
  };

  const filtersButton = ui.button('Add filters', () => {
    tv.dockManager.dock(filtersDiv, 'left', null, 'Filters', 0.2);
    if (!queryBuilder)
      initializeQueryBuilder(entityLevel);
  });
  tv.setRibbonPanels([[filtersButton]]);
}

export async function getPropertiesForEntity(entityLevel: Scope): Promise<DG.Property[]> {
  const entityTypeArr = searchTypeMapping.filter((it) => it.level === entityLevel);
  if (!entityTypeArr.length)
    return [];
  const entityType = entityTypeArr[0].propEntityType;
  if (!molTrackEntitiesProps) {
    molTrackEntitiesProps = {};
    const props = await MolTrackDockerService.fetchSchema();
    for (const prop of props) {
      const propOptions: {[key: string]: any} = {
        name: prop.name,
        type: prop.value_type,
      };
      const dgProp = DG.Property.fromOptions(propOptions);
      if (!molTrackEntitiesProps[prop.entity_type])
        molTrackEntitiesProps[prop.entity_type] = [dgProp];
      else
        molTrackEntitiesProps[prop.entity_type].push(dgProp);
    }
  }
  return molTrackEntitiesProps[entityType];
}

export function convertQueryBuilderConditionToMolTrackQuery(
  cond: ComplexCondition,
  level: Scope,
  output: string[],
): MolTrackSearchQuery {
  // If query builder condition contains only one condition in array,
  // create MolTrackSearchQuery where filter is a MolTrackSimpleCondition
  if (cond.conditions.length === 1 && isSimpleCondition(cond.conditions[0])) {
    const filter = convertSimpleCondition(cond.conditions[0], level, true);

    return {
      level: level,
      output: output,
      filter: filter,
    };
  }

  const filter = convertComplexCondition(cond, level);

  return {
    level: level,
    output: output,
    filter: filter,
  };
}

function convertComplexCondition(cond: ComplexCondition, level: Scope): MolTrackComplexCondition {
  if (cond.conditions.length === 0) {
    return {//@ts-ignore
      operator: Operators.Logical.and.toUpperCase(),
      conditions: [],
    };
  }

  const convertedConditions: MolTrackFilter[] = cond.conditions.map((condition) => {
    if (isSimpleCondition(condition))
      return convertSimpleCondition(condition, level, true);
    else
      return convertComplexCondition(condition as ComplexCondition, level);
  });

  return { //@ts-ignore
    operator: cond.logicalOperator.toUpperCase(),
    conditions: convertedConditions,
  };
}

function convertSimpleCondition(cond: any, level: Scope, isDynamicProp?: boolean): MolTrackFilter {
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
    value: cond.value,
    threshold: cond.threshold || null,
  };
}

function isSimpleCondition(condition: any): boolean {
  return 'field' in condition && 'operator' in condition && 'value' in condition;
}

