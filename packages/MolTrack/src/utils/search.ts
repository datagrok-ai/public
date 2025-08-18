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
  MolTrackEntityType,
  MolTrackSearchQuery,
  MolTrackFilter,
  MolTrackComplexCondition,
} from './types';

export let molTrackEntitiesProps: {[key: string]: DG.Property[]} | null = null;

export async function createSearchPanel(tv: DG.TableView, entityType: MolTrackEntityType) {
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

  const initializeQueryBuilder = async (entityType: MolTrackEntityType) => {
    ui.setUpdateIndicator(filtersDiv, true, 'Loading filters...');
    const filterFields: DG.Property[] = await getPropertiesForEntity(entityType);
    queryBuilder = new QueryBuilder(filterFields, undefined, QueryBuilderLayout.Narrow);
    const runSearchButton = ui.bigButton('Search', async () => {
      ui.setUpdateIndicator(tv.grid.root, true, 'Searching...');
     // const molTrackQuery = convertQueryBuilderConditionToMolTrackQuery(queryBuilder!.condition);
      // const resultDf = await MolTrackDockerService.search(queryBuilder);
      // tv.dataFrame = resultDf;
      ui.setUpdateIndicator(tv.grid.root, false);
    });
    ui.setUpdateIndicator(filtersDiv, false);
    filtersDiv.append(queryBuilder.root);
    filtersDiv.append(ui.div(runSearchButton, { style: { paddingLeft: '4px' } }));
  };

  const filtersButton = ui.button('Add filters', () => {
    tv.dockManager.dock(filtersDiv, 'left', null, 'Filters', 0.2);
    if (!queryBuilder)
      initializeQueryBuilder(entityType);
  });
  tv.setRibbonPanels([[filtersButton]]);
}

export async function getPropertiesForEntity(entityType: MolTrackEntityType): Promise<DG.Property[]> {
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
  level: string,
  output: string[],
): MolTrackSearchQuery {
  const filter = convertComplexCondition(cond);

  return {
    level: level,
    output: output,
    filter: filter,
  };
}

function convertComplexCondition(cond: ComplexCondition): MolTrackComplexCondition {
  if (cond.conditions.length === 0) {
    return {
      operator: Operators.Logical.and,
      conditions: [],
    };
  }

  const convertedConditions: MolTrackFilter[] = cond.conditions.map((condition) => {
    if (isSimpleCondition(condition)) {
      return convertSimpleCondition(condition);
    } else {
      return convertComplexCondition(condition as ComplexCondition);
    }
  });

  return {
    operator: cond.logicalOperator,
    conditions: convertedConditions,
  };
}

function convertSimpleCondition(cond: any): MolTrackFilter {
  // Handle chemical similarity search special case
  if (cond.operator === 'is_similar' && cond.value && typeof cond.value === 'object' && 'molecule' in cond.value) {
    return {
      field: cond.field,
      operator: Operators.IS_SIMILAR,
      value: cond.value.molecule,
      threshold: cond.value.threshold || null,
    };
  }

  // Handle regular simple conditions - operators are similar, no mapping needed
  return {
    field: cond.field,
    operator: cond.operator,
    value: cond.value,
    threshold: cond.threshold || null,
  };
}

function isSimpleCondition(condition: any): boolean {
  return 'field' in condition && 'operator' in condition && 'value' in condition;
}

