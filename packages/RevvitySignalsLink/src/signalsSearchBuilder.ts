import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { SignalsSearchQuery, QueryOperator, SignalsSearchParams } from './signalsSearchQuery';
import * as grok from 'datagrok-api/grok';

export enum OPERATORS {
  MATCH = '$match',
  IN = '$in',
  GT = '$gt',
  LT = '$lt',
  GTE = '$gte',
  LTE = '$lte',
  RANGE = '$range',
  PREFIX = '$prefix',
  EXISTS = '$exists',
  SIMPLE = '$simple',
  INTERSECT = '$intersect',
  CHEMSEARCH = '$chemsearch',
  AND = '$and',
  OR = '$or',
  NOT = '$not',
}

export const OPERATOR_DISPLAY_NAMES: Record<string, OPERATORS> = {
  'MATCH': OPERATORS.MATCH,
  'IN': OPERATORS.IN,
  'GREATER THAN': OPERATORS.GT,
  'LESS THAN': OPERATORS.LT,
  'GREATER THAN OR EQUAL': OPERATORS.GTE,
  'LESS THAN OR EQUAL': OPERATORS.LTE,
  'RANGE': OPERATORS.RANGE,
  'PREFIX': OPERATORS.PREFIX,
  'EXISTS': OPERATORS.EXISTS,
  'SIMPLE': OPERATORS.SIMPLE,
  'INTERSECT': OPERATORS.INTERSECT,
  'CHEMICAL SEARCH': OPERATORS.CHEMSEARCH,
  'AND': OPERATORS.AND,
  'OR': OPERATORS.OR,
  'NOT': OPERATORS.NOT,
};

export enum SourceType {
  SN = 'SN',
  CONNECTED = 'CONNECTED',
  IVT = 'IVT',
  CHEMICALS = 'CHEMICALS',
}

export enum IncludeType {
  CREATED_BY = 'createdBy',
  EDITED_BY = 'editedBy',
  OWNER = 'owner',
}


// Helper: Create a default operator object
export function createDefaultOperator(): QueryOperator {
  return { $match: { field: '', value: '' } } as QueryOperator;
}

// Recursive UI builder for a QueryOperator
export function buildOperatorUI(
  operator: QueryOperator,
  onValueChange: () => void,
  onStructureChange: () => void
): HTMLElement {
  // Detect operator type
  const opKey = Object.keys(operator)[0] as keyof QueryOperator;
  const opValue = (operator as any)[opKey];

  // If the current operator is logical, show group UI, else show field-based UI
  if (opKey === OPERATORS.AND || opKey === OPERATORS.OR || opKey === OPERATORS.NOT) {
    const group = opValue as QueryOperator[];
    // Show a dropdown for all operators
    const operatorInput = ui.input.choice('Operator', {
      items: Object.keys(OPERATOR_DISPLAY_NAMES), value: Object.keys(OPERATOR_DISPLAY_NAMES).find(key => OPERATOR_DISPLAY_NAMES[key] === opKey),
      onValueChanged: () => {
        const newOp = OPERATOR_DISPLAY_NAMES[operatorInput.value!];
        let newOperator: QueryOperator;
        if (newOp === OPERATORS.AND || newOp === OPERATORS.OR || newOp === OPERATORS.NOT) {
          newOperator = { [String(newOp)]: group } as unknown as QueryOperator;
        } else {
          newOperator = { [String(newOp)]: { field: '', value: '' } } as unknown as QueryOperator;
        }
        Object.keys(operator).forEach((k) => delete (operator as any)[k]);
        Object.assign(operator, newOperator);
        onStructureChange();
      }
    });
    return ui.divV([
      operatorInput.root,
      ...group.map((childOp, idx) =>
        ui.divH([
          buildOperatorUI(childOp, onValueChange, onStructureChange),
          ui.icons.delete(() => {
            group.splice(idx, 1);
            onStructureChange();
          }, 'Remove condition'),
        ])
      ),
      ui.icons.add(() => {
        group.push(createDefaultOperator());
        onStructureChange();
      }, 'Add condition'),
    ], 'revvity-signals-search-logical-operator');
  }

  // Operator selector
  const operatorInput = ui.input.choice('Operator', {
    items: Object.keys(OPERATOR_DISPLAY_NAMES), value: Object.keys(OPERATOR_DISPLAY_NAMES).find(key => OPERATOR_DISPLAY_NAMES[key] === opKey),
    onValueChanged: () => {
      // Replace operator type
      const newOp = OPERATOR_DISPLAY_NAMES[operatorInput.value!];
      let newOperator: QueryOperator = createDefaultOperator();
      if (newOp === OPERATORS.AND)
        newOperator = { [OPERATORS.AND]: [createDefaultOperator()] };
      else if (newOp === OPERATORS.OR)
        newOperator = { [OPERATORS.OR]: [createDefaultOperator()] };
      else if (newOp === OPERATORS.NOT)
        newOperator = { [OPERATORS.NOT]: [createDefaultOperator()] };
      else
        newOperator = { [String(newOp)]: { field: 'type', value: '' } } as unknown as QueryOperator;
      Object.keys(operator).forEach((k) => delete (operator as any)[k]);
      Object.assign(operator, newOperator);
      onStructureChange();
    }
  });

  const choiceInput = (fieldName: string, choices: string[], required?: boolean, updateOpValueOnCreate?: boolean) => {
    const value = opValue[fieldName] ?? '';
    const input = ui.input.choice(fieldName.charAt(0).toUpperCase() + fieldName.slice(1), {
      items: choices,
      value: value && choices.includes(value) ? value : required ? choices[0] : null,
      nullable: !(!!required),
      onValueChanged: () => {
        if (input.value || required)
          opValue[fieldName] = input.value;
        else
          delete opValue[fieldName];
        onValueChange();
      }
    });
    if (updateOpValueOnCreate)
      input.fireChanged();
    return input;
  }


  const multiValuesInput = () => {
    const values = opValue.values || [];
    const valuesInput = ui.input.string('Values (comma separated)', {
      value: values.join(','),
      nullable: false,
      onValueChanged: () => {
        opValue.values = valuesInput.value.split(',').map((v: string) => v.trim());
        onValueChange();
      }
    });
    return valuesInput;
  }

  const stringInput = (fieldName: string, required: boolean) => {
    const valueStringInput = ui.input.string(fieldName.charAt(0).toUpperCase() + fieldName.slice(1), {
      value: opValue[fieldName] ?? '',
      nullable: !(!!required),
      onValueChanged: () => {
        opValue[fieldName] = valueStringInput.value;
        onValueChange();
      }
    });
    return valueStringInput;
  }

  const boolInInput = () => {
    const valueBoolInput = ui.input.bool('Search in tags', {
      value: opValue['in'] ?? false,
      nullable: true,
      onValueChanged: () => {
        if (valueBoolInput.value)
          opValue['in'] = 'tags';
        else
          delete opValue['in'];
        onValueChange();
      }
    });
    return valueBoolInput;
  }

  const floatValueInput = (fieldName: string) => {
    const valueFloatInput = ui.input.float(fieldName.charAt(0).toUpperCase() + fieldName.slice(1), {
      value: typeof opValue[fieldName] === 'number' ? opValue[fieldName] : undefined,
      nullable: false,
      onValueChanged: () => {
        opValue[fieldName] = typeof valueFloatInput.value === 'number' ? valueFloatInput.value : undefined;
        onValueChange();
      }
    });
    return valueFloatInput;
  }

  const minimumInput = () => {
    const minimum = typeof opValue.minimum === 'number' ? opValue.minimum : undefined;
    const minimumInput = ui.input.int('Minimum', {
      value: minimum,
      nullable: false,
      onValueChanged: () => {
        opValue.minimum = typeof minimumInput.value === 'number' ? minimumInput.value : undefined;
        onValueChange();
      }
    });
    return minimumInput;
  }

  const moleculeInput = () => {
    const molecule = opValue.molecule || '';
    const moleculeInput = ui.input.molecule('Molecule', {
      value: molecule,
      nullable: false,
      onValueChanged: () => {
        opValue.molecule = moleculeInput.value; onValueChange();
      }
    });
    return moleculeInput;
  }


  let inputs: DG.InputBase[] = [];
  if (opKey === OPERATORS.MATCH || opKey === OPERATORS.PREFIX)
    inputs = [stringInput('field', true), boolInInput(), stringInput('value', true),
      choiceInput('mode', opKey === OPERATORS.MATCH ? ['keyword'] : ['keyword'], opKey === OPERATORS.PREFIX, opKey === OPERATORS.PREFIX)];
  else if (opKey === OPERATORS.IN)
    inputs = [stringInput('field', true), boolInInput(), multiValuesInput()];
  else if (opKey === OPERATORS.INTERSECT)
    inputs = [stringInput('field', true), boolInInput(), multiValuesInput(), minimumInput()];
  else if (opKey === OPERATORS.GT || opKey === OPERATORS.LT || opKey === OPERATORS.GTE || opKey === OPERATORS.LTE)
    inputs = [stringInput('field', true), boolInInput(), floatValueInput('value')];
  else if (opKey === OPERATORS.RANGE)
    inputs = [stringInput('field', true), boolInInput(), floatValueInput('from'), floatValueInput('to')];
  else if (opKey === OPERATORS.EXISTS)
    inputs = [stringInput('field', true), boolInInput()];
  else if (opKey === OPERATORS.SIMPLE)
    inputs = [stringInput('query', true)];
  else if (opKey === OPERATORS.CHEMSEARCH) {
    inputs = [moleculeInput(), stringInput('mime', true), stringInput('options', false)]
  } else {
    //@ts-ignore
    inputs = [ui.label('Unsupported operator')];
  }

  return ui.form([
    operatorInput,
    ...inputs,
  ], 'revvity-signals-search-query-operator');
}


// Main UI builder function
export function signalsSearchBuilderUI(): HTMLElement {
  const requestParams: SignalsSearchParams = {};
  let rootQuery: SignalsSearchQuery = {
    query: createDefaultOperator(),
    options: { sort: { modifiedAt: 'desc' } },
  };

  // Inputs for request parameters
  const pageOffsetInput = ui.input.int('Page Offset', { value: requestParams['page[offset]'] ?? 0, 
    onValueChanged: () => { requestParams['page[offset]'] = pageOffsetInput.value ?? 0; } });

  const pageLimitInput = ui.input.int('Page Limit', { value: requestParams['page[limit]'] ?? 20, 
    onValueChanged: () => { requestParams['page[limit]'] = pageLimitInput.value ?? 20; } });

  const includeInput = ui.input.multiChoice('Include', {
    items: Object.values(IncludeType),
    value:  requestParams.include ?? [IncludeType.OWNER],
    onValueChanged: () => {requestParams.include = includeInput.value ?? undefined;}});

  const sortInput = ui.input.string('Sort', { value: requestParams.sort ?? '',
    onValueChanged: () => { requestParams.sort = sortInput.value; } });

  const sourceInput = ui.input.choice('Source', {
    items: Object.values(SourceType),
    value: requestParams.source ?? SourceType.SN,
    onValueChanged: () => { requestParams.source = sourceInput.value ?? undefined; }
  });

  const fieldsEntityInput = ui.input.string('Fields[entity]', { value: requestParams['fields[entity]'] ?? '', 
    onValueChanged: () => { requestParams['fields[entity]'] = fieldsEntityInput.value ?? undefined; } });

  const stopAfterItemsInput = ui.input.int('Stop After Items', { value: requestParams.stopAfterItems ?? undefined, 
    onValueChanged: () => { requestParams.stopAfterItems = stopAfterItemsInput.value ?? undefined; } });

  // Use a plain textarea for preview
  const jsonPreview = document.createElement('textarea');
  jsonPreview.readOnly = true;
  jsonPreview.classList.add('revvity-signals-search-query-preview');

  function updatePreview() {
    jsonPreview.value = JSON.stringify(rootQuery, null, 2);
  }

  function rebuildUI() {
    updatePreview();
    builderDiv.replaceWith(builderDiv = buildOperatorUI(rootQuery.query, updatePreview, rebuildUI));
  }

  let builderDiv = buildOperatorUI(rootQuery.query, updatePreview, rebuildUI);

  const updateSorting = () => {
    rootQuery.options = rootQuery.options || {};
    rootQuery.options.sort = { [String(sortFieldInput.value)]: sortOrderInput.value as 'asc' | 'desc' };
    updatePreview();
  }

  // Sort options (now free text for field)
  const sortFieldInput = ui.input.string('Sort by', {
    onValueChanged: () => {
      updateSorting();
    }
  });
  const sortOrderInput = ui.input.choice('Order', {
    items: ['desc', 'asc'], value: 'desc',
    onValueChanged: () => {
      updateSorting();
    }
  });

  const submitBtn = ui.button('Submit Query', async () => {
    ui.empty(reGridDiv);
    const df = await grok.functions.call('RevvitySignalsLink:searchEntities', {
      query: JSON.stringify(rootQuery),
      params: JSON.stringify(requestParams)
    });
    //const df = await searchEntities(JSON.stringify(rootQuery), JSON.stringify(requestParams));
    const grid = df.plot.grid();
    grid.root.classList.add('revvity-signals-search-res-grid');
    reGridDiv.append(grid.root);
  });

  const acc = ui.accordion('Query builder');
  acc.addPane('Condition', () =>
    ui.divV([
      builderDiv,
      ui.divH([sortFieldInput.root, sortOrderInput.root])
    ], {style: {width: '100%'}}), true
  );
  acc.addPane('Parameters', () =>
    ui.form([
      pageOffsetInput,
      pageLimitInput,
      includeInput,
      sortInput,
      sourceInput,
      fieldsEntityInput,
      stopAfterItemsInput,
    ]),
  );
  acc.addPane('Query preview', () => ui.div([ui.label('Preview JSON'), jsonPreview], {style: {width: '100%'}}), true);


  const builderWithSubmitDiv = ui.divV([
    acc.root,
    submitBtn,
  ],'revvity-signals-search');

  const reGridDiv = ui.div();

  const builderWithResults = ui.splitV([
    builderWithSubmitDiv,
    reGridDiv
  ], {style: {width: '100%'}}, true)

  return builderWithResults;
} 