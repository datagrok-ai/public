import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { SignalsSearchQuery, QueryOperator } from './signalsSearchQuery';

const OPERATORS = [
  '$match', '$in', '$gt', '$lt', '$gte', '$lte', '$range', '$prefix', '$exists', '$simple', '$intersect', '$chemsearch', '$and', '$or', '$not',
];

// Helper: Create a default operator object
export function createDefaultOperator(): QueryOperator {
  return { $match: { field: 'type', value: '' } } as QueryOperator;
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

  // Logical operators
  if (opKey === '$and' || opKey === '$or' || opKey === '$not') {
    const group = opValue as QueryOperator[];
    return ui.divV([
      ui.label(opKey),
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
    ], { style: { border: '1px solid #ccc', margin: '4px', padding: '4px' } });
  }

  // Field-based operators
  let field = opValue.field || '';
  let value = opValue.value;
  let values = opValue.values || [];
  let from = typeof opValue.from === 'number' ? opValue.from : undefined;
  let to = typeof opValue.to === 'number' ? opValue.to : undefined;
  let minimum = typeof opValue.minimum === 'number' ? opValue.minimum : undefined;
  let molecule = opValue.molecule || '';
  let mime = opValue.mime || '';
  let options = opValue.options || '';
  let query = opValue.query || '';
  let inValue = opValue.in || '';

  function appendAsterisk(input: DG.InputBase) {
    const label = input.captionLabel;
    if (label && !label.querySelector('.required-field')) {
      label.classList.add('revvity-signals-search-required-field');
      const asterisk = document.createElement('span');
      asterisk.textContent = '*';
      asterisk.className = 'required-field';
      asterisk.style.color = 'red';
      asterisk.style.margin = '0 0 0 2px'; // Only left margin
      label.appendChild(asterisk);
    }
  }

  // Operator selector
  const operatorInput = ui.input.choice('Operator', {items: OPERATORS, value: opKey});
  operatorInput.onChanged.subscribe(() => {
    // Replace operator type
    const newOp = operatorInput.value;
    let newOperator: QueryOperator = createDefaultOperator();
    if (newOp === '$and')
      newOperator = { "$and": [createDefaultOperator()] };
    else if (newOp === '$or')
      newOperator = { "$or": [createDefaultOperator()] };
    else if (newOp === '$not')
      newOperator = { "$not": [createDefaultOperator()] };
    else
      newOperator = { [String(newOp)]: { field: 'type', value: '' } } as unknown as QueryOperator;
    Object.keys(operator).forEach((k) => delete (operator as any)[k]);
    Object.assign(operator, newOperator);
    onStructureChange();
  });

  // Field selector (now a free text input)
  const fieldInput = ui.input.string('Field', {value: field});
  appendAsterisk(fieldInput);
  fieldInput.onChanged.subscribe(() => {
    opValue.field = fieldInput.value;
    onValueChange();
  });

  // 'in' parameter input (for all applicable operators)
  const inInput = ui.input.string('In', {value: inValue}); // optional, no asterisk
  inInput.onChanged.subscribe(() => {
    opValue.in = inInput.value || undefined;
    onValueChange();
  });

  // Value input (adapt to operator)
  let valueInput: any | HTMLElement;
  if (opKey === '$match' || opKey === '$prefix') {
    const valueStringInput = ui.input.string('Value', {value: value});
    appendAsterisk(valueStringInput);
    valueStringInput.onChanged.subscribe(() => {
      opValue.value = valueStringInput.value;
      onValueChange();
    });
    valueInput = valueStringInput.root;
  } else if (opKey === '$in') {
    const valuesInput = ui.input.string('Values (comma separated)', {value: values.join(',')});
    appendAsterisk(valuesInput);
    valuesInput.onChanged.subscribe(() => {
      opValue.values = valuesInput.value.split(',').map((v: string) => v.trim());
      onValueChange();
    });
    valueInput = valuesInput.root;
  } else if (opKey === '$intersect') {
    const valuesInput = ui.input.string('Values (comma separated)', {value: values.join(',')});
    appendAsterisk(valuesInput);
    valuesInput.onChanged.subscribe(() => {
      opValue.values = valuesInput.value.split(',').map((v: string) => v.trim());
      onValueChange();
    });
    const minimumInput = ui.input.int('Minimum', {value: minimum});
    minimumInput.onChanged.subscribe(() => {
      opValue.minimum = typeof minimumInput.value === 'number' ? minimumInput.value : undefined;
      onValueChange();
    });
    valueInput = ui.divV([valuesInput.root, minimumInput.root]);
  } else if (opKey === '$gt' || opKey === '$lt' || opKey === '$gte' || opKey === '$lte') {
    const valueFloatInput = ui.input.float('Value', {value: typeof value === 'number' ? value : undefined});
    appendAsterisk(valueFloatInput);
    valueFloatInput.onChanged.subscribe(() => {
      opValue.value = typeof valueFloatInput.value === 'number' ? valueFloatInput.value : undefined;
      onValueChange();
    });
    valueInput = valueFloatInput.root;
  } else if (opKey === '$range') {
    const fromInput = ui.input.float('From', {value: from});
    const toInput = ui.input.float('To', {value: to});
    appendAsterisk(fromInput);
    appendAsterisk(toInput);
    fromInput.onChanged.subscribe(() => { opValue.from = typeof fromInput.value === 'number' ? fromInput.value : undefined; onValueChange(); });
    toInput.onChanged.subscribe(() => { opValue.to = typeof toInput.value === 'number' ? toInput.value : undefined; onValueChange(); });
    valueInput = ui.divV([fromInput.root, toInput.root]);
  } else if (opKey === '$exists') {
    valueInput = undefined;
  } else if (opKey === '$simple') {
    const simpleInput = ui.input.string('Query', {value: query});
    appendAsterisk(simpleInput);
    simpleInput.onChanged.subscribe(() => {
      opValue.query = simpleInput.value;
      onValueChange();
    });
    valueInput = simpleInput.root;
  } else if (opKey === '$chemsearch') {
    const moleculeInput = ui.input.molecule('Molecule', {value: molecule});
    const mimeInput = ui.input.string('Mime', {value: mime});
    appendAsterisk(moleculeInput);
    appendAsterisk(mimeInput);
    const optionsInput = ui.input.string('Options', {value: options});
    moleculeInput.onChanged.subscribe(() => { opValue.molecule = moleculeInput.value; onValueChange(); });
    mimeInput.onChanged.subscribe(() => { opValue.mime = mimeInput.value; onValueChange(); });
    optionsInput.onChanged.subscribe(() => { opValue.options = optionsInput.value; onValueChange(); });
    valueInput = ui.divV([moleculeInput.root, mimeInput.root, optionsInput.root]);
  } else {
    valueInput = ui.label('Unsupported operator');
  }

  // Compose the order: field, in, value
  const children = [fieldInput.root, inInput.root];
  if (valueInput) children.push(valueInput);

  return ui.divV([
    operatorInput.root,
    ...children,
  ], { style: { border: '1px solid #eee', margin: '4px', padding: '4px' } });
}

// Main UI builder function
export function signalsSearchBuilderUI(onSubmit: (query: SignalsSearchQuery) => void): HTMLElement {
  let rootQuery: SignalsSearchQuery = {
    query: createDefaultOperator(),
    options: { sort: { modifiedAt: 'desc' } },
  };

  // Use a plain textarea for preview
  const jsonPreview = document.createElement('textarea');
  jsonPreview.readOnly = true;
  jsonPreview.style.height = '200px';
  jsonPreview.style.width = '100%';

  function updatePreview() {
    jsonPreview.value = JSON.stringify(rootQuery, null, 2);
  }

  function rebuildUI() {
    updatePreview();
    builderDiv.replaceWith(builderDiv = buildOperatorUI(rootQuery.query, updatePreview, rebuildUI));
  }

  let builderDiv = buildOperatorUI(rootQuery.query, updatePreview, rebuildUI);

  // Sort options (now free text for field)
  const sortFieldInput = ui.input.string('Sort by');
  const sortOrderInput = ui.input.choice('Order', {items: ['desc', 'asc'], value: 'desc'});
  sortFieldInput.onChanged.subscribe(() => {
    rootQuery.options = rootQuery.options || {};
    rootQuery.options.sort = { [String(sortFieldInput.value)]: sortOrderInput.value as 'asc' | 'desc' };
    updatePreview();
  });
  sortOrderInput.onChanged.subscribe(() => {
    rootQuery.options = rootQuery.options || {};
    rootQuery.options.sort = { [String(sortFieldInput.value)]: sortOrderInput.value as 'asc' | 'desc' };
    updatePreview();
  });

  // Submit button
  const submitBtn = ui.button('Submit Query', () => {
    onSubmit(rootQuery);
  });

  return ui.divV([
    builderDiv,
    ui.divH([sortFieldInput.root, sortOrderInput.root]),
    ui.div([ui.label('Preview JSON'), jsonPreview]),
    submitBtn,
  ], 'revvity-signals-search');
} 