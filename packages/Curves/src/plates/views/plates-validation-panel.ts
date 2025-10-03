/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Plate} from '../../plate/plate';
import {PlateTemplate} from '../plates-crud';

const REQUIRED_ANALYSIS_FIELDS = [
  {name: 'Activity', required: true},
  {name: 'Concentration', required: true},
  {name: 'SampleID', required: true},
];

const MOCKED_REQUIRED_TEMPLATE_FIELDS = ['Target', 'Assay Format'];

function createDynamicMappingRow(
  sourceColumns: string[],
  onMap: (source: string, target: string) => void,
  onCancel: () => void
): HTMLElement {
  let propName: string | null = null;
  let sourceCol: string | null = null;

  const tryApplyMapping = () => {
    if (propName && sourceCol)
      onMap(sourceCol, propName);
  };

  const propInput = ui.input.string('', {
    placeholder: 'Property name...',
    onValueChanged: (value) => {
      propName = value;
      tryApplyMapping();
    }
  });

  const colChoice = ui.input.choice('', {
    value: null,
    items: [null, ...sourceColumns],
    nullable: true,
    onValueChanged: (value: string | null) => {
      sourceCol = value;
      tryApplyMapping();
    }
  });

  const cancelBtn = ui.iconFA('times', onCancel, 'Cancel');
  cancelBtn.classList.add('mapping-cancel-icon');

  const rightCell = ui.divH([colChoice.root, cancelBtn], 'dynamic-row-right-cell');
  const newRow = ui.divH([propInput.root, rightCell], 'plate-validation-table-row');

  propInput.input.focus();
  return newRow;
}


export function renderValidationResults(
  tableElement: HTMLElement,
  plate: Plate,
  template: PlateTemplate,
  onMap: (sourceColumn: string, targetProperty: string) => void,
  reconciliationMap: Map<string, string>,
  onUndo: (targetProperty: string) => void
): { element: HTMLElement, conflictCount: number } {
  const sourceColumns = plate.data.columns.names();

  const templateProps = template.wellProperties
    .filter((p) => p && p.name)
    .map((p) => ({
      name: p.name!,
      isRequired: MOCKED_REQUIRED_TEMPLATE_FIELDS.includes(p.name!),
    }));

  const analysisProps = REQUIRED_ANALYSIS_FIELDS.map((p) => ({
    name: p.name,
    isRequired: p.required,
  }));

  const allPropsMap = new Map<string, { name: string, isRequired: boolean }>();
  [...templateProps, ...analysisProps].forEach((p) => {
    const existing = allPropsMap.get(p.name);
    if (existing)
      existing.isRequired = existing.isRequired || p.isRequired;
    else
      allPropsMap.set(p.name, {...p});
  });

  reconciliationMap.forEach((_, targetCol) => {
    if (!allPropsMap.has(targetCol))
      allPropsMap.set(targetCol, {name: targetCol, isRequired: false});
  });

  const allTargetProps = Array.from(allPropsMap.values()).sort((a, b) => {
    if (a.isRequired !== b.isRequired) return a.isRequired ? -1 : 1;
    return a.name.localeCompare(b.name);
  });

  ui.empty(tableElement);
  if (sourceColumns.length === 0) {
    tableElement.appendChild(ui.divText('Import a CSV file to map columns.', 'info-message'));
    return {element: tableElement, conflictCount: 0};
  }

  tableElement.className = 'plate-validation-table';

  let conflictCount = 0;

  allTargetProps.forEach((prop) => {
    const mappedSource = reconciliationMap.get(prop.name);

    const propNameEl = ui.divH([ui.span([prop.name])]);
    if (prop.isRequired) {
      const asterisk = ui.element('sup');
      asterisk.className = 'required-asterisk';
      asterisk.innerText = '*';
      propNameEl.appendChild(asterisk);
      if (!mappedSource)
        conflictCount++;
    }

    const choiceControl = ui.input.choice('', {
      value: mappedSource || null,
      items: [null, ...sourceColumns],
      nullable: true,
      onValueChanged: (v: string | null) => {
        if (v)
          onMap(v, prop.name);
        else if (mappedSource)
          onUndo(prop.name);
      },
    });

    const rightCell = ui.divH([choiceControl.root], 'mapping-input-container');
    if (mappedSource) {
      const undoIcon = ui.iconFA('times', () => onUndo(prop.name), 'Undo mapping');
      undoIcon.classList.add('mapping-undo-icon');
      rightCell.appendChild(undoIcon);
    }

    const row = ui.divH([propNameEl, rightCell], 'plate-validation-table-row');
    tableElement.appendChild(row);
  });

  const skeletonRow = ui.divH([], 'plate-validation-add-row');

  const addBtn = ui.button(ui.iconFA('plus'), () => {
    const newRow = createDynamicMappingRow(sourceColumns, onMap, () => newRow.remove());
    tableElement.insertBefore(newRow, skeletonRow);
  }, 'Add new property mapping');
  addBtn.classList.add('curves-icon-button', 'curves-add-button');
  skeletonRow.appendChild(addBtn);
  tableElement.appendChild(skeletonRow);

  return {element: tableElement, conflictCount};
}
