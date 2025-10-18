/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Plate} from '../../plate/plate';
import {PlateTemplate} from '../plates-crud';


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
  cancelBtn.classList.add('assay-plates--mapping-cancel-icon');

  const rightCell = ui.divH([colChoice.root, cancelBtn], 'assay-plates--dynamic-row-right-cell');
  const newRow = ui.divH([propInput.root, rightCell], 'assay-plates--validation-table-row');

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

  const requiredPropIds = new Set(template.required_props.map((tuple) => tuple[0]));

  const templateProps = template.wellProperties
    .filter((p) => p && p.name)
    .map((p) => ({
      name: p.name!,
      isRequired: requiredPropIds.has(p.id!),
    }));


  const allPropsMap = new Map<string, { name: string, isRequired: boolean }>();
  [...templateProps].forEach((p) => {
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
    tableElement.appendChild(ui.divText('Import a CSV file to map columns.', 'assay-plates--info-message'));
    return {element: tableElement, conflictCount: 0};
  }

  tableElement.className = 'assay-plates--validation-table';

  let conflictCount = 0;

  allTargetProps.forEach((prop) => {
    const mappedSource = reconciliationMap.get(prop.name);

    const propNameEl = ui.divH([ui.span([prop.name])]);
    if (prop.isRequired) {
      const asterisk = ui.span([' *'], 'assay-plates--required-asterisk');
      asterisk.style.color = 'red';
      asterisk.style.marginLeft = '2px';
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

    const rightCell = ui.divH([choiceControl.root], 'assay-plates--mapping-input-container');
    if (mappedSource) {
      const undoIcon = ui.iconFA('times', () => onUndo(prop.name), 'Undo mapping');
      undoIcon.classList.add('assay-plates--mapping-undo-icon');
      rightCell.appendChild(undoIcon);
    }

    const row = ui.divH([propNameEl, rightCell], 'assay-plates--validation-table-row');
    tableElement.appendChild(row);
  });

  const skeletonRow = ui.divH([], 'assay-plates--validation-add-row');

  const addBtn = ui.button(ui.iconFA('plus'), () => {
    const newRow = createDynamicMappingRow(sourceColumns, onMap, () => newRow.remove());
    tableElement.insertBefore(newRow, skeletonRow);
  }, 'Add new property mapping');
  addBtn.classList.add('assay-plates--icon-button', 'assay-plates--add-button');
  skeletonRow.appendChild(addBtn);
  tableElement.appendChild(skeletonRow);

  return {element: tableElement, conflictCount};
}
