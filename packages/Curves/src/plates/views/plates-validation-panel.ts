/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Plate} from '../../plate/plate';
import {PlateTemplate} from '../plates-crud';

type ReconciliationDragObject = {
  sourceField: string,
  sourceType: DG.TYPE
};

/**
 * Renders a validation table comparing the plate's data schema to the template's requirements.
 */
export function renderValidationResults(
  tableElement: HTMLElement,
  plate: Plate,
  template: PlateTemplate,
  onMap: (source: string, target: string) => void
): { element: HTMLElement, conflictCount: number } {
  const templateWellProps = template.wellProperties
    .filter((p) => p && p.name)
    .map((p) => ({name: p.name!.toLowerCase(), type: p.type as DG.TYPE, originalName: p.name!}));
  const plateColumns = plate.data.columns.toList().map((c) => ({name: c.name.toLowerCase(), type: c.type as DG.TYPE, originalName: c.name}));
  const validationData: any[] = [];
  let conflictCount = 0;

  const mappedPlateCols = new Set<string>();
  for (const requiredProp of templateWellProps) {
    const foundCol = plateColumns.find((p) => p.name === requiredProp.name);
    if (foundCol)
      mappedPlateCols.add(foundCol.name);
  }
  const unmappedPlateCols = plateColumns.filter((p) => !mappedPlateCols.has(p.name));

  for (const requiredProp of templateWellProps) {
    const foundCol = plateColumns.find((p) => p.name === requiredProp.name);
    if (foundCol) {
      const typeMatch = foundCol.type === requiredProp.type;
      if (!typeMatch) conflictCount++;
      validationData.push({
        templateField: requiredProp.originalName, templateType: requiredProp.type,
        plateField: foundCol.originalName, status: typeMatch ? 'Match' : 'Type Mismatch', isConflict: !typeMatch,
      });
    } else {
      conflictCount++;
      validationData.push({
        templateField: requiredProp.originalName, templateType: requiredProp.type, plateField: '',
        status: 'Not Mapped', isConflict: true,
      });
    }
  }

  ui.empty(tableElement);
  if (validationData.length === 0) {
    tableElement.appendChild(ui.divText('No well properties defined in template.'));
    return {element: tableElement, conflictCount: 0};
  }

  tableElement.className = 'plate-validation-table';

  const header = ui.divH([
    ui.divText('Property', {style: {fontWeight: 'bold'}}),
    ui.divText('Column', {style: {fontWeight: 'bold'}}),
  ], {style: {padding: '8px 12px', borderBottom: '1px solid var(--grey-2)', display: 'grid', gridTemplateColumns: '1fr 2fr'}});
  tableElement.appendChild(header);

  validationData.forEach((item) => {
    const propertyCell = ui.divV([
      ui.divText(item.templateField),
      ui.divText(item.templateType, {style: {color: 'var(--grey-5)', fontSize: '11px'}})
    ]);
    ui.tooltip.bind(propertyCell, () => `Expected type: ${item.templateType}`);
    propertyCell.style.cursor = 'grab';
    ui.makeDraggable(propertyCell, {
      getDragObject: (): ReconciliationDragObject => ({sourceField: item.templateField, sourceType: item.templateType}),
      getDragCaption: () => `Map: ${item.templateField}`,
    });

    let columnCell: HTMLElement;
    if (item.status === 'Not Mapped') {
      const choices = ['', ...unmappedPlateCols.map((p) => p.originalName)];
      const choiceInput = ui.input.choice('', {items: choices, value: ''});
      choiceInput.onInput.subscribe(() => {
        if (choiceInput.value)
          onMap(item.templateField, choiceInput.value);
      });
      ui.tooltip.bind(choiceInput.root, `Property "${item.templateField}" is not mapped.`);
      columnCell = choiceInput.root;
    } else {
      const statusIcon = item.status === 'Match' ? ui.iconFA('check', null, 'Match') : ui.iconFA('exchange-alt', null, 'Type Mismatch');
      statusIcon.style.color = item.status === 'Match' ? 'var(--green-3)' : 'var(--orange-3)';
      columnCell = ui.divH([ui.divText(item.plateField), statusIcon], {style: {gap: '8px', alignItems: 'center'}});
    }

    ui.makeDroppable(columnCell, {
      acceptDrop: (dragObject: ReconciliationDragObject) => dragObject && typeof dragObject.sourceField === 'string',
      doDrop: (dragObject: ReconciliationDragObject) => onMap(dragObject.sourceField, item.plateField),
    });

    const row = ui.divH([propertyCell, columnCell], {style: {padding: '4px 12px', display: 'grid', gridTemplateColumns: '1fr 2fr'}});
    row.classList.add('d4-table-row-hover');
    Array.from(row.children).forEach((c) => (c as HTMLElement).style.alignSelf = 'center');

    tableElement.appendChild(row);
  });

  return {element: tableElement, conflictCount};
}
