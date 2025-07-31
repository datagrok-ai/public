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
 * Shows a resolution dialog for a specific validation item.
 * @param item The data for the clicked row.
 * @param extraColumns A list of unmapped columns from the plate.
 * @param onMap The callback to execute when a resolution is applied.
 */
function showResolutionDialog(item: any, extraColumns: string[], onMap: (source: string, target: string) => void) {
  const dialog = ui.dialog('Property Resolution');

  const resolutionInput = ui.div();
  if (item.status === 'Missing' || item.status === 'Type Mismatch') {
    if (extraColumns.length > 0) {
      const availableColumns = ui.input.choice('Map from', {items: extraColumns});
      ui.tooltip.bind(availableColumns.root, 'Select a column from your file to map to this template property');
      resolutionInput.appendChild(availableColumns.root);
      dialog.onOK(() => {
        if (availableColumns.value)
          onMap(item.templateField, availableColumns.value);
      });
    } else {
      resolutionInput.appendChild(ui.divText('No unmapped columns available to map from.'));
    }
  } else {
    resolutionInput.appendChild(ui.divText('No actions available for this status.'));
  }

  dialog.add(ui.divV([
    ui.h2('Resolve Property Mismatch'),
    ui.tableFromMap({
      'Template requires': `${item.templateField} (type: ${item.templateType})`,
      'Plate provides': item.plateField || 'N/A',
      'Status': item.status,
    }),
    ui.divV([
      ui.h3('Resolution UX', {style: {marginTop: '10px'}}),
      resolutionInput
    ], {style: {borderTop: '1px solid var(--grey-2)', marginTop: '10px'}})
  ]))
    .show();
}

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
  const extraColumns: string[] = [];
  let conflictCount = 0;

  const ICONS = {OK: '✔', ERROR: '✖', INFO: 'ⓘ'};

  for (const requiredProp of templateWellProps) {
    const foundCol = plateColumns.find((p) => p.name === requiredProp.name);
    if (foundCol) {
      const typeMatch = foundCol.type === requiredProp.type;
      if (!typeMatch) conflictCount++;
      validationData.push({
        templateField: requiredProp.originalName, templateType: requiredProp.type,
        plateField: foundCol.originalName, status: typeMatch ? 'Match' : 'Type Mismatch',
        icon: typeMatch ? ICONS.OK : ICONS.ERROR, color: typeMatch ? 'var(--green-3)' : 'var(--red-3)',
        isConflict: !typeMatch,
      });
    } else {
      conflictCount++;
      validationData.push({
        templateField: requiredProp.originalName, templateType: requiredProp.type, plateField: '',
        status: 'Missing', icon: ICONS.ERROR, color: 'var(--red-3)', isConflict: true,
      });
    }
  }
  for (const plateCol of plateColumns) {
    if (!templateWellProps.some((p) => p.name === plateCol.name)) {
      extraColumns.push(plateCol.originalName);
      validationData.push({
        templateField: '', templateType: null, plateField: plateCol.originalName,
        status: 'Extra', icon: ICONS.INFO, color: 'var(--blue-3)', isConflict: false,
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
    ui.divText(template.name, {style: {fontWeight: 'bold'}}),
    ui.divText(plate.barcode ?? 'Plate', {style: {fontWeight: 'bold'}}),
    ui.divText('Status', {style: {fontWeight: 'bold'}}),
  ], {style: {padding: '8px 12px', borderBottom: '1px solid var(--grey-2)', display: 'grid', gridTemplateColumns: '1fr 1fr 100px'}});
  tableElement.appendChild(header);

  validationData.forEach((item) => {
    const templateCell = ui.divText(item.templateField || '');
    if (item.templateType)
      ui.tooltip.bind(templateCell, () => `Expected type: ${item.templateType}`);

    const plateCell = ui.divText(item.plateField || '');
    const statusCell = ui.span([ui.iconFA(item.icon), ui.divText(item.status)], {style: {color: item.color, gap: '8px'}});

    const row = ui.divH([templateCell, plateCell, statusCell], {style: {padding: '8px 12px', display: 'grid', gridTemplateColumns: '1fr 1fr 100px', cursor: 'pointer'}});
    Array.from(row.children).forEach((c) => (c as HTMLElement).style.alignSelf = 'center');

    if (item.templateField) {
      templateCell.style.cursor = 'grab';
      ui.makeDraggable(templateCell, {
        getDragObject: (): ReconciliationDragObject => ({sourceField: item.templateField, sourceType: item.templateType}),
        getDragCaption: () => `Map: ${item.templateField}`,
        onDragStart: (me: MouseEvent) => {
          tableElement.classList.add('is-dragging');
          return true;
        },
        onDragEnd: () => tableElement.classList.remove('is-dragging'),
      });
    }

    if (item.plateField) {
      const dragOverHandler = (e: DragEvent) => { e.preventDefault(); plateCell.classList.add('d4-drag-over'); };
      const dragLeaveHandler = () => plateCell.classList.remove('d4-drag-over');
      plateCell.addEventListener('dragenter', dragOverHandler);
      plateCell.addEventListener('dragover', dragOverHandler);
      plateCell.addEventListener('dragleave', dragLeaveHandler);
      plateCell.addEventListener('drop', dragLeaveHandler);
      ui.makeDroppable(plateCell, {
        acceptDrop: (dragObject: ReconciliationDragObject) => dragObject && typeof dragObject.sourceField === 'string',
        doDrop: (dragObject: ReconciliationDragObject) => onMap(dragObject.sourceField, item.plateField),
      });
    }


    tableElement.appendChild(row);
  });

  return {element: tableElement, conflictCount};
}
