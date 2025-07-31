/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Plate} from '../../plate/plate';
import {PlateTemplate} from '../plates-crud';

/**
 * Shows a resolution dialog for a specific validation item.
 */
function showResolutionDialog(item: any) {
  const dialogContent = ui.divV([
    ui.h3('Resolve Property Mismatch'),
    ui.divText(`Template requires: ${item.templateField || 'N/A'} (type: ${item.templateType || 'N/A'})`),
    ui.divText(`Plate provides: ${item.plateField || 'N/A'}`),
    ui.divText(`Status: ${item.status}`),
    ui.divText('RESOLUTION UX', {
      style: {
        marginTop: '20px',
        padding: '24px',
        border: '2px dashed var(--grey-3)',
        textAlign: 'center',
        color: 'var(--grey-6)',
      },
    }),
  ]);

  ui.dialog('Property Resolution')
    .add(dialogContent)
    .show();
}

/**
 * Renders a validation table comparing the plate's data schema to the template's requirements.
 * @param plate The Plate object parsed from the CSV.
 * @param template The currently selected PlateTemplate.
 * @returns An object containing the validation element and the number of conflicts.
 */
export function renderValidationResults(plate: Plate, template: PlateTemplate): { element: HTMLElement, conflictCount: number } {
  const templateWellProps = template.wellProperties.map((p) => ({name: p.name!.toLowerCase(), type: p.type as DG.TYPE, originalName: p.name!}));
  const plateColumns = plate.data.columns.toList().map((c) => ({name: c.name.toLowerCase(), type: c.type as DG.TYPE, originalName: c.name}));
  const validationData: any[] = [];
  let conflictCount = 0;

  const ICONS = {
    OK: '✔',
    ERROR: '✖',
    INFO: 'ⓘ',
  };

  // Check for required properties
  for (const requiredProp of templateWellProps) {
    const foundCol = plateColumns.find((p) => p.name === requiredProp.name);
    if (foundCol) {
      const typeMatch = foundCol.type === requiredProp.type;
      if (!typeMatch) conflictCount++;
      validationData.push({
        templateField: requiredProp.originalName,
        templateType: requiredProp.type,
        plateField: foundCol.originalName,
        status: typeMatch ? 'Match' : 'Type Mismatch',
        icon: typeMatch ? ICONS.OK : ICONS.ERROR,
        color: typeMatch ? 'var(--green-3)' : 'var(--red-3)',
        isConflict: !typeMatch,
      });
    } else {
      conflictCount++;
      validationData.push({
        templateField: requiredProp.originalName,
        templateType: requiredProp.type,
        plateField: '', // Empty string instead of '—'
        status: 'Missing',
        icon: ICONS.ERROR,
        color: 'var(--red-3)',
        isConflict: true,
      });
    }
  }

  // Check for extra columns in the plate
  for (const plateCol of plateColumns) {
    if (!templateWellProps.some((p) => p.name === plateCol.name)) {
      validationData.push({
        templateField: '', // Empty string
        templateType: null,
        plateField: plateCol.originalName,
        status: 'Extra',
        icon: ICONS.INFO,
        color: 'var(--blue-3)',
        isConflict: false, // Extra columns are not considered blocking conflicts
      });
    }
  }

  if (validationData.length === 0) {
    return {
      element: ui.divText('No well properties defined in template.', {style: {padding: '10px', color: 'var(--grey-5)'}}),
      conflictCount: 0,
    };
  }

  const table = ui.divV([], {style: {border: '1px solid var(--grey-2)', width: '100%'}});

  // Header - just bold text, no background
  const header = ui.divH([
    ui.divText(template.name, {style: {fontWeight: 'bold'}}),
    ui.divText(plate.barcode ?? 'Plate', {style: {fontWeight: 'bold'}}),
    ui.divText('Status', {style: {fontWeight: 'bold'}}),
  ], {style: {padding: '8px 12px', borderBottom: '1px solid var(--grey-2)'}});
  header.style.display = 'grid';
  header.style.gridTemplateColumns = '1fr 1fr 100px';

  table.appendChild(header);

  // Rows
  validationData.forEach((item) => {
    const templateCell = ui.divText(item.templateField);
    if (item.templateType)
      ui.tooltip.bind(templateCell, () => `Expected type: ${item.templateType}`);

    const plateCell = ui.divText(item.plateField);
    const statusCell = ui.span([ui.iconFA(item.icon), ui.divText(item.status)], {style: {color: item.color, gap: '8px'}});

    const row = ui.divH([templateCell, plateCell, statusCell], {style: {padding: '8px 12px', borderBottom: '1px solid var(--grey-1)', cursor: 'pointer'}});
    row.style.display = 'grid';
    row.style.gridTemplateColumns = '1fr 1fr 100px';

    // Center all cell content vertically
    Array.from(row.children).forEach((c) => (c as HTMLElement).style.alignSelf = 'center');

    row.addEventListener('mouseenter', () => row.style.backgroundColor = 'var(--grey-1)');
    row.addEventListener('mouseleave', () => row.style.backgroundColor = 'transparent');
    row.addEventListener('click', () => showResolutionDialog(item));

    table.appendChild(row);
  });

  return {element: table, conflictCount};
}
