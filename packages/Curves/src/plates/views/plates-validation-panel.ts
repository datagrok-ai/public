import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Plate} from '../../plate/plate';
import {PlateTemplate} from '../plates-crud';

/**
 * Renders a validation table comparing the plate's data schema to the template's requirements.
 * @param plate The Plate object parsed from the CSV.
 * @param template The currently selected PlateTemplate.
 * @returns An HTMLElement containing the validation results.
 */
export function renderValidationResults(plate: Plate, template: PlateTemplate): HTMLElement {
  // ... (the internal logic of this function remains exactly the same) ...
  const templateWellProps = template.wellProperties.map((p) => ({name: p.name!.toLowerCase(), type: p.type as DG.TYPE}));
  const plateColumns = plate.data.columns.toList().map((c) => ({name: c.name.toLowerCase(), type: c.type as DG.TYPE, originalName: c.name}));
  const validationData = [];
  const OK_ICON = '✔';
  const ERROR_ICON = '✖';
  const INFO_ICON = 'ⓘ';

  for (const requiredProp of templateWellProps) {
    const foundCol = plateColumns.find((p) => p.name === requiredProp.name);
    if (foundCol) {
      const typeMatch = foundCol.type === requiredProp.type;
      validationData.push({
        'Template Field': requiredProp.name,
        'Expected Type': requiredProp.type,
        'Found Column': foundCol.originalName,
        'Status': typeMatch ?
          ui.divText(`${OK_ICON} Types Match`, {style: {color: DG.Color.toRgb(DG.Color.green)}}) :
          ui.divText(`${ERROR_ICON} Type Mismatch`, {style: {color: DG.Color.toRgb(DG.Color.red)}})
      });
    } else {
      validationData.push({
        'Template Field': requiredProp.name, 'Expected Type': requiredProp.type, 'Found Column': '—',
        'Status': ui.divText(`${ERROR_ICON} Missing`, {style: {color: DG.Color.toRgb(DG.Color.red)}})
      });
    }
  }

  for (const plateCol of plateColumns) {
    const isRequired = templateWellProps.some((p) => p.name === plateCol.name);
    if (!isRequired) {
      validationData.push({
        'Template Field': '—', 'Expected Type': '—', 'Found Column': plateCol.originalName,
        'Status': ui.divText(`${INFO_ICON} Extra Column`, {style: {color: DG.Color.toRgb(DG.Color.blue)}})
      });
    }
  }

  if (validationData.length === 0)
    return ui.divV([ui.h2('Well Properties'), ui.divText('No well properties are defined in the selected template to validate against.')]);

  // MODIFIED: Added a title
  return ui.divV([
    ui.h2('Well Properties'),
    DG.HtmlTable.create(validationData, (item: any) =>
      [item['Template Field'], item['Expected Type'], item['Found Column'], item.Status]
    ).root
  ]);
}
