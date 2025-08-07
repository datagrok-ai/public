/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Plate} from '../../plate/plate';
import {PlateTemplate} from '../plates-crud';
import {fromEvent} from 'rxjs';
import {debounceTime} from 'rxjs/operators';

/**
 * Creates a custom input that allows both free text entry and selection from a list of suggestions.
 */
function createMappingInput(
  initialValue: string,
  suggestions: string[],
  onCommit: (newValue: string) => void
): DG.InputBase<string> {
  const textInput = ui.input.string('', {value: initialValue});

  // Commit the value when the user presses Enter or the input loses focus
  fromEvent(textInput.input, 'change').subscribe(() => {
    if (textInput.value)
      onCommit(textInput.value);
  });

  // Show a suggestions popup menu when the user types
  fromEvent(textInput.input, 'input').pipe(debounceTime(200)).subscribe(() => {
    const currentVal = textInput.value.toLowerCase();
    const filtered = suggestions.filter((s) => s.toLowerCase().includes(currentVal));

    if (filtered.length > 0) {
      const menu = DG.Menu.popup();
      for (const prop of filtered) {
        menu.item(prop, () => {
          textInput.value = prop;
          onCommit(prop);
        });
      }
      menu.show();
    }
  });

  return textInput;
}


/**
 * Renders a validation and mapping table, driven by the columns of the uploaded plate data.
 */
export function renderValidationResults(
  tableElement: HTMLElement,
  plate: Plate,
  template: PlateTemplate,
  onMap: (currentColName: string, templatePropName: string) => void,
  reconciliationMap: Map<string, string>,
  onUndo: (mappedField: string) => void
): { element: HTMLElement, conflictCount: number } {
  const templateWellProps = template.wellProperties
    .filter((p) => p && p.name)
    .map((p) => ({name: p.name!.toLowerCase(), type: p.type as DG.TYPE, originalName: p.name!}));

  const plateColumns = plate.data.columns.toList().map((c) => ({
    name: c.name.toLowerCase(),
    type: c.type as DG.TYPE,
    originalName: c.name
  }));

  const mappingData: any[] = [];
  let conflictCount = 0;

  for (const pCol of plateColumns) {
    const matchedProp = templateWellProps.find((tProp) => tProp.name === pCol.name);

    if (matchedProp) {
      const typeMatch = matchedProp.type === pCol.type;
      if (!typeMatch) conflictCount++;
      mappingData.push({
        plateField: pCol.originalName,
        templateField: matchedProp.originalName,
        status: typeMatch ? 'Match' : 'Type Mismatch',
      });
    } else {
      mappingData.push({
        plateField: pCol.originalName,
        templateField: '', 
        status: 'Not Mapped',
      });
    }
  }

  ui.empty(tableElement);
  if (plateColumns.length === 0) {
    tableElement.appendChild(ui.divText('No columns found in the imported file.'));
    return {element: tableElement, conflictCount: 0};
  }

  tableElement.className = 'plate-validation-table';

  const header = ui.divH([
    ui.divText('Column', {style: {fontWeight: 'bold'}}),
    ui.divText('Maps To', {style: {fontWeight: 'bold'}}),
  ], 'plate-validation-table-header');
  tableElement.appendChild(header);

  mappingData.forEach((item) => {
    const currentName = item.plateField;
    const originalName = reconciliationMap.get(currentName);
    const isMapped = !!originalName;

    const columnCell = ui.divText(isMapped ? originalName : currentName);
    let rightCell: HTMLElement;

    if (isMapped) {
      const host = ui.divText(currentName);
      const icon = ui.iconFA('times', () => onUndo(currentName), 'Undo this mapping');
      rightCell = ui.divH([host, icon], 'locked-mapping-input');
    } else {
      const availableProps = templateWellProps.map((p) => p.originalName);
      const mappingInput = createMappingInput(item.templateField, availableProps, (newValue) => {
        onMap(currentName, newValue);
      });
      rightCell = mappingInput.root;
    }

    const row = ui.divH([columnCell, rightCell], 'plate-validation-table-row');
    row.classList.add('d4-table-row-hover');
    tableElement.appendChild(row);
  });

  // (Conflict count logic can be refined here if needed)

  return {element: tableElement, conflictCount};
}
