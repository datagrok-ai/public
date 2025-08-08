/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Plate} from '../../plate/plate';
import {PlateTemplate} from '../plates-crud';
import {fromEvent, Subscription} from 'rxjs';
import {debounceTime} from 'rxjs/operators';

/**
 * Creates a custom input that shows a custom-styled dropdown with suggestions on focus.
 */
function createMappingInput(
  initialValue: string,
  template: PlateTemplate,
  onCommit: (newValue: string) => void
): DG.InputBase<string> {
  const textInput = ui.input.string('', {value: initialValue});
  const popup = ui.div([], 'custom-suggestion-popup');
  let clickOutsideSubscription: Subscription | null = null;

  const requiredFields = ['activity', 'concentration', 'sampleid'];
  const templateWellProps = template.wellProperties
    .filter((p) => p && p.name)
    .map((p) => ({
      name: p.name!,
      required: requiredFields.includes(p.name!.toLowerCase()),
    }));

  const hideSuggestions = () => {
    if (popup.parentElement)
      popup.remove();
    clickOutsideSubscription?.unsubscribe();
    clickOutsideSubscription = null;
  };

  const showSuggestions = () => {
    if (popup.parentElement)
      hideSuggestions();

    popup.innerHTML = '';
    const currentVal = textInput.value.toLowerCase();
    const filtered = templateWellProps.filter((p) => p.name.toLowerCase().includes(currentVal));

    if (filtered.length === 0) {
      hideSuggestions();
      return;
    }

    for (const prop of filtered) {
      const nameEl = ui.span([prop.name], 'custom-suggestion-name');
      const annotationText = `from ${template.name}${prop.required ? ' (required)' : ''}`;
      const annotationEl = ui.span([annotationText], 'custom-suggestion-annotation');
      const itemEl = ui.divH([nameEl, annotationEl], 'custom-suggestion-item');

      itemEl.addEventListener('mousedown', (e) => {
        e.preventDefault(); // Prevent input from losing focus
        textInput.value = prop.name;
        onCommit(prop.name);
        hideSuggestions();
      });
      popup.appendChild(itemEl);
    }

    const rect = textInput.input.getBoundingClientRect();
    popup.style.left = `${rect.left}px`;
    popup.style.top = `${rect.bottom + 2}px`;
    // FIXED: Use min-width to allow the popup to expand, instead of a fixed width.
    popup.style.minWidth = `${rect.width}px`;
    document.body.appendChild(popup);

    // Use a timeout to prevent the same click that opened the popup from closing it.
    setTimeout(() => {
      clickOutsideSubscription = fromEvent(document, 'click').subscribe((e: Event) => {
        if (!textInput.root.contains(e.target as Node))
          hideSuggestions();
      });
    }, 0);
  };

  fromEvent(textInput.input, 'focus').subscribe(showSuggestions);
  fromEvent(textInput.input, 'input').pipe(debounceTime(150)).subscribe(showSuggestions);
  fromEvent(textInput.input, 'change').subscribe(() => {
    if (textInput.value)
      onCommit(textInput.value);
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
    originalName: c.name,
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
        templateField: '', // Start with empty mapping
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
      const mappingInput = createMappingInput(item.templateField, template, (newValue) => {
        onMap(currentName, newValue);
      });
      rightCell = mappingInput.root;
    }

    const row = ui.divH([columnCell, rightCell], 'plate-validation-table-row');
    row.classList.add('d4-table-row-hover');
    tableElement.appendChild(row);
  });

  return {element: tableElement, conflictCount};
}
