/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Plate} from '../../plate/plate';
import {PlateTemplate} from '../plates-crud';
import {fromEvent, Subscription} from 'rxjs';
import {debounceTime} from 'rxjs/operators';

// --- Configuration for autosuggestions ---
const REQUIRED_ANALYSIS_FIELDS = [
  {name: 'Activity', required: true},
  {name: 'Concentration', required: true},
  {name: 'SampleID', required: true},
];

const MOCKED_REQUIRED_TEMPLATE_FIELDS = ['Target', 'Assay Format'];

/**
 * Creates a custom input that shows a dropdown with icon-based suggestions.
 */
function createMappingInput(
  initialValue: string,
  template: PlateTemplate,
  usedTemplateProperties: Set<string>,
  onCommit: (newValue: string) => void
): DG.InputBase<string> {
  const textInput = ui.input.string('', {value: initialValue});
  const popup = ui.div([], 'custom-suggestion-popup');
  let clickOutsideSubscription: Subscription | null = null;

  const hideSuggestions = () => {
    if (popup.parentElement)
      popup.remove();
    clickOutsideSubscription?.unsubscribe();
    clickOutsideSubscription = null;
  };

  const showSuggestions = () => {
    if (popup.parentElement)
      hideSuggestions();

    popup.innerHTML = ''; // Clear previous content

    // --- 1. Create and add the Legend ---
    const legend = ui.divH([
      ui.divH([ui.iconFA('file-alt'), ui.span(['Template'])]),
      ui.divH([ui.iconFA('chart-line'), ui.span(['Analysis'])]),
      ui.divH([ui.iconFA('lock'), ui.span(['Required'])]),
    ], 'custom-suggestion-legend');
    popup.appendChild(legend);
    popup.appendChild(ui.div([], 'dg-separator'));

    // --- 2. Get suggestions from different sources ---
    const templateSuggestions = template.wellProperties
      .filter((p) => p && p.name && !usedTemplateProperties.has(p.name))
      .map((p) => ({
        name: p.name!,
        isTemplate: true,
        isAnalysis: false,
        isRequired: MOCKED_REQUIRED_TEMPLATE_FIELDS.includes(p.name!),
      }));

    const analysisSuggestions = REQUIRED_ANALYSIS_FIELDS
      .filter((ap) => !usedTemplateProperties.has(ap.name))
      .map((ap) => ({
        name: ap.name,
        isTemplate: false,
        isAnalysis: true,
        isRequired: ap.required,
      }));

    // --- 3. Merge and de-duplicate suggestions ---
    const suggestionMap = new Map<string, {name: string, isTemplate: boolean, isAnalysis: boolean, isRequired: boolean}>();
    [...templateSuggestions, ...analysisSuggestions].forEach((s) => {
      const existing = suggestionMap.get(s.name);
      if (existing) {
        existing.isTemplate = existing.isTemplate || s.isTemplate;
        existing.isAnalysis = existing.isAnalysis || s.isAnalysis;
        existing.isRequired = existing.isRequired || s.isRequired;
      } else {
        suggestionMap.set(s.name, {...s});
      }
    });
    const allAvailableSuggestions = Array.from(suggestionMap.values());

    // --- 4. Filter by user input & Sort ---
    const currentVal = textInput.value.toLowerCase();
    const filtered = allAvailableSuggestions
      .filter((s) => s.name.toLowerCase().includes(currentVal))
      .sort((a, b) => {
        const scoreA = a.isRequired ? 100 : 0;
        const scoreB = b.isRequired ? 100 : 0;
        if (scoreA !== scoreB) return scoreB - scoreA; // Required items first
        return a.name.localeCompare(b.name); // Alphabetical tie-break
      });

    // --- 5. Render the sorted list with icons ---
    if (filtered.length === 0) {
      hideSuggestions();
      return;
    }

    for (const prop of filtered) {
      const nameEl = ui.span([prop.name], 'custom-suggestion-name');
      const iconsHost = ui.divH([], 'custom-suggestion-icons');

      if (prop.isRequired)
        iconsHost.appendChild(ui.iconFA('lock'));
      if (prop.isTemplate)
        iconsHost.appendChild(ui.iconFA('file-alt'));
      if (prop.isAnalysis)
        iconsHost.appendChild(ui.iconFA('chart-line'));

      const itemEl = ui.divH([nameEl, iconsHost], 'custom-suggestion-item');

      itemEl.addEventListener('mousedown', (e) => {
        e.preventDefault();
        textInput.value = prop.name;
        onCommit(prop.name);
        hideSuggestions();
      });
      popup.appendChild(itemEl);
    }

    const rect = textInput.input.getBoundingClientRect();
    popup.style.left = `${rect.left}px`;
    popup.style.top = `${rect.bottom + 2}px`;
    popup.style.minWidth = `${rect.width}px`;
    document.body.appendChild(popup);

    setTimeout(() => {
      clickOutsideSubscription = fromEvent(document, 'click').subscribe((e: Event) => {
        if (!popup.contains(e.target as Node) && !textInput.root.contains(e.target as Node))
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

  // Track which template properties have been used, either by direct match or by user mapping.
  const usedTemplateProperties = new Set<string>();
  reconciliationMap.forEach((_, mappedName) => usedTemplateProperties.add(mappedName));

  for (const pCol of plateColumns) {
    const matchedProp = templateWellProps.find((tProp) => tProp.name === pCol.name);

    if (matchedProp) {
      usedTemplateProperties.add(matchedProp.originalName);
      const typeMatch = matchedProp.type === pCol.type;
      if (!typeMatch) conflictCount++;
      mappingData.push({
        plateField: pCol.originalName,
        templateField: matchedProp.originalName,
        status: typeMatch ? 'Match' : 'Type Mismatch',
      });
    } else {
      conflictCount++;
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
    const isMappedByRecon = reconciliationMap.has(currentName);
    const originalNameIfRecon = reconciliationMap.get(currentName);
    const isDirectMatch = item.status === 'Match' && item.templateField && !isMappedByRecon;

    const columnCell = ui.divText(isMappedByRecon ? originalNameIfRecon! : currentName);
    let rightCell: HTMLElement;

    const onCommit = (newValue: string) => {
      onMap(currentName, newValue);
      // Add the newly mapped property to the used set so other inputs can update
      usedTemplateProperties.add(newValue);
    };

    if (isMappedByRecon) {
      const host = ui.divText(currentName);
      const icon = ui.iconFA('times', () => onUndo(currentName), 'Undo this mapping');
      rightCell = ui.divH([host, icon], 'locked-mapping-input');
    } else if (isDirectMatch) {
      rightCell = ui.divText(item.templateField, 'd4-read-only-input');
    } else {
      const mappingInput = createMappingInput('', template, usedTemplateProperties, onCommit);
      rightCell = mappingInput.root;
    }

    const row = ui.divH([columnCell, rightCell], 'plate-validation-table-row');
    row.classList.add('d4-table-row-hover');
    tableElement.appendChild(row);
  });

  return {element: tableElement, conflictCount};
}
