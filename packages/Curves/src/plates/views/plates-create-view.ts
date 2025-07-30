/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {Plate} from '../../plate/plate';
import {
  PlateTemplate,
  plateTemplates,
  PlateType,
  plateTypes,
  savePlate,
  savePlateAsTemplate,
  createNewPlateForTemplate,
  initPlates
} from '../plates-crud';
import {PlateWidget} from '../../plate/plate-widget';
import {parsePlateFromCsv} from '../../plate/csv-plates';
import {renderValidationResults} from './plates-validation-panel';

// Define a type for our more complex state
type PlateFile = {plate: Plate, file: DG.FileInfo};
type TemplateState = {
  plates: PlateFile[],
  activePlateIdx: number
};

export function createPlatesView(): DG.View {
  const view = DG.View.create();
  view.name = 'Create Plate';

  // The state "store" now supports multiple plates per template
  const templateState = new Map<number, TemplateState>();

  // --- UI Hosts ---
  const platePropertiesHost = ui.divV([]);
  const validationHost = ui.divV([]);
  const fileTabsHost = ui.divH([], {style: {'gap': '4px', 'flexWrap': 'wrap'}});

  let plateType = plateTypes[0];
  let plateTemplate = plateTemplates[0];
  if (!plateTemplate) {
    grok.shell.warning('No plate templates found. Please create a template first.');
    throw new Error('No plate templates found. Please create a template first.');
  }

  const plateWidget = PlateWidget.fromPlate(new Plate(plateType.rows, plateType.cols));
  plateWidget.editable = true;

  const setTemplate = async (template: PlateTemplate) => {
    plateTemplate = template;
    const state = templateState.get(template.id);
    const activePlateFile = state ? state.plates[state.activePlateIdx] : null;

    ui.empty(platePropertiesHost);
    ui.empty(validationHost);
    ui.empty(fileTabsHost);

    renderFileTabs(template, state);

    if (activePlateFile) { // If there's an active plate, restore its state
      plateWidget.plate = activePlateFile.plate;
      validationHost.appendChild(renderValidationResults(activePlateFile.plate, template));
      const form = ui.input.form(activePlateFile.plate.details, template.plateProperties.map((p) => DG.Property.js(p.name!, p.type! as DG.TYPE)));
      platePropertiesHost.appendChild(form);
    } else {
      plateWidget.plate = await createNewPlateForTemplate(plateType, template);
      const form = ui.input.form({}, template.plateProperties.map((p) => DG.Property.js(p.name!, p.type! as DG.TYPE)));
      platePropertiesHost.appendChild(form);
    }
    plateWidget.refresh();
  };

  const renderFileTabs = (template: PlateTemplate, state: TemplateState | undefined) => {
    state?.plates.forEach((plateFile, idx) => {
      const tabLabel = plateFile.plate.barcode ? `${plateFile.file.name} (${plateFile.plate.barcode})` : plateFile.file.name;
      const tab = ui.divText(tabLabel, 'd4-tag-editor');
      tab.style.padding = '4px 8px';
      tab.style.borderRadius = '4px';
      tab.style.cursor = 'pointer';
      tab.style.border = '1px solid var(--grey-2)';
      tab.style.display = 'flex';
      tab.style.alignItems = 'center';

      if (idx === state.activePlateIdx) {
        tab.style.backgroundColor = 'var(--blue-1)';
        tab.style.color = 'var(--blue-4)';
      }

      tab.onclick = () => {
        state.activePlateIdx = idx;
        setTemplate(template);
      };

      const clearIcon = ui.iconFA('times', () => {
        state.plates.splice(idx, 1);
        if (state.activePlateIdx >= idx)
          state.activePlateIdx = Math.max(0, state.activePlateIdx - 1);
        if (state.plates.length === 0)
          templateState.delete(template.id);
        setTemplate(template);
      }, 'Clear imported data');
      clearIcon.style.marginLeft = '8px';

      tab.appendChild(clearIcon);
      fileTabsHost.appendChild(tab);
    });

    // Render the "Import" button
    const importInput = ui.input.file('', {
      onValueChanged: async (file: DG.FileInfo) => {
        if (!file) return;

        const parsedPlates = await parsePlateFromCsv(await file.readAsString());
        const currentState = templateState.get(template.id) ?? {plates: [], activePlateIdx: -1};

        for (const plate of parsedPlates)
          currentState.plates.push({plate, file});


        currentState.activePlateIdx = currentState.plates.length - 1;
        templateState.set(template.id, currentState);
        await setTemplate(template);
      }
    });

    importInput.root.querySelector('button')?.replaceChildren(ui.iconFA('plus'));
    importInput.root.querySelector('label')?.remove();
    importInput.root.title = 'Import new CSV file';
    fileTabsHost.appendChild(importInput.root);
  };

  // --- Top Level Selectors ---
  const plateTypeSelector = ui.input.choice('Plate type', {
    value: plateType.name, items: plateTypes.map((pt) => pt.name),
    onValueChanged: (v) => { plateType = plateTypes.find((pt) => pt.name === v)!; setTemplate(plateTemplate); }
  });
  const plateTemplateSelector = ui.input.choice('Template', {
    value: plateTemplate.name, items: plateTemplates.map((pt) => pt.name),
    onValueChanged: (v) => setTemplate(plateTemplates.find((pt) => pt.name === v)!)
  });

  // --- Assemble the View Layout ---
  const leftPanel = ui.divV([
    ui.h2('Plate Properties'),
    platePropertiesHost,
    validationHost,
  ], {style: {'minWidth': '320px', 'flexGrow': '0'}});

  const rightPanel = ui.divV([
    fileTabsHost,
    plateWidget.root,
  ], {style: {'flexGrow': '1'}});

  view.root.appendChild(ui.divV([
    ui.form([plateTypeSelector, plateTemplateSelector]),
    ui.divH([leftPanel, rightPanel], {style: {gap: '20px'}})
  ]));

  // Initial setup call
  setTemplate(plateTemplates[0]);

  // --- Ribbon Panels ---
  const getPlate = () => {
    const state = templateState.get(plateTemplate.id);
    if (state && state.plates.length > 0) {
      const activePlate = state.plates[state.activePlateIdx].plate;
      return activePlate;
    }
    return plateWidget.plate;
  };

  view.setRibbonPanels([[
    ui.bigButton('CREATE', async () => {
      const plateToSave = getPlate();
      if (!plateToSave) {
        grok.shell.warning('No active plate to save.');
        return;
      }
      await savePlate(plateToSave);
      grok.shell.info(`Plate created: ${plateToSave.id}`);
    }),
    ui.button('SAVE TEMPLATE', async () => {
      const plateToSave = getPlate();
      if (!plateToSave) {
        grok.shell.warning('No active plate to save as template.');
        return;
      }
      await savePlateAsTemplate(plateToSave, plateTemplate);
      await initPlates(true);
      grok.shell.info(`Plate template saved: ${plateToSave.id}`);
    })
  ]]);

  return view;
}
