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
  initPlates,
} from '../plates-crud';
import {PlateWidget} from '../../plate/plate-widget';
// Assuming your new CSV parser will live here
import {renderValidationResults} from './plates-validation-panel';
import {parsePlateFromCsv} from './../../plate/csv-plates';

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
  const validationHost = ui.divV([], {style: {marginTop: '10px'}});
  const fileTabsHost = ui.divH([], {style: {'gap': '4px', 'flexWrap': 'wrap', 'alignItems': 'center'}});
  const wellPropsHeaderHost = ui.div();

  let plateType = plateTypes[0];
  let plateTemplate = plateTemplates[0];
  if (!plateTemplate) {
    grok.shell.error('No plate templates found. Please create a template first.');
    // Return an empty view to avoid further errors
    return view;
  }

  const plateWidget = PlateWidget.fromPlate(new Plate(plateType.rows, plateType.cols));
  plateWidget.editable = true;

  // --- File Tabs UI Renderer ---
  // Moved this function definition before `setTemplate` to fix the scope issue.
  const renderFileTabs = (template: PlateTemplate, state: TemplateState | undefined) => {
    ui.empty(fileTabsHost);

    // Render existing file tabs
    state?.plates.forEach((plateFile, idx) => {
      const validation = renderValidationResults(plateFile.plate, template);
      const hasConflicts = validation.conflictCount > 0;
      const tabLabel = plateFile.plate.barcode ?? `Plate ${idx + 1}`;

      const tab = ui.divText(tabLabel, 'd4-tag-editor');
      ui.tooltip.bind(tab, () => plateFile.file.name);

      Object.assign(tab.style, {
        padding: '4px 8px', borderRadius: '4px', cursor: 'pointer',
        border: '1px solid var(--grey-2)', display: 'flex', alignItems: 'center',
        position: 'relative', transition: 'background-color 0.1s, border-color 0.1s, box-shadow 0.1s'
      });

      if (idx === state.activePlateIdx) {
        tab.style.borderColor = 'var(--blue-2)';
        tab.style.boxShadow = '0 0 0 1px var(--blue-2)';
        tab.style.backgroundColor = 'var(--blue-0)';
      } else {
        tab.addEventListener('mouseenter', () => tab.style.backgroundColor = 'var(--grey-1)');
        tab.addEventListener('mouseleave', () => tab.style.backgroundColor = '');
      }

      tab.onclick = () => {
        if (state) {
          state.activePlateIdx = idx;
          // We can call setTemplate here since it will be defined by the time this event is fired.
          setTemplate(template);
        }
      };

      if (hasConflicts) {
        const conflictDot = ui.div('');
        Object.assign(conflictDot.style, {
          position: 'absolute', top: '-3px', right: '-3px', width: '8px', height: '8px',
          backgroundColor: 'var(--red-3)', borderRadius: '50%',
          border: '1px solid var(--grey-0)'
        });
        tab.appendChild(conflictDot);
      }

      const clearIcon = ui.iconFA('times', () => {
        if (state) {
          state.plates.splice(idx, 1);
          if (state.activePlateIdx >= idx)
            state.activePlateIdx = Math.max(0, state.activePlateIdx - 1);

          if (state.plates.length === 0)
            templateState.delete(template.id);

          setTemplate(template);
        }
      }, 'Clear imported data');
      clearIcon.style.marginLeft = '8px';

      tab.appendChild(clearIcon);
      fileTabsHost.appendChild(tab);
    });

    // Render the "Import" button
    const importInput = ui.input.file('', {
      onValueChanged: async (file: DG.FileInfo) => {
        if (!file) return;

        try {
          const parsedPlates = await parsePlateFromCsv(await file.readAsString());
          const currentState = templateState.get(template.id) ?? {plates: [], activePlateIdx: -1};

          for (const plate of parsedPlates)
            currentState.plates.push({plate, file});

          currentState.activePlateIdx = currentState.plates.length - 1;
          templateState.set(template.id, currentState);
          await setTemplate(template);
        } catch (e: any) {
          grok.shell.error(`Failed to parse CSV: ${e.message}`);
        }
      }
    });

    const buttonElement = importInput.root.querySelector('button');
    if (buttonElement instanceof HTMLElement) {
      ui.empty(buttonElement);
      buttonElement.appendChild(ui.iconFA('plus'));
    }
    importInput.root.querySelector('label')?.remove();
    ui.tooltip.bind(importInput.root, 'Import new CSV file');
    fileTabsHost.appendChild(importInput.root);
  };

  // --- Main Controller Function ---
  const setTemplate = async (template: PlateTemplate) => {
    plateTemplate = template;
    const state = templateState.get(template.id);
    const activePlateFile = state ? state.plates[state.activePlateIdx] : null;

    // Clear all dynamic UI hosts
    ui.empty(platePropertiesHost);
    ui.empty(validationHost);
    ui.empty(wellPropsHeaderHost);

    // Render the row of file tabs
    renderFileTabs(template, state);

    let conflictCount = 0;
    if (activePlateFile) {
      plateWidget.plate = activePlateFile.plate;
      const validation = renderValidationResults(activePlateFile.plate, template);
      validationHost.appendChild(validation.element);
      conflictCount = validation.conflictCount;
      const form = ui.input.form(activePlateFile.plate.details, template.plateProperties.map((p) => DG.Property.js(p.name!, p.type! as DG.TYPE)));
      platePropertiesHost.appendChild(form);
    } else {
      plateWidget.plate = await createNewPlateForTemplate(plateType, template);
      const form = ui.input.form({}, template.plateProperties.map((p) => DG.Property.js(p.name!, p.type! as DG.TYPE)));
      platePropertiesHost.appendChild(form);
      validationHost.appendChild(ui.divText('Import a file to see validation results.', {style: {color: 'var(--grey-5)', padding: '10px'}}));
    }

    // --- (4) Render Well Properties Header with Fix Indicator ---
    const appliedFixes = 3; // Placeholder as requested

    // FIX: Create a "chip" manually using a styled ui.div, since ui.chip doesn't exist.
    const indicator = ui.divText(`${conflictCount} pending, ${appliedFixes} applied`, {
      style: {
        backgroundColor: 'var(--grey-1)',
        border: '1px solid var(--grey-2)',
        borderRadius: '12px',
        padding: '2px 8px',
        fontSize: '11px',
        cursor: 'pointer',
        color: 'var(--grey-6)',
        fontWeight: '500'
      }
    });
    indicator.addEventListener('click', () => {
      ui.dialog('Reconciliation Summary')
        .add(ui.divV([
          ui.h2('Pending Fixes'),
          ui.divText(`- ${conflictCount} properties have mismatches or are missing.`),
          ui.h2('Applied Fixes'),
          ui.divText('- 3 properties have been manually mapped.'), // Placeholder
        ]))
        .show();
    });
    ui.tooltip.bind(indicator, 'View reconciliation summary');

    const header = ui.h2('Well Properties');
    wellPropsHeaderHost.appendChild(ui.divH([header, indicator], {style: {justifyContent: 'space-between', alignItems: 'center'}}));

    plateWidget.refresh();
  };

  // --- (2 & 3) Improved Top Level Selectors ---
  const plateTypeSelector = ui.input.choice('Plate Type', {
    value: plateType.name, items: plateTypes.map((pt) => pt.name),
    onValueChanged: (v) => { plateType = plateTypes.find((pt) => pt.name === v)!; setTemplate(plateTemplate); }
  });
  const plateTemplateSelector = ui.input.choice('Template', {
    value: plateTemplate.name, items: plateTemplates.map((pt) => pt.name),
    onValueChanged: (v) => setTemplate(plateTemplates.find((pt) => pt.name === v)!)
  });

  // Add styling to make inputs stand out
  plateTypeSelector.input.style.border = '1px solid var(--grey-3)';
  plateTypeSelector.input.style.backgroundColor = 'var(--grey-1)';
  plateTemplateSelector.input.style.border = '1px solid var(--grey-3)';
  plateTemplateSelector.input.style.backgroundColor = 'var(--grey-1)';

  const topPanel = ui.divH([
    plateTypeSelector.root,
    plateTemplateSelector.root,
  ], {style: {gap: '24px', alignItems: 'center', marginBottom: '10px'}});

  // --- Assemble the View Layout ---
  const leftPanel = ui.divV([
    ui.h2('Plate Properties'),
    platePropertiesHost,
    wellPropsHeaderHost, // The new header with indicator goes here
    validationHost,
  ], {style: {minWidth: '320px', maxWidth: '400px', flexGrow: '0', gap: '10px'}});

  const rightPanel = ui.divV([
    fileTabsHost,
    plateWidget.root,
  ], {style: {flexGrow: '1', gap: '10px'}});

  view.root.appendChild(ui.divV([
    topPanel,
    ui.divH([leftPanel, rightPanel], {style: {gap: '20px'}}),
  ]));

  // --- Initial setup call ---
  setTemplate(plateTemplates[0]);

  // --- Ribbon Panels ---

  const getPlate = () => {
    const state = templateState.get(plateTemplate.id);
    if (state && state.plates.length > 0)
      return state.plates[state.activePlateIdx].plate;

    return plateWidget.plate;
  };

  view.setRibbonPanels([[
    ui.bigButton('CREATE', async () => {
      const plateToSave = getPlate();
      if (!plateToSave) {
        grok.shell.warning('No active plate to save.');
        return;
      }

      // Check for conflicts before saving
      const validation = renderValidationResults(plateToSave, plateTemplate);
      if (validation.conflictCount > 0) {
        // As requested, skipping the confirm dialog for now.
        // You can re-enable this later if grok.shell.confirm becomes available or you implement a custom dialog.
        // if (!await grok.shell.confirm('There are unresolved validation issues. Are you sure you want to create the plate?'))
        //    return;
        grok.shell.warning('Saving plate with unresolved validation issues.');
      }

      // FIX: The `savePlate` function expects only the plate object and an optional `options` object.
      // We are not passing any special options here, so we just pass the plate.
      // The template association seems to be handled differently (e.g., via `savePlateAsTemplate`).
      plateToSave.plateTemplateId = plateTemplate.id; // Store the template ID on the plate object itself.
      await savePlate(plateToSave);

      grok.shell.info(`Plate created: ${plateToSave.barcode}`);
    }),
    ui.button('SAVE TEMPLATE', async () => {
      const plateToSave = getPlate();
      if (!plateToSave) {
        grok.shell.warning('No active plate to save as template.');
        return;
      }
      await savePlateAsTemplate(plateToSave, plateTemplate);
      await initPlates(true); // Force reload of templates
      grok.shell.info(`Plate template updated: ${plateTemplate.name}`);
    })
  ]]);


  return view;
}
