/* eslint-disable prefer-const */
/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {Plate} from '../../plate/plate';
import {
  PlateTemplate, plateTemplates, PlateType, plateTypes,
  savePlate, savePlateAsTemplate, createNewPlateForTemplate, initPlates,
} from '../plates-crud';
import {PlateWidget} from '../../plate/plate-widget';
import {renderValidationResults} from './plates-validation-panel';
import {parsePlateFromCsv} from '../../plate/csv-plates';
import {toExcelPosition} from '../../plate/utils';

type PlateFile = {
    plate: Plate;
    file: DG.FileInfo;
    reconciliationMap: Map<string, string>; // Maps NEW name -> ORIGINAL name
};
type TemplateState = {
    plates: PlateFile[],
    activePlateIdx: number
};

export function createPlatesView(): DG.View {
  const view = DG.View.create();
  view.name = 'Create Plate';

  const templateState = new Map<number, TemplateState>();

  const platePropertiesHost = ui.divV([]);
  const validationHost = ui.divV([]);
  const wellPropsHeaderHost = ui.div();

  const fileTabs = ui.tabControl(null, false);
  fileTabs.root.classList.add('plate-file-tabs');
  fileTabs.root.style.marginLeft = '140px';


  let plateType = plateTypes[0];
  let plateTemplate = plateTemplates[0];
  if (!plateTemplate) {
    grok.shell.error('No plate templates found.');
    return view;
  }

  const plateWidget = PlateWidget.fromPlate(new Plate(plateType.rows, plateType.cols));
  plateWidget.editable = true;

  let setTemplate: (template: PlateTemplate) => Promise<void>;

  const handleMapping = (sourceField: string, targetField: string) => {
    const state = templateState.get(plateTemplate.id);
    if (!state || !state.plates[state.activePlateIdx]) return;

    const activeFile = state.plates[state.activePlateIdx];
    const plate = activeFile.plate;
    const targetColumn = plate.data.col(targetField);
    if (!targetColumn) return;

    let originalName = targetField;
    for (const [newName, oldName] of activeFile.reconciliationMap.entries()) {
      if (newName === targetField) {
        originalName = oldName;
        break;
      }
    }

    activeFile.reconciliationMap.set(sourceField, originalName);
    if (activeFile.reconciliationMap.has(targetField))
      activeFile.reconciliationMap.delete(targetField);

    targetColumn.name = sourceField;
    grok.shell.info(`Mapped '${targetField}' to '${sourceField}'.`);
    setTemplate(plateTemplate);
  };

  const handleUndoMapping = (mappedField: string) => {
    const state = templateState.get(plateTemplate.id);
    if (!state || !state.plates[state.activePlateIdx]) {
      console.warn('No active plate state found for undo mapping');
      return;
    }

    const activeFile = state.plates[state.activePlateIdx];
    const originalName = activeFile.reconciliationMap.get(mappedField);

    if (originalName) {
      const columnToRevert = activeFile.plate.data.col(mappedField);
      if (columnToRevert) {
        columnToRevert.name = originalName;
        activeFile.reconciliationMap.delete(mappedField);
        grok.shell.info(`Reverted mapping for '${mappedField}'.`);
        setTemplate(plateTemplate);
      } else {
        console.warn(`Column '${mappedField}' not found for reverting`);
      }
    } else {
      console.warn(`No original mapping found for '${mappedField}'`);
    }
  };

  const createReconciliationSummary = (activeFile: PlateFile | null, validationResults: {conflictCount: number}) => {
    const appliedCount = activeFile ? activeFile.reconciliationMap.size : 0;
    const pendingCount = validationResults.conflictCount;

    const pendingIcon = ui.iconFA('exclamation-circle');
    const appliedIcon = ui.iconFA('link');

    const indicator = ui.divH([
      ui.divH([pendingIcon, ui.divText(`${pendingCount}`)], {style: {gap: '4px', alignItems: 'center'}}),
      ui.divH([appliedIcon, ui.divText(`${appliedCount}`)], {style: {gap: '4px', alignItems: 'center'}}),
    ], {
      style: {
        backgroundColor: pendingCount > 0 ? 'var(--orange-1)' : 'var(--green-1)',
        border: `1px solid ${pendingCount > 0 ? 'var(--orange-3)' : 'var(--green-3)'}`,
        borderRadius: '12px',
        padding: '4px 12px',
        fontSize: '12px',
        cursor: 'pointer',
        color: pendingCount > 0 ? 'var(--orange-6)' : 'var(--green-6)',
        fontWeight: '500',
        userSelect: 'none',
        gap: '12px',
      },
    });

    ui.tooltip.bind(indicator, 'Click to view mapping details and undo changes');

    indicator.addEventListener('click', () => {
      const summaryDialog = ui.dialog('Field Mapping Summary');
      summaryDialog.root.style.minWidth = '400px';
      const content = ui.divV([], {style: {gap: '16px'}});
      if (appliedCount > 0) {
        content.appendChild(ui.h3('Applied Mappings', {style: {color: 'var(--green-6)', margin: '0'}}));
        const mappingsList = ui.divV([], {style: {gap: '8px'}});
        Array.from(activeFile!.reconciliationMap.entries()).forEach(([newName, oldName]) => {
          const mappingRow = ui.divH([
            ui.divV([
              ui.divText(oldName, {style: {fontSize: '11px', color: 'var(--grey-5)'}}),
              (() => {
                const arrow = ui.iconFA('arrow-down');
                arrow.style.fontSize = '10px';
                arrow.style.color = 'var(--grey-4)';
                return arrow;
              })(),
              ui.divText(newName, {style: {fontWeight: '500'}}),
            ], {style: {alignItems: 'center', gap: '2px'}}),
            (() => {
              const closeIcon = ui.iconFA('times', () => {
                handleUndoMapping(newName);
                summaryDialog.close();
              }, 'Undo this mapping');
              Object.assign(closeIcon.style, {
                cursor: 'pointer', color: 'var(--red-3)', padding: '4px',
                borderRadius: '2px', fontSize: '12px',
              });
              return closeIcon;
            })(),
          ], {
            style: {
              justifyContent: 'space-between', alignItems: 'center', padding: '8px 12px',
              backgroundColor: 'var(--green-0)', border: '1px solid var(--green-2)', borderRadius: '4px',
            },
          });
          mappingsList.appendChild(mappingRow);
        });
        content.appendChild(mappingsList);
      }
      if (pendingCount > 0) {
        if (appliedCount > 0)
          content.appendChild(ui.divText('', {style: {borderTop: '1px solid var(--grey-2)', margin: '8px 0'}}));
        content.appendChild(ui.h3('Pending Conflicts', {style: {color: 'var(--orange-6)', margin: '0'}}));
        content.appendChild(ui.divText(
          `${pendingCount} field${pendingCount > 1 ? 's' : ''} still need${pendingCount === 1 ? 's' : ''} to be mapped. Use drag & drop in the validation table to resolve.`,
          {style: {color: 'var(--grey-6)', fontSize: '14px'}}
        ));
      }
      if (appliedCount === 0 && pendingCount === 0)
        content.appendChild(ui.divText('All fields are properly mapped!', {style: {color: 'var(--green-6)'}}));
      summaryDialog.add(content);
      summaryDialog.show();
    });
    return indicator;
  };

  const renderFileTabs = (template: PlateTemplate, state: TemplateState | undefined) => {
    fileTabs.clear();

    if (!state || state.plates.length === 0) {
      fileTabs.root.style.display = 'none';
      return;
    }
    fileTabs.root.style.display = 'flex';

    state.plates.forEach((plateFile, idx) => {
      const dummyTable = ui.div();
      const validation = renderValidationResults(dummyTable, plateFile.plate, template, handleMapping);
      const hasConflicts = validation.conflictCount > 0;
      const tabLabelText = plateFile.plate.barcode ?? `Plate ${idx + 1}`;

      const pane = fileTabs.addPane(tabLabelText, () => ui.div());
      const header = pane.header;

      ui.empty(header);
      header.style.display = 'flex';
      header.style.alignItems = 'center';

      const labelElement = ui.divText(tabLabelText);

      header.onclick = (e) => {
        e.preventDefault();
        e.stopPropagation();
        if (state) {
          state.activePlateIdx = idx;
          setTemplate(template);
        }
      };

      const controlsContainer = ui.divH([], {style: {gap: '8px', alignItems: 'center', marginLeft: '8px'}});

      if (hasConflicts) {
        const conflictDot = ui.div('', {
          style: {
            width: '8px', height: '8px',
            backgroundColor: 'var(--red-3)', borderRadius: '50%',
            border: '1px solid var(--white)',
            boxShadow: '0 1px 2px rgba(0,0,0,0.1)'
          }
        });
        ui.tooltip.bind(conflictDot, `${validation.conflictCount} unresolved fields`);
        controlsContainer.appendChild(conflictDot);
      }

      const clearIcon = ui.iconFA('times', (e: MouseEvent) => {
        e.stopPropagation();
        if (state) {
          state.plates.splice(idx, 1);
          if (state.activePlateIdx >= idx) state.activePlateIdx = Math.max(0, state.activePlateIdx - 1);
          if (state.plates.length === 0) templateState.delete(template.id);
          setTemplate(template);
        }
      }, 'Remove this plate');
      clearIcon.style.opacity = '0.6';
      controlsContainer.appendChild(clearIcon);

      header.appendChild(labelElement);
      header.appendChild(controlsContainer);
      ui.tooltip.bind(header, () => plateFile.file.name);
    });

    if (state.activePlateIdx !== -1 && fileTabs.panes.length > state.activePlateIdx)
      fileTabs.currentPane = fileTabs.panes[state.activePlateIdx];
  };

  setTemplate = async (template: PlateTemplate) => {
    try {
      plateTemplate = template;
      const state = templateState.get(template.id);
      const activeFile = state ? state.plates[state.activePlateIdx] : null;
      ui.empty(platePropertiesHost);
      ui.empty(validationHost);
      ui.empty(wellPropsHeaderHost);
      renderFileTabs(template, state);
      let validationResults = {conflictCount: 0};
      if (activeFile && activeFile.plate && activeFile.plate.data) {
        try {
          plateWidget.plate = activeFile.plate;
          validationResults = renderValidationResults(validationHost, activeFile.plate, template, handleMapping);
          const plateProperties = template.plateProperties
            .filter((p) => p && p.name && p.type)
            .map((p) => DG.Property.js(p.name!, p.type! as DG.TYPE));
          if (plateProperties.length > 0) {
            const form = ui.input.form(activeFile.plate.details || {}, plateProperties);
            platePropertiesHost.appendChild(form);
          }
        } catch (validationError) {
          console.error('Error during validation rendering:', validationError);
          validationHost.appendChild(ui.divText(
            'Error loading validation results. Please try reimporting the file.',
            {style: {color: 'var(--red-5)', padding: '20px', textAlign: 'center'}}
          ));
        }
      } else {
        try {
          plateWidget.plate = await createNewPlateForTemplate(plateType, template);
          const templateProperties = template.plateProperties
            .filter((p) => p && p.name && p.type)
            .map((p) => DG.Property.js(p.name!, p.type! as DG.TYPE));
          if (templateProperties.length > 0) {
            const form = ui.input.form({}, templateProperties);
            platePropertiesHost.appendChild(form);
          }
          validationHost.appendChild(ui.divText(
            'Import a CSV file to see validation results and start mapping fields.',
            {style: {color: 'var(--grey-5)', padding: '20px', textAlign: 'center', fontStyle: 'italic'}}
          ));
        } catch (plateCreationError) {
          console.error('Error creating new plate:', plateCreationError);
        }
      }
      const indicator = createReconciliationSummary(activeFile, validationResults);
      const header = ui.h2('Well Properties');
      wellPropsHeaderHost.appendChild(ui.divH([header, indicator], {style: {justifyContent: 'space-between', alignItems: 'center'}}));
      try {
        plateWidget.refresh();
        plateWidget.updateRoleSummary();
      } catch (refreshError) {
        console.error('Error refreshing plate widget:', refreshError);
      }
    } catch (error) {
      console.error('Critical error in setTemplate:', error);
      grok.shell.error(`Error updating template: ${error}`);
    }
  };

  const selectionStatusDiv = ui.divText('Drag on the plate to select wells.');
  selectionStatusDiv.style.marginTop = '10px';
  const roles = ['Control', 'Buffer', 'Assay Reagent', 'Sample'];
  const roleInput = ui.input.multiChoice('Roles', {items: roles, value: []});

  const assignButton = ui.button('Assign Roles', async () => {
    const selectedRoles = roleInput.value;
    const plate = plateWidget.plate;
    const selection = plate.data.selection;

    if (selectedRoles?.length === 0 || selection.trueCount === 0) return;

    const roleCol = plate.data.columns.getOrCreate('Role', DG.COLUMN_TYPE.STRING);
    const selectedIndexes = selection.getSelectedIndexes();

    for (const i of selectedIndexes)
      roleCol.set(i, selectedRoles?.join(',') ?? '');

    plateWidget._colorColumn = roleCol;
    plateWidget.grid.invalidate();
    plateWidget.updateRoleSummary();

    selection.setAll(false, true);
  });
  assignButton.disabled = true;

  const roleAssignmentPanel = ui.divV([
    ui.h2('Role Assignment'),
    selectionStatusDiv,
    roleInput.root,
    assignButton,
  ], {style: {borderTop: '1px solid var(--grey-2)', marginTop: '20px', paddingTop: '10px'}});

  const updateAssignButtonState = () => {
    const selectionCount = plateWidget.plate.data.selection.trueCount;
    const rolesCount = roleInput.value?.length ?? 0;
    assignButton.disabled = selectionCount === 0 || rolesCount === 0;
  };

  plateWidget.plate.data.selection.onChanged.subscribe(() => {
    const selectionCount = plateWidget.plate.data.selection.trueCount;

    if (selectionCount > 0) {
      selectionStatusDiv.textContent = `${selectionCount} wells selected`;
      const selectedPositions = Array.from(plateWidget.plate.data.selection.getSelectedIndexes()).map((idx) => {
        const [row, col] = plateWidget.plate.rowIndexToExcel(idx);
        return toExcelPosition(row, col);
      });
      ui.tooltip.bind(selectionStatusDiv, () => ui.divText(selectedPositions.join(', ')));
    } else {
      selectionStatusDiv.textContent = 'Drag on the plate to select wells.';
      ui.tooltip.bind(selectionStatusDiv, () => null);
    }
    updateAssignButtonState();
  });

  roleInput.onInput.subscribe(() => updateAssignButtonState());

  const importInput = ui.input.file('', {
    onValueChanged: async (file: DG.FileInfo) => {
      if (!file) return;
      try {
        const parsedPlates = await parsePlateFromCsv(await file.readAsString());
        const currentState = templateState.get(plateTemplate.id) ?? {plates: [], activePlateIdx: -1};
        for (const plate of parsedPlates)
          currentState.plates.push({plate, file, reconciliationMap: new Map<string, string>()});
        currentState.activePlateIdx = currentState.plates.length - 1;
        templateState.set(plateTemplate.id, currentState);
        await setTemplate(plateTemplate);
      } catch (e: any) { grok.shell.error(`Failed to parse CSV: ${e.message}`); }
    },
  });
  const buttonElement = importInput.root.querySelector('button');
  if (buttonElement) {
    ui.empty(buttonElement);
    buttonElement.appendChild(ui.iconFA('plus'));
    buttonElement.style.padding = '6px 8px';
    buttonElement.style.borderRadius = '6px';
  }
  importInput.root.querySelector('label')?.remove();
  ui.tooltip.bind(importInput.root, 'Import new CSV file');

  const plateTypeSelector = ui.input.choice('Plate Type', {value: plateType.name, items: plateTypes.map((pt) => pt.name), onValueChanged: (v) => { plateType = plateTypes.find((pt) => pt.name === v)!; setTemplate(plateTemplate); }});
  const plateTemplateSelector = ui.input.choice('Template', {value: plateTemplate.name, items: plateTemplates.map((pt) => pt.name), onValueChanged: (v) => setTemplate(plateTemplates.find((pt) => pt.name === v)!)});

  const topPanel = ui.divH([plateTypeSelector.root, plateTemplateSelector.root, ui.div([], {style: {flexGrow: '1'}}), importInput.root], {style: {gap: '24px', alignItems: 'center', marginBottom: '10px'}});

  const leftPanel = ui.divV([ui.h2('Plate Properties'), platePropertiesHost, wellPropsHeaderHost, validationHost, roleAssignmentPanel], {style: {minWidth: '320px', maxWidth: '400px', flexGrow: '0', gap: '10px'}});

  const rightPanel = ui.divV([fileTabs.root, plateWidget.root], {style: {flexGrow: '1', gap: '0px'}});

  view.root.appendChild(ui.divV([topPanel, ui.divH([leftPanel, rightPanel], {style: {gap: '20px'}})]));

  setTemplate(plateTemplates[0]);

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
      const dummyTableElement = ui.div();
      const validation = renderValidationResults(dummyTableElement, plateToSave, plateTemplate, handleMapping);
      if (validation.conflictCount > 0)
        grok.shell.warning('Saving plate with unresolved validation issues.');
      plateToSave.plateTemplateId = plateTemplate.id;
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
      await initPlates(true);
      grok.shell.info(`Plate template updated: ${plateTemplate.name}`);
    }),
  ]]);
  return view;
}
