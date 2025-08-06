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
import {parsePlates} from '../../plate/csv-plates';
import {toExcelPosition} from '../../plate/utils';
import {PlateDrcAnalysis} from '../../plate/plate-drc-analysis';

import './plates-create-view.css';

type PlateFile = {
    plate: Plate;
    file: DG.FileInfo;
    reconciliationMap: Map<string, string>;
};
type TemplateState = {
    plates: PlateFile[],
    activePlateIdx: number
};

export function createPlatesView(): DG.View {
  const view = DG.View.create();
  view.name = 'Create Plate';
  view.root.classList.add('create-plate-view');


  const templateState = new Map<number, TemplateState>();

  const platePropertiesHost = ui.divV([]);
  const validationHost = ui.divV([]);
  const wellPropsHeaderHost = ui.div();
  const drcAnalysisHost = ui.divV([], 'create-plate-view__drc-host');

  drcAnalysisHost.style.width = '100%';
  drcAnalysisHost.style.flexGrow = '1';

  const tabHeaderContainer = ui.divH([], 'plate-file-tabs-container');

  let plateType = plateTypes[0];
  let plateTemplate = plateTemplates[0];
  if (!plateTemplate) {
    grok.shell.error('No plate templates found.');
    return view;
  }

  let sourceDataFrame: DG.DataFrame | null = null;
  let plateIdentifierColumn: string | null = 'Destination Plate Barcode';
  let replicateIdentifierColumn: string | null = 'Technical Duplicate ID';
  let plateNumberColumn: string | null = 'Plate Number'; // New state variable
  const plateIdentifierHost = ui.div([]);
  const replicateIdentifierHost = ui.div([]);
  const plateNumberHost = ui.div([]); // New UI host


  let reprocessPlates: () => Promise<void>;

  const updatePlateIdentifierControl = () => {
    ui.empty(plateIdentifierHost);
    if (!sourceDataFrame) return; // Do nothing if no file is loaded

    const choiceInput = ui.input.choice('Plate Index', {
      value: plateIdentifierColumn,
      items: [null, ...sourceDataFrame.columns.names()], // Allow 'None' to treat as single plate
      onValueChanged: (newColumn: string | null) => {
        plateIdentifierColumn = newColumn;
        reprocessPlates(); // Trigger re-parsing and UI update
      },
    });
    ui.tooltip.bind(choiceInput.root, 'Select the column that identifies individual plates in the file.');
    plateIdentifierHost.appendChild(choiceInput.root);
  };

  const updateReplicateIdentifierControl = () => {
    ui.empty(replicateIdentifierHost);
    if (!sourceDataFrame) return;

    const choiceInput = ui.input.choice('Replicate Index', {
      value: replicateIdentifierColumn,
      items: [null, ...sourceDataFrame.columns.names()], // Allow 'None'
      onValueChanged: (newColumn: string | null) => {
        replicateIdentifierColumn = newColumn;
        reprocessPlates(); // Trigger re-parsing
      },
    });
    ui.tooltip.bind(choiceInput.root, 'Optional: Select the column that identifies technical replicates.');
    replicateIdentifierHost.appendChild(choiceInput.root);
  };
  const updatePlateNumberControl = () => {
    ui.empty(plateNumberHost);
    if (!sourceDataFrame) return;

    const choiceInput = ui.input.choice('Plate Number Index', {
      value: plateNumberColumn,
      items: [null, ...sourceDataFrame.columns.names()],
      onValueChanged: (newColumn: string | null) => {
        plateNumberColumn = newColumn;
        reprocessPlates();
      },
    });
    ui.tooltip.bind(choiceInput.root, 'Optional: Select the column that identifies the plate number within a barcode.');
    plateNumberHost.appendChild(choiceInput.root);
  };

  reprocessPlates = async () => {
    if (!sourceDataFrame) return;

    // Pass BOTH column names to the updated parsing function
    const parsedPlates = parsePlates(sourceDataFrame, plateIdentifierColumn, replicateIdentifierColumn, plateNumberColumn);


    const currentState = templateState.get(plateTemplate.id) ?? {plates: [], activePlateIdx: -1};
    currentState.plates = [];
    for (const plate of parsedPlates) {
      const dummyFile = {name: `Plate: ${plate.barcode}`, fullPath: ''} as DG.FileInfo;
      currentState.plates.push({plate, file: dummyFile, reconciliationMap: new Map<string, string>()});
    }
    currentState.activePlateIdx = currentState.plates.length > 0 ? 0 : -1;
    templateState.set(plateTemplate.id, currentState);
    await setTemplate(plateTemplate);
  };


  const plateWidget = PlateWidget.fromPlate(new Plate(plateType.rows, plateType.cols));
  plateWidget.editable = true;
  plateWidget.root.style.width = '100%';
  plateWidget.root.style.minWidth = '0';
  plateWidget.root.style.flexGrow = '0';


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
    const appliedIcon = ui.iconFA('exchange-alt');

    const indicator = ui.divH([
      ui.divH([pendingIcon, ui.divText(`${pendingCount}`)], 'reco-summary__item'),
      ui.divH([appliedIcon, ui.divText(`${appliedCount}`)], 'reco-summary__item'),
    ], 'reco-summary');

    indicator.classList.toggle('reco-summary--pending', pendingCount > 0);

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
    ui.empty(tabHeaderContainer);


    if (!state || state.plates.length === 0) {
      tabHeaderContainer.style.visibility = 'hidden';
      return;
    }

    tabHeaderContainer.style.visibility = 'visible';

    tabHeaderContainer.style.display = 'flex';

    state.plates.forEach((plateFile, idx) => {
      const dummyTable = ui.div();
      const validation = renderValidationResults(dummyTable, plateFile.plate, template, handleMapping);
      const hasConflicts = validation.conflictCount > 0;
      const tabLabelText = plateFile.plate.barcode ?? `Plate ${idx + 1}`;

      const labelElement = ui.divText(tabLabelText);
      const controlsContainer = ui.divH([], 'plate-file-tab__controls');

      if (hasConflicts) {
        const conflictDot = ui.div('', 'plate-file-tab__conflict-dot');
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
      clearIcon.classList.add('plate-file-tab__clear-icon');
      controlsContainer.appendChild(clearIcon);

      const header = ui.divH([labelElement, controlsContainer], 'plate-file-tab__header');
      header.classList.toggle('active', idx === state.activePlateIdx);

      header.onclick = (e) => {
        e.preventDefault();
        e.stopPropagation();
        if (state) {
          state.activePlateIdx = idx;
          setTemplate(template);
        }
      };

      ui.tooltip.bind(header, () => plateFile.file.name);
      tabHeaderContainer.appendChild(header);
    });
  };

  setTemplate = async (template: PlateTemplate) => {
    try {
      plateTemplate = template;
      const state = templateState.get(template.id);
      const activeFile = state ? state.plates[state.activePlateIdx] : null;
      ui.empty(platePropertiesHost);
      ui.empty(validationHost);
      ui.empty(wellPropsHeaderHost);
      // --- MODIFIED: Clear the DRC host on every template change ---
      ui.empty(drcAnalysisHost);
      renderFileTabs(template, state);
      let validationResults = {conflictCount: 0};
      if (activeFile && activeFile.plate && activeFile.plate.data) {
        try {
          plateWidget.plate = activeFile.plate;

          const drcCurvesElement = PlateDrcAnalysis.createCurvesGrid(activeFile.plate, plateWidget);

          if (drcCurvesElement) {
            drcAnalysisHost.appendChild(drcCurvesElement);
          } else {
            const info = ui.divText('To see dose-response curves, ensure your CSV has "SampleID", "Concentration", and "Activity" columns.', 'info-message');
            info.style.padding = '10px';
            drcAnalysisHost.appendChild(info);
          }

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
          validationHost.appendChild(ui.divText('Error loading validation results.', 'error-message'));
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
            'Import a CSV file to see validation results and start mapping fields.', 'info-message'
          ));
        } catch (plateCreationError) {
          console.error('Error creating new plate:', plateCreationError);
        }
      }
      const indicator = createReconciliationSummary(activeFile, validationResults);
      const header = ui.h2('Well Properties');
      const templateNameSpan = ui.span([template.name], 'well-props-header__template-name');
      wellPropsHeaderHost.appendChild(ui.divH([header, templateNameSpan, indicator], 'space-between-center'));
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
  ], 'left-panel-section');

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

  const fileInput = ui.input.file('', {
    onValueChanged: async (file: DG.FileInfo) => {
      if (!file) return;
      try {
        sourceDataFrame = DG.DataFrame.fromCsv(await file.readAsString());

        updatePlateIdentifierControl();
        updateReplicateIdentifierControl();
        updatePlateNumberControl(); // Update the new control

        await reprocessPlates();
      } catch (e: any) {
        grok.shell.error(`Failed to parse CSV: ${e.message}`);
      }
    },
  });


  fileInput.root.classList.add('plate-import-button');
  const fileInputButton = fileInput.root.querySelector('button');
  if (fileInputButton) {
    ui.empty(fileInputButton);
    fileInputButton.appendChild(ui.iconFA('upload'));
  }
  fileInput.root.querySelector('label')?.remove();
  ui.tooltip.bind(fileInput.root, 'Import a plate from a CSV file');

  const importContainer = ui.divH([
    ui.divText('Import Plate File'),
    fileInput.root,
    plateIdentifierHost,
    replicateIdentifierHost,
    plateNumberHost, // Add the new host here
  ], 'plate-import-container'); // We'll add this class in the CSS step


  const plateTypeSelector = ui.input.choice('Plate Type', {value: plateType.name, items: plateTypes.map((pt) => pt.name), onValueChanged: (v) => {
    plateType = plateTypes.find((pt) => pt.name === v)!; setTemplate(plateTemplate);
  }});
  const plateTemplateSelector = ui.input.choice('Template', {value: plateTemplate.name, items: plateTemplates.map((pt) => pt.name), onValueChanged: (v) => setTemplate(plateTemplates.find((pt) => pt.name === v)!)});

  const templatePanel = ui.divV([
    ui.h2('Template'),
    plateTypeSelector.root,
    plateTemplateSelector.root,
  ], 'left-panel-section');

  const leftPanel = ui.divV([
    templatePanel,
    ui.divV([ui.h2('Plate Properties'), platePropertiesHost], 'left-panel-section'),
    ui.divV([wellPropsHeaderHost, validationHost], 'left-panel-section'),
    roleAssignmentPanel,
  ], 'create-plate-view__left-panel');

  const rightPanel = ui.divV([
    importContainer,
    tabHeaderContainer,
    plateWidget.root,
    drcAnalysisHost,
  ], 'create-plate-view__right-panel');

  rightPanel.style.width = '100%';
  rightPanel.style.height = '100%';
  rightPanel.style.display = 'flex';
  rightPanel.style.flexDirection = 'column';


  const mainLayout = ui.divH([leftPanel, rightPanel], 'create-plate-view__main-layout');
  view.root.appendChild(mainLayout);

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
