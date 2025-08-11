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

type PlateFile = {plate: Plate; file: DG.FileInfo; reconciliationMap: Map<string, string>;
};
type TemplateState = {plates: PlateFile[], activePlateIdx: number
};

function createCollapsiblePanel(header: HTMLElement, content: HTMLElement, expanded: boolean = true): HTMLElement {
  const icon = ui.iconFA(expanded ? 'chevron-down' : 'chevron-right', () => {
    const isExpanded = content.style.display !== 'none';
    content.style.display = isExpanded ? 'none' : 'block';
    icon.classList.toggle('fa-chevron-down', !isExpanded);
    icon.classList.toggle('fa-chevron-right', isExpanded);
  });
  icon.style.marginRight = '8px';
  icon.style.color = 'var(--grey-4)';
  icon.style.cursor = 'pointer';

  // Make the entire header clickable
  header.style.cursor = 'pointer';
  header.onclick = () => icon.click();

  const headerContainer = ui.divH([icon, header]);
  headerContainer.style.alignItems = 'center';

  content.style.display = expanded ? 'block' : 'none';
  content.style.paddingLeft = '24px'; // Indent content to align under the header text

  return ui.divV([headerContainer, content], 'left-panel-section');
}
/**
 * Shows a dialog to manage and apply field mappings to multiple plates.
 * @param allPlates - Array of all currently loaded PlateFile objects.
 * @param sourceMappings - The mappings from the active plate, used as the template to apply.
 * @param onSync - Callback to execute when synchronizing mappings with the selected plates.
 * @param onUndo - Callback to undo a mapping definition from the source.
 */
function showMultiPlateMappingDialog(
  allPlates: PlateFile[],
  sourceMappings: Map<string, string>,
  onSync: (mappings: Map<string, string>, selectedIndexes: number[]) => void,
  onUndo: (mappedField: string) => void
) {
  const dialog = ui.dialog('Apply Field Mappings');
  dialog.root.style.minWidth = '500px';
  const content = ui.divV([], {style: {gap: '20px'}});

  // --- Section 1: Mappings to Apply ---
  const mappingsSection = ui.divV([], {style: {gap: '8px'}});
  mappingsSection.appendChild(ui.h3('Mappings to Apply'));

  if (sourceMappings.size > 0) {
    const mappingsList = ui.divV([], 'mapping-summary-list');
    sourceMappings.forEach((oldName, newName) => {
      const row = ui.divH([
        ui.divText(oldName, 'mapping-summary-old-name'),
        ui.iconFA('long-arrow-alt-right', null, 'maps to'),
        ui.divText(newName, 'mapping-summary-new-name'),
        ui.iconFA('times', () => {
          onUndo(newName);
          dialog.close();
          grok.shell.info(`Mapping for '${newName}' has been reset. Open "Manage Mappings" again to proceed.`);
        }, 'Undo this mapping definition')
      ], 'mapping-summary-row');
      mappingsList.appendChild(row);
    });
    mappingsSection.appendChild(mappingsList);
  } else {
    mappingsSection.appendChild(ui.divText('No mappings defined yet. Use the validation panel to map columns on the active plate.', 'info-message'));
  }

  // --- Section 2: Target Plates ---
  const platesSection = ui.divV([], {style: {gap: '8px'}});
  platesSection.appendChild(ui.h3('Apply to Plates'));

  const plateCheckboxes: DG.InputBase<boolean>[] = [];
  const plateSelectionHost = ui.divV([], 'plate-selection-host');

  allPlates.forEach((plateFile, idx) => {
    // A plate is considered "mapped" if all source mappings are present in its reconciliation map.
    const isCurrentlyMapped = [...sourceMappings.keys()].every((key) => plateFile.reconciliationMap.has(key));
    const cb = ui.input.bool(plateFile.plate.barcode ?? `Plate ${idx + 1}`, {value: isCurrentlyMapped});
    plateCheckboxes.push(cb);
    plateSelectionHost.appendChild(cb.root);
  });

  const allPlatesCheckbox = ui.input.bool('Select All', {value: false, onValueChanged: (v) => {
    plateCheckboxes.forEach((cb) => cb.value = v);
  }});
  allPlatesCheckbox.value = plateCheckboxes.every((cb) => cb.value);

  platesSection.appendChild(allPlatesCheckbox.root);
  platesSection.appendChild(plateSelectionHost);

  // --- Dialog Assembly ---
  content.appendChild(mappingsSection);
  if (plateCheckboxes.length > 0)
    content.appendChild(platesSection);


  dialog.add(content);
  dialog.onOK(() => {
    const selectedIndexes = plateCheckboxes
      .map((cb, idx) => cb.value ? idx : -1)
      .filter((idx) => idx !== -1);

    if (sourceMappings.size > 0)
      onSync(sourceMappings, selectedIndexes);
  });

  dialog.show();
}

function createAnalysisSkeleton(): HTMLElement {
  const curveSvg = `
  <svg width="200" height="150" viewBox="0 0 200 150" xmlns="http://www.w3.org/2000/svg">
    <style>
      .axis { stroke: #e0e0e0; stroke-width: 2; }
      .curve { fill: none; stroke: #d0d0d0; stroke-width: 3; stroke-dasharray: 4 4; }
      .grid-line { stroke: #f0f0f0; stroke-width: 1; }
    </style>
    <line x1="20" y1="20" x2="20" y2="130" class="axis" />
    <line x1="20" y1="130" x2="180" y2="130" class="axis" />
    <line x1="60" y1="130" x2="60" y2="20" class="grid-line" />
    <line x1="100" y1="130" x2="100" y2="20" class="grid-line" />
    <line x1="140" y1="130" x2="140" y2="20" class="grid-line" />
    <line x1="20" y1="95" x2="180" y2="95" class="grid-line" />
    <line x1="20" y1="60" x2="180" y2="60" class="grid-line" />
    <path d="M 30 25 Q 80 30, 100 80 T 170 125" class="curve"/>
  </svg>`;
  const svgDiv = ui.div([]);
  svgDiv.innerHTML = curveSvg;

  const message = ui.divText('To see dose-response curves, ensure your CSV has "SampleID", "Concentration", and "Activity" columns.');

  const skeleton = ui.divV([
    svgDiv,
    message,
  ], 'drc-skeleton');
  return skeleton;
}

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
  let plateNumberColumn: string | null = 'Plate Number';
  const plateIdentifierHost = ui.div([]);
  const replicateIdentifierHost = ui.div([]);
  const plateNumberHost = ui.div([]);


  let reprocessPlates: () => Promise<void>;

  const updatePlateIdentifierControl = () => {
    ui.empty(plateIdentifierHost);
    if (!sourceDataFrame) return;

    const choiceInput = ui.input.choice('Plate Index', {
      value: plateIdentifierColumn,
      items: [null, ...sourceDataFrame.columns.names()],
      onValueChanged: (newColumn: string | null) => {
        plateIdentifierColumn = newColumn;
        reprocessPlates();
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
      items: [null, ...sourceDataFrame.columns.names()],
      onValueChanged: (newColumn: string | null) => {
        replicateIdentifierColumn = newColumn;
        reprocessPlates();
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

  /** Applies mapping to the active plate ONLY */
  const handleMapping = (currentColName: string, templatePropName: string) => {
    const state = templateState.get(plateTemplate.id);
    if (!state || !state.plates[state.activePlateIdx]) return;

    const activeFile = state.plates[state.activePlateIdx];
    const plate = activeFile.plate;

    const columnToRename = plate.data.col(currentColName);
    if (!columnToRename) {
      console.warn(`Column '${currentColName}' not found for renaming.`);
      return;
    }

    columnToRename.name = templatePropName;
    activeFile.reconciliationMap.set(templatePropName, currentColName);
    grok.shell.info(`Mapped column '${currentColName}' to '${templatePropName}'.`);
    setTemplate(plateTemplate);
  };

  /** Undoes mapping for the active plate ONLY */
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
      } else { console.warn(`Column '${mappedField}' not found for reverting`); }
    } else { console.warn(`No original mapping found for '${mappedField}'`); }
  };

  /** NEW: Synchronizes mappings for all plates based on dialog selection. */
  const syncMappingsWithSelection = (sourceMappings: Map<string, string>, selectedIndexes: number[]) => {
    const state = templateState.get(plateTemplate.id);
    if (!state) return;

    const selectedIndexesSet = new Set(selectedIndexes);
    let platesModified = 0;

    state.plates.forEach((plateFile, index) => {
      const shouldBeMapped = selectedIndexesSet.has(index);
      const isCurrentlyMapped = [...sourceMappings.keys()].every((key) => plateFile.reconciliationMap.has(key));

      if (shouldBeMapped && !isCurrentlyMapped) {
        // APPLY MAPPING
        sourceMappings.forEach((oldName, newName) => {
          const columnToRename = plateFile.plate.data.col(oldName);
          if (columnToRename && oldName !== newName) {
            columnToRename.name = newName;
            plateFile.reconciliationMap.set(newName, oldName);
          }
        });
        platesModified++;
      } else if (!shouldBeMapped && isCurrentlyMapped) {
        // UNDO MAPPING
        sourceMappings.forEach((oldName, newName) => {
          const columnToRevert = plateFile.plate.data.col(newName);
          // Ensure we are reverting the correct mapping
          if (columnToRevert && plateFile.reconciliationMap.get(newName) === oldName) {
            columnToRevert.name = oldName;
            plateFile.reconciliationMap.delete(newName);
          }
        });
        platesModified++;
      }
    });

    if (platesModified > 0)
      grok.shell.info(`Synchronized mappings for ${platesModified} plate(s).`);

    setTemplate(plateTemplate); // Refresh UI
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
      const validation = renderValidationResults(dummyTable, plateFile.plate, template, handleMapping, plateFile.reconciliationMap, handleUndoMapping);
      const hasConflicts = validation.conflictCount > 0;
      const hasMappings = plateFile.reconciliationMap.size > 0;
      const tabLabelText = plateFile.plate.barcode ?? `Plate ${idx + 1}`;

      const headerItems = [];

      // Add mapping indicator on the left
      if (hasMappings) {
        const mappedIcon = ui.iconFA('exchange-alt');
        const pill = ui.div(mappedIcon, 'plate-file-tab__mapped-pill');
        ui.tooltip.bind(pill, `${plateFile.reconciliationMap.size} field(s) mapped`);
        headerItems.push(pill);
      }

      const labelElement = ui.divText(tabLabelText, 'plate-file-tab__label');
      headerItems.push(labelElement);

      const controlsContainer = ui.divH([], 'plate-file-tab__controls');
      headerItems.push(controlsContainer);

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

      const header = ui.divH(headerItems, 'plate-file-tab__header');
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
      ui.empty(drcAnalysisHost);
      renderFileTabs(template, state);

      if (activeFile && activeFile.plate && activeFile.plate.data) {
        try {
          plateWidget.plate = activeFile.plate;
          const drcCurvesElement = PlateDrcAnalysis.createCurvesGrid(activeFile.plate, plateWidget);
          if (drcCurvesElement)
            drcAnalysisHost.appendChild(drcCurvesElement);
          else
            drcAnalysisHost.appendChild(createAnalysisSkeleton());

          renderValidationResults(validationHost, activeFile.plate, template, handleMapping, activeFile.reconciliationMap, handleUndoMapping);
          const plateProperties = template.plateProperties
            .filter((p) => p && p.name && p.type)
            .map((p) => DG.Property.js(p.name!, p.type! as DG.TYPE));
          if (plateProperties.length > 0)
            platePropertiesHost.appendChild(ui.input.form(activeFile.plate.details || {}, plateProperties));
        } catch (validationError) {
          console.error('Error during validation rendering:', validationError);
          validationHost.appendChild(ui.divText('Error loading validation results.', 'error-message'));
        }
      } else {
        plateWidget.plate = await createNewPlateForTemplate(plateType, template);
        const templateProperties = template.plateProperties
          .filter((p) => p && p.name && p.type)
          .map((p) => DG.Property.js(p.name!, p.type! as DG.TYPE));
        if (templateProperties.length > 0)
          platePropertiesHost.appendChild(ui.input.form({}, templateProperties));
        validationHost.appendChild(ui.divText('Import a CSV file to see validation results.', 'info-message'));
        drcAnalysisHost.appendChild(createAnalysisSkeleton());
      }

      // --- UPDATED: Well Properties Header with Icon Button ---
      const manageMappingsButton = ui.button(ui.divH([ui.iconFA('exchange-alt'), ui.span(['Manage Mappings'])]), () => {
        if (state && activeFile && activeFile.reconciliationMap.size > 0)
          showMultiPlateMappingDialog(state.plates, activeFile.reconciliationMap, syncMappingsWithSelection, handleUndoMapping);
        else
          grok.shell.warning('Please import a plate and define at least one mapping on the active plate first.');
      }, 'Synchronize field mappings across multiple plates');
      manageMappingsButton.classList.add('manage-mappings-button');


      if (state && state.plates.length > 0) {
        const badges = ui.divH([]);
        let totalConflicts = 0;
        state.plates.forEach((pf) => {
          const validation = renderValidationResults(ui.div(), pf.plate, template, ()=>{}, pf.reconciliationMap, ()=>{});
          totalConflicts += validation.conflictCount;
        });

        if (totalConflicts > 0) {
          const conflictBadge = ui.span([`${totalConflicts}`], 'ui-badge-red');
          ui.tooltip.bind(conflictBadge, `${totalConflicts} unresolved fields across all plates`);
          badges.appendChild(conflictBadge);
        }

        const totalMappings = state.plates.reduce((sum, pf) => sum + pf.reconciliationMap.size, 0);
        if (totalMappings > 0) {
          const mappingBadge = ui.span([`${totalMappings}`], 'ui-badge-blue');
          ui.tooltip.bind(mappingBadge, `${totalMappings} total mappings applied across all plates`);
          badges.appendChild(mappingBadge);
        }
        if (badges.children.length > 0)
          manageMappingsButton.prepend(badges);
      }

      const header = ui.h2('Well Properties');
      const templateNameSpan = ui.span([template.name], 'well-props-header__template-name');
      const headerContainer = ui.divH([header, templateNameSpan, ui.divH([manageMappingsButton], {style: {marginLeft: 'auto'}})], 'space-between-center');
      wellPropsHeaderHost.appendChild(headerContainer);
      // --- End of updated header section ---

      plateWidget.refresh();
      plateWidget.updateRoleSummary();
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
        updatePlateNumberControl();
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
    plateNumberHost,
  ], 'plate-import-container');


  const plateTypeSelector = ui.input.choice('Plate Type', {value: plateType.name, items: plateTypes.map((pt) => pt.name), onValueChanged: (v) => {
    plateType = plateTypes.find((pt) => pt.name === v)!; setTemplate(plateTemplate);
  }});
  const plateTemplateSelector = ui.input.choice('Template', {value: plateTemplate.name, items: plateTemplates.map((pt) => pt.name), onValueChanged: (v) => setTemplate(plateTemplates.find((pt) => pt.name === v)!)});

  // --- Left Panel with Collapsible Sections ---
  const templateIcon = ui.iconFA('file-alt', null, 'Template-defined properties');
  templateIcon.classList.add('legend-icon', 'legend-icon-template');
  const templateHeader = ui.divH([templateIcon, ui.h2('Template')], {style: {alignItems: 'center', gap: '8px', flexGrow: '1'}});
  const templateContent = ui.divV([
    plateTypeSelector.root,
    plateTemplateSelector.root,
  ], 'left-panel-section-content');
  const templatePanel = createCollapsiblePanel(templateHeader, templateContent, true);

  const platePropsHeader = ui.h2('Plate Properties', {style: {flexGrow: '1'}});
  const platePropsContent = platePropertiesHost;
  platePropsContent.classList.add('left-panel-section-content');
  const platePropsPanel = createCollapsiblePanel(platePropsHeader, platePropsContent, true);

  const wellPropsContent = validationHost;
  wellPropsContent.classList.add('left-panel-section-content');
  // wellPropsHeaderHost is the complex header that is populated dynamically
  const wellPropsPanel = createCollapsiblePanel(wellPropsHeaderHost, wellPropsContent, true);

  const roleAssignmentHeader = ui.h2('Role Assignment', {style: {flexGrow: '1'}});
  const roleAssignmentContent = ui.divV([
    selectionStatusDiv,
    roleInput.root,
    assignButton,
  ], 'left-panel-section-content');
  const roleAssignmentPanel = createCollapsiblePanel(roleAssignmentHeader, roleAssignmentContent, true);

  const leftPanel = ui.divV([
    templatePanel,
    platePropsPanel,
    wellPropsPanel,
    roleAssignmentPanel,
  ], 'create-plate-view__left-panel');
  // --- End of Left Panel ---

  const rightPanelTop = ui.divV([
    importContainer,
    tabHeaderContainer,
    plateWidget.root,
  ]);
  rightPanelTop.style.flexShrink = '0';

  const analysisIcon = ui.iconFA('chart-line');
  analysisIcon.classList.add('legend-icon', 'legend-icon-analysis');
  const analysisTitle = ui.divH([analysisIcon, ui.h2('Analyses')], {style: {alignItems: 'center', gap: '8px'}});
  const drcContainer = ui.divV([analysisTitle, drcAnalysisHost], 'drc-container');

  const rightPanel = ui.divV([
    rightPanelTop,
    drcContainer,
  ], 'create-plate-view__right-panel');

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
      const validation = renderValidationResults(dummyTableElement, plateToSave, plateTemplate, handleMapping, new Map(), () => {});
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
