/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {PlateWidget} from '../../plate/plate-widget';
import {Plate} from '../../plate/plate';
import {plateTemplates, plateTypes, savePlate, savePlateAsTemplate, initPlates} from '../plates-crud';
import {PlateStateManager} from './shared/plate-state-manager';
import {TemplatePanel} from './components/template-panel/template-panel';
import {PlateDrcAnalysis} from '../../plate/plate-drc-analysis';
import {renderValidationResults} from './plates-validation-panel';
import {Subscription} from 'rxjs';

import './components/template-panel/template-panel.css';
import './components/plate-grid-manager/plate-grid-manager.css';
import './components/plate-analysis-panel/plate-analysis-panel.css';
import {PlateGridManager} from './components/plate-grid-manager/plate-grid-manager';

export function createPlatesView(): DG.View {
  const view = DG.View.create();
  view.name = 'Create Plate';
  view.root.classList.add('create-plate-view');

  const plateType = plateTypes[0];
  const plateTemplate = plateTemplates[0];

  if (!plateTemplate) {
    grok.shell.error('No plate templates found.');
    return view;
  }

  const stateManager = new PlateStateManager(plateTemplate, plateType);
  const subscriptions: Subscription[] = [];

  // Create plate widget
  const plateWidget = new PlateWidget();
  plateWidget.editable = true;
  plateWidget.plate = new Plate(plateType.rows, plateType.cols);

  // Use Datagrok's native TabControl for Wells/Analysis
  const tabControl = ui.tabControl();

  // Function to create the DRC analysis content
  const createAnalysisContent = () => {
    const activePlate = stateManager.activePlate;
    if (activePlate) {
      const mappedData = stateManager.getActivePlateMappedData();
      console.log('[DEBUG] Mapped data columns:', mappedData?.columns.names());
      if (mappedData) {
        const hasActivity = mappedData.columns.contains('Activity');
        const hasConcentration = mappedData.columns.contains('Concentration');
        const hasSampleID = mappedData.columns.contains('SampleID');
        console.log('[DEBUG] DRC requirements check:', {
          hasActivity,
          hasConcentration,
          hasSampleID,
          mappings: activePlate.reconciliationMap,
        });

        if (hasActivity && hasConcentration && hasSampleID) {
          const mappedPlate = new Plate(activePlate.plate.rows, activePlate.plate.cols);
          mappedPlate.data = mappedData;
          mappedPlate.barcode = activePlate.plate.barcode;
          mappedPlate.id = activePlate.plate.id;

          const drcGrid = PlateDrcAnalysis.createCurvesGrid(
            mappedPlate, plateWidget, {}, 'csv'
          );
          if (drcGrid) {
            drcGrid.style.width = '100%';
            drcGrid.style.height = '100%';
            return drcGrid;
          } else {
            console.log('[DEBUG] DRC grid creation failed despite having required columns');
            return createAnalysisSkeleton();
          }
        } else {
          console.log('[DEBUG] Missing required columns for DRC analysis');
          const message = ui.divV([
            createAnalysisSkeleton(),
            ui.divText(
              `Missing required columns: ${[
                !hasActivity && 'Activity',
                !hasConcentration && 'Concentration',
                !hasSampleID && 'SampleID'
              ].filter(Boolean).join(', ')}`,
              'warning-message'
            )
          ]);
          return message;
        }
      } else {
        return createAnalysisSkeleton();
      }
    } else {
      return createAnalysisSkeleton();
    }
  };

  // Add the Wells and Analysis panes
  tabControl.addPane('Wells', () => {
    // Ensure the plate widget takes full width and height
    plateWidget.root.style.width = '100%';
    plateWidget.root.style.height = '100%';
    return plateWidget.root;
  });
  tabControl.addPane('Analysis', () => createAnalysisContent());

  // Fixes for the Datagrok tabs layout
  tabControl.root.style.width = '100%';
  tabControl.root.style.flexGrow = '1';
  tabControl.root.style.display = 'flex';
  tabControl.root.style.flexDirection = 'column';

  const plateGridManager = new PlateGridManager(stateManager);

  const handleManageMappings = () => {
    const state = stateManager.currentState;
    const activePlate = stateManager.activePlate;

    if (state && activePlate && activePlate.reconciliationMap.size > 0) {
      plateGridManager.showMultiPlateMappingDialog({
        allPlates: state.plates,
        sourceMappings: activePlate.reconciliationMap,
        onSync: (mappings, selectedIndexes) => {
          stateManager.syncMappings(mappings, selectedIndexes);
        },
        onUndo: (mappedField) => {
          if (state.activePlateIdx >= 0)
            stateManager.undoMapping(state.activePlateIdx, mappedField);
        },
      });
    } else {
      grok.shell.warning('Please import a plate and define at least one mapping on the active plate first.');
    }
  };

  const templatePanel = new TemplatePanel(
    stateManager,
    plateWidget,
    () => handleManageMappings()
  );

  // Subscribe to state changes and update the plate widget
  stateManager.onStateChange.subscribe(async (event) => {
    console.log('[DEBUG] State change event:', event);
    const activePlate = stateManager.activePlate;

    if (activePlate) {
      console.log('[DEBUG] Setting active plate:', activePlate.plate.barcode);
      plateWidget.plate = activePlate.plate;
      plateWidget.refresh();
    } else {
      console.log('[DEBUG] No active plate, showing default');
      plateWidget.plate = await stateManager.getOrCreateDefaultPlate();
      plateWidget.refresh();
    }

    plateWidget.updateRoleSummary();

    if (plateWidget.grid)
      plateWidget.grid.invalidate();

    // If analysis view is active, refresh it when mappings change
    if (tabControl.currentPane.name === 'Analysis' && (event.type === 'mapping-changed' || event.type === 'plate-selected')) {
      console.log('[DEBUG] Refreshing DRC analysis view due to mapping change');
      // Re-create the analysis content to force a full refresh
      const newContent = createAnalysisContent();
      ui.empty(tabControl.getPane('Analysis').content);
      tabControl.getPane('Analysis').content.appendChild(newContent);
    }
  });

  // Create the right panel with proper flex distribution
  const rightPanel = ui.divV([
    plateGridManager.root,
    tabControl.root,
  ], 'create-plate-view__right-panel');

  const mainLayout = ui.divH(
    [templatePanel.root, rightPanel],
    'create-plate-view__main-layout'
  );

  view.root.appendChild(mainLayout);

  // Initialize with wells view
  tabControl.currentPane = tabControl.getPane('Wells');

  // Initialize with the first template
  stateManager.setTemplate(plateTemplates[0]);

  view.setRibbonPanels([[
    ui.bigButton('CREATE', async () => {
      const plateToSave = await stateManager.getOrCreateDefaultPlate();
      if (!plateToSave) {
        grok.shell.warning('No active plate to save.');
        return;
      }

      const activePlate = stateManager.activePlate;
      if (activePlate) {
        const dummyTableElement = ui.div();
        const validation = renderValidationResults(
          dummyTableElement,
          plateToSave,
          stateManager.currentTemplate,
          () => {},
          new Map(),
          () => {}
        );

        if (validation.conflictCount > 0)
          grok.shell.warning('Saving plate with unresolved validation issues.');
      }

      plateToSave.plateTemplateId = stateManager.currentTemplate.id;
      await savePlate(plateToSave);
      grok.shell.info(`Plate created: ${plateToSave.barcode}`);
    }),

    ui.button('SAVE TEMPLATE', async () => {
      const plateToSave = await stateManager.getOrCreateDefaultPlate();
      if (!plateToSave) {
        grok.shell.warning('No active plate to save as template.');
        return;
      }

      await savePlateAsTemplate(plateToSave, stateManager.currentTemplate);
      await initPlates(true);
      grok.shell.info(`Plate template updated: ${stateManager.currentTemplate.name}`);
    }),
  ]]);

  // --- NEW: Global File Import Interceptor ---
  const handleDroppedFile = async (file: File) => {
    const pi = DG.TaskBarProgressIndicator.create(`Importing ${file.name}...`);
    try {
      const csvString = await file.text();
      const df = DG.DataFrame.fromCsv(csvString);
      await stateManager.loadDataFrame(df);
      grok.shell.info(`File "${file.name}" imported successfully.`);
    } catch (error: any) {
      grok.shell.error(`Failed to import file: ${error.message}`);
    } finally {
      pi.close();
    }
  };

  const fileImportSub = grok.events.onFileImportRequest.subscribe((args) => {
    // 1. Context Check: Only handle if our view is active.
    if (grok.shell.v.name !== view.name)
      return;

    const file = args.args.file;
    if (!file) return;

    // 2. File Type Check: Ensure it's a CSV file.
    if (!file.name.toLowerCase().endsWith('.csv')) {
      // Optional: Inform user if they drop a wrong file type
      grok.shell.warning(`Only .csv files can be dropped here. To open other files, navigate away from the 'Create Plate' view.`);
      args.preventDefault(); // Prevent opening other files even if they are valid Datagrok tables
      return;
    }

    // 3. Prevent Default Behavior for this file.
    args.preventDefault();

    // 4. Route to our custom handler.
    handleDroppedFile(file);
  });
  subscriptions.push(fileImportSub);


  // --- Updated View Detach Logic ---
  const originalDetach = view.detach.bind(view);
  view.detach = () => {
    subscriptions.forEach((sub) => sub.unsubscribe()); // <-- Clean up the subscription
    stateManager.destroy();
    templatePanel.destroy();
    plateGridManager.destroy();
    originalDetach();
  };

  return view;
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

  const message = ui.divText(
    'To see dose-response curves, ensure your CSV has "SampleID", "Concentration", and "Activity" columns.'
  );

  const skeleton = ui.divV([svgDiv, message], 'drc-skeleton');

  skeleton.style.display = 'flex';
  skeleton.style.alignItems = 'center';
  skeleton.style.justifyContent = 'center';
  skeleton.style.width = '100%';
  skeleton.style.flexGrow = '1';
  skeleton.style.textAlign = 'center';
  skeleton.style.border = '1px dashed var(--grey-2)';
  skeleton.style.borderRadius = '8px';
  skeleton.style.backgroundColor = 'var(--grey-0)';

  return skeleton;
}
