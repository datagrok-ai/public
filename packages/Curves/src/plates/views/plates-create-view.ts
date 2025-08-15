import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {PlateWidget} from '../../plate/plate-widget';
import {Plate} from '../../plate/plate';
import {plateTemplates, plateTypes, savePlate, savePlateAsTemplate, initPlates} from '../plates-crud';
import {PlateStateManager} from './shared/plate-state-manager';
import {TemplatePanel} from './components/template-panel/template-panel';
import {PlateGridManager} from './components/plate-grid-manager/plate-grid-manager';
import {PlateDrcAnalysis} from '../../plate/plate-drc-analysis';
import {renderValidationResults} from './plates-validation-panel';

import './components/template-panel/template-panel.css';
import './components/plate-grid-manager/plate-grid-manager.css';
import './components/plate-analysis-panel/plate-analysis-panel.css';

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

  // Create plate widget
  const plateWidget = new PlateWidget();
  plateWidget.editable = true;
  plateWidget.plate = new Plate(plateType.rows, plateType.cols);

  // DEBUG: Log initial plate setup
  console.log('[DEBUG] Initial plate widget setup with plate:', plateWidget.plate);

  // Create the DRC analysis container
  const drcAnalysisContainer: HTMLElement | null = null;

  // Create toggle buttons for Wells/Analysis view
  const wellsButton = ui.button('Wells', () => toggleView('wells'));
  const analysisButton = ui.button('Analysis', () => toggleView('analysis'));

  wellsButton.classList.add('view-toggle-button', 'active');
  analysisButton.classList.add('view-toggle-button');

  const toggleButtonsContainer = ui.divH([wellsButton, analysisButton], 'view-toggle-container');

  // Create the bottom panel that will switch between plate widget and analysis
  const bottomPanel = ui.div([], 'bottom-panel-container');
  bottomPanel.style.width = '100%';
  bottomPanel.style.flexGrow = '1';
  bottomPanel.style.display = 'flex';
  bottomPanel.style.minHeight = '0';

  function toggleView(viewType: 'wells' | 'analysis') {
    console.log(`[DEBUG] Toggling to view: ${viewType}`);

    if (viewType === 'wells') {
      wellsButton.classList.add('active');
      analysisButton.classList.remove('active');

      ui.empty(bottomPanel);

      // For plate widget, we need to ensure it takes full width
      plateWidget.root.style.width = '100%';
      plateWidget.root.style.height = '100%';
      plateWidget.root.style.display = 'flex';
      plateWidget.root.style.flexDirection = 'column';

      // Make sure the internal grid takes full width
      if (plateWidget.grid && plateWidget.grid.root)
        plateWidget.grid.root.style.width = '100%';


      // Fix the tabs container to take full width
      if (plateWidget.tabsContainer) {
        plateWidget.tabsContainer.style.width = '100%';
        plateWidget.tabsContainer.style.display = 'flex';
        plateWidget.tabsContainer.style.flexDirection = 'column';
      }

      bottomPanel.appendChild(plateWidget.root);

      // Force refresh after DOM update
      setTimeout(() => {
        plateWidget.refresh();
        plateWidget.grid?.invalidate();
      }, 0);
    } else {
      wellsButton.classList.remove('active');
      analysisButton.classList.add('active');

      ui.empty(bottomPanel);

      // Create or update DRC analysis
      const activePlate = stateManager.activePlate;
      if (activePlate) {
        // Get the mapped data which applies the reconciliation map
        const mappedData = stateManager.getActivePlateMappedData();
        console.log('[DEBUG] Mapped data columns:', mappedData?.columns.names());

        if (mappedData) {
          // Check if all required columns exist after mapping
          const hasActivity = mappedData.columns.contains('Activity');
          const hasConcentration = mappedData.columns.contains('Concentration');
          const hasSampleID = mappedData.columns.contains('SampleID');

          console.log('[DEBUG] DRC requirements check:', {
            hasActivity,
            hasConcentration,
            hasSampleID,
            mappings: activePlate.reconciliationMap
          });

          if (hasActivity && hasConcentration && hasSampleID) {
            // Create a new plate with the mapped data for DRC analysis
            const mappedPlate = new Plate(activePlate.plate.rows, activePlate.plate.cols);
            mappedPlate.data = mappedData;
            mappedPlate.barcode = activePlate.plate.barcode;
            mappedPlate.id = activePlate.plate.id;

            // Use the mapped plate for DRC analysis with the correct column names
            const drcGrid = PlateDrcAnalysis.createCurvesGrid(
              mappedPlate,
              plateWidget,
              {
                roleName: 'SampleID',
                concentrationName: 'Concentration',
                valueName: 'Activity',
                normalize: true,
                controlColumns: ['High Control', 'Low Control'],
                autoFilterOutliers: true,
                categorizeFormula: '${rSquared} > 0.8 && ${Hill} > 0.25 && ${Max} > 80 && ${Max} < 120',
                statisticsColumns: ['Min', 'Max', 'rSquared', 'Hill', 'IC50', 'AUC'],
              },
              'csv'
            );

            if (drcGrid) {
              drcGrid.style.width = '100%';
              drcGrid.style.height = '100%';
              bottomPanel.appendChild(drcGrid);
            } else {
              console.log('[DEBUG] DRC grid creation failed despite having required columns');
              bottomPanel.appendChild(createAnalysisSkeleton());
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
            bottomPanel.appendChild(message);
          }
        } else {
          bottomPanel.appendChild(createAnalysisSkeleton());
        }
      } else {
        bottomPanel.appendChild(createAnalysisSkeleton());
      }
    }
  }

  const plateGridManager = new PlateGridManager(stateManager);

  // Fix plate grid manager to take full width
  plateGridManager.root.style.width = '100%';
  plateGridManager.root.style.display = 'flex';
  plateGridManager.root.style.flexDirection = 'column';

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

    // Force grid invalidation to ensure visual update
    if (plateWidget.grid)
      plateWidget.grid.invalidate();


    // If analysis view is active, refresh it when mappings change
    if (analysisButton.classList.contains('active') &&
        (event.type === 'mapping-changed' || event.type === 'plate-selected')) {
      console.log('[DEBUG] Refreshing DRC analysis view due to mapping change');
      toggleView('analysis');
    }
  });

  // Create the right panel with proper flex distribution
  const rightPanel = ui.divV([
    plateGridManager.root,
    toggleButtonsContainer,
    bottomPanel,
  ], 'create-plate-view__right-panel');

  const mainLayout = ui.divH(
    [templatePanel.root, rightPanel],
    'create-plate-view__main-layout'
  );

  view.root.appendChild(mainLayout);

  // Initialize with wells view
  toggleView('wells');

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

  const originalDetach = view.detach.bind(view);
  view.detach = () => {
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
