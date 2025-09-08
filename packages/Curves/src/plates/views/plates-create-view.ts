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
// ADD THIS IMPORT
import {PlateDoseRatioAnalysis} from '../../plate/dose-ratio-analysis';

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

  const plateWidget = new PlateWidget();
  plateWidget.editable = true;
  plateWidget.plate = new Plate(plateType.rows, plateType.cols);
  const tabControl = ui.tabControl();

  tabControl.root.classList.remove('ui-box');
  tabControl.root.style.cssText = `
    width: 100%;
    flex-grow: 1;
    display: flex;
    flex-direction: column;
  `;


  const createDrcAnalysisContent = (): HTMLElement => {
    try {
      console.log('[DEBUG] createDrcAnalysisContent called');

      const activePlate = stateManager.activePlate;
      const activeIndex = stateManager.currentState?.activePlateIdx ?? -1;

      console.log('[DEBUG] activePlate:', activePlate?.plate.barcode);
      console.log('[DEBUG] activeIndex:', activeIndex);
      console.log('[DEBUG] currentState:', stateManager.currentState);

      if (!activePlate || activeIndex < 0) {
        console.log('[DEBUG] Returning skeleton - no active plate or bad index');
        return createAnalysisSkeleton('Dose Response', ['SampleID', 'Concentration', 'Activity']);
      }

      console.log('[DEBUG] About to call createAnalysisView...');


      // Get current mappings for this analysis (for now, just use empty map until you implement the state manager methods)
      // Get current mappings for this analysis (for now, just use empty map)
      const currentMappings = stateManager.getAnalysisMapping(activeIndex, 'drc');


      const handleMap = (target: string, source: string) => {
        stateManager.remapAnalysisProperty(activeIndex, 'drc', target, source);
      };

      const handleUndo = (target: string) => {
        stateManager.undoAnalysisMapping(activeIndex, 'drc', target);
      };

      return PlateDrcAnalysis.createAnalysisViewWithMapping(
        activePlate.plate,
        currentMappings,
        handleMap,
        handleUndo,
        plateWidget
      );
    } catch (e: any) {
      console.error('Error creating Dose Response view:', e);
      return ui.divText(`Error displaying Dose Response analysis: ${e.message}`, 'error-message');
    }
  };

  const createDoseRatioContent = (): HTMLElement => {
    console.log('[DEBUG] PlateDoseRatioAnalysis methods:', Object.getOwnPropertyNames(PlateDoseRatioAnalysis));
    console.log('[DEBUG] createAnalysisView exists:', typeof PlateDoseRatioAnalysis.createAnalysisView);
    try {
      const activePlate = stateManager.activePlate;
      const activeIndex = stateManager.currentState?.activePlateIdx ?? -1;

      if (!activePlate || activeIndex < 0)
        return createAnalysisSkeleton('Dose Ratio', ['Agonist_Concentration_M', 'Antagonist_Concentration_M', 'Percent_Inhibition']);


      // Get current mappings for this analysis (for now, just use empty map)
      const currentMappings = stateManager.getAnalysisMapping(activeIndex, 'doseRatio');


      const handleMap = (target: string, source: string) => {
        stateManager.remapAnalysisProperty(activeIndex, 'doseRatio', target, source);
      };

      const handleUndo = (target: string) => {
      // For now, just use the regular mapping
        stateManager.undoMapping(activeIndex, target);
      };

      return PlateDoseRatioAnalysis.createAnalysisView(
        activePlate.plate,
        currentMappings,
        handleMap,
        handleUndo
      );
    } catch (e: any) {
      console.error('Error creating Dose Ratio view:', e);
      return ui.divText(`Error displaying Dose Ratio analysis: ${e.message}`, 'error-message');
    }
  };

  tabControl.addPane('Wells', () => {
    plateWidget.root.style.width = '100%';
    plateWidget.root.style.height = '100%';
    return plateWidget.root;
  });
  tabControl.addPane('Dose Response', () => createDrcAnalysisContent());
  tabControl.addPane('Dose Ratio', () => createDoseRatioContent());
  ui.tools.waitForElementInDom(tabControl.root).then(() => {
    tabControl.panes.forEach((pane) => {
      pane.content.style.cssText = `
        width: 100%;
        height: 100%;
        max-width: none; /* Override any max-width constraints */
        flex: 1;
        display: flex;
        flex-direction: column;
      `;
    });

    // Ensure the PlateWidget's root itself also conforms. This might be redundant
    // but ensures consistency.
    plateWidget.root.style.width = '100%';
    plateWidget.root.style.height = '100%';
  });


  const tabContentHost = tabControl.root.querySelector('.d4-tab-host');
  if (tabContentHost) {
    (tabContentHost as HTMLElement).style.flexGrow = '1';
    (tabContentHost as HTMLElement).style.width = '100%';
    (tabContentHost as HTMLElement).style.display = 'flex';
  }
  // tabControl.root.style.width = '100%';
  // tabControl.root.style.flexGrow = '1';
  // tabControl.root.style.display = 'flex';
  // tabControl.root.style.flexDirection = 'column';

  const plateGridManager = new PlateGridManager(stateManager);

  const handleManageMappings = () => {
    const state = stateManager.currentState;
    const activePlate = stateManager.activePlate;

    if (state && activePlate) {
    // 1. Create the mappings map from the plate's aliases (the new source of truth)
      const currentMappings = new Map<string, string>();
      for (const layer of activePlate.plate.getLayerNames()) {
        const aliases = activePlate.plate.getAliases(layer);
        for (const alias of aliases)
          currentMappings.set(alias, layer);
      }

      // 2. Check if there are any mappings to sync
      if (currentMappings.size > 0) {
        plateGridManager.showMultiPlateMappingDialog({
          allPlates: state.plates,
          sourceMappings: currentMappings, // 3. Pass the new map
          onSync: (mappings, selectedIndexes) => {
            stateManager.syncMappings(mappings, selectedIndexes);
          },
          onUndo: (mappedField) => {
            if (state.activePlateIdx >= 0)
              stateManager.undoMapping(state.activePlateIdx, mappedField);
          },
        });
      } else {
        grok.shell.warning('Please define at least one mapping on the active plate first.');
      }
    } else {
      grok.shell.warning('Please import a plate first.');
    }
  }; ;

  const templatePanel = new TemplatePanel(
    stateManager,
    plateWidget,
    () => handleManageMappings()
  );

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

    // Update this condition to include analysis mapping changes
    if (event.type === 'mapping-changed' ||
    event.type === 'plate-selected' ||
    event.type === 'analysis-mapping-changed' ||
    event.type === 'plate-added') { // Add this line
      console.log('[DEBUG] Refreshing analysis views due to state change');

      const drcPane = tabControl.getPane('Dose Response');
      if (drcPane) {
        ui.empty(drcPane.content);
        drcPane.content.appendChild(createDrcAnalysisContent());
      }

      const doseRatioPane = tabControl.getPane('Dose Ratio');
      if (doseRatioPane) {
        ui.empty(doseRatioPane.content);
        doseRatioPane.content.appendChild(createDoseRatioContent());
      }
    }
  });

  const rightPanel = ui.divV([
    plateGridManager.root,
    tabControl.root,
  ], 'create-plate-view__right-panel');

  const mainLayout = ui.divH(
    [templatePanel.root, rightPanel],
    'create-plate-view__main-layout'
  );

  view.root.appendChild(mainLayout);

  tabControl.currentPane = tabControl.getPane('Wells');

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
    if (grok.shell.v.name !== view.name)
      return;

    const file = args.args.file;
    if (!file) return;

    if (!file.name.toLowerCase().endsWith('.csv')) {
      grok.shell.warning(`Only .csv files can be dropped here. To open other files, navigate away from the 'Create Plate' view.`);
      args.preventDefault(); // Prevent opening other files even if they are valid Datagrok tables
      return;
    }

    args.preventDefault();

    handleDroppedFile(file);
  });
  subscriptions.push(fileImportSub);


  const originalDetach = view.detach.bind(view);
  view.detach = () => {
    subscriptions.forEach((sub) => sub.unsubscribe());
    stateManager.destroy();
    templatePanel.destroy();
    plateGridManager.destroy();
    originalDetach();
  };

  return view;
}

function createAnalysisSkeleton(analysisType: string, requiredColumns: string[] = []): HTMLElement {
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

  const messageText = `To see ${analysisType.toLowerCase()} curves, please map the following required properties: ${requiredColumns.join(', ')}.`;
  const message = ui.divText(messageText);

  const skeleton = ui.divV([svgDiv, message], 'drc-skeleton');
  skeleton.style.cssText = `
    display: flex; flex-direction: column; align-items: center; justify-content: center;
    width: 100%; flex-grow: 1; text-align: center;
    border: 1px dashed var(--grey-2); border-radius: 8px; background-color: var(--grey-0);`;

  return skeleton;
}
