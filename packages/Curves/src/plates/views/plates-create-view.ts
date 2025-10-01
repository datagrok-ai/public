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
import {PlateDoseRatioAnalysis} from '../../plate/dose-ratio-analysis';
import {MAPPING_SCOPES} from './shared/scopes';
import {PlateQpcrAnalysis} from '../../plate/relative-qpcr-analysis';
import {AnalysisManager} from '../../plate/analyses/analysis-manager';

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
  plateWidget.root.style.cssText = `
    width: 100%;
    height: 100%;
    display: flex;
    flex-direction: column;
  `;
  // plateWidget.roleSummaryDiv.style.display = 'none';

  const tabControl = ui.tabControl();
  tabControl.root.classList.remove('ui-box');
  tabControl.root.style.cssText = `
    width: 100%;
    flex-grow: 1;
    display: flex;
    flex-direction: column;
  `;

  // const createDrcAnalysisContent = (): HTMLElement => {
  //   try {
  //     const activePlate = stateManager.activePlate;
  //     const activeIndex = stateManager.currentState?.activePlateIdx ?? -1;


  //     if (!activePlate || activeIndex < 0)
  //       return createAnalysisSkeleton('Dose Response', ['SampleID', 'Concentration', 'Activity']);

  //     const currentMappings = stateManager.getScopedMapping(activeIndex, MAPPING_SCOPES.DRC);

  //     const handleMap = (target: string, source: string) => {
  //       stateManager.remapScopedProperty(activeIndex, MAPPING_SCOPES.DRC, target, source);
  //     };

  //     const handleUndo = (target: string) => {
  //       stateManager.undoScopedMapping(activeIndex, MAPPING_SCOPES.DRC, target);
  //     };

  //     return PlateDrcAnalysis.createAnalysisViewWithMapping(
  //       activePlate.plate,
  //       currentMappings,
  //       handleMap,
  //       handleUndo,
  //       plateWidget
  //     );
  //   } catch (e: any) {
  //     console.error('Error creating Dose Response view:', e);
  //     return ui.divText(`Error displaying Dose Response analysis: ${e.message}`, 'error-message');
  //   }
  // };

  // const createDoseRatioContent = (): HTMLElement => {
  //   try {
  //     const activePlate = stateManager.activePlate;
  //     const activeIndex = stateManager.currentState?.activePlateIdx ?? -1;

  //     if (!activePlate || activeIndex < 0)
  //       return createAnalysisSkeleton('Dose Ratio', ['Agonist_Concentration_M', 'Antagonist_Concentration_M', 'Percent_Inhibition']);

  //     const currentMappings = stateManager.getScopedMapping(activeIndex, MAPPING_SCOPES.DOSE_RATIO);

  //     const handleMap = (target: string, source: string) => {
  //       stateManager.remapScopedProperty(activeIndex, MAPPING_SCOPES.DOSE_RATIO, target, source);
  //     };

  //     const handleUndo = (target: string) => {
  //       stateManager.undoScopedMapping(activeIndex, MAPPING_SCOPES.DOSE_RATIO, target);
  //     };

  //     return PlateDoseRatioAnalysis.createAnalysisView(
  //       activePlate.plate,
  //       currentMappings,
  //       handleMap,
  //       handleUndo
  //     );
  //   } catch (e: any) {
  //     console.error('Error creating Dose Ratio view:', e);
  //     return ui.divText(`Error displaying Dose Ratio analysis: ${e.message}`, 'error-message');
  //   }
  // };

  // const createQpcrContent = (): HTMLElement => {
  //   try {
  //     const activePlate = stateManager.activePlate;
  //     const activeIndex = stateManager.currentState?.activePlateIdx ?? -1;

  //     if (!activePlate || activeIndex < 0)
  //       return createAnalysisSkeleton('qPCR Analysis', ['Ct']);

  //     const currentMappings = stateManager.getScopedMapping(activeIndex, MAPPING_SCOPES.QPCR);

  //     const handleMap = (target: string, source: string) => {
  //       stateManager.remapScopedProperty(activeIndex, MAPPING_SCOPES.QPCR, target, source);
  //     };

  //     const handleUndo = (target: string) => {
  //       stateManager.undoScopedMapping(activeIndex, MAPPING_SCOPES.QPCR, target);
  //     };

  //     // --- START: NEW TRIGGER FUNCTION ---
  //     const handleDataChanged = () => {
  //       stateManager.notifyPlateDataChanged();
  //     };
  //     // --- END: NEW TRIGGER FUNCTION ---

  //     return PlateQpcrAnalysis.createAnalysisView(
  //       activePlate.plate,
  //       currentMappings,
  //       handleMap,
  //       handleUndo,
  //       handleDataChanged // <-- PASS THE NEW TRIGGER HERE
  //     );
  //   } catch (e: any) {
  //     console.error('Error creating qPCR view:', e);
  //     return ui.divText(`Error displaying qPCR analysis: ${e.message}`, 'error-message');
  //   }
  // };
  // ; ;

  tabControl.addPane('Plate View', () => plateWidget.root);
  const analysisManager = AnalysisManager.instance;
  for (const analysis of analysisManager.analyses) {
    const createTabContent = () => {
      try {
        const activePlate = stateManager.activePlate;
        const activeIndex = stateManager.currentState?.activePlateIdx ?? -1;

        if (!activePlate || activeIndex < 0)
          return createAnalysisSkeleton(analysis.friendlyName, analysis.getRequiredFields().map((f) => f.name));

        // The analysis's unique `name` (e.g., 'DRC') is used as the scope for its mappings.
        const currentMappings = stateManager.getScopedMapping(activeIndex, analysis.name);

        const handleMap = (target: string, source: string) => {
          stateManager.remapScopedProperty(activeIndex, analysis.name, target, source);
        };
        const handleUndo = (target: string) => {
          stateManager.undoScopedMapping(activeIndex, analysis.name, target);
        };

        // Each analysis class is now responsible for creating its own view.
        return analysis.createView(activePlate.plate, plateWidget, currentMappings, handleMap, handleUndo);
      } catch (e: any) {
        console.error(`Error creating ${analysis.friendlyName} view:`, e);
        return ui.divText(`Error displaying ${analysis.friendlyName} analysis: ${e.message}`, 'error-message');
      }
    };

    tabControl.addPane(analysis.friendlyName, createTabContent);
  }
  // tabControl.addPane('Dose Response', () => createDrcAnalysisContent());
  // tabControl.addPane('Dose Ratio', () => createDoseRatioContent());
  // tabControl.addPane('qPCR Analysis', () => createQpcrContent()); // 2. ADD THE NEW TAB


  ui.tools.waitForElementInDom(tabControl.root).then(() => {
    tabControl.panes.forEach((pane) => {
      pane.content.style.cssText = `
        width: 100%;
        height: 100%;
        max-width: none;
        flex: 1;
        display: flex;
        flex-direction: column;
      `;
    });
  });

  const tabContentHost = tabControl.root.querySelector('.d4-tab-host');
  if (tabContentHost) {
    (tabContentHost as HTMLElement).style.flexGrow = '1';
    (tabContentHost as HTMLElement).style.width = '100%';
    (tabContentHost as HTMLElement).style.display = 'flex';
  }

  const plateGridManager = new PlateGridManager(stateManager);

  const templatePanel = new TemplatePanel(
    stateManager,
    plateWidget,
    () => { console.warn('[DEBUG] TemplatePanel is not implemented'); }
  );

  stateManager.onStateChange.subscribe(async (event) => {
    // This part updates the main plate widget itself
    const activePlate = stateManager.activePlate;
    if (activePlate)
      plateWidget.plate = activePlate.plate;
    else
      plateWidget.plate = await stateManager.getOrCreateDefaultPlate();

    plateWidget.refresh();
    if (plateWidget.grid)
      plateWidget.grid.invalidate();

    // This part handles refreshing the content of the currently active analysis tab
    const eventsThatRequireRefresh = [
      'mapping-changed', 'plate-selected', 'analysis-mapping-changed',
      'plate-data-changed', 'plate-added'
    ];

    if (eventsThatRequireRefresh.includes(event.type)) {
      const currentPane = tabControl.currentPane;

      // Only refresh if the active tab is an analysis tab (not 'Plate View')
      if (currentPane && currentPane.name !== 'Plate View') {
        // 1. Find the correct analysis object using the tab's friendly name.
        const analysis = analysisManager.byFriendlyName(currentPane.name);

        if (analysis) {
          // 2. Re-create the content for that analysis using the latest state.
          // This logic is the same as in the initial `createTabContent` factory.
          const activeIndex = stateManager.currentState?.activePlateIdx ?? -1;
          const plateToAnalyze = stateManager.activePlate;
          const handleRerender = () => stateManager.notifyPlateDataChanged();


          let newContent: HTMLElement;
          if (!plateToAnalyze || activeIndex < 0) {
            newContent = createAnalysisSkeleton(analysis.friendlyName, analysis.getRequiredFields().map((f) => f.name));
          } else {
            const currentMappings = stateManager.getScopedMapping(activeIndex, analysis.name);
            const handleMap = (target: string, source: string) => stateManager.remapScopedProperty(activeIndex, analysis.name, target, source);
            const handleUndo = (target: string) => stateManager.undoScopedMapping(activeIndex, analysis.name, target);

            newContent = analysis.createView(plateToAnalyze.plate, plateWidget, currentMappings, handleMap, handleUndo, handleRerender);
          }

          // 3. Replace the old content with the newly generated content.
          ui.empty(currentPane.content);
          currentPane.content.appendChild(newContent);
        }
      }
    }
  }); ;

  const rightPanel = ui.divV([
    plateGridManager.root,
    tabControl.root,
  ], 'create-plate-view__right-panel');

  const mainLayout = ui.divH(
    [templatePanel.root, rightPanel],
    'create-plate-view__main-layout'
  );

  view.root.appendChild(mainLayout);

  // START of CHANGE: Set the new 'Plate View' as the default tab.
  tabControl.currentPane = tabControl.getPane('Plate View');
  // END of CHANGE

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
      args.preventDefault();
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
