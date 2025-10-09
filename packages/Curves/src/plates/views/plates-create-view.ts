/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {PlateWidget} from '../../plate/plate-widget';
import {Plate} from '../../plate/plate';
import {plateTemplates, plateTypes, savePlate, savePlateAsTemplate, initPlates} from '../plates-crud';
import {PlateStateManager} from './shared/plate-state-manager';
import {TemplatePanel} from './components/template-panel/template-panel';
import {renderValidationResults} from './plates-validation-panel';
import {Subscription} from 'rxjs';
import './components/plate-grid-manager/plate-grid-manager.css';
import {PlateGridManager} from './components/plate-grid-manager/plate-grid-manager';
import {AnalysisManager} from '../../plate/analyses/analysis-manager';
import './plates-create-view.css';

export function createPlatesView(): DG.View {
  const view = DG.View.create();
  view.name = 'Create Plate';

  view.root.classList.add('assay_plates__create-plate-view');

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
  tabControl.root.classList.add('assay_plates__create-plate-view__tab-control'); // Added class


  tabControl.addPane('Plate View', () => plateWidget.root);
  const analysisManager = AnalysisManager.instance;
  for (const analysis of analysisManager.analyses) {
    const createTabContent = () => {
      try {
        const activePlate = stateManager.activePlate;
        const activeIndex = stateManager.currentState?.activePlateIdx ?? -1;

        if (!activePlate || activeIndex < 0)
          return createAnalysisSkeleton(analysis.friendlyName, analysis.getRequiredFields().map((f) => f.name));

        const currentMappings = stateManager.getScopedMapping(activeIndex, analysis.name);

        const handleMap = (target: string, source: string) => {
          stateManager.remapScopedProperty(activeIndex, analysis.name, target, source);
        };
        const handleUndo = (target: string) => {
          stateManager.undoScopedMapping(activeIndex, analysis.name, target);
        };

        return analysis.createView(activePlate.plate, plateWidget, currentMappings, handleMap, handleUndo);
      } catch (e: any) {
        console.error(`Error creating ${analysis.friendlyName} view:`, e);
        return ui.divText(`Error displaying ${analysis.friendlyName} analysis: ${e.message}`, 'error-message');
      }
    };

    tabControl.addPane(analysis.friendlyName, createTabContent);
  }

  // The logic to apply styles via waitForElementInDom has been removed and moved to plates-create-view.css
  // The direct styling of tabContentHost has also been removed and moved to plates-create-view.css

  const plateGridManager = new PlateGridManager(stateManager);

  const templatePanel = new TemplatePanel(
    stateManager,
    plateWidget,
    () => { console.warn('[DEBUG] TemplatePanel is not implemented'); }
  );

  stateManager.onStateChange.subscribe(async (event) => {
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

      if (currentPane && currentPane.name !== 'Plate View') {
        const analysis = analysisManager.byFriendlyName(currentPane.name);

        if (analysis) {
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

          ui.empty(currentPane.content);
          currentPane.content.appendChild(newContent);
        }
      }
    }
  }); ;

  const rightPanel = ui.divV([
    plateGridManager.root,
    tabControl.root,
  ], 'assay_plates__create-plate-view__right-panel');

  const mainLayout = ui.divH(
    [templatePanel.root, rightPanel],
    'assay_plates__create-plate-view__main-layout'
  );

  view.root.appendChild(mainLayout);

  tabControl.currentPane = tabControl.getPane('Plate View');

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

  const skeleton = ui.divV([svgDiv, message], 'assay_plates__drc-skeleton');

  return skeleton;
}
