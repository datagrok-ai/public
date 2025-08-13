import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {PlateWidget} from '../../plate/plate-widget';
import {Plate} from '../../plate/plate';
import {plateTemplates, plateTypes, savePlate, savePlateAsTemplate, initPlates} from '../plates-crud';
import {PlateStateManager} from './shared/plate-state-manager';
import {TemplatePanel} from './components/template-panel/template-panel';
import {PlateTabsManager} from './components/plate-tabs-manager/plate-tabs-manager';
import {PlateDrcAnalysis} from '../../plate/plate-drc-analysis';
import {renderValidationResults} from './plates-validation-panel';

import './components/template-panel/template-panel.css';
import './components/plate-tabs-manager/plate-tabs-manager.css';
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

  const plateWidget = PlateWidget.fromPlate(new Plate(plateType.rows, plateType.cols));
  plateWidget.editable = true;
  plateWidget.root.style.width = '100%';
  plateWidget.root.style.minWidth = '0';
  plateWidget.root.style.flexGrow = '0';
  const drcAnalysisHost = ui.divV([], 'create-plate-view__drc-host');
  drcAnalysisHost.style.width = '100%';
  drcAnalysisHost.style.flexGrow = '1';

  const templatePanel = new TemplatePanel(
    stateManager,
    plateWidget,
    () => handleManageMappings()
  );

  const plateTabsManager = new PlateTabsManager(stateManager);

  const handleManageMappings = () => {
    const state = stateManager.currentState;
    const activePlate = stateManager.activePlate;

    if (state && activePlate && activePlate.reconciliationMap.size > 0) {
      plateTabsManager.showMultiPlateMappingDialog({
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
  ;

  stateManager.onStateChange.subscribe(async (event) => {
    const activePlate = stateManager.activePlate;

    if (activePlate) {
      plateWidget.plate = activePlate.plate;
      ui.empty(drcAnalysisHost);

      const mappedDf = stateManager.getActivePlateMappedData();
      let drcCurvesElement: HTMLElement | null = null;

      if (mappedDf) {
        // Create a temporary plate object with the mapped data for analysis
        const tempPlate = new Plate(activePlate.plate.rows, activePlate.plate.cols);
        tempPlate.data = mappedDf;
        tempPlate.barcode = activePlate.plate.barcode;

        drcCurvesElement = PlateDrcAnalysis.createCurvesGrid(
          tempPlate,
          plateWidget
        );
      }

      if (drcCurvesElement)
        drcAnalysisHost.appendChild(drcCurvesElement);
      else
        drcAnalysisHost.appendChild(createAnalysisSkeleton());
    } else {
      // Show default plate if no file is loaded
      plateWidget.plate = await stateManager.getOrCreateDefaultPlate();
      ui.empty(drcAnalysisHost);
      drcAnalysisHost.appendChild(createAnalysisSkeleton());
    }

    plateWidget.refresh();
    plateWidget.updateRoleSummary();
  });

  const rightPanelTop = ui.divV([
    plateTabsManager.root,
    plateWidget.root,
  ]);
  rightPanelTop.style.flexShrink = '0';

  const analysisIcon = ui.iconFA('chart-line');
  analysisIcon.classList.add('legend-icon', 'legend-icon-analysis');
  const analysisTitle = ui.divH(
    [analysisIcon, ui.h2('Analyses')],
    {style: {alignItems: 'center', gap: '8px'}}
  );
  const drcContainer = ui.divV([analysisTitle, drcAnalysisHost], 'drc-container');


  const rightPanel = ui.divV([
    rightPanelTop,
    drcContainer,
  ], 'create-plate-view__right-panel');

  const mainLayout = ui.divH(
    [templatePanel.root, rightPanel],
    'create-plate-view__main-layout'
  );

  view.root.appendChild(mainLayout);

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
    plateTabsManager.destroy();
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
