/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {Plate} from '../../plate/plate';
import {plateTemplates, plateTypes, savePlate, savePlateAsTemplate, initPlates} from '../plates-crud';
import {PlateStateManager} from './shared/plate-state-manager';
import {TemplatePanel} from './components/template-panel/template-panel';
import {Subscription} from 'rxjs';
import './components/plate-grid-manager/plate-grid-manager.css';
import {PlateGridManager} from './components/plate-grid-manager/plate-grid-manager';
import {AnalysisManager} from '../../plate/analyses/analysis-manager';
import './plates-create-view.css';
import {PlateWidget} from '../../plate/plate-widget/plate-widget';
import {PlateSelectionController} from '../../plate/plate-widget/plate-selection-controller';
import {MAPPING_SCOPES} from './shared/scopes';

type InteractionMode = 'view' | 'outlier-marking';

export function createPlatesView(): DG.View {
  const view = DG.View.create();
  view.name = 'Create Plate';
  view.root.classList.add('assay-plates--create-plate-view');

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


  let currentMode: InteractionMode = 'view';
  let selectionController: PlateSelectionController | null = null;

  const setInteractionMode = (mode: InteractionMode) => {
    currentMode = mode;

    selectionController?.enable();

    if (mode === 'outlier-marking') {
      selectionController?.setMode('passthrough');
      plateWidget.grid.root.classList.add('outlier-mode');
      plateWidget.plate.data.selection.setAll(false, true);
    } else {
      selectionController?.setMode('select');
      plateWidget.grid.root.classList.remove('outlier-mode');
    }

    plateWidget.grid.invalidate();
  };

  ui.tools.waitForElementInDom(plateWidget.grid.root).then(() => {
    selectionController = new PlateSelectionController(plateWidget, plateWidget.grid.overlay!);

    const isOutlierTab = plateWidget.tabs.currentPane?.name.includes('Outliers');
    setInteractionMode(isOutlierTab ? 'outlier-marking' : 'view');

    const plateTabChangeSub = plateWidget.tabs.onTabChanged.subscribe(() => {
      const tabName = plateWidget.tabs.currentPane?.name || '';
      const isOutlierTab = tabName.includes('Outliers');
      const isSummaryTab = tabName === 'Summary';
      if (isOutlierTab)
        setInteractionMode('outlier-marking');
      else if (isSummaryTab)
        setInteractionMode('view');
      else
        selectionController?.disable();
    });
    subscriptions.push(plateTabChangeSub);

    const wellClickSub = plateWidget.onWellClick.subscribe((clickEvent) => {
      if (!plateWidget.editable)
        return;

      if (currentMode === 'outlier-marking') {
        const currentState = plateWidget.plate.isOutlier(clickEvent.row, clickEvent.col);
        plateWidget.plate.markOutlier(clickEvent.row, clickEvent.col, !currentState);
        plateWidget.plate.data.selection.setAll(false, true);
      } else if (currentMode === 'view') {
        const selection = plateWidget.plate.data.selection;
        selection.set(clickEvent.dataIndex, !selection.get(clickEvent.dataIndex), true);
      }
    });

    subscriptions.push(wellClickSub);

    const selectionSub = selectionController.onSelectionComplete.subscribe((selectionEvent) => {
      if (currentMode === 'outlier-marking') {
        if (selectionEvent.isDrag) {
          const allAreOutliers = selectionEvent.wells.every((well) =>
            plateWidget.plate.isOutlier(well.row, well.col)
          );
          const shouldMark = !allAreOutliers;


          selectionEvent.wells.forEach((well) => {
            plateWidget.plate.markOutlier(well.row, well.col, shouldMark);
          });
        }
    selectionController!.clearSelection();
    plateWidget.grid.invalidate();
      } else if (currentMode === 'view') {
      }
    });
    subscriptions.push(selectionSub);
  });

  const tabControl = ui.tabControl();
  tabControl.root.classList.remove('ui-box');
  tabControl.root.classList.add('assay-plates--create-plate-view__tab-control');
  tabControl.addPane('Plate View', () => plateWidget.root);

  const analysisManager = AnalysisManager.instance;
  for (const analysis of analysisManager.analyses) {
    const createTabContent = () => {
      try {
        const activePlate = stateManager.activePlate;
        const activeIndex = stateManager.currentState?.activePlateIdx ?? -1;

        if (!activePlate || activeIndex < 0)
          return createAnalysisSkeleton(analysis.friendlyName, analysis.getRequiredFields().map((f) => f.name));

        const currentMappings = stateManager.getMappings(activeIndex, analysis.name);

        const handleMap = (target: string, source: string) => {
          stateManager.setMapping(activeIndex, analysis.name, target, source);
        };
        const handleUndo = (target: string) => {
          stateManager.removeMapping(activeIndex, analysis.name, target);
        };

        return analysis.createView(activePlate.plate, plateWidget, currentMappings, handleMap, handleUndo);
      } catch (e: any) {
        return ui.divText(`Error displaying ${analysis.friendlyName} analysis: ${e.message}`, 'error-message');
      }
    };

    tabControl.addPane(analysis.friendlyName, createTabContent);
  }

  const plateGridManager = new PlateGridManager(stateManager);

  const templatePanel = new TemplatePanel(
    stateManager
  );


  templatePanel.root.style.minWidth = '350px';
  templatePanel.root.style.maxWidth = '600px';

  stateManager.onStateChange$.subscribe(async (event) => {
    const activePlate = stateManager.activePlate;
    if (activePlate)
      plateWidget.plate = activePlate.plate;
    else
      plateWidget.plate = await stateManager.getOrCreateDefaultPlate();

    plateWidget.refresh();
    if (plateWidget.grid)
      plateWidget.grid.invalidate();

    // Refresh analysis tabs when needed
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
            const currentMappings = stateManager.getMappings(activeIndex, analysis.name);
            const handleMap = (target: string, source: string) => stateManager.setMapping(activeIndex, analysis.name, target, source);
            const handleUndo = (target: string) => stateManager.removeMapping(activeIndex, analysis.name, target);

            newContent = analysis.createView(plateToAnalyze.plate, plateWidget, currentMappings, handleMap, handleUndo, handleRerender);
          }

          ui.empty(currentPane.content);
          currentPane.content.appendChild(newContent);
        }
      }
    }
  });

  const rightSideSplitter = ui.splitV([
    plateGridManager.root,
    tabControl.root,
  ], {
    style: {
      flex: '1',
      minHeight: '0',
      width: '100%', // Explicitly set width
      minWidth: '0' // Allow shrinking
    }
  }, true);

  // Force the tab control to be flexible
  tabControl.root.style.flex = '1';
  tabControl.root.style.minHeight = '0';
  tabControl.root.style.minWidth = '0';

  // Force the plate grid manager to be flexible
  plateGridManager.root.style.flex = '1';
  plateGridManager.root.style.minHeight = '0';
  plateGridManager.root.style.minWidth = '0';

  const rightPanelWrapper = ui.divV([
    rightSideSplitter
  ], 'assay-plates--create-plate-view__right-panel');

  const mainLayout = ui.splitH(
    [templatePanel.root, rightPanelWrapper],
    {style: {width: '100%', height: '100%', minHeight: '0'}},
    true
  );
  mainLayout.classList.add('assay-plates--create-plate-view__main-layout');
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

      const template = stateManager.currentTemplate;
      const activeIndex = stateManager.currentState?.activePlateIdx ?? -1;

      if (activeIndex < 0) {
        grok.shell.warning('No active plate is selected.');
        return;
      }

      const wellPropertyMappings = stateManager.getMappings(activeIndex, MAPPING_SCOPES.TEMPLATE);
      const plateDetails = plateToSave.details || {};
      const requiredPropIds = new Set(template.required_props.map((tuple) => tuple[0]));

      const missingPlateProps: string[] = template.plateProperties
        .filter((p) => p && p.name && requiredPropIds.has(p.id!))
        .filter((prop) => {
          const value = plateDetails[prop.name!];
          return value === null || value === undefined || value === '';
        })
        .map((prop) => prop.name!);

      const missingWellProps: string[] = template.wellProperties
        .filter((p) => p && p.name && requiredPropIds.has(p.id!))
        .filter((prop) => !wellPropertyMappings.has(prop.name!))
        .map((prop) => prop.name!);

      const totalErrors = missingPlateProps.length + missingWellProps.length;

      if (totalErrors > 0) {
        const plateErrors = missingPlateProps.map((prop) => {
          const li = ui.element('li');
          li.appendChild(ui.span([prop], 'ui-label'));
          return li;
        });
        const wellErrors = missingWellProps.map((prop) => {
          const li = ui.element('li');
          li.appendChild(ui.span([prop], 'ui-label'));
          return li;
        });

        const errorContent = ui.divV([]);

        if (missingPlateProps.length > 0) {
          errorContent.appendChild(ui.h3('Missing Plate Properties:'));
          const plateErrorList = ui.element('ul', 'ui-list');
          plateErrors.forEach((li) => plateErrorList.appendChild(li));
          errorContent.appendChild(plateErrorList);
        }

        if (missingWellProps.length > 0) {
          errorContent.appendChild(ui.h3('Missing Well Mappings:'));
          const wellErrorList = ui.element('ul', 'ui-list');
          wellErrors.forEach((li) => wellErrorList.appendChild(li));
          errorContent.appendChild(wellErrorList);
        }

        ui.dialog('Validation Failed')
          .add(ui.divText('Please resolve the following required items before creating the plate:'))
          .add(errorContent)
          .onOK(() => {})
          .show();

        return;
      }

      plateToSave.plateTemplateId = stateManager.currentTemplate.id;
      await savePlate(plateToSave);
      const plateInfo = ui.divV([
        ui.divText(`Barcode: ${plateToSave.barcode}`),
        ui.divText(`Template: ${template.name}`),
      ]);
      grok.shell.info(ui.divV([ui.h3('Plate Created Successfully'), plateInfo]));
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
    selectionController?.destroy();
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

  const skeleton = ui.divV([svgDiv, message], 'assay-plates--drc-skeleton');

  return skeleton;
}
