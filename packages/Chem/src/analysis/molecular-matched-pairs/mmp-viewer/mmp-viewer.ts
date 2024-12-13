import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getRdKitModule} from '../../../utils/chem-common-rdkit';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {Subject} from 'rxjs/internal/Subject';
import $ from 'cash-dom';
import {ILineSeries, MouseOverLineEvent, ScatterPlotLinesRenderer}
  from '@datagrok-libraries/utils/src/render-lines-on-sp';

import {MMPA} from '../mmp-analysis/mmpa';
import {CLIFFS_TAB_TOOLTIP, FRAGMENTS_GRID_TOOLTIP, FRAGMENTS_TAB_TOOLTIP,
  MATHED_MOLECULAR_PAIRS_TOOLTIP_CLIFFS, MATHED_MOLECULAR_PAIRS_TOOLTIP_FRAGS,
  MATHED_MOLECULAR_PAIRS_TOOLTIP_TRANS, MMP_CONTEXT_PANE_CLASS, MMP_NAMES, SHOW_FRAGS_MODE,
  TrellisAxis, TrellisSortByProp, TrellisSortType} from './mmp-constants';

import {PaletteCodes, getPalette} from './palette';
import {getMmpTrellisPlot} from './mmp-frag-vs-frag';
import {getMmpScatterPlot, runMmpChemSpace, fillPairInfo} from './mmp-cliffs';
import {getGenerations} from './mmp-generations';

import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {FormsViewer} from '@datagrok-libraries/utils/src/viewers/forms-viewer';
import {getEmbeddingColsNames} from
  '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/reduce-dimensionality';
import {getMmpFilters, MmpFilters} from './mmp-filters';
import {getSigFigs} from '../../../utils/chem-common';
import {createLines} from './mmp-lines';
import {MmpPairedGrids} from './mmp-grids';
import {chemDescriptor} from '../../../package';
import {getZoomCoordinates} from '../../../utils/ui-utils';

export type MmpInput = {
  table: DG.DataFrame,
  molecules: DG.Column,
  activities: DG.ColumnList,
  fragmentCutoff: number
};

export type TrellisSorting = {
  [TrellisAxis.From]: SortType,
  [TrellisAxis.To]: SortType
}

export type SortType = {
  property: TrellisSortByProp,
  type?: TrellisSortType
}

export type SortData = {
  frequency: number,
  mw?: number
}

export class MatchedMolecularPairsViewer extends DG.JsViewer {
  static TYPE: string = 'MMP';

  //properties
  molecules: string | null = null;
  activities: string[] | null = null;
  fragmentCutoff: number | null;
  totalData: string;
  totalDataUpdated: boolean = false;
  onPropertyChangedObs : Subject<DG.Property | null> = new Subject<DG.Property | null>();

  moleculesCol: DG.Column | null = null;
  activitiesCols: DG.ColumnList | null = null;
  mmpa: MMPA | null = null;

  parentTable: DG.DataFrame | null = null;

  //mmpRules: MmpRules | null = null;
  mmpView: DG.View | null = null;
  colorPalette: PaletteCodes | null = null;

  pairedGrids: MmpPairedGrids | null = null;

  //cliffs tab objects
  diffs: Array<Float32Array> | null = null;
  linesIdxs: Uint32Array | null = null;
  cutoffMasks: Array<DG.BitSet> | null = null;
  totalCutoffMask: DG.BitSet | null = null;
  linesMask: BitArray | null = null;
  linesRenderer: ScatterPlotLinesRenderer | null = null;
  lines: ILineSeries | null = null;
  linesActivityCorrespondance: Uint32Array | null = null;
  //generations tab objects
  generationsGrid: DG.Grid | null = null;
  corrGrid: DG.Grid | null = null;
  generationsGridDiv = ui.div(ui.divText('Generations in progress...'));
  generationsSp: DG.ScatterPlotViewer | null = null;

  rdkitModule: RDModule | null = null;
  currentTab = '';
  lastSelectedPair: number | null = null;
  propPanelViewer: FormsViewer | null = null;
  mmpFilters: MmpFilters | null = null;

  calculatedOnGPU: boolean | null = null;
  parentTableFilterBackup: DG.BitSet | null = null;
  cliffsFiltered = false;
  mWCalulationsReady = false;

  fragSortingInfo: {[key: string]: SortData} = {};
  sp: DG.ScatterPlotViewer | null = null;
  spAxesNames: string[] = [];
  tabs: DG.TabControl | null = null;
  mutationObserver: MutationObserver | null = null;

  constructor() {
    super();
    DG.debounce(this.onPropertyChangedObs, 50).subscribe(this.onPropertyChangedDebounced.bind(this));
    //properties
    this.molecules = this.string('molecules');
    this.activities = this.stringList('activities');
    this.fragmentCutoff = this.float('fragmentCutoff');

    this.totalData = this.string('totalData', 'null', {userEditable: false, includeInLayout: false});
  }

  onPropertyChangedDebounced() {
    if (!this.dataFrame)
      return;
    if (this.totalDataUpdated) {
      this.moleculesCol = this.dataFrame.col(this.molecules!);
      this.activitiesCols = DG.DataFrame.fromColumns(this.dataFrame.columns.byNames(this.activities!)).columns;

      this.render();

      return;
    }

    this.moleculesCol = this.dataFrame.col(this.molecules!);
    this.activitiesCols = DG.DataFrame.fromColumns(this.dataFrame.columns.byNames(this.activities!)).columns;
    if (this.molecules && this.activities && this.fragmentCutoff) {
      this.render();
      return;
    }
  }

  onPropertyChanged(property: DG.Property | null): void {
    console.log(property!.name);
    super.onPropertyChanged(property);
    if (property?.name === 'totalData')
      this.totalDataUpdated = true;

    this.onPropertyChangedObs.next(property);
  }

  async render() {
    $(this.root).empty();
    if (this.dataFrame) {
      const loader = ui.div(ui.loader());
      loader.classList.add('mmpa-loader');
      this.root.appendChild(loader);
      const progressMMP = DG.TaskBarProgressIndicator.create(`Running MMP analysis...`);

      try {
        await this.runMMP(
          {table: this.dataFrame,
            molecules: this.moleculesCol!,
            activities: this.activitiesCols!,
            fragmentCutoff: this.fragmentCutoff!,
          });
      } catch (e: any) {
        const errMsg = e instanceof Error ? e.message : e.toString();
        grok.log.error(e);
        grok.shell.error(errMsg);
      } finally {
        $(this.root).empty();
        if (this.mmpView)
          this.root.appendChild(this.mmpView!.root);
        else
          this.close();
        progressMMP.close();
      }
    }
  }

  setupFilters(mmpFilters: MmpFilters, linesActivityCorrespondance: Uint32Array, tp: DG.Viewer): void {
    for (let i = 0; i < mmpFilters.activitySliderInputs.length; i++) {
      mmpFilters.activityActiveInputs[i].onChanged.subscribe(() => {
        this.refilterCliffs(mmpFilters.activitySliderInputs.map((si) => si.value),
          mmpFilters.activityActiveInputs.map((ai) => ai.value), true);
      });

      mmpFilters.activitySliderInputs[i].onChanged.subscribe(() => {
        mmpFilters.activityValuesDivs[i].innerText = mmpFilters.activitySliderInputs[i].value === 0 ? '0' :
          getSigFigs(mmpFilters.activitySliderInputs[i].value, 4).toString();
        this.refilterCliffs(mmpFilters.activitySliderInputs.map((si) => si.value),
          mmpFilters.activityActiveInputs.map((ai) => ai.value), true);
      });

      mmpFilters.activityColorInputs[i].value = this.colorPalette!.hex[i];
      mmpFilters.activityColorInputs[i].onChanged.subscribe(() => {
        const progressRendering = DG.TaskBarProgressIndicator.create(`Changing colors...`);

        //refresh lines
        for (let i =0; i < mmpFilters.activityColorInputs.length; i++) {
          this.colorPalette!.hex[i] = mmpFilters.activityColorInputs[i].value;
          this.colorPalette!.numerical[i] = DG.Color.fromHtml(mmpFilters.activityColorInputs[i].value);
          this.colorPalette!.rgb[i] = DG.Color.toRgb(this.colorPalette!.numerical[i]);
          this.colorPalette!.rgbCut[i] = this.colorPalette!.rgb[i].replace('rgb(', '').replace(')', '');
        }

        const colors = this.lines!.colors;
        for (let i = 0; i < colors!.length; i++)
          colors![i] = this.colorPalette!.rgbCut[linesActivityCorrespondance[i]];

        const lines: ILineSeries = {
          from: this.lines!.from,
          to: this.lines!.to,
          drawArrows: true,
          colors: colors,
          arrowSize: 10,
          skipMultiLineCalculation: true,
          width: 0.5,
        };

        this.lines = lines;
        this.linesRenderer!.linesToRender = lines;

        //refresh trellis plot
        const schemes = new Array<any>(mmpFilters.activityColorInputs.length);
        for (let i = 0; i < mmpFilters.activityColorInputs.length; i++)
          schemes[i] = [this.colorPalette!.numerical[i]];

        tp.setOptions({'innerViewerLook': {'colorSchemes': schemes}});
        progressRendering.close();
      });

      this.cutoffMasks![i] = DG.BitSet.create(this.parentTable!.rowCount);
      this.cutoffMasks![i].setAll(true);
    }
  }

  setupCliffsTab(sp: DG.Viewer, mmpFilters: MmpFilters, linesEditor: ScatterPlotLinesRenderer): HTMLElement {
    sp.root.style.width = '100%';

    this.totalCutoffMask!.setAll(true);
    this.linesMask!.setAll(true);

    this.linesRenderer = linesEditor;

    this.refilterCliffs(mmpFilters.activitySliderInputs.map((si) => si.value),
      mmpFilters.activityActiveInputs.map((ai) => ai.value), false);

    return sp.root;
  }

  getTabs(tp: DG.Viewer, mmpFilters: MmpFilters, cliffs: HTMLElement): DG.TabControl {
    const tabs = ui.tabControl(null, false);

    //const mmPairsDiv = ui.div('', {style: {width: '100%', height: '100%'}});
    const mmPairsRoot1 = this.createGridDiv(MMP_NAMES.PAIRS_GRID,
      this.pairedGrids!.mmpGridTrans, MATHED_MOLECULAR_PAIRS_TOOLTIP_TRANS, this.pairedGrids!.mmpGridTransMessage);

    const mmPairsRoot2 = this.createGridDiv(MMP_NAMES.PAIRS_GRID,
      this.pairedGrids!.mmpGridFrag, MATHED_MOLECULAR_PAIRS_TOOLTIP_FRAGS, this.pairedGrids!.mmpGridFragMessage);

    const mmPairsRoot3 = this.createGridDiv(MMP_NAMES.PAIRS_GRID, this.pairedGrids!.pairsGridCliffsTab,
      MATHED_MOLECULAR_PAIRS_TOOLTIP_CLIFFS, this.pairedGrids!.pairsGridCliffsTabMessage);
    mmPairsRoot3.classList.add('mmp-pairs-grid-cliffs-tab', 'cliffs-opened');

    const showFragsChoice = ui.input.choice('', {items: [SHOW_FRAGS_MODE.All, SHOW_FRAGS_MODE.Current],
      nullable: false, value: SHOW_FRAGS_MODE.All,
      onValueChanged: (value) => {
        this.pairedGrids!.fragsShowAllMode = value === SHOW_FRAGS_MODE.All;
        if (value === SHOW_FRAGS_MODE.All) {
          this.pairedGrids!.fpMaskByMolecule.setAll(true);
          this.pairedGrids!.fpGrid.dataFrame.filter.setAll(true);
        } else
          this.pairedGrids!.refilterFragmentPairsByMolecule(true);
      }});
    showFragsChoice.root.classList.add('chem-mmp-fragments-grid-mode-choice');

    const fpGrid = this.createGridDiv(MMP_NAMES.FRAGMENTS_GRID,
      this.pairedGrids!.fpGrid, FRAGMENTS_GRID_TOOLTIP, this.pairedGrids!.fpGridMessage, showFragsChoice.root);
    fpGrid.prepend(
      ui.divText('No substitutions found for current molecule. Select another molecule or switch to \'All\' mode.',
        'chem-mmpa-no-fragments-error'));
    this.subs.push(this.pairedGrids!.showErrorEvent.subscribe((showError: boolean) => {
      showError ? fpGrid.classList.add('chem-mmp-no-fragments') : fpGrid.classList.remove('chem-mmp-no-fragments');
    }));


    const gridsDiv = ui.splitV([
      fpGrid,
      mmPairsRoot1,
    ], {}, true);


    const trellisHeader = ui.h1('Fragment vs Fragment', 'chem-mmpa-transformation-tab-header');

    let dockNode: DG.DockNode | null = null;
    const filterIcon = ui.icons.filter(() => {
      if (!dockNode?.parent) {
        dockNode = grok.shell.tv.dockManager
          .dock(this.pairedGrids!.filters.root, DG.DOCK_TYPE.RIGHT, null, 'Fragment filters', 0.2);
      }
    }, 'Open fragments filters');
    filterIcon.classList.add('chem-mmpa-fragments-filters-icon');

    const trellisSortState: TrellisSorting = {
      [TrellisAxis.From]: {property: TrellisSortByProp.Frequency, type: TrellisSortType.Desc},
      [TrellisAxis.To]: {property: TrellisSortByProp.None, type: TrellisSortType.Asc},
    };


    const summaryColsButton = ui.div('', 'mmp-trellis-summary-column');
    ui.tooltip.bind(summaryColsButton, 'Select columns to show in trellis plot');

    tp.root.prepend(trellisHeader);
    const tpDiv = ui.splitV([
      ui.box(
        ui.divH([trellisHeader, filterIcon, summaryColsButton,
          this.helpButton('chem-mmpa-grid-help-icon', FRAGMENTS_TAB_TOOLTIP)]),
        {style: {maxHeight: '30px'}},
      ),
      tp.root,
    ], {style: {width: '100%', height: '100%'}});

    tp.onEvent('d4-viewer-rendered').subscribe(() => {
      this.createSortIcon(trellisSortState, TrellisAxis.From, tp, 'chem-mmpa-fragments-sort-icon-x-axis');
      this.createSortIcon(trellisSortState, TrellisAxis.To, tp, 'chem-mmpa-fragments-sort-icon-y-axis');
      const tpButtons = Array.from(tp.root.getElementsByTagName('button'));
      if (tpButtons.length) {
        //workaround to disable aggregation functions selections
        tpButtons[0].onmousedown = () => {
          if (!this.mutationObserver) {
            this.mutationObserver = new MutationObserver((mutationsList) => {
              for (let i = 0; i < mutationsList.length; i++) {
                for (const node of Array.from(mutationsList[i].addedNodes)) {
                  const dlgHeaders = Array.from((node as HTMLElement).getElementsByClassName('d4-dialog-title'))
                    .filter((el) => (el as HTMLElement).innerText === 'Edit columns aggregations');
                  if (dlgHeaders.length) {
                    dlgHeaders[0].parentElement?.parentElement?.classList.add('mmp-trellis-summary-col-dlg');
                    this.mutationObserver?.disconnect();
                    this.mutationObserver = null;
                    return;
                  }
                }
              }
            });
            this.mutationObserver.observe(document.body, {attributes: true, childList: true});
          }
        };
        ui.empty(summaryColsButton);
        summaryColsButton.append(tpButtons[0]);
      }
    });

    this.sortTrellis(TrellisAxis.From, trellisSortState[TrellisAxis.From], tp);

    const fragmentsDiv = ui.splitV([
      tpDiv,
      mmPairsRoot2,
    ], {}, true);

    let cliffsOpened = true;
    const cliffsHeader = ui.h1('2D Molecules Map', 'chem-mmpa-transformation-tab-header'); ;
    const cliffsNumButton = ui.button(`Close pairs`, () => {
      cliffsOpened = !cliffsOpened;
      if (cliffsOpened) {
        cliffsNumButton.innerText = 'Close pairs';
        mmPairsRoot3.classList.replace('cliffs-closed', 'cliffs-opened');
        //workaround for case when pairs grid is opened for the first time with minimal height
        if (parseFloat(mmPairsRoot3.style.height.replace('px', '')) < 5)
          mmPairsRoot3.style.height = '200px';
      } else {
        cliffsNumButton.innerText = 'Open pairs';
        mmPairsRoot3.classList.replace('cliffs-opened', 'cliffs-closed');
      }
    });
    cliffsNumButton.classList.add('chem-mmp-open-cliffs-button');
    const spDiv = ui.splitV([
      ui.box(
        ui.divH([cliffsHeader, cliffsNumButton, this.helpButton('chem-mmpa-grid-help-icon', CLIFFS_TAB_TOOLTIP)]),
        {style: {maxHeight: '30px'}},
      ),
      cliffs,
    ]);

    const cliffsDiv = ui.splitV([
      spDiv,
      mmPairsRoot3,
    ], {style: {width: '100%', height: '100%'}}, true);


    this.pairedGrids!.pairsGridCliffsTab.dataFrame.onCurrentRowChanged.subscribe(() => {
      const currentRowIdx = this.pairedGrids!.pairsGridCliffsTab.dataFrame.currentRowIdx;
      if (currentRowIdx !== -1 && this.currentTab === MMP_NAMES.TAB_CLIFFS) {
        const fromIdx = this.pairedGrids!.pairsGridCliffsTab.dataFrame.get(MMP_NAMES.PAIRNUM_FROM, currentRowIdx);
        const toIdx = this.pairedGrids!.pairsGridCliffsTab.dataFrame.get(MMP_NAMES.PAIRNUM_TO, currentRowIdx);
        //this.lineIdxs contain idxs if pairs from pairs dataset
        let currentLineIdx: number | null = null;
        for (let i = 0; i < this.linesIdxs!.length; i++) {
          if (this.linesIdxs![i] === currentRowIdx && this.linesRenderer!.visibility?.getBit(i)) {
            currentLineIdx = i;
            break;
          }
        }
        if (currentLineIdx) {
          this.linesRenderer!.currentLineId = currentLineIdx;
          const {zoomLeft, zoomRight, zoomTop, zoomBottom} = getZoomCoordinates(
              this.sp!.viewport.width,
              this.sp!.viewport.height,
              this.sp!.dataFrame.get(this.spAxesNames[0], fromIdx),
              this.sp!.dataFrame.get(this.spAxesNames[1], fromIdx),
              this.sp!.dataFrame.get(this.spAxesNames[0], toIdx),
              this.sp!.dataFrame.get(this.spAxesNames[1], toIdx),
          );
            this.sp!.zoom(zoomLeft,
              zoomTop,
              zoomRight,
              zoomBottom);
        }
      }
    });


    //tabs generation
    tabs.addPane(MMP_NAMES.TAB_TRANSFORMATIONS, () => {
      //grok.shell.o = mmPairsRoot;
      this.pairedGrids!.enableFilters = true;
      return gridsDiv;
    });
    const fragmentsPane = tabs.addPane(MMP_NAMES.TAB_FRAGMENTS, () => {
      return fragmentsDiv;
      //return tp.root;
    });
    fragmentsPane.content.classList.add('mmpa-fragments-tab');
    const cliffsTab = tabs.addPane(MMP_NAMES.TAB_CLIFFS, () => {
      return cliffsDiv;
    });

    cliffsTab.content.classList.add('mmpa-cliffs-tab');
    tabs.addPane(MMP_NAMES.TAB_GENERATION, () => {
      return this.generationsGridDiv;
    });

    let refilter = true;
    tabs.onTabChanged.subscribe(() => {
      if (this.currentTab === MMP_NAMES.TAB_FRAGMENTS && tabs.currentPane.name !== MMP_NAMES.TAB_FRAGMENTS) {
        grok.shell.tv.dockManager.close(this.pairedGrids!.filters.root);
        this.pairedGrids!.fpMaskFragmentsTab.copyFrom(this.pairedGrids!.fpGrid.dataFrame.filter);
      }
      if (this.currentTab === MMP_NAMES.TAB_CLIFFS && tabs.currentPane.name !== MMP_NAMES.TAB_CLIFFS) {
        if (this.parentTableFilterBackup)
          this.parentTable!.filter.copyFrom(this.parentTableFilterBackup);
      }
      this.currentTab = tabs.currentPane.name;
      this.pairedGrids!.currentTab = tabs.currentPane.name as MMP_NAMES;
      if (tabs.currentPane.name == MMP_NAMES.TAB_TRANSFORMATIONS) {
        this.pairedGrids!.enableFilters = true;
        //setting masks on fragments grid and pairs grid
        this.pairedGrids!.mmpGridTrans.dataFrame.filter.copyFrom(this.pairedGrids!.mmpMaskTrans);
        this.pairedGrids!.fpGrid!.dataFrame.filter.copyFrom(this.pairedGrids!.fpMaskByMolecule!);
      } else if (tabs.currentPane.name == MMP_NAMES.TAB_FRAGMENTS) {
        this.pairedGrids!.refreshMaskFragmentPairsFilter();
        this.pairedGrids!.fpGrid.dataFrame.filter.copyFrom(this.pairedGrids!.fpMaskFragmentsTab);
      } else if (tabs.currentPane.name == MMP_NAMES.TAB_CLIFFS) {
        this.sp!.root.append(mmpFilters.filtersDiv);

        if (refilter)
          grok.shell.warning('Cutoff filters were applied for all activities');

        this.refilterCliffs(mmpFilters.activitySliderInputs.map((si) => si.value),
          mmpFilters.activityActiveInputs.map((ai) => ai.value), refilter);
        refilter = false;
        if (this.lastSelectedPair) {
          setTimeout(() => {
            grok.shell.windows.showContextPanel = true;
            grok.shell.o = fillPairInfo(this.mmpa!, this.lastSelectedPair!, this.linesIdxs!,
              this.linesActivityCorrespondance![this.lastSelectedPair!],
              this.pairedGrids!.mmpGridTrans.dataFrame, this.diffs!, this.parentTable!, this.rdkitModule!);
          }, 500);
        }
      } else if (tabs.currentPane.name == MMP_NAMES.TAB_GENERATION) {
        if (this.generationsGrid) {
          if (this.generationsSp) {
            const header = ui.h1('Observed vs Predicted', 'chem-mmpa-generation-tab-cp-header');
            this.generationsSp.root.prepend(header);
            const spCorrDiv = ui.splitV([
              ui.box(
                ui.divH([header]),
                {style: {maxHeight: '30px'}},
              ),
              this.generationsSp.root,
            ], {style: {width: '100%', height: '100%'}});
            spCorrDiv.classList.add(MMP_CONTEXT_PANE_CLASS);
            grok.shell.windows.showContextPanel = true;
            grok.shell.o = spCorrDiv;
          }
        }
      }
    });

    const decript1 = 'View all fragment substitutions found in the dataset';
    const decript2 = 'Analyze activity changes across fragment substitutions using a trellis plot';
    const decript3 = 'Molecule pairs analysis on 2d scatter plot';
    const decript4 = 'Generation of molecules based on obtained rules';

    ui.tooltip.bind(tabs.getPane(MMP_NAMES.TAB_TRANSFORMATIONS).header, decript1);
    ui.tooltip.bind(tabs.getPane(MMP_NAMES.TAB_FRAGMENTS).header, decript2);
    ui.tooltip.bind(tabs.getPane(MMP_NAMES.TAB_CLIFFS).header, decript3);
    ui.tooltip.bind(tabs.getPane(MMP_NAMES.TAB_GENERATION).header, decript4);

    return tabs;
  }

  createGridDiv(name: string, grid: DG.Grid, helpTooltip: string, messageBox: HTMLElement, extraEl?: HTMLElement) {
    const header = ui.h1(name, 'chem-mmpa-transformation-tab-header');
    const helpBUtton = this.helpButton('chem-mmpa-grid-help-icon', helpTooltip);
    grid.root.prepend(header);
    return ui.splitV([
      ui.box(
        ui.divH([ui.divH([header, extraEl ?? ui.div(),
          this.addToWorkspaceButton(grid.dataFrame, name, 'chem-mmpa-add-to-workspace-button'), helpBUtton]),
        messageBox], {style: {justifyContent: 'space-between'}}),
        {style: {maxHeight: '30px'}},
      ),
      grid.root,
    ], {style: {width: '100%', height: '100%'}});
  };

  addToWorkspaceButton(table: DG.DataFrame, name: string, className: string) {
    const button = ui.icons.add(() => {
      const clonedTable = table.clone();
      clonedTable.name = name;
      const tv = grok.shell.addTableView(clonedTable);
      //without setTimeout tableView is not set as current view
      setTimeout(() => grok.shell.v = tv, 10);
    }, 'Add table to workspace');
    button.classList.add(className);
    return button;
  };

  helpButton(className: string, tooltip: string) {
    const button = ui.icons.help(() => {
      grok.shell.windows.showHelp = true;
      grok.shell.windows.help.showHelp('/help/datagrok/solutions/domains/chem/chem#matched-molecular-pairs');
    });
    ui.tooltip.bind(button, () => ui.divText(tooltip, {style: {width: '300px'}}));
    button.classList.add(className);
    return button;
  };

  createSortIcon(trellisSortState: TrellisSorting, axis: TrellisAxis, tp: DG.Viewer, className: string) {
    const tpAxisSelectorArr = Array.from(tp.root.getElementsByClassName('d4-column-selector-column'))
      .filter((it) => (it as HTMLElement).innerText === axis);
    if (tpAxisSelectorArr.length) {
      const tpAxisSelector = tpAxisSelectorArr[0] as HTMLElement;
      const sortIconExists = tpAxisSelector.parentElement!.parentElement!.getElementsByClassName(className).length;
      if (!sortIconExists) {
        const sortPropChoice = ui.input.choice('Sort by', {
          value: trellisSortState[axis].property, items: Object.keys(TrellisSortByProp), nullable: false,
          onValueChanged: () => {
            if (trellisSortState[axis].property !== sortPropChoice.value) {
              if (sortPropChoice.value === TrellisSortByProp.None)
                sortTypeChoice.root.classList.add('chem-mmp-trellis-sort-type-invisible');
              else
                sortTypeChoice.root.classList.remove('chem-mmp-trellis-sort-type-invisible');
              trellisSortState[axis].property = sortPropChoice.value as TrellisSortByProp;
              this.sortTrellis(axis, trellisSortState[axis], tp);
            }
          },
        });
        const sortTypeChoice = ui.input.choice('Order', {
          value: trellisSortState[axis].type, items: Object.keys(TrellisSortType), nullable: false,
          onValueChanged: () => {
            if (trellisSortState[axis].type !== sortTypeChoice.value) {
              trellisSortState[axis].type = sortTypeChoice.value as TrellisSortType;
              this.sortTrellis(axis, trellisSortState[axis], tp);
            }
          },
        });
        if (trellisSortState[axis].property === TrellisSortByProp.None)
          sortTypeChoice.root.classList.add('chem-mmp-trellis-sort-type-invisible');
        const sortIcon = ui.iconFA('sort-alt', () => {
          ui.showPopup(ui.inputs([
            sortPropChoice, sortTypeChoice,
          ], 'chem-mmp-trellis-plot-sort-div'), sortIcon);
        }, `Sort ${axis} axis`);
        sortIcon.classList.add(className);
          tpAxisSelector.parentElement!.parentElement!.append(sortIcon);
          tpAxisSelector.onclick = (e: MouseEvent) => {
            e.stopImmediatePropagation();
            e.preventDefault();
          };
      }
    }
  }

  sortTrellis(axis: TrellisAxis, sorting: SortType, tp: DG.Viewer) {
    const filterBackup = this.pairedGrids!.fpGrid.dataFrame.filter.clone();
    switch (sorting.property) {
    case TrellisSortByProp.None:
      this.sortTrellisByProp(axis, tp);
      break;
    case TrellisSortByProp.Frequency:
      this.sortTrellisByProp(axis, tp, sorting.type, ([_1, val1], [_2, val2]) => val1.frequency - val2.frequency,
        ([_1, val1], [_2, val2]) => val2.frequency - val1.frequency);
      break;
    case TrellisSortByProp.MW:
      if (!this.mWCalulationsReady) {
        grok.shell.warning('MW calculations for fragments in progress. Please try again later');
        return;
      }
      this.sortTrellisByProp(axis, tp, sorting.type, ([_1, val1], [_2, val2]) => val1.mw! - val2.mw!,
        ([_1, val1], [_2, val2]) => val2.mw! - val1.mw!);
    }
    this.pairedGrids!.fpGrid.dataFrame.filter.copyFrom(filterBackup);
  }

  sortTrellisByProp(axis: TrellisAxis, tp: DG.Viewer, type?: TrellisSortType,
    ascSortFunc?: (arg1: [string, SortData], arg2: [string, SortData]) => number,
    descSortFunc?: (arg1: [string, SortData], arg2: [string, SortData]) => number) {
    const fragCol = this.pairedGrids!.fpGrid.dataFrame.col(axis);
    const axisName = tp.props.yColumnNames[0] === axis ? 'yColumnNames' :
      tp.props.xColumnNames[0] === axis ? 'xColumnNames' : null;

    if (fragCol) {
      let cats: string[] = [];
      switch (type) {
      case undefined:
        cats = fragCol.categories.map((it) => it).sort();
        break;
      case TrellisSortType.Asc:
        cats = Object.entries(this.fragSortingInfo).sort(ascSortFunc).map((it) => it[0]);
        break;
      case TrellisSortType.Desc:
        cats = Object.entries(this.fragSortingInfo).sort(descSortFunc).map((it) => it[0]);
        break;
      }
      if (cats.length)
        fragCol.setCategoryOrder(cats);
    }
    if (axisName)
      tp.props[axisName] = tp.props[axisName];
  }

  fillAll(mmpInput: MmpInput, palette: PaletteCodes,
    mmpa: MMPA, diffs: Array<Float32Array>,
    linesIdxs: Uint32Array, pairedGrids: MmpPairedGrids,
    tp: DG.Viewer, sp: DG.ScatterPlotViewer, spAxesNames: string[], mmpFilters: MmpFilters,
    linesEditor: ScatterPlotLinesRenderer, lines: ILineSeries, linesActivityCorrespondance: Uint32Array,
    rdkitModule: RDModule): void {
    this.rdkitModule = rdkitModule;

    this.parentTable = mmpInput.table;
    this.parentTableFilterBackup = DG.BitSet.create(this.parentTable.rowCount).copyFrom(this.parentTable.filter);

    this.subs.push(DG.debounce(this.parentTable!.onFilterChanged, 1000).subscribe(() => {
      if (!this.cliffsFiltered) {
          this.parentTableFilterBackup!.copyFrom(this.parentTable!.filter);
          if (this.currentTab === MMP_NAMES.TAB_CLIFFS) {
            setTimeout(() => {
              this.cliffsFiltered = true;
              this.parentTable!.filter.and(this.totalCutoffMask!);
            }, 10);
          }
      } else
        this.cliffsFiltered = false;
    }));

    this.colorPalette = palette;
    this.mmpa = mmpa;

    this.mmpView = DG.View.create();
    this.mmpView!.name = MMP_NAMES.VIEW_NAME;
    this.mmpView!.box = true;

    //main grids
    this.pairedGrids = pairedGrids;
    this.pairedGrids.setupGrids();

    //Cliffs tab setup
    this.sp = sp;
    this.spAxesNames = spAxesNames;
    this.diffs = diffs;
    this.lines = lines;
    this.linesActivityCorrespondance = linesActivityCorrespondance;
    this.linesIdxs = linesIdxs;
    this.mmpFilters = mmpFilters;
    this.cutoffMasks = new Array<DG.BitSet>(mmpFilters.activitySliderInputs.length);
    this.totalCutoffMask = DG.BitSet.create(this.parentTable!.rowCount);
    this.linesMask = new BitArray(linesIdxs.length);
    this.setupFilters(mmpFilters, linesActivityCorrespondance, tp);
    const cliffs = this.setupCliffsTab(sp, mmpFilters, linesEditor);

    //tabs
    this.tabs = this.getTabs(tp, mmpFilters, cliffs);

      this.linesRenderer!.lineClicked.subscribe((event: MouseOverLineEvent) => {
        this.linesRenderer!.currentLineId = event.id;
        if (event.id !== -1) {
          setTimeout(() => {
            grok.shell.windows.showContextPanel = true;
            grok.shell.o = fillPairInfo(this.mmpa!, event.id, linesIdxs, linesActivityCorrespondance[event.id],
              pairedGrids.mmpGridTrans.dataFrame, diffs, mmpInput.table, rdkitModule, this.propPanelViewer!);
            this.lastSelectedPair = event.id;
            this.propPanelViewer!.fitHeaderToLabelWidth(100);
          }, 500);
        }
      });

      this.mmpView.append(this.tabs);

      const propertiesColumnsNames = this.parentTable!.columns.names()
        .filter((name) => !name.startsWith('~'));
      this.propPanelViewer = new FormsViewer();
      this.propPanelViewer.dataframe = this.parentTable!;
      this.propPanelViewer.columns = propertiesColumnsNames;
      this.propPanelViewer.inputClicked.subscribe(() => {
        setTimeout(() => {
          grok.shell.windows.showContextPanel = true;
          grok.shell.o = fillPairInfo(this.mmpa!,
            this.lastSelectedPair!, linesIdxs, linesActivityCorrespondance[this.lastSelectedPair!],
            pairedGrids.mmpGridTrans.dataFrame, diffs, mmpInput.table, rdkitModule, this.propPanelViewer!);
          this.propPanelViewer!.fitHeaderToLabelWidth(100);
        }, 500);
      });

      this.calculatedOnGPU = mmpa.gpu;
      this.prepareMwForSorting();
  }

  async prepareMwForSorting() {
    const frags = Object.keys(this.fragSortingInfo);
    chemDescriptor(DG.Column.fromStrings('smiles', frags), 'MolWt').then((res: DG.Column) => {
      let errorCount = 0;
      frags.forEach((key, idx) => {
        let resMW = res.get(idx);
        if (!resMW || typeof resMW === 'string') {
          errorCount++;
          resMW = undefined;
        }
        this.fragSortingInfo[key].mw = resMW ?? 0;
      });
      this.mWCalulationsReady = true;
      if (errorCount > 0)
        grok.shell.warning(`Molecular weight hasn't been calculated for ${errorCount} fragments`);
    });
  }

  async runMMP(mmpInput: MmpInput) {
    //console.profile('MMP');
    const moleculesArray = mmpInput.molecules.toList();
    const variates = mmpInput.activities.length;
    const activitiesArrays = new Array<Float32Array>(variates);
    const activitiesNames = new Array<string>(variates);
    for (let i = 0; i < variates; i++) {
      activitiesArrays[i] = mmpInput.activities.byIndex(i).asDoubleList();
      activitiesNames[i] = mmpInput.activities.byIndex(i).name;
    }

    let mmpa: MMPA;

    try {
      if (this.totalDataUpdated) {
        mmpa = await MMPA.fromData(
          mmpInput.molecules.name, this.totalData, moleculesArray,
          activitiesArrays, activitiesNames, this.fragSortingInfo);
        this.totalDataUpdated = false;
      } else {
        mmpa = await MMPA.init(
          mmpInput.molecules.name, moleculesArray, mmpInput.fragmentCutoff,
          activitiesArrays, activitiesNames, this.fragSortingInfo);
      }
    } catch (err: any) {
      const errMsg = err instanceof Error ? err.message : err.toString();
      grok.log.error(err);
      grok.shell.error(errMsg);
      throw new Error(errMsg);
    }

    const palette = getPalette(mmpInput.activities.length);
    const activityMeanNames = Array<string>(mmpa.initData.activitiesCount);
    for (let i = 0; i < mmpa.initData.activitiesCount; i++) {
      const name = MMP_NAMES.MEANDIFF + ' ' + mmpa.initData.activitiesNames[i];
      activityMeanNames[i] = name;
    }

    //Transformations tab
    //const {activityMeanNames, fpGrid, mmpGrid} = getMmpPairsGrids(mmpa);
    const pairedGrids = new MmpPairedGrids(this.subs, mmpInput, mmpa, activityMeanNames);

    //Fragments tab
    const tp = getMmpTrellisPlot(pairedGrids.fpGrid, activityMeanNames, palette);
    tp.onEvent('d4-trellis-plot-inner-viewer-clicked').subscribe((cats) => {
      //grok.shell.info(cats);
      this.pairedGrids?.refilterMatchedPairsByFragments(cats);
    });

    //Cliffs tab
    const {linesIdxs, lines, linesActivityCorrespondance} = createLines(mmpa, palette);
    const embedColsNames = getEmbeddingColsNames(mmpInput.table).map((it) => `~${it}`);

    const mmpFilters = getMmpFilters(mmpInput, mmpa.allCasesBased.maxActs);
    console.log(`created mmpa filters`);

    const sp = getMmpScatterPlot(mmpInput, embedColsNames, mmpInput.molecules.name);

    //running internal chemspace
    const module = getRdKitModule();
    const [linesEditor, chemSpaceParams] = runMmpChemSpace(mmpInput, sp, lines, linesIdxs, linesActivityCorrespondance,
      pairedGrids.mmpGridTrans.dataFrame, mmpa, module, embedColsNames);

    const progressBarSpace = DG.TaskBarProgressIndicator.create(`Running Chemical space...`);
    mmpa.chemSpace(chemSpaceParams).then((res) => {
      const embeddings = res.coordinates;
      for (const col of embeddings)
        mmpInput.table.columns.replace(col.name, col);

      this.totalData = mmpa.toJSON();
      progressBarSpace.close();
    });

    getGenerations(mmpa, pairedGrids.fpGrid).then(([genGrid, corrGrid]) => {
      this.generationsGrid = genGrid;
      this.corrGrid = corrGrid;

      this.corrGrid.dataFrame.meta.formulaLines.addLine({
        title: 'Identity',
        formula: '${Observed} = ${Predicted}',
        color: '#717581',
        width: 1,
        visible: true,
      });

      this.generationsSp = DG.Viewer.scatterPlot(this.corrGrid?.dataFrame!, {
        x: 'Observed',
        y: 'Predicted',
        zoomAndFilter: 'no action',
        color: 'Activity',
        showXSelector: true,
        showXAxis: true,
        showYSelector: true,
        showYAxis: true,
        showColorSelector: true,
        showSizeSelector: true,
        markerDefaultSize: mmpa.initData.molecules.length > 10000 ? 1 : 2,
        markerType: 'circle',
        showRegressionLine: true,
      });

      ui.empty(this.generationsGridDiv);
      this.generationsGridDiv.append(this.createGridDiv('Generated Molecules', this.generationsGrid!, '', ui.div()));
    }).catch((error: any) => {
      const errorStr = `Generations haven't been completed due to error: ${error}`;
      ui.empty(this.generationsGridDiv);
      this.generationsGridDiv.append(ui.divText(errorStr));
      grok.shell.error(errorStr);
    });

    this.fillAll(mmpInput, palette, mmpa, mmpa.allCasesBased.diffs, linesIdxs, pairedGrids,
      tp, sp, embedColsNames, mmpFilters, linesEditor, lines, linesActivityCorrespondance, module);

    this.totalData = mmpa.toJSON();
    //console.profileEnd('MMP');
  }

  refilterCliffs(cutoffs: number[], isActiveVar: boolean[], refilter: boolean): void {
    if (refilter) {
      for (let i = 0; i < this.cutoffMasks!.length; i ++)
      this.cutoffMasks![i].setAll(false);

      this.totalCutoffMask!.setAll(false);
      this.linesMask!.setAll(false);
      this.pairedGrids?.pairsMaskCliffsTab.setAll(false);

      //setting activity associated masks
      for (let i = 0; i < this.lines!.from.length; i++) {
        const activityNumber = this.linesActivityCorrespondance![i];
        //line is idx of pair in pairs dataset
        const line = this.linesIdxs![i];
        //if 'isActive' checkbox for variable is unset
        //setting line to invisible and continue, otherwise look at cutoff
        if (!isActiveVar[activityNumber]) {
          this.linesMask!.setBit(i, false, false);
          continue;
        }
        if (this.diffs![activityNumber][line] >= cutoffs[activityNumber]) {
          this.cutoffMasks![activityNumber].set(this.lines!.from[i], true, false);
          this.cutoffMasks![activityNumber].set(this.lines!.to[i], true, false);
          this.linesMask!.setBit(i, true, false);
          this.pairedGrids?.pairsMaskCliffsTab.set(line, true);
        }
      }

      //setting resulting masks
      //TODO: refine
      for (let j = 0; j < this.totalCutoffMask!.length; j++) {
        let res = false;
        for (let i = 0; i < this.cutoffMasks!.length; i++)
          res = res || this.cutoffMasks![i].get(j);
        this.totalCutoffMask!.set(j, res, false);
      }
    }

    this.cliffsFiltered = true;
    this.parentTable!.filter.copyFrom(this.parentTableFilterBackup!).and(this.totalCutoffMask!);
    this.pairedGrids?.pairsGridCliffsTab.dataFrame.filter.copyFrom(this.pairedGrids.pairsMaskCliffsTab);
    this.linesRenderer!.linesVisibility = this.linesMask!;
  }

  detach(): void {
    if ((grok.shell.o as HTMLElement).classList?.contains(MMP_CONTEXT_PANE_CLASS))
      grok.shell.o = ui.div();
    super.detach();
  }
}
