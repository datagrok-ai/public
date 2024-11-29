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
import {CLIFFS_TAB_TOOLTIP, FRAGMENTS_GRID_TOOLTIP, FRAGMENTS_TAB_TOOLTIP, MATHED_MOLECULAR_PAIRS_TOOLTIP_FRAGS,
  MATHED_MOLECULAR_PAIRS_TOOLTIP_TRANS, MMP_NAMES, TrellisAxis, TrellisSortByProp,
  TrellisSortType} from './mmp-constants';

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
import {getMolProperty} from '../../../package';

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
  type: TrellisSortType
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

    const addToWorkspaceButton = (table: DG.DataFrame, name: string, className: string) => {
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

    const helpButton = (className: string, tooltip: string) => {
      const button = ui.icons.help(() => {});
      ui.tooltip.bind(button, () => ui.divText(tooltip, {style: {width: '200px'}}));
      button.classList.add(className);
      return button;
    };

    const createGridDiv = (name: string, grid: DG.Grid, helpTooltip: string) => {
      const header = ui.h1(name, 'chem-mmpa-transformation-tab-header');
      const helpBUtton = helpButton('chem-mmpa-grid-help-icon', helpTooltip);
      grid.root.prepend(header);
      return ui.splitV([
        ui.box(
          ui.divH([ui.divH([header, helpBUtton]),
            addToWorkspaceButton(grid.dataFrame, name, 'chem-mmpa-add-to-workspace-button')],
          {style: {justifyContent: 'space-between'}}),
          {style: {maxHeight: '30px'}},
        ),
        grid.root,
      ], {style: {width: '100%', height: '100%'}});
    };

    //const mmPairsDiv = ui.div('', {style: {width: '100%', height: '100%'}});
    const mmPairsRoot1 = createGridDiv('Matched Molecular Pairs',
      this.pairedGrids!.mmpGridTrans, MATHED_MOLECULAR_PAIRS_TOOLTIP_TRANS);
    const mmPairsRoot2 = createGridDiv('Matched Molecular Pairs',
      this.pairedGrids!.mmpGridFrag, MATHED_MOLECULAR_PAIRS_TOOLTIP_FRAGS);
    const fpGrid = createGridDiv('Fragment Pairs',
      this.pairedGrids!.fpGrid, FRAGMENTS_GRID_TOOLTIP);
    fpGrid.prepend(ui.divText('No substitutions found for current molecule. Please select another molecule.',
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
      [TrellisAxis.From]: {property: TrellisSortByProp.Frequency, type: TrellisSortType.None},
      [TrellisAxis.To]: {property: TrellisSortByProp.Frequency, type: TrellisSortType.None},
    };

    const sortAxisChoice = ui.input.choice('Axis', {
      value: TrellisAxis.From, items: Object.keys(TrellisAxis), nullable: false, onValueChanged: () => {
        sortPropChoice.value = trellisSortState[sortAxisChoice.value as TrellisAxis].property;
        sortTypeChoice.value = trellisSortState[sortAxisChoice.value as TrellisAxis].type;
      },
    });
    const sortPropChoice = ui.input.choice('Property', {
      value: trellisSortState[TrellisAxis.From].property, items: Object.keys(TrellisSortByProp), nullable: false,
      onValueChanged: () => {
        const axis = sortAxisChoice.value as TrellisAxis;
        if (trellisSortState[axis].property !== sortPropChoice.value) {
          trellisSortState[axis].property = sortPropChoice.value as TrellisSortByProp;
          this.sortTrellis(axis, trellisSortState[axis], tp);
        }
      },
    });
    const sortTypeChoice = ui.input.choice('Type', {
      value: trellisSortState[TrellisAxis.From].type, items: Object.keys(TrellisSortType), nullable: false,
      onValueChanged: () => {
        const axis = sortAxisChoice.value as TrellisAxis;
        if (trellisSortState[axis].type !== sortTypeChoice.value) {
          trellisSortState[axis].type = sortTypeChoice.value as TrellisSortType;
          this.sortTrellis(axis, trellisSortState[axis], tp);
        }
      },
    });

    let initialSortPerformed = false;
    const sortIcon = ui.iconFA('sort-alt', () => {
      ui.showPopup(ui.inputs([
        sortAxisChoice, sortPropChoice, sortTypeChoice,
      ], 'chem-mmp-trellis-plot-sort-div'), sortIcon);
    }, 'Sort trellis plot axes');
    sortIcon.classList.add('chem-mmpa-fragments-filters-icon');

    tp.root.prepend(trellisHeader);
    const tpDiv = ui.splitV([
      ui.box(
        ui.divH([trellisHeader, filterIcon, sortIcon, helpButton('chem-mmpa-grid-help-icon', FRAGMENTS_TAB_TOOLTIP)]),
        {style: {maxHeight: '30px'}},
      ),
      tp.root,
    ], {style: {width: '100%', height: '100%'}});

    const fragmentsDiv = ui.splitV([
      tpDiv,
      mmPairsRoot2,
    ], {}, true);

    const cliffsHeader = ui.h1('2D Molecules Map', 'chem-mmpa-transformation-tab-header'); ;
    const cliffsDiv = ui.splitV([
      ui.box(
        ui.divH([cliffsHeader, helpButton('chem-mmpa-grid-help-icon', CLIFFS_TAB_TOOLTIP)]),
        {style: {maxHeight: '30px'}},
      ),
      cliffs,
    ], {style: {width: '100%', height: '100%'}});


    //tabs generation
    tabs.addPane(MMP_NAMES.TAB_TRANSFORMATIONS, () => {
      //grok.shell.o = mmPairsRoot;
      this.pairedGrids!.enableFilters = true;
      this.pairedGrids!.refilterFragmentPairsByMolecule(false);
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
      return createGridDiv('Generated Molecules', this.generationsGrid!, '');
    });

    let refilter = true;
    tabs.onTabChanged.subscribe(() => {
      if (tabs.currentPane.name !== MMP_NAMES.TAB_FRAGMENTS) {
        grok.shell.tv.dockManager.close(this.pairedGrids!.filters.root);
        this.pairedGrids!.fpMaskFragmentsTab.copyFrom(this.pairedGrids!.fpGrid.dataFrame.filter);
      }
      if (tabs.currentPane.name !== MMP_NAMES.TAB_CLIFFS) {
        if (this.parentTableFilterBackup)
          this.parentTable!.filter.copyFrom(this.parentTableFilterBackup);
      }
      this.currentTab = tabs.currentPane.name;
      if (tabs.currentPane.name == MMP_NAMES.TAB_TRANSFORMATIONS) {
        this.pairedGrids!.enableFilters = true;
        this.pairedGrids!.refilterFragmentPairsByMolecule(false);
        //gridsDiv.append(mmPairsRoot);
        //grok.shell.o = ui.div();
        //grok.shell.o = mmPairsRoot;
      } else if (tabs.currentPane.name == MMP_NAMES.TAB_FRAGMENTS) {
        this.pairedGrids!.refreshMaskFragmentPairsFilter();
        this.pairedGrids!.fpGrid.dataFrame.filter.copyFrom(this.pairedGrids!.fpMaskFragmentsTab);
        if (!initialSortPerformed) {
          initialSortPerformed = true;
          sortTypeChoice.value = TrellisSortType.Desc;
        }
        // this.pairedGrids!.enableFilters = false;
        // this.pairedGrids!.mmpMaskByFragment.setAll(false);
        // this.pairedGrids!.mmpGrid.dataFrame.filter.copyFrom(this.pairedGrids!.mmpMaskByFragment);
        // grok.shell.o = mmPairsRoot;
      } else if (tabs.currentPane.name == MMP_NAMES.TAB_CLIFFS) {
        tabs.currentPane.content.append(mmpFilters.filtersDiv);

        if (refilter)
          grok.shell.warning('Cutoff filters were applied for all activities');

        this.refilterCliffs(mmpFilters.activitySliderInputs.map((si) => si.value),
          mmpFilters.activityActiveInputs.map((ai) => ai.value), refilter);
        refilter = false;
        if (this.lastSelectedPair) {
          setTimeout(() => {
            grok.shell.o = fillPairInfo(this.lastSelectedPair!, this.linesIdxs!,
              this.linesActivityCorrespondance![this.lastSelectedPair!],
              this.pairedGrids!.mmpGridTrans.dataFrame, this.diffs!, this.parentTable!, this.rdkitModule!);
          }, 500);
        }
      } else if (tabs.currentPane.name == MMP_NAMES.TAB_GENERATION) {
        const spCorr = DG.Viewer.scatterPlot(this.generationsGrid?.dataFrame!, {
          x: '~sss',
          y: 'Prediction',
          zoomAndFilter: 'no action',
          color: 'Activity',
          showXSelector: true,
          showXAxis: true,
          showYSelector: true,
          showYAxis: true,
          showColorSelector: true,
          showSizeSelector: true,
          markerDefaultSize: 7,
        });
        const header = ui.h1('Observed vs Predicted', 'chem-mmpa-generation-tab-cp-header');
        spCorr.root.prepend(header);
        const spCorrDiv = ui.splitV([
          ui.box(
            ui.divH([header]),
            {style: {maxHeight: '30px'}},
          ),
          spCorr.root,
        ], {style: {width: '100%', height: '100%'}});
        grok.shell.o = spCorrDiv;
      }
    });

    const decript1 = 'Shows all fragmental substitutions for a given molecule';
    const decript2 = 'Analysis of fragments versus explored value';
    const decript3 = 'Cliffs analysis';
    const decript4 = 'Generation of molecules based on obtained rules';

    ui.tooltip.bind(tabs.getPane(MMP_NAMES.TAB_TRANSFORMATIONS).header, decript1);
    ui.tooltip.bind(tabs.getPane(MMP_NAMES.TAB_FRAGMENTS).header, decript2);
    ui.tooltip.bind(tabs.getPane(MMP_NAMES.TAB_CLIFFS).header, decript3);
    ui.tooltip.bind(tabs.getPane(MMP_NAMES.TAB_GENERATION).header, decript4);

    return tabs;
  }

  sortTrellis(axis: TrellisAxis, sorting: SortType, tp: DG.Viewer) {
    const filterBackup = this.pairedGrids!.fpGrid.dataFrame.filter.clone();
    switch (sorting.property) {
    case TrellisSortByProp.Frequency:
      this.sortTrellisByProp(axis, sorting.type, tp, ([_1, val1], [_2, val2]) => val1.frequency - val2.frequency,
        ([_1, val1], [_2, val2]) => val2.frequency - val1.frequency);
      break;
    case TrellisSortByProp.MW:
      if (!this.mWCalulationsReady) {
        grok.shell.warning('MW calculations for fragments in progress. Please try again later');
        return;
      }
      this.sortTrellisByProp(axis, sorting.type, tp, ([_1, val1], [_2, val2]) => val1.mw! - val2.mw!,
        ([_1, val1], [_2, val2]) => val2.mw! - val1.mw!);
    }
    this.pairedGrids!.fpGrid.dataFrame.filter.copyFrom(filterBackup);
  }

  sortTrellisByProp(axis: TrellisAxis, type: TrellisSortType, tp: DG.Viewer,
    ascSortFunc: (arg1: [string, SortData], arg2: [string, SortData]) => number,
    descSortFunc: (arg1: [string, SortData], arg2: [string, SortData]) => number) {
    const fragCol = this.pairedGrids!.fpGrid.dataFrame.col(axis);
    const axisName = tp.props.yColumnNames[0] === axis ? 'yColumnNames' :
      tp.props.xColumnNames[0] === axis ? 'xColumnNames' : null;

    if (fragCol) {
      let cats: string[] = [];
      switch (type) {
      case TrellisSortType.None:
        cats = fragCol.categories.map((it) => it).sort();
        break;
      case TrellisSortType.Asc:
        cats = Object.entries(this.fragSortingInfo).sort(ascSortFunc).map((it) => it[0]);
        break;
      case TrellisSortType.Desc:
        cats = Object.entries(this.fragSortingInfo).sort(descSortFunc).map((it) => it[0]);
        break;
      }
      if (cats.length) {
        fragCol.setCategoryOrder(cats);
        if (axis === TrellisAxis.From)
          this.pairedGrids!.fpCatsFrom = cats;
        else
          this.pairedGrids!.fpCatsTo = cats;
      }
    }
    if (axisName)
      tp.props[axisName] = tp.props[axisName];
  }

  fillAll(mmpInput: MmpInput, palette: PaletteCodes,
    mmpa: MMPA, diffs: Array<Float32Array>,
    linesIdxs: Uint32Array, pairedGrids: MmpPairedGrids, generationsGrid: DG.Grid,
    tp: DG.Viewer, sp: DG.Viewer, mmpFilters: MmpFilters,
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
    this.diffs = diffs;
    this.generationsGrid = generationsGrid;
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
    const tabs = this.getTabs(tp, mmpFilters, cliffs);

      this.linesRenderer!.lineClicked.subscribe((event: MouseOverLineEvent) => {
        this.linesRenderer!.currentLineId = event.id;
        if (event.id !== -1) {
          setTimeout(() => {
            grok.shell.o = fillPairInfo(event.id, linesIdxs, linesActivityCorrespondance[event.id],
              pairedGrids.mmpGridTrans.dataFrame, diffs, mmpInput.table, rdkitModule, this.propPanelViewer!);
            this.lastSelectedPair = event.id;
            this.propPanelViewer!.fitHeaderToLabelWidth(100);
          }, 500);
        }
      });

      this.mmpView.append(tabs);

      const propertiesColumnsNames = this.parentTable!.columns.names()
        .filter((name) => !name.startsWith('~'));
      this.propPanelViewer = new FormsViewer();
      this.propPanelViewer.dataframe = this.parentTable!;
      this.propPanelViewer.columns = propertiesColumnsNames;
      this.propPanelViewer.inputClicked.subscribe(() => {
        setTimeout(() => {
          grok.shell.o = fillPairInfo(
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
    const smilesArray = frags.map((it) =>
      it.replace(new RegExp('\[\*:1\]|\[R1\]|\[R:1\]|\[1\*\]', 'gm'), '').replaceAll('()', ''));
    getMolProperty(DG.Column.fromStrings('smiles', smilesArray), 'MW').then((res: DG.Column) => {
      frags.forEach((key, idx) => this.fragSortingInfo[key].mw = res.get(idx) ?? 0);
      this.mWCalulationsReady = true;
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

    const generationsGrid: DG.Grid = await getGenerations(mmpa, pairedGrids.fpGrid);

    this.fillAll(mmpInput, palette, mmpa, mmpa.allCasesBased.diffs, linesIdxs, pairedGrids,
      generationsGrid, tp, sp, mmpFilters, linesEditor, lines, linesActivityCorrespondance, module);

    this.totalData = mmpa.toJSON();
    //console.profileEnd('MMP');
  }

  refilterCliffs(cutoffs: number[], isActiveVar: boolean[], refilter: boolean): void {
    if (refilter) {
      for (let i = 0; i < this.cutoffMasks!.length; i ++)
      this.cutoffMasks![i].setAll(false);

      this.totalCutoffMask!.setAll(false);
      this.linesMask!.setAll(false);

      //setting activity associated masks
      for (let i = 0; i < this.lines!.from.length; i++) {
        const activityNumber = this.linesActivityCorrespondance![i];
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
    this.linesRenderer!.linesVisibility = this.linesMask!;
  }

  detach(): void {
    super.detach();
  }
}
