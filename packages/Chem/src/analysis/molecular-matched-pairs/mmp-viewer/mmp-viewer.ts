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
import {MMP_NAMES} from './mmp-constants';

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

export type MmpInput = {
  table: DG.DataFrame,
  molecules: DG.Column,
  activities: DG.ColumnList,
  fragmentCutoff: number
};


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

  setupCliffsTab(sp: DG.Viewer, mmpFilters: MmpFilters, linesEditor: ScatterPlotLinesRenderer): HTMLDivElement {
    sp.root.style.width = '100%';

    this.totalCutoffMask!.setAll(true);
    this.linesMask!.setAll(true);

    this.linesRenderer = linesEditor;

    this.refilterCliffs(mmpFilters.activitySliderInputs.map((si) => si.value),
      mmpFilters.activityActiveInputs.map((ai) => ai.value), false);

    return ui.box(sp.root);
  }

  getTabs(tp: DG.Viewer, mmpFilters: MmpFilters, cliffs: HTMLDivElement): DG.TabControl {
    const tabs = ui.tabControl(null, false);

    const addToWorkspaceButton = (table: DG.DataFrame, name: string, className: string) => {
      const button = ui.iconFA('arrow-square-down', () => {
        const clonedTable = table.clone();
        clonedTable.name = name;
        const tv = grok.shell.addTableView(clonedTable);
        //without setTimeout tableView is not set as current view
        setTimeout(() => grok.shell.v = tv, 10);
      }, 'Add table to workspace');
      button.classList.add(className);
      return button;
    };

    const createGridDiv = (name: string, grid: DG.Grid) => {
      const header = ui.h1(name, 'chem-mmpa-transformation-tab-header');
      grid.root.prepend(header);
      return ui.splitV([
        ui.box(
          ui.divH([header, addToWorkspaceButton(grid.dataFrame, name, 'chem-mmpa-add-to-workspace-button')]),
          {style: {maxHeight: '30px'}},
        ),
        grid.root,
      ], {style: {width: '100%', height: '100%'}});
    };

    //const mmPairsDiv = ui.div('', {style: {width: '100%', height: '100%'}});
    const mmPairsRoot1 = createGridDiv('Matched Molecular Pairs', this.pairedGrids!.mmpGridTrans);
    const mmPairsRoot2 = createGridDiv('Matched Molecular Pairs', this.pairedGrids!.mmpGridFrag);
    const gridsDiv = ui.splitV([
      createGridDiv('Fragment Pairs', this.pairedGrids!.fpGrid),
      mmPairsRoot1,
    ], {}, true);


    const header = ui.h1('Fragment vs fragment', 'chem-mmpa-transformation-tab-header');
    const trellisTv = DG.TableView.create(this.pairedGrids!.fpGrid.dataFrame, false);
    let filters: DG.FilterGroup | null = null;

    let dockNode: DG.DockNode | null = null;
    const filterIcon = ui.icons.filter(() => {
      if (!filters)
        filters = trellisTv.getFiltersGroup();
      if (!dockNode?.parent)
        dockNode = grok.shell.tv.dockManager.dock(filters.root, DG.DOCK_TYPE.RIGHT, null, 'Fragment filters', 0.2);
    }, 'Open fragments filters');
    filterIcon.classList.add('chem-mmpa-fragments-filters-icon');

    tp.root.prepend(header);
    const tpDiv = ui.splitV([
      ui.box(
        ui.divH([header, filterIcon]),
        {style: {maxHeight: '30px'}},
      ),
      tp.root,
    ], {style: {width: '100%', height: '100%'}});

    const fragmentsDiv = ui.splitV([
      tpDiv,
      mmPairsRoot2,
    ], {}, true);

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
      return cliffs;
    });
    cliffsTab.content.classList.add('mmpa-cliffs-tab');
    const genTab = tabs.addPane(MMP_NAMES.TAB_GENERATION, () => {
      return this.generationsGrid!.root;
    });
    genTab.header.append(addToWorkspaceButton(this.generationsGrid!.dataFrame,
      'Generation', 'chem-mmpa-add-generation-to-workspace-button'));

    let refilter = true;
    tabs.onTabChanged.subscribe(() => {
      this.currentTab = tabs.currentPane.name;
      if (tabs.currentPane.name == MMP_NAMES.TAB_TRANSFORMATIONS) {
        this.pairedGrids!.enableFilters = true;
        this.pairedGrids!.refilterFragmentPairsByMolecule(false);
        //gridsDiv.append(mmPairsRoot);
        //grok.shell.o = ui.div();
        //grok.shell.o = mmPairsRoot;
      } else if (tabs.currentPane.name == MMP_NAMES.TAB_FRAGMENTS) {
        this.pairedGrids!.refreshMaskFragmentPairsFilter();
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
    const decript4 = 'Genneration of molecules based on obtained rules';

    ui.tooltip.bind(tabs.getPane(MMP_NAMES.TAB_TRANSFORMATIONS).header, decript1);
    ui.tooltip.bind(tabs.getPane(MMP_NAMES.TAB_FRAGMENTS).header, decript2);
    ui.tooltip.bind(tabs.getPane(MMP_NAMES.TAB_CLIFFS).header, decript3);
    ui.tooltip.bind(tabs.getPane(MMP_NAMES.TAB_GENERATION).header, decript4);

    return tabs;
  }

  fillAll(mmpInput: MmpInput, palette: PaletteCodes,
    mmpa: MMPA, diffs: Array<Float32Array>,
    linesIdxs: Uint32Array, pairedGrids: MmpPairedGrids, generationsGrid: DG.Grid,
    tp: DG.Viewer, sp: DG.Viewer, mmpFilters: MmpFilters,
    linesEditor: ScatterPlotLinesRenderer, lines: ILineSeries, linesActivityCorrespondance: Uint32Array,
    rdkitModule: RDModule): void {
    this.rdkitModule = rdkitModule;

    this.parentTable = mmpInput.table;

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
          mmpInput.molecules.name, this.totalData, moleculesArray, activitiesArrays, activitiesNames);
        this.totalDataUpdated = false;
      } else {
        mmpa = await MMPA.init(
          mmpInput.molecules.name, moleculesArray, mmpInput.fragmentCutoff, activitiesArrays, activitiesNames);
      }
    } catch (err: any) {
      const errMsg = err instanceof Error ? err.message : err.toString();
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

    const mmpFilters = getMmpFilters(mmpInput, mmpa.allCasesBased.maxActs,
      pairedGrids.fpGrid.dataFrame.col(MMP_NAMES.PAIRS)!.stats.max);
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

    this.parentTable!.filter.copyFrom(this.totalCutoffMask!);
    this.linesRenderer!.linesVisibility = this.linesMask!;
  }

  detach(): void {
    super.detach();
  }
}
