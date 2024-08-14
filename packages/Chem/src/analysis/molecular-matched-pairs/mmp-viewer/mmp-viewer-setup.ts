import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {MmpRules, MmpInput} from '../mmp-constants';
import {PaletteCodes} from '../mmp-mol-rendering';
import {MmpFilters} from '../mmp-filters';
import {ILineSeries, MouseOverLineEvent, ScatterPlotLinesRenderer}
  from '@datagrok-libraries/utils/src/render-lines-on-sp';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {FormsViewer} from '@datagrok-libraries/utils/src/viewers/forms-viewer';
import {fillPairInfo} from '../mmp-cliffs';
import {debounceTime} from 'rxjs/operators';
import {MMP_NAMES} from '../mmp-constants';
import {getSigFigs} from '../../../utils/chem-common';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {MatchedMolecularPairsViewer} from './mmp-viewer';

declare global {
  interface MatchedMolecularPairsViewer {
    setupTransformationTab(): void;
    setupFilters(mmpFilters: MmpFilters, linesActivityCorrespondance: Uint32Array, tp: DG.Viewer): void;
    setupCliffsTab(sp: DG.Viewer, mmpFilters: MmpFilters, linesEditor: ScatterPlotLinesRenderer): HTMLDivElement;
    getTabs(tp: DG.Viewer, mmpFilters: MmpFilters, cliffs: HTMLDivElement): DG.TabControl;
    fillAll(mmpInput: MmpInput, palette: PaletteCodes,
      rules: MmpRules, diffs: Array<Float32Array>,
      linesIdxs: Uint32Array, transFragmentsGrid: DG.Grid, transPairsGrid: DG.Grid, generationsGrid: DG.Grid,
      tp: DG.Viewer, sp: DG.Viewer, mmpFilters: MmpFilters,
      linesEditor: ScatterPlotLinesRenderer, lines: ILineSeries, linesActivityCorrespondance: Uint32Array,
      rdkitModule: RDModule, gpuUsed: boolean): void;
  }
}

MatchedMolecularPairsViewer.prototype.setupTransformationTab = function(): void {
  this.transFragmentsMask!.setAll(true);
  this.parentTable!.onCurrentRowChanged.pipe(debounceTime(1000)).subscribe(() => {
    if (this.parentTable!.currentRowIdx !== -1) {
      this.refilterAllFragments(true);
      this.refreshPair(this.rdkitModule!);
    }
  });

  this.refilterAllFragments(true);

  this.transFragmentsGrid!.table.onCurrentRowChanged.subscribe(() => {
    if (this.transFragmentsGrid!.table.currentRowIdx !== -1)
      this.refreshPair(this.rdkitModule!);
  });

  this.transPairsGrid!.table.onCurrentRowChanged.subscribe(() => {
    if (this.transPairsGrid!.table.currentRowIdx !== -1)
      this.pinPair(this.transPairsGrid!.table.currentRowIdx);
  });

  this.mmpView!.name = MMP_NAMES.VIEW_NAME;
  this.mmpView!.box = true;

  this.transPairsMask!.setAll(false);
  this.refreshPair(this.rdkitModule!);
};

MatchedMolecularPairsViewer.prototype.setupFilters =
function(mmpFilters: MmpFilters, linesActivityCorrespondance: Uint32Array, tp: DG.Viewer): void {
  for (let i = 0; i < mmpFilters.activitySliderInputs.length; i ++) {
    mmpFilters.activityActiveInputs[i].onChanged(() => {
      this.refilterCliffs(mmpFilters.activitySliderInputs.map((si) => si.value),
        mmpFilters.activityActiveInputs.map((ai) => ai.value), true);
    });

    mmpFilters.activitySliderInputs[i].onChanged(() => {
      mmpFilters.activityValuesDivs[i].innerText = mmpFilters.activitySliderInputs[i].value === 0 ? '0' :
        getSigFigs(mmpFilters.activitySliderInputs[i].value, 4).toString();
      this.refilterCliffs(mmpFilters.activitySliderInputs.map((si) => si.value),
        mmpFilters.activityActiveInputs.map((ai) => ai.value), true);
    });

    mmpFilters.activityColorInputs[i].value = this.colorPalette!.hex[i];
    mmpFilters.activityColorInputs[i].onChanged(() => {
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

  mmpFilters.pairsSliderInput.onChanged(() => {
    mmpFilters.pairsValueDiv.innerText = mmpFilters.pairsSliderInput.value.toString();
    // const value = mmpFilters.pairsSliderInput.value;
  });
};

MatchedMolecularPairsViewer.prototype.setupCliffsTab =
function(sp: DG.Viewer, mmpFilters: MmpFilters, linesEditor: ScatterPlotLinesRenderer): HTMLDivElement {
  sp.root.style.width = '100%';

  this.totalCutoffMask!.setAll(true);
  this.linesMask!.setAll(true);

  this.linesRenderer = linesEditor;

  this.refilterCliffs(mmpFilters.activitySliderInputs.map((si) => si.value),
    mmpFilters.activityActiveInputs.map((ai) => ai.value), false);

  return ui.box(sp.root);
};

MatchedMolecularPairsViewer.prototype.getTabs =
function(tp: DG.Viewer, mmpFilters: MmpFilters, cliffs: HTMLDivElement): DG.TabControl {
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

  tabs.addPane(MMP_NAMES.TAB_TRANSFORMATIONS, () => {
    const createGridDiv = (name: string, grid: DG.Grid) => {
      const header = ui.h1(name, 'chem-mmpa-transformation-tab-header');
      grid.root.prepend(header);
      return ui.splitV([
        ui.box(
          ui.divH([header, addToWorkspaceButton(grid.dataFrame, name, 'chem-mmpa-add-to-workspace-button')]),
          {style: {maxHeight: '30px'}},
        ),
        grid.root,
      ]);
    };

    return ui.splitV([
      createGridDiv('Fragments', this.transFragmentsGrid!),
      createGridDiv('Pairs', this.transPairsGrid!),
    ], {}, true);
  });
  tabs.addPane(MMP_NAMES.TAB_FRAGMENTS, () => {
    return tp.root;
  });
  tabs.addPane(MMP_NAMES.TAB_CLIFFS, () => {
    return cliffs;
  });
  const genTab = tabs.addPane(MMP_NAMES.TAB_GENERATION, () => {
    return this.generationsGrid!.root;
  });
  genTab.header.append(addToWorkspaceButton(this.generationsGrid!.dataFrame,
    'Generation', 'chem-mmpa-add-generation-to-workspace-button'));

  let refilter = true;
  tabs.onTabChanged.subscribe(() => {
    this.currentTab = tabs.currentPane.name;
    if (tabs.currentPane.name == MMP_NAMES.TAB_TRANSFORMATIONS) {
      this.enableFilters = true;
      this.refilterAllFragments(false);
    } else if (tabs.currentPane.name == MMP_NAMES.TAB_FRAGMENTS) {
      tabs.currentPane.content.append(mmpFilters.filtersDiv);
      this.refreshFilterAllFragments();
      this.enableFilters = false;
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
            this.transPairsGrid!.dataFrame, this.diffs!, this.parentTable!, this.rdkitModule!);
        }, 500);
      }
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
};

MatchedMolecularPairsViewer.prototype.fillAll =
function(mmpInput: MmpInput, palette: PaletteCodes,
  rules: MmpRules, diffs: Array<Float32Array>,
  linesIdxs: Uint32Array, transFragmentsGrid: DG.Grid, transPairsGrid: DG.Grid, generationsGrid: DG.Grid,
  tp: DG.Viewer, sp: DG.Viewer, mmpFilters: MmpFilters,
  linesEditor: ScatterPlotLinesRenderer, lines: ILineSeries, linesActivityCorrespondance: Uint32Array,
  rdkitModule: RDModule, gpuUsed: boolean) {
  this.rdkitModule = rdkitModule;

  this.parentTable = mmpInput.table;
  this.parentCol = mmpInput.molecules;
  this.colorPalette = palette;
  this.mmpRules = rules;

  this.diffs = diffs;
  this.transFragmentsGrid = transFragmentsGrid;
  this.transPairsGrid = transPairsGrid;
  this.generationsGrid = generationsGrid;
  this.lines = lines;
  this.linesActivityCorrespondance = linesActivityCorrespondance;
  this.linesIdxs = linesIdxs;
  this.calculatedOnGPU = gpuUsed;

  //transformations tab setup
  this.transFragmentsMask = DG.BitSet.create(this.transFragmentsGrid.dataFrame.rowCount);
  this.transPairsMask = DG.BitSet.create(this.transPairsGrid.dataFrame.rowCount);
  this.mmpView = DG.View.create();
  this.setupTransformationTab();

  //Cliffs tab setup
  this.mmpFilters = mmpFilters;
  this.cutoffMasks = new Array<DG.BitSet>(mmpFilters.activitySliderInputs.length);
  this.totalCutoffMask = DG.BitSet.create(this.parentTable.rowCount);
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
          transPairsGrid.dataFrame, diffs, mmpInput.table, rdkitModule, this.propPanelViewer!);
        this.lastSelectedPair = event.id;
        this.propPanelViewer!.fitHeaderToLabelWidth(100);
      }, 500);
    }
  });

  this.mmpView.append(tabs);

  const propertiesColumnsNames = this.parentTable.columns.names()
    .filter((name) => !name.startsWith('~'));
  this.propPanelViewer = new FormsViewer();
  this.propPanelViewer.dataframe = this.parentTable;
  this.propPanelViewer.columns = propertiesColumnsNames;
  this.propPanelViewer.inputClicked.subscribe(() => {
    setTimeout(() => {
      grok.shell.o = fillPairInfo(
        this.lastSelectedPair!, linesIdxs, linesActivityCorrespondance[this.lastSelectedPair!],
        transPairsGrid.dataFrame, diffs, mmpInput.table, rdkitModule, this.propPanelViewer!);
      this.propPanelViewer!.fitHeaderToLabelWidth(100);
    }, 500);
  });
};
