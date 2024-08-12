import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getRdKitModule} from '../../utils/chem-common-rdkit';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

import {MMP_NAMES, MMP_CONSTRICTIONS, MMP_ERRORS} from './mmp-constants';
import {MmpRules, MmpInput} from './mmp-constants';
import {getInverseSubstructuresAndAlign, PaletteCodes, getPalette} from './mmp-mol-rendering';
import {getMmpFrags, getMmpRules} from './mmp-fragments';
import {getMmpActivityPairsAndTransforms} from './mmp-pairs-transforms';
import {getMmpTrellisPlot} from './mmp-frag-vs-frag';
import {fillPairInfo, getMmpScatterPlot, runMmpChemSpace} from './mmp-cliffs';
import {getGenerations} from './mmp-generations';

import {debounceTime} from 'rxjs/operators';
import {getSigFigs} from '../../utils/chem-common';
import {drawMoleculeLabels} from '../../rendering/molecule-label';

import {ILineSeries, MouseOverLineEvent, ScatterPlotLinesRenderer}
  from '@datagrok-libraries/utils/src/render-lines-on-sp';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {FormsViewer} from '@datagrok-libraries/utils/src/viewers/forms-viewer';
import {getEmbeddingColsNames} from
  '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/reduce-dimensionality';
import $ from 'cash-dom';
import {getGPUDevice} from '@datagrok-libraries/math/src/webGPU/getGPUDevice';
import {getMmpFilters, MmpFilters} from './mmp-filters';

export class MatchedMolecularPairsViewer extends DG.JsViewer {
  static TYPE: string = 'MMP';

  molecules: string | null = null;
  moleculesCol: DG.Column | null = null;
  activities: string[] | null = null;
  activitiesCols: DG.ColumnList | null = null;
  fragmentCutoff: number | null;

  parentTable: DG.DataFrame | null = null;
  parentCol: DG.Column| null = null;
  mmpRules: MmpRules | null = null;
  mmpView: DG.View | null = null;
  enableFilters: boolean = true;
  colorPalette: PaletteCodes | null = null;
  //transformations tab objects
  transFragmentsGrid: DG.Grid | null = null;
  transFragmentsMask: DG.BitSet | null = null;
  transPairsGrid: DG.Grid | null = null;
  transPairsMask: DG.BitSet | null = null;
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
  //rdkit
  rdkitModule: RDModule | null = null;
  currentTab = '';
  lastSelectedPair: number | null = null;
  propPanelViewer: FormsViewer | null = null;
  mmpFilters: MmpFilters | null = null;
  calculatedOnGPU: boolean | null = null;

  constructor() {
    super();

    //properties
    this.molecules = this.string('molecules');
    this.activities = this.stringList('activities');
    this.fragmentCutoff = this.float('fragmentCutoff');
  }

  onTableAttached() {

  }

  onPropertyChanged(property: DG.Property | null): void {
    console.log(property!.name);
    super.onPropertyChanged(property);
    if (property?.name === 'molecules')
      this.moleculesCol = this.dataFrame.col(property!.get(this));
    if (property?.name === 'activities')
      this.activitiesCols = DG.DataFrame.fromColumns(this.dataFrame.columns.byNames(this.activities!)).columns;
    if (this.molecules && this.activities && this.fragmentCutoff) {
      this.render();
      return;
    }
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

      } finally {
        $(this.root).empty();
        if (this.mmpView)
          this.root.appendChild(this.mmpView!.root);
        else
          this.close();
        progressMMP.close();
      }

      // let tableInfo = this.parentTable!.getTableInfo();
      // await grok.dapi.tables.uploadDataFrame(this.parentTable!);
      // await grok.dapi.tables.save(tableInfo);
      // console.log(`${grok.dapi.root}/entities/${tableInfo.id}/token`);
      // const response = await fetch(`${grok.dapi.root}/entities/${tableInfo.id}/token`, {method: "POST"});
      //or without POST
      // const respText = await response.text();
    }
  }

  private setupTransformationTab(): void {
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
  }

  private setupFilters(mmpFilters: MmpFilters, linesActivityCorrespondance: Uint32Array, tp: DG.Viewer) {
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
      //TODO: implement filtering
    });
  }

  private setupCliffsTab(sp: DG.Viewer, mmpFilters: MmpFilters,
    linesEditor: ScatterPlotLinesRenderer): HTMLDivElement {
    sp.root.style.width = '100%';

    this.totalCutoffMask!.setAll(true);
    this.linesMask!.setAll(true);

    this.linesRenderer = linesEditor;

    this.refilterCliffs(mmpFilters.activitySliderInputs.map((si) => si.value),
      mmpFilters.activityActiveInputs.map((ai) => ai.value), false);

    return ui.box(sp.root);
  }

  private getTabs(tp: DG.Viewer, mmpFilters: MmpFilters, cliffs: HTMLDivElement): DG.TabControl {
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
  }

  private fillAll(mmpInput: MmpInput, palette: PaletteCodes,
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
  }

  async runMMP(mmpInput: MmpInput) {
    console.profile('MMP');
    //rdkit module

    const module = getRdKitModule();
    const moleculesArray = mmpInput.molecules.toList();

    const gpuCheck = await getGPUDevice();
    const gpu: boolean = !gpuCheck ? false : true;

    if (!gpu && mmpInput.molecules.length > MMP_CONSTRICTIONS.CPU) {
      grok.shell.error(MMP_ERRORS.FRAGMENTS_CPU);
      throw new Error(MMP_ERRORS.FRAGMENTS_CPU);
    } else if (mmpInput.molecules.length > MMP_CONSTRICTIONS.GPU) {
      grok.shell.error(MMP_ERRORS.FRAGMENTS_GPU);
      throw new Error(MMP_ERRORS.FRAGMENTS_GPU);
    }

    //initial calculations
    const fragsOut = await getMmpFrags(moleculesArray);
    const [mmpRules, allCasesNumber] = await getMmpRules(fragsOut, mmpInput.fragmentCutoff, gpu);

    const palette = getPalette(mmpInput.activities.length);
    //Transformations tab
    const {maxActs, diffs, meanDiffs, activityMeanNames,
      linesIdxs, transFragmentsGrid, transPairsGrid, lines, linesActivityCorrespondance} =
      getMmpActivityPairsAndTransforms(mmpInput, mmpRules, allCasesNumber, palette);

    //Fragments tab
    const tp = getMmpTrellisPlot(transFragmentsGrid, activityMeanNames, palette);

    //Cliffs tab
    const embedColsNames = getEmbeddingColsNames(mmpInput.table).map((it) => `~${it}`);

    const mmpFilters = getMmpFilters(mmpInput, maxActs, transFragmentsGrid.dataFrame.col(MMP_NAMES.PAIRS)!.stats.max);
    console.log(`created mmpa filters`);

    const sp = getMmpScatterPlot(mmpInput, embedColsNames);

    drawMoleculeLabels(mmpInput.table, mmpInput.molecules, sp as DG.ScatterPlotViewer, 20, 7, 100, 110);

    //running internal chemspace
    const linesEditor = runMmpChemSpace(mmpInput, sp, lines, linesIdxs, linesActivityCorrespondance,
      transPairsGrid.dataFrame, diffs, module, embedColsNames);


    const generationsGrid: DG.Grid =
      await getGenerations(mmpInput, moleculesArray, fragsOut, meanDiffs, transFragmentsGrid, activityMeanNames, gpu);


    this.fillAll(mmpInput, palette, mmpRules, diffs, linesIdxs, transFragmentsGrid, transPairsGrid, generationsGrid,
      tp, sp, mmpFilters, linesEditor, lines, linesActivityCorrespondance, module, gpu);

    console.profileEnd('MMP');
  }

  findSpecificRule(diffFromSubstrCol: DG.Column): [idxPairs: number, cases: number[]] {
    const idxParent = this.parentTable!.currentRowIdx;
    let idxPairs = -1;
    const cases: number[] = [];
    const idx = this.transFragmentsGrid!.table.currentRowIdx;
    if (idx !== -1) {
      const ruleSmi1 = this.transFragmentsGrid!.table.getCol(MMP_NAMES.FROM).get(idx);
      const ruleSmi2 = this.transFragmentsGrid!.table.getCol(MMP_NAMES.TO).get(idx);
      const ruleSmiNum1 = this.mmpRules!.smilesFrags.indexOf(ruleSmi1);
      const ruleSmiNum2 = this.mmpRules!.smilesFrags.indexOf(ruleSmi2);

      let counter = 0;

      for (let i = 0; i < this.mmpRules!.rules.length; i++) {
        const first = this.mmpRules!.rules[i].smilesRule1;
        const second = this.mmpRules!.rules[i].smilesRule2;
        for (let j = 0; j < this.mmpRules!.rules[i].pairs.length; j++) {
          if (ruleSmiNum1 == first && ruleSmiNum2 == second) {
            this.transPairsMask!.set(counter, true, false);
            if (diffFromSubstrCol.get(counter) === null)
              cases.push(counter);
            if (this.mmpRules!.rules[i].pairs[j].firstStructure == idxParent)
              idxPairs = counter;
          }
          counter++;
        }
      }
    }
    return [idxPairs, cases];
  }

  async recoverHighlights(
    cases: number[], diffFrom : DG.Column, diffTo: DG.Column,
    diffFromSubstrCol: DG.Column, diffToSubstrCol: DG.Column, rdkitModule: RDModule) : Promise<void> {
    const pairsFrom = Array<string>(cases.length);
    const pairsTo = Array<string>(cases.length);
    for (let i = 0; i < cases.length; i++) {
      pairsFrom[i] = diffFrom.get(cases[i]);
      pairsTo[i] = diffTo.get(cases[i]);
    }
    const {inverse1, inverse2, fromAligned, toAligned} =
      await getInverseSubstructuresAndAlign(pairsFrom, pairsTo, rdkitModule);
    for (let i = 0; i < cases.length; i++) {
      diffFrom.set(cases[i], fromAligned[i]);
      diffTo.set(cases[i], toAligned[i]);
      diffFromSubstrCol.set(cases[i], inverse1[i]);
      diffToSubstrCol.set(cases[i], inverse2[i]);
    }
  }

  /**
  * Prepares all the entities to show for selected pair.
  * @param {RDModule} rdkitModule RDkit module instance
  */
  async refreshPair(rdkitModule: RDModule) : Promise<void> {
    const progressBarPairs = DG.TaskBarProgressIndicator.create(`Refreshing pairs...`);

    this.transPairsMask!.setAll(false);
    const diffFromSubstrCol = this.transPairsGrid!.dataFrame.getCol(MMP_NAMES.STRUCT_DIFF_FROM_NAME);
    const diffToSubstrCol = this.transPairsGrid!.dataFrame.getCol(MMP_NAMES.STRUCT_DIFF_TO_NAME);
    const diffFrom = this.transPairsGrid!.dataFrame.getCol(MMP_NAMES.FROM);
    const diffTo = this.transPairsGrid!.dataFrame.getCol(MMP_NAMES.TO);

    const [idxPairs, cases] = this.findSpecificRule(diffFromSubstrCol);
    await this.recoverHighlights(cases, diffFrom, diffTo, diffFromSubstrCol, diffToSubstrCol, rdkitModule);

    this.transPairsGrid!.dataFrame.filter.copyFrom(this.transPairsMask!);

    this.transPairsGrid!.setOptions({
      pinnedRowColumnNames: [],
      pinnedRowValues: [],
    });
    if (idxPairs >= 0) {
      this.transPairsGrid!.setOptions({
        pinnedRowValues: [idxPairs.toString()],
        pinnedRowColumnNames: [MMP_NAMES.PAIRNUM],
      });
    }

    progressBarPairs.close();
  }

  refreshFilterAllFragments() {
    const consistsBitSet: DG.BitSet = DG.BitSet.create(this.transFragmentsGrid!.dataFrame.rowCount);
    consistsBitSet.setAll(true);
    this.transFragmentsGrid!.dataFrame.filter.copyFrom(consistsBitSet);
  }

  /**
  * Gets visible all fragment pairs according to selected molecule in the table.
  * @param {boolean} rowChanged flag if row was changed
  */
  refilterAllFragments(rowChanged: boolean) {
    let idxTrue = -1;
    if (rowChanged) {
      const idx = this.parentTable!.currentRowIdx;
      this.transFragmentsMask!.setAll(false);

      for (let i = 0; i < this.mmpRules!.rules.length; i++) {
        for (let j = 0; j < this.mmpRules!.rules[i].pairs.length; j++) {
          const fs = this.mmpRules!.rules[i].pairs[j].firstStructure;
          if (idx == fs) {
            if (idxTrue == -1)
              idxTrue = i;
            this.transFragmentsMask!.set(i, true, false);
            break;
          }
        }
      }
    }

    if (this.enableFilters) {
      this.transFragmentsGrid!.dataFrame.filter.copyFrom(this.transFragmentsMask!);
      this.transFragmentsGrid!.invalidate();
      if (rowChanged)
        this.transFragmentsGrid!.table.currentRowIdx = idxTrue;
    }
  }

  /**
  * Makes the dataset filtering according to user specified cutoffs and 'isActive' checkbox for each activity.
  * @param {number[]} cutoffs for each activity the corresponding cuttof is specified by user
  * @param {boolean} refilter flag if calculation is required before filtering
  */
  refilterCliffs(cutoffs: number[], isActiveVar: boolean[], refilter: boolean) {
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

  pinPair(idx: number) {
    const columns = this.transPairsGrid?.dataFrame.columns;
    const idxFrom: number = columns!.byName(MMP_NAMES.PAIRNUM_FROM).get(idx);
    const idxToTo: number = columns!.byName(MMP_NAMES.PAIRNUM_TO).get(idx);
    const molFrom = this.parentTable!.columns.byName(this.molecules!).get(idxFrom);
    const molTo = this.parentTable!.columns.byName(this.molecules!).get(idxToTo);
    const grid = grok.shell.tv.grid;
    const indexesPairs = this.transPairsMask?.getSelectedIndexes();
    const indexesAllFrom = indexesPairs?.map((ip) => columns!.byName(MMP_NAMES.PAIRNUM_FROM).get(ip));
    const indexesAllTo = indexesPairs?.map((ip) => columns!.byName(MMP_NAMES.PAIRNUM_TO).get(ip));
    indexesAllFrom?.forEach((i) => grid.dataFrame.selection.set(i, true));
    indexesAllTo?.forEach((i) => grid.dataFrame.selection.set(i, true));

    grid.setOptions({
      pinnedRowValues: [molFrom, molTo],
      pinnedRowColumnNames: [this.molecules!, this.molecules!],
    });
  }
}
